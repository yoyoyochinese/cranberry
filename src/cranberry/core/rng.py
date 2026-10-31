"""Keyed random-number streams (G11, PLAN §2.11, §7.3).

Every randomness consumer gets its own ``numpy.random.Generator`` derived from the run's root seed
and a **stable hash of its name** (normally the block instance name, e.g. ``"lock1.pid"``):

    child SeedSequence = SeedSequence(entropy=root.entropy, spawn_key=root.spawn_key + key(name))
    key(name)          = 4 × uint32 words of SHA-256(name as UTF-8)       (:func:`stream_key`)
    generator          = Generator(PCG64(child))

Properties this guarantees (and the tests check):

* reproducible: same (root, name) → bit-identical draws, in any process, on any machine;
* independent of creation order and of which other streams exist — switching one noise source off
  does not change any other stream (needed for the metric-(ii) decomposition reruns, G13);
* chunk-length invariant as long as consumers draw sequentially from their own stream (G11).

Blocks must use only these streams, never the global ``numpy.random`` state (PLAN §7.3). The builder
(WP-F1, ``core.pipeline.build``) passes ``stream(root, block.name)`` to ``Block.init``. A block that
needs several independent sub-streams derives them with :func:`substream` (order-independent too).
The root entropy and each stream's key are what ``runs/<id>/seeds.json`` records (WP-F2).

Random numbers inside the compiled loop (WP-F1): kernels cannot call a ``Generator``, so a block
that needs per-step random numbers owns a :class:`DrawBuffer`. Before each piece of the time loop
(``Block.prepare``, ``core.rates``) the buffer is refilled in Python with **exactly** the number of
draws the block consumes in that piece, drawn sequentially from the block's own stream; the
kernel (and the pure-Python step) read the buffer one entry at a time. Because the stream is
consumed sequentially and never over-drawn, the realisation is bit-identical for every chunk or
piece length (G11, VT-F2), and identical between the compiled and the reference loop (VT-44).
"""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

import numpy as np

from cranberry.core.rates import PIECE_SAMPLES, max_steps_in, n_steps_in

if TYPE_CHECKING:
    from cranberry.core.layout import SlotLayout
    from cranberry.core.state import State

__all__ = [
    "stream_key",
    "child_seed_sequence",
    "stream",
    "substream",
    "root_seed_sequence",
    "DRAW_DISTS",
    "DrawBuffer",
]

_KEY_WORDS = 4


def stream_key(name: str) -> tuple[int, ...]:
    """Stable spawn-key words for ``name``: the first 4 little-endian uint32 words of SHA-256(name).

    Python's built-in ``hash`` is salted per process and must not be used for this.
    """
    if not isinstance(name, str) or not name:
        raise ValueError(f"stream name must be a non-empty str, got {name!r}")
    digest = hashlib.sha256(name.encode("utf-8")).digest()
    return tuple(int.from_bytes(digest[4 * i : 4 * i + 4], "little") for i in range(_KEY_WORDS))


def root_seed_sequence(root: int | np.random.SeedSequence) -> np.random.SeedSequence:
    """Normalise a root seed (non-negative int, or an existing SeedSequence) to a SeedSequence."""
    if isinstance(root, np.random.SeedSequence):
        return root
    if isinstance(root, (bool, np.bool_)) or not isinstance(root, (int, np.integer)):
        raise TypeError(f"root seed must be a non-negative int or SeedSequence, got {root!r}")
    if root < 0:
        raise ValueError(f"root seed must be >= 0, got {root}")
    return np.random.SeedSequence(int(root))


def child_seed_sequence(root: int | np.random.SeedSequence, name: str) -> np.random.SeedSequence:
    """The SeedSequence of stream ``name`` under ``root`` (spawn key = root key + :func:`stream_key`)."""
    ss = root_seed_sequence(root)
    return np.random.SeedSequence(entropy=ss.entropy, spawn_key=tuple(ss.spawn_key) + stream_key(name))


def stream(root: int | np.random.SeedSequence, name: str) -> np.random.Generator:
    """Independent, reproducible Generator (PCG64) for consumer ``name`` under root seed ``root``."""
    return np.random.Generator(np.random.PCG64(child_seed_sequence(root, name)))


def substream(parent: np.random.Generator, label: str) -> np.random.Generator:
    """Order-independent sub-stream ``label`` of a Generator produced by :func:`stream`.

    Equivalent to ``stream(root, name)`` with the key extended by ``stream_key(label)``; it does not
    consume draws from ``parent``.
    """
    seed_seq = getattr(parent.bit_generator, "seed_seq", None)
    if not isinstance(seed_seq, np.random.SeedSequence):
        raise TypeError("substream needs a Generator whose bit generator was seeded by a SeedSequence")
    child = np.random.SeedSequence(
        entropy=seed_seq.entropy, spawn_key=tuple(seed_seq.spawn_key) + stream_key(label)
    )
    return np.random.Generator(np.random.PCG64(child))


#: Distributions a :class:`DrawBuffer` can pre-draw: ``Generator`` method names (float64 output).
DRAW_DISTS: tuple[str, ...] = ("standard_normal", "random")


class DrawBuffer:
    """Pre-drawn random numbers for one block, readable from kernels (G11; see module docstring).

    Slots (registered by :meth:`init`): ``rng.<owner>[.<label>].buf`` float64 of shape
    ``(capacity * per_step,)`` and ``rng.<owner>[.<label>].idx`` int64 (next unread entry).

    Kernel use (offsets from :meth:`bind`)::

        w = S[buf_off + I[idx_off]]; I[idx_off] += 1        # one draw

    Python reference use: :meth:`take`.

    Args:
        owner: block name (slot names and error messages).
        per_step: draws consumed per step of the owning block (>= 1).
        dist: one of :data:`DRAW_DISTS`.
        label: optional sub-name when a block owns several buffers.
    """

    def __init__(self, owner: str, per_step: int, dist: str, label: str | None = None):
        if dist not in DRAW_DISTS:
            raise ValueError(f"dist must be one of {DRAW_DISTS}, got {dist!r}")
        if isinstance(per_step, bool) or int(per_step) != per_step or per_step < 1:
            raise ValueError(f"per_step must be an int >= 1, got {per_step!r}")
        self.owner = owner
        self.per_step = int(per_step)
        self.dist = dist
        prefix = f"rng.{owner}" if label is None else f"rng.{owner}.{label}"
        self.buf_name = f"{prefix}.buf"
        self.idx_name = f"{prefix}.idx"
        self.rng: np.random.Generator | None = None
        self.rate_div = 1
        self.capacity = 0
        self.buf_off = -1
        self.idx_off = -1

    @property
    def slot_names(self) -> tuple[str, str]:
        return (self.buf_name, self.idx_name)

    def init(self, layout: SlotLayout, rng: np.random.Generator, rate_div: int) -> None:
        """Register the slots for a block stepping every ``rate_div`` samples; keep ``rng``."""
        self.rng = rng
        self.rate_div = int(rate_div)
        self.capacity = max_steps_in(PIECE_SAMPLES, self.rate_div)
        layout.add(self.buf_name, (self.capacity * self.per_step,), np.float64, "1")
        layout.add(self.idx_name, (), np.int64, "1")

    def bind(self, layout: SlotLayout) -> None:
        self.buf_off = layout.slot(self.buf_name).offset
        self.idx_off = layout.slot(self.idx_name).offset

    def refill(self, state: State, n0: int, n_samples: int) -> None:
        """Draw exactly the entries consumed in base samples [n0, n0 + n_samples) and reset the index."""
        if self.rng is None:
            raise RuntimeError(f"DrawBuffer {self.buf_name!r} used before init")
        if n_samples > PIECE_SAMPLES:
            raise ValueError(f"piece of {n_samples} samples exceeds PIECE_SAMPLES={PIECE_SAMPLES}")
        count = n_steps_in(n0, n_samples, self.rate_div) * self.per_step
        buf = state.view(self.buf_name)
        getattr(self.rng, self.dist)(out=buf[:count])
        state.view(self.idx_name)[...] = 0

    def take(self, state: State) -> float:
        """Python reference: return the next entry and advance the index."""
        idx = state.view(self.idx_name)
        i = int(idx)
        idx[...] = i + 1
        return float(state.view(self.buf_name)[i])
