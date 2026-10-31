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
"""

from __future__ import annotations

import hashlib

import numpy as np

__all__ = ["stream_key", "child_seed_sequence", "stream", "substream", "root_seed_sequence"]

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
