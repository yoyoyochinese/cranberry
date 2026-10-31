"""Compilation layer: compiled inner loop and the pure-Python reference loop (G01, G02). Owner: WP-F1.

* :func:`compile_loop` composes the ``kernel`` of every leaf (in pipeline order, honouring
  ``rate_div``) into one ``@njit`` inner loop over samples. The loop body is generated as flat
  source (one guarded call per leaf kernel, the kernels bound as globals) and compiled eagerly for
  C-contiguous arrays, so the leaf kernels inline into a single machine-code loop.
* :func:`run_reference` runs the same pipeline with the pure-Python ``__call__`` of each block.
* :func:`run_compiled` drives a compiled loop with the same piece / ``prepare`` protocol as
  :func:`run_reference`.
* VT-44: both must agree to 1e−12 (relative) on every slot after N steps.

``LoopFn`` signature: ``loop(S, Z, I, P, n_steps) -> None``: advances the state ``n_steps`` samples
in place, including the runner-owned slots ``n`` and ``t``, and returns nothing. A raw ``LoopFn``
call does not run ``prepare`` hooks; use :func:`run_compiled` when any leaf has one.

Per-sample order (both loops): with n = state.n, every leaf with ``n % rate_div == 0`` steps in
execution order (t = n / f_sim is visible to the blocks); then n ← n + 1 and t ← n / f_sim.

Pieces (``core.rates``): the drivers advance in pieces of at most ``PIECE_SAMPLES`` samples and call
``leaf.prepare(state, params, n0, n_samples)`` (when present) for every leaf before each piece.

Kernel reuse (G02): numba's on-disk cache does not apply to closures, so compiled loops are kept in
an in-process cache keyed by the *structure*: leaf classes, kinds, names, rate_divs, structural
params (runtime-param values stripped), the slot / param table of the layout and f_sim. A sweep over
runtime parameters or seeds with the same structure reuses the compiled loop (compile once, then
only ``Params.set`` / a new ``State``). This relies on the kernel rule of ``core.block``: kernels
capture only offsets, structural params, f_sim and rate_div.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import numba
import numpy as np

from cranberry.core.rates import PIECE_SAMPLES, should_step

if TYPE_CHECKING:
    from cranberry.core.block import Block
    from cranberry.core.layout import SlotLayout
    from cranberry.core.state import Params, State

__all__ = [
    "LoopFn",
    "NotCompilableError",
    "compile_loop",
    "structure_key",
    "clear_loop_cache",
    "loop_cache_size",
    "run_reference",
    "run_compiled",
    "is_compilable",
]

#: ``loop(S, Z, I, P, n_steps) -> None``
LoopFn = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int], None]

_LOOP_SIG = "void(float64[::1], complex128[::1], int64[::1], float64[::1], int64)"
_LOOP_CACHE: dict[str, LoopFn] = {}


class NotCompilableError(TypeError):
    """A leaf has no compiled kernel (``kernel is None``)."""


def _leaves(block: Block) -> list[Block]:
    from cranberry.core.pipeline import iter_leaves

    return list(iter_leaves(block))


def _layout_of(block: Block, layout: SlotLayout | None) -> SlotLayout:
    lay = layout if layout is not None else getattr(block, "layout", None)
    if lay is None:
        raise ValueError("compile_loop needs the layout: build the pipeline with core.pipeline.build or pass layout=")
    if not lay.frozen:
        raise ValueError("compile_loop needs a frozen (built) layout")
    return lay


def is_compilable(block: Block) -> bool:
    """True if every leaf of the bound pipeline has a kernel."""
    return all(leaf.kernel is not None for leaf in _leaves(block))


def _canon(obj: Any) -> Any:
    """Hashable, repr-stable canonical form of spec params (arrays hashed by content)."""
    if isinstance(obj, np.ndarray):
        return ("ndarray", obj.dtype.str, obj.shape, hashlib.sha256(np.ascontiguousarray(obj).tobytes()).hexdigest())
    if isinstance(obj, dict) or hasattr(obj, "items"):
        return ("map", tuple(sorted((str(k), _canon(v)) for k, v in obj.items())))
    if isinstance(obj, (list, tuple)):
        return ("seq", tuple(_canon(v) for v in obj))
    if isinstance(obj, float):
        return ("f", obj.hex())
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return _canon(obj.item())
    if callable(obj):
        return ("fn", getattr(obj, "__module__", "?"), getattr(obj, "__qualname__", "?"), id(obj))
    return (type(obj).__name__, repr(obj))


def structure_key(block: Block, layout: SlotLayout) -> str:
    """Hash of everything a compiled loop depends on (module docstring), excluding runtime values."""
    leaves = []
    for leaf in _leaves(block):
        params = dict(getattr(leaf, "params", {}) or {})
        for k in leaf.runtime_params:
            params.pop(k, None)
        leaves.append((type(leaf).__module__, type(leaf).__qualname__, leaf.kind, leaf.name, leaf.rate_div,
                       _canon(params)))
    slots = tuple((n, s.array, s.shape, s.ring_len, s.offset, s.head) for n in layout.slot_names()
                  for s in [layout.slot(n)])
    pars = tuple((n, layout.param(n).offset, layout.param(n).shape) for n in layout.param_names())
    key = repr((tuple(leaves), slots, pars, layout.f_sim_hz.hex(), tuple(sorted(layout.sizes.items()))))
    return hashlib.sha256(key.encode()).hexdigest()


def _generate(leaves: list[Block], n_off: int, t_off: int, f_sim: float) -> LoopFn:
    lines = ["def loop(S, Z, I, P, n_steps):", "    for _ in range(n_steps):", f"        n = I[{n_off}]"]
    glb: dict[str, Any] = {}
    for i, leaf in enumerate(leaves):
        glb[f"k{i}"] = leaf.kernel
        call = f"k{i}(S, Z, I, P)"
        if leaf.rate_div == 1:
            lines.append(f"        {call}")
        else:
            lines.append(f"        if n % {int(leaf.rate_div)} == 0:")
            lines.append(f"            {call}")
    lines += ["        n += 1", f"        I[{n_off}] = n", f"        S[{t_off}] = n / {float(f_sim)!r}"]
    src = "\n".join(lines) + "\n"
    exec(compile(src, "<cranberry.core.compile loop>", "exec"), glb)  # noqa: S102 (generated from offsets only)
    return numba.njit(_LOOP_SIG)(glb["loop"])


def compile_loop(block: Block, layout: SlotLayout | None = None, *, cache: bool = True) -> LoopFn:
    """Compile the bound pipeline ``block`` into a numba inner loop.

    Args:
        block: a built pipeline root (``core.pipeline.build``); every leaf must have a kernel.
        layout: its frozen layout (default: ``block.layout``, set by ``build``).
        cache: reuse a loop compiled earlier for the same structure (module docstring).

    Raises:
        NotCompilableError: some leaves have ``kernel is None`` (run them with :func:`run_reference`).
    """
    lay = _layout_of(block, layout)
    leaves = _leaves(block)
    missing = [leaf.name for leaf in leaves if leaf.kernel is None]
    if missing:
        raise NotCompilableError(f"leaves without a compiled kernel: {missing}; use run_reference")
    key = structure_key(block, lay) if cache else None
    if key is not None and key in _LOOP_CACHE:
        return _LOOP_CACHE[key]
    loop = _generate(leaves, lay.slot("n").offset, lay.slot("t").offset, lay.f_sim_hz)
    if key is not None:
        _LOOP_CACHE[key] = loop
    return loop


def clear_loop_cache() -> None:
    """Drop all cached compiled loops."""
    _LOOP_CACHE.clear()


def loop_cache_size() -> int:
    return len(_LOOP_CACHE)


def _pieces(n0: int, n_steps: int):
    if n_steps < 0:
        raise ValueError(f"n_steps must be >= 0, got {n_steps}")
    done = 0
    while done < n_steps:
        m = min(PIECE_SAMPLES, n_steps - done)
        yield n0 + done, m
        done += m


def _prepare(leaves: list[Block], state: State, params: Params, n0: int, m: int) -> None:
    for leaf in leaves:
        prep = getattr(leaf, "prepare", None)
        if prep is not None:
            prep(state, params, n0, m)


def run_reference(block: Block, state: State, params: Params, n_steps: int) -> State:
    """Advance ``state`` by ``n_steps`` samples with the pure-Python steps (updates n, t); returns ``state``."""
    leaves = _leaves(block)
    lay = state.layout
    n_off, t_off, f_sim = lay.slot("n").offset, lay.slot("t").offset, lay.f_sim_hz
    S, I = state.S, state.I
    root_div = block.rate_div
    for n0, m in _pieces(state.n, n_steps):
        _prepare(leaves, state, params, n0, m)
        for _ in range(m):
            n = int(I[n_off])
            if should_step(n, root_div):
                block(state, params)
            n += 1
            I[n_off] = n
            S[t_off] = n / f_sim
    return state


def run_compiled(block: Block, state: State, params: Params, n_steps: int, loop: LoopFn | None = None) -> State:
    """Advance ``state`` by ``n_steps`` samples with the compiled loop (``compile_loop(block)`` if not
    given), running ``prepare`` hooks per piece exactly as :func:`run_reference`; returns ``state``."""
    if loop is None:
        loop = compile_loop(block, state.layout)
    leaves = _leaves(block)
    has_prepare = any(getattr(leaf, "prepare", None) is not None for leaf in leaves)
    S, Z, I, P = state.S, state.Z, state.I, params.P
    for n0, m in _pieces(state.n, n_steps):
        if has_prepare:
            _prepare(leaves, state, params, n0, m)
        loop(S, Z, I, P, m)
    return state
