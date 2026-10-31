"""Compilation layer: compiled inner loop and the pure-Python reference loop (G01, G02). Owner: WP-F1.

* :func:`compile_loop` composes the ``kernel`` of every leaf (in pipeline order, honouring
  ``rate_div``) into one ``@njit`` inner loop over samples. Composition is a kernel calling the
  captured leaf kernels in order (numba compiles calls to jit functions captured by a closure).
* :func:`run_reference` runs the same pipeline with the pure-Python ``__call__`` of each block.
* VT-44: both must agree to 1e−12 (relative) on every slot after N steps.

``LoopFn`` signature: ``loop(S, Z, I, P, n_steps) -> None``: advances the state ``n_steps`` samples
in place, including the runner-owned slots ``n`` and ``t`` (t = n / f_sim), and returns nothing.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from cranberry.core.block import Block
    from cranberry.core.state import Params, State

__all__ = ["LoopFn", "compile_loop", "run_reference"]

#: ``loop(S, Z, I, P, n_steps) -> None``
LoopFn = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int], None]


def compile_loop(block: Block) -> LoopFn:
    """Compile the bound pipeline ``block`` into a numba inner loop. Implemented by WP-F1."""
    raise NotImplementedError("core.compile.compile_loop is implemented by WP-F1")


def run_reference(block: Block, state: State, params: Params, n_steps: int) -> State:
    """Advance ``state`` by ``n_steps`` samples with the pure-Python steps (updates n, t). WP-F1."""
    raise NotImplementedError("core.compile.run_reference is implemented by WP-F1")
