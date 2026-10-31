"""Multirate stepping and time-step helpers (G10, PLAN §2.7, §6.4). Owner: WP-F1.

Hold semantics
    A block with ``rate_div = r`` is stepped only at base samples with ``n % r == 0`` (n = state.n,
    the sample index at f_sim) and its written slots hold their values in between (zero-order hold:
    nothing else writes them). Its own time step is therefore r / f_sim (:func:`block_dt_s`); blocks
    discretise themselves with that step. Phase alignment is fixed: all decimated blocks step at
    n = 0. Composites (Seq / Parallel / Stack) always have ``rate_div = 1``; decimation is a property
    of leaves.

Hold delay (rule, Phase 5 F1-3)
    A block held at ``rate_div = r`` adds **(r − 1)/(2 f_sim)** of delay to what it passes on (the
    phase-averaged group delay of the discrete sample-and-hold: an input seen at n = r·⌊n/r⌋ is used
    for n … n + r − 1, mean lag (r − 1)/2 base samples), not the r/(2 f_sim) of a continuous-time
    zero-order hold of width r/f_sim. The remaining half base sample 1/(2 f_sim) of the physical D/A
    hold appears **only** when the downstream continuous-time block is discretised with ``zoh``
    (whose DC group delay is the analogue one + h/2); with ``bilinear`` / ``tustin`` it does not. So a
    preset subtracting explicitly modelled delays from a loop budget (D61-2, ``tau_s = τ_L − 1/f_sim −
    explicit``) counts a held instrument stage as (r − 1)/(2 f_sim) plus 1/(2 f_sim) only if the next
    block is ``zoh``-discretised (e.g. Ref D: DSC1 ``dac`` at r = 2 into the ``zoh`` 914 filter, VT-24
    test 10). ``tests/unit/core/test_rates_draws_recorder.py::test_hold_delay_rule`` fixes both numbers.

Pieces
    The time loop advances in *pieces* of at most :data:`PIECE_SAMPLES` base samples (the runner's
    ``chunk_samples`` is split into pieces). Before each piece, the optional leaf hook
    ``prepare(state, params, n0, n_samples)`` runs in Python (open-loop work: pre-drawing random
    numbers, G11). Buffers filled by ``prepare`` therefore need room for at most
    ``max_steps_in(PIECE_SAMPLES, r)`` steps of their block. Pieces are an implementation detail:
    results do not depend on piece or chunk boundaries (VT-F2).

f_sim convergence hooks (PLAN §6.4: "f_sim ×2, ×4, compare")
    :func:`with_f_sim` returns a copy of a ``RunSpec`` with f_sim multiplied by an integer factor.
    Every physical parameter is in SI units, so blocks re-discretise themselves; blocks whose
    update rate is a device property (e.g. the DSC1 PID at 500 kHz) are listed in ``keep_rate`` so
    their ``rate_div`` is multiplied too and their physical update rate is unchanged.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable, Iterable
from typing import Any

__all__ = [
    "PIECE_SAMPLES",
    "should_step",
    "block_dt_s",
    "n_steps_in",
    "max_steps_in",
    "scale_rate_divs",
    "with_f_sim",
]

#: Maximum number of base samples advanced between two ``prepare`` calls (see module docstring).
PIECE_SAMPLES: int = 4096


def should_step(n: int, rate_div: int) -> bool:
    """True if a block with decimation ``rate_div`` steps at base sample ``n``."""
    return n % rate_div == 0


def block_dt_s(f_sim_hz: float, rate_div: int) -> float:
    """Time step of a block with decimation ``rate_div``: rate_div / f_sim [s]."""
    return rate_div / f_sim_hz


def n_steps_in(n0: int, n_samples: int, rate_div: int) -> int:
    """Number of steps a block with decimation ``rate_div`` takes in base samples [n0, n0 + n_samples)."""
    if n0 < 0 or n_samples < 0 or rate_div < 1:
        raise ValueError(f"need n0 >= 0, n_samples >= 0, rate_div >= 1; got {n0}, {n_samples}, {rate_div}")
    r = rate_div
    return (n0 + n_samples + r - 1) // r - (n0 + r - 1) // r


def max_steps_in(n_samples: int, rate_div: int) -> int:
    """Upper bound of :func:`n_steps_in` over all start samples: ceil(n_samples / rate_div)."""
    return -(-n_samples // rate_div)


def scale_rate_divs(pipeline_spec: Any, factor: int, keep_rate: Iterable[str] | Callable[[Any], bool] = ()) -> Any:
    """Copy of a pipeline spec where leaves selected by ``keep_rate`` have ``rate_div *= factor``.

    Args:
        pipeline_spec: a ``core.pipeline`` spec tree (BlockSpec / SeqSpec / ParallelSpec / StackSpec).
        factor: integer f_sim multiplier (>= 1).
        keep_rate: leaf names, or a predicate on ``BlockSpec``, of blocks whose physical update
            rate must stay fixed when f_sim is multiplied by ``factor``.
    """
    from cranberry.core.block import BlockSpec
    from cranberry.core.pipeline import ParallelSpec, SeqSpec, StackSpec

    if isinstance(factor, bool) or int(factor) != factor or factor < 1:
        raise ValueError(f"factor must be an int >= 1, got {factor!r}")
    if callable(keep_rate):
        pred = keep_rate
    else:
        names = set(keep_rate)
        pred = lambda b: b.name in names  # noqa: E731
    seen: set[str] = set()

    def walk(s: Any) -> Any:
        if isinstance(s, BlockSpec):
            if pred(s):
                seen.add(s.name)
                return dataclasses.replace(s, rate_div=s.rate_div * int(factor))
            return s
        if isinstance(s, (SeqSpec, ParallelSpec)):
            return dataclasses.replace(s, children=tuple(walk(c) for c in s.children))
        if isinstance(s, StackSpec):
            return dataclasses.replace(
                s, stages=tuple(walk(c) for c in s.stages), actuators=tuple(walk(c) for c in s.actuators)
            )
        raise TypeError(f"not a pipeline spec: {s!r}")

    out = walk(pipeline_spec)
    if not callable(keep_rate):
        missing = set(keep_rate) - seen
        if missing:
            raise KeyError(f"keep_rate names not found in the pipeline: {sorted(missing)}")
    return out


def with_f_sim(spec: Any, factor: int, keep_rate: Iterable[str] | Callable[[Any], bool] = ()) -> Any:
    """Copy of a ``RunSpec`` with ``run.f_sim_hz *= factor`` and :func:`scale_rate_divs` applied.

    ``duration_s``, ``chunk_samples`` and ``seed`` are unchanged (results never depend on the chunk
    length). Note that the random realisation changes with f_sim (more draws per second), so f_sim
    convergence is judged on statistics or on deterministic configurations.

    This operates on the *resolved* spec: expressions that depend on ``run.f_sim_hz`` (e.g. a loop
    delay ``tau_s = τ_L − 1/f_sim``) are NOT re-evaluated. Use ``runner.convergence(spec, ["f_sim"])``
    or ``runner.sweep.expand`` for that (they start from ``RunSpec.unresolved``, D64).
    """
    run = dataclasses.replace(spec.run, f_sim_hz=spec.run.f_sim_hz * int(factor))
    return dataclasses.replace(spec, run=run, pipeline=scale_rate_divs(spec.pipeline, factor, keep_rate))
