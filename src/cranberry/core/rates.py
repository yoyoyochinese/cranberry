"""Multirate stepping (G10, PLAN §2.7). Owner: WP-F1.

A block with ``rate_div = r`` is stepped only at samples with ``n % r == 0`` (n = state.n, the base
sample index at f_sim) and its written slots hold their values in between (zero-order hold). Its
own time step is therefore r / f_sim. Phase alignment is fixed: all decimated blocks step at n = 0.
"""

from __future__ import annotations

__all__ = ["should_step", "block_dt_s"]


def should_step(n: int, rate_div: int) -> bool:
    """True if a block with decimation ``rate_div`` steps at base sample ``n``."""
    return n % rate_div == 0


def block_dt_s(f_sim_hz: float, rate_div: int) -> float:
    """Time step of a block with decimation ``rate_div``: rate_div / f_sim [s]."""
    return rate_div / f_sim_hz
