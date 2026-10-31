"""Delay buffers: integer-sample ring + fractional part (protected list; PLAN §2.6, G09). Owner: WP-F1.

The delay buffer is the only legal feedback path in a pipeline: a controller writes ``ctrl.<x>.u``,
a ``delay_buffer`` block delays it by τ, and the actuator reads the delayed slot. Loop delay is the
sum of its parts (controller pure delay from the device preset + ZOH h/2 when not discretised exactly
+ optical propagation n_g·L/c from actuator to detector + optional extra; PLAN §2.6).

τ = (k + μ)/f_sim with integer k ≥ 0 and fractional μ ∈ [0, 1); the fractional part is applied by
``method``:

* ``"linear"``  : linear interpolation between samples k and k+1;
* ``"thiran1"`` : first-order Thiran all-pass (unit magnitude, approximately flat group delay).

The same ring implementation serves the emission-time histories of ``laser.*`` and ``drv.*`` used by
``path.beams`` (beam_assembler, WP-F6).

Block kind ``delay_buffer`` (factory :func:`make_delay_buffer`), params (all required, D34):

    in      : str    slot name to delay (must be a scalar float64 slot)
    out     : str    slot name written (``delay.<name>`` ring + output; naming fixed by WP-F1)
    tau_s   : float  delay [s], >= 0
    method  : "linear" | "thiran1"
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from cranberry.core.block import Block
from cranberry.core.pipeline import DELAY_BUFFER_KIND
from cranberry.core.registry import register_block

__all__ = ["DELAY_METHODS", "DelayLine", "make_delay_buffer"]

DELAY_METHODS: tuple[str, ...] = ("linear", "thiran1")


class DelayLine:
    """Ring buffer with fractional read-out, usable from Python and (via its arrays) from kernels.

    Args:
        tau_s: delay [s] (>= 0).
        f_sim_hz: sample rate of the writer [Hz].
        method: one of :data:`DELAY_METHODS`.
    """

    def __init__(self, tau_s: float, f_sim_hz: float, method: str):
        raise NotImplementedError("core.delay.DelayLine is implemented by WP-F1")

    def push(self, x: float) -> None:
        """Write the newest sample."""
        raise NotImplementedError("WP-F1")

    def read(self) -> float:
        """Value delayed by tau_s relative to the newest sample."""
        raise NotImplementedError("WP-F1")


@register_block(DELAY_BUFFER_KIND)
def make_delay_buffer(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``delay_buffer`` block (see module docstring). Implemented by WP-F1."""
    raise NotImplementedError("block 'delay_buffer' (core.delay.make_delay_buffer) is implemented by WP-F1")
