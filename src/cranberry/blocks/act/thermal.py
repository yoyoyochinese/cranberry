"""``thermal_actuator``: TEC / spool temperature actuator that writes temperature (G08; T16 §3.4, T08, T10).

Model (lumped first-order lag, T16 §3.4)::

    ΔT_h(s) = K_T / (1 + s τ_th) · u(s),     ΔT_min ≤ ΔT_h ≤ ΔT_max

discretised exactly for an input held over each step: x ← a x + (1 − a) K_T u, a = exp(−h/τ_th)
(h = rate_div / f_sim; τ_th is a runtime value, a is computed per step). The clip acts on the
temperature state itself (toggle ``saturation``; ``act.<name>.sat``).

Output: ΔT_h [K] → ``act.<name>.dT`` and ``drv.<segment>.<driver>`` (normally driver ``tec``), copied to
every subsegment of ``sub`` (temperature is intensive). Phase and group delay follow in
``segment_phase`` from ``coef[segment][tec][beam]`` (c_T·L per subsegment, per axis): the TEC never
writes phase. "Phase written directly" (T16's σ c_T L_h ΔT model) is reproduced by a coefficient
preset (equal axis coefficients, zero group-delay coefficient), not by a code path (G08, D52).
A transport delay τ_tr is a ``delay_buffer`` in front of the block. The TEC usually runs decimated
(``rate_div``, 10–100 Hz, G10); then ``ring_len`` must be 0.

Params: ``in`` (V or A), ``segment``, ``driver``, ``n_sub``, ``sub``, ``ring_len``; runtime
``k_t_k_per_v`` (signed K_T [K/V]), ``tau_th_s``, ``dt_min_k``, ``dt_max_k``. Toggle ``saturation``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import DrvTarget, KitBlock, drv_write, drv_write_nb
from cranberry.core.block import require
from cranberry.core.registry import register_block

__all__ = ["ThermalActuator", "make_thermal_actuator"]


def _make_step(write):
    def step(S, I, P, i_in, i_out, i_sat, p_k, p_tau, p_lo, p_hi, h, sat_on,  # noqa: N803
             d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale):
        a = math.exp(-h / P[p_tau])
        x = a * S[i_out] + (1.0 - a) * P[p_k] * S[i_in]
        sat = 0
        if sat_on:
            if x > P[p_hi]:
                x, sat = P[p_hi], 1
            elif x < P[p_lo]:
                x, sat = P[p_lo], -1
        S[i_out] = x
        I[i_sat] = sat
        write(S, I, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale, x)

    return step


class ThermalActuator(KitBlock):
    KIND = "thermal_actuator"
    STEP = staticmethod(_make_step(drv_write))
    STEP_NB = numba.njit(_make_step(drv_write_nb))
    toggle_names = ("saturation",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", block=name)
        self.drv = DrvTarget(params, name, "K", split=False)
        self.a = {"dT": f"act.{name}.dT", "sat": f"act.{name}.sat"}
        super().__init__(params, name, reads=(params["in"],), writes=(*self.a.values(), self.drv.slot),
                         runtime=("k_t_k_per_v", "tau_th_s", "dt_min_k", "dt_max_k"))

    def init(self, layout, rng):
        self.drv.check_rate(self.rate_div)
        super().init(layout, rng)

    def slots(self):
        return [(self.a["dT"], (), np.float64, "K", 0, 0.0), (self.a["sat"], (), np.int64, "1", 0, 0),
                self.drv.slot_def()]

    def make_args(self, layout):
        return (self.o(layout, self.params["in"]), self.o(layout, self.a["dT"]), self.o(layout, self.a["sat"]),
                self.p(layout, "k_t_k_per_v"), self.p(layout, "tau_th_s"), self.p(layout, "dt_min_k"),
                self.p(layout, "dt_max_k"), self.h, self.toggles["saturation"], *self.drv.args(layout))


@register_block("thermal_actuator")
def make_thermal_actuator(params: Mapping[str, Any], *, name: str | None = None) -> ThermalActuator:
    """Factory of the ``thermal_actuator`` block (module docstring)."""
    return ThermalActuator(params, name)
