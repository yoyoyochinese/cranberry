"""``gated_hold``: pass an error signal only inside a gate window, hold / zero / noise outside (Phase 6 WP6-G, D85-2;
P2 controller hold, P3 detector gating, P8 intermittent loop; notes/lit-madsen2022.md §4.2 N1, lit-nardelli2026.md
§4.2-2, lit-luo2025.md §4.2).

Time loop (CTRL), placed between the error signal (e.g. ``v.adc1``) and the controller. Step (h = rate_div/f_sim)::

    open   = gate > 0.5                                   (gate: e.g. gate.<name>.lock of a ``gate_schedule``)
    k      = steps since the last change of ``open`` (0 at the step of the change)
    active = open and k·h ≥ settle_s                      (settle_s blanking after each gate edge: the error is
                                                           not trusted until the probe has reached the detector,
                                                           n_g L/c, and the edge has passed)
    out    = in                                  if active (and ``held`` ← in)
           = 0                                   mode zero       (an integrator-only PID then holds its output
                                                                  exactly; with K_p ≠ 0 the P term drops to 0)
           = held                                mode hold       (last active input: the controller keeps acting on
                                                                  the last error, e.g. the integrator keeps ramping)
           = noise_v_rms · w, w ~ N(0, 1)        mode noise_only (the detector's electronic noise alone, no fringe:
                                                                  what the ADC sees with the probe off)

``held`` starts at 0. Params (D34, all required): ``in``, ``gate``, ``out``, ``unit`` (of ``out``), ``mode``
(``zero`` | ``hold`` | ``noise_only``); runtime ``settle_s`` (≥ 0) and, for ``noise_only``, ``noise_v_rms`` (per-step
standard deviation, white to the block's Nyquist frequency; a density converts as √(S_v f_b/2)). Toggle ``gating``
(off: out = in, identity). ``noise_only`` draws one standard normal per step (consumed also while active, G11).
Slots: ``out``; ``ctrl.<name>.x`` = [held, k, open of the previous step].
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["GatedHold", "make_gated_hold", "HOLD_MODES"]

KIND = "gated_hold"
HOLD_MODES: dict[str, int] = {"zero": 0, "hold": 1, "noise_only": 2}


def _step(S, I, P, i_in, i_gate, o_out, st, ps, pn, h, mode, gating, buf, idx):  # noqa: N803
    x = S[i_in]
    w = 0.0
    if mode == 2:
        w = S[buf + I[idx]]
        I[idx] += 1
    if not gating:
        S[o_out] = x
        return
    is_open = S[i_gate] > 0.5
    k = 0.0 if is_open != (S[st + 2] > 0.5) else S[st + 1] + 1.0
    S[st + 1] = k
    S[st + 2] = 1.0 if is_open else 0.0
    if is_open and k >= P[ps] / h - 1e-9:
        S[st] = x
        S[o_out] = x
    elif mode == 0:
        S[o_out] = 0.0
    elif mode == 1:
        S[o_out] = S[st]
    else:
        S[o_out] = P[pn] * w


class GatedHold(KitBlock):
    """The ``gated_hold`` block (module docstring)."""

    KIND = KIND
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("gating",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "in", "gate", "out", "unit", "mode", block=name)
        if params["mode"] not in HOLD_MODES:
            raise ValueError(f"block {name!r}: mode must be one of {sorted(HOLD_MODES)}")
        self.mode = HOLD_MODES[params["mode"]]
        runtime = ("settle_s", "noise_v_rms") if self.mode == 2 else ("settle_s",)
        self.state_slot = f"ctrl.{name}.x"
        super().__init__(params, name, reads=(params["in"], params["gate"]), writes=(params["out"], self.state_slot),
                         runtime=runtime)
        if not float(params["settle_s"]) >= 0.0:
            raise ValueError(f"block {name!r}: settle_s must be >= 0")
        self.draws = DrawBuffer(self.name, 1, "standard_normal")
        if self.mode == 2:
            self.writes = (*self.writes, *self.draws.slot_names)
        self.slot_units = {params["out"]: str(params["unit"])}

    def slots(self):
        return [(self.params["out"], (), np.float64, str(self.params["unit"]), 0, 0.0),
                (self.state_slot, (3,), np.float64, "1", 0, 0.0)]

    def init(self, layout, rng):
        super().init(layout, rng)
        if self.mode == 2:
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self.mode == 2:
            self.draws.refill(state, n0, n_samples)

    def make_args(self, layout):
        if self.mode == 2:
            self.draws.bind(layout)
        return (self.o(layout, self.params["in"]), self.o(layout, self.params["gate"]),
                self.o(layout, self.params["out"]), self.o(layout, self.state_slot), self.p(layout, "settle_s"),
                self.p(layout, "noise_v_rms") if self.mode == 2 else -1, float(self.h), self.mode,
                self.toggles["gating"], self.draws.buf_off, self.draws.idx_off)


@register_block(KIND)
def make_gated_hold(params: Mapping[str, Any], *, name: str | None = None) -> GatedHold:
    """Factory of the ``gated_hold`` block (module docstring)."""
    return GatedHold(params, name)
