"""``gate_schedule``: periodic lock / measurement window waveform (Phase 6 WP6-G, D85-2, D85-10; P2 Borealis windows,
P3 chopper, P8 intermittent loop; notes/lit-madsen2022.md §4.2, lit-nardelli2026.md §4.2, lit-luo2025.md §4.2).

Time loop (MOD), no inputs (reads ``t``). Writes, every step:

    gate.<name>.lock    probe transmission in [leak, 1]: 1 in the lock window, ``leak`` when off, linear edges
    gate.<name>.meas    measurement (detector / GBS) window, 0 or 1 (float slot, usable as a metric ``mask_slot``)
    gate.<name>.phase   optional (``write_phase: true``): 2π·x [rad], x = position in the current period

Schedule (cycles u = t/T + phase, T = ``period_s``; position x = u − n − δ_n/T in period n)::

    lock window   x ∈ [0, D_lock)                 D_lock = lock_fraction (also the probe duty of a chopper)
    edges         rising over edge_s from x = 0, falling over edge_s from x = D_lock (probe transmission ramps;
                  a chopper blade crossing the beam, an AOM rise)
    lock          = leak + (1 − leak)·g(x),  g ∈ [0, 1] the edge-shaped window
    meas          complement : 1 for τ = xT ≥ D_lock·T + edge_s (the whole probe-off time after the falling edge,
                               the Borealis measurement window)
                  anti_sync  : 1 for τ ∈ [D_lock·T + settle_s, D_lock·T + settle_s + duty_det·T) (P3: the detector
                               chopper in antiphase; settle_s = (1 − D_lock − duty_det)T/2 centres it in the probe-off
                               time, P3's 50 % / 25 % at 180°)

Timing jitter (``timing_jitter_s`` = σ_t): every period's schedule is shifted as a whole by δ_n ~ N(0, σ_t²), white
from period to period (chopper wheel / trigger jitter; P3 p.142 gives no model, so this is the simplest one with a
knob). δ_{n+1} is drawn when period n starts, so an early start (δ < 0) is honoured. The first period is unjittered.
One standard normal is drawn at every step (used only at period starts; G11 consume-and-discard).
Rounding: a period start is taken within 1e−9 cycles (so a period of an integer number of samples stays periodic in
samples); a window edge that falls exactly on a sample is resolved only to rounding (choose D_lock·T·f_sim off-integer
when the sample count matters).

Not modelled here: the probe flight to the detector (put a ``delay_buffer`` on ``gate.<name>.lock`` before the probe's
error-signal block, or blank the controller with ``gated_hold.settle_s``), the stray light of the gated probe (offline:
``window_stray_light``), the backscatter tail limit (``gated_loop_stability_check`` ``round_trip``).

Uses (WP6-G): (a) switch a probe off: ``liu_probe_interferometer`` ``in_gate: gate.<name>.lock`` scales the probe power
(fringe amplitude ∝ √lock, probe shot noise ∝ lock); (b) gate a detector / controller: ``gated_hold`` ``gate`` reads
``gate.<name>.lock`` (or ``.meas``); (c) metric mask: ``metric_i`` / ``metric_ii`` ``mask_slot: gate.<name>.meas``.

Params (D34, all required): ``detector_gate`` (``complement`` | ``anti_sync``), ``write_phase`` (bool); runtime
``period_s``, ``lock_fraction`` (0 < D_lock < 1; duty 1 = toggle ``gating`` off), ``phase`` (cycles, [0, 1)),
``edge_s``, ``timing_jitter_s`` (≤ period_s/10), ``leak`` (power extinction when off, [0, 1]); ``anti_sync`` also
``duty_det``, ``settle_s``. Toggles: ``gating`` (off: lock = meas = 1 always, i.e. no windows: the continuous lock and
continuous measurement of the ungated presets), ``jitter``, ``edges``, ``leak`` (off: that effect is 0).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["GateSchedule", "make_gate_schedule", "DETECTOR_GATES"]

KIND = "gate_schedule"
DETECTOR_GATES: dict[str, int] = {"complement": 0, "anti_sync": 1}
TWO_PI = 2.0 * math.pi
#: tolerance [cycles] of the period boundary: t/T is rounded, and a period of an integer number of samples must start
#: on its sample (else every few periods the window loses its first sample)
U_TOL = 1e-9


def _step(S, I, P, i_t, o_lock, o_meas, o_phase, st, pt, pd, pph, pe, pj, pl, pdd, pse, buf, idx,  # noqa: N803
          det_mode, gating, jit_on, edges_on, leak_on):
    w = 0.0
    if jit_on:
        w = S[buf + I[idx]]
        I[idx] += 1
    if not gating:
        S[o_lock] = 1.0
        S[o_meas] = 1.0
        if o_phase >= 0:
            S[o_phase] = 0.0
        return
    T = P[pt]
    u = S[i_t] / T + P[pph]
    n_cur, d_cur, d_next = S[st], S[st + 1], S[st + 2]
    if u >= n_cur + 1.0 + d_next / T - U_TOL:
        n_cur += 1.0
        d_cur = d_next
        d_next = P[pj] * w if jit_on else 0.0
        S[st], S[st + 1], S[st + 2] = n_cur, d_cur, d_next
    x = u - n_cur - d_cur / T
    if x < 0.0:
        x = 0.0
    tau = x * T
    toff = P[pd] * T
    edge = P[pe] if edges_on else 0.0
    if tau < toff:
        g = 1.0 if (edge <= 0.0 or tau >= edge) else tau / edge
    else:
        r = tau - toff
        g = 0.0 if (edge <= 0.0 or r >= edge) else 1.0 - r / edge
    leak = P[pl] if leak_on else 0.0
    S[o_lock] = leak + (1.0 - leak) * g
    if det_mode == 0:
        S[o_meas] = 1.0 if tau >= toff + edge else 0.0
    else:
        s0 = toff + P[pse]
        S[o_meas] = 1.0 if (tau >= s0 and tau < s0 + P[pdd] * T) else 0.0
    if o_phase >= 0:
        S[o_phase] = TWO_PI * x


class GateSchedule(KitBlock):
    """The ``gate_schedule`` block (module docstring)."""

    KIND = KIND
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("gating", "jitter", "edges", "leak")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "detector_gate", "write_phase", block=name)
        if params["detector_gate"] not in DETECTOR_GATES:
            raise ValueError(f"block {name!r}: detector_gate must be one of {sorted(DETECTOR_GATES)}")
        self.det_mode = DETECTOR_GATES[params["detector_gate"]]
        runtime = ["period_s", "lock_fraction", "phase", "edge_s", "timing_jitter_s", "leak"]
        if self.det_mode == 1:
            runtime += ["duty_det", "settle_s"]
        self.lock_slot, self.meas_slot = f"gate.{name}.lock", f"gate.{name}.meas"
        self.phase_slot = f"gate.{name}.phase" if bool(params["write_phase"]) else None
        self.state_slot = f"gate.{name}.state"
        writes = [self.lock_slot, self.meas_slot, self.state_slot] + ([self.phase_slot] if self.phase_slot else [])
        super().__init__(params, name, reads=("t",), writes=writes, runtime=runtime)
        self._validate()
        self.draws = DrawBuffer(self.name, 1, "standard_normal")
        if self.toggles["jitter"]:
            self.writes = (*self.writes, *self.draws.slot_names)
        self.slot_units = {self.lock_slot: "1", self.meas_slot: "1"}

    def _validate(self) -> None:
        p, n = self.params, self.name
        T, d = float(p["period_s"]), float(p["lock_fraction"])
        if not T > 0.0:
            raise ValueError(f"block {n!r}: period_s must be > 0")
        if not 0.0 < d < 1.0:
            raise ValueError(f"block {n!r}: lock_fraction must be in (0, 1) (duty 1: toggle gating off)")
        if not 0.0 <= float(p["phase"]) < 1.0:
            raise ValueError(f"block {n!r}: phase is in cycles, [0, 1)")
        if not float(p["edge_s"]) >= 0.0:
            raise ValueError(f"block {n!r}: edge_s must be >= 0")
        if not 0.0 <= float(p["timing_jitter_s"]) <= 0.1 * T:
            raise ValueError(f"block {n!r}: timing_jitter_s must be in [0, period_s/10]")
        if not 0.0 <= float(p["leak"]) <= 1.0:
            raise ValueError(f"block {n!r}: leak is a power fraction in [0, 1]")
        if self.det_mode == 1:
            dd, se = float(p["duty_det"]), float(p["settle_s"])
            if not (dd > 0.0 and se >= 0.0 and d + se / T + dd <= 1.0 + 1e-12):
                raise ValueError(f"block {n!r}: anti_sync needs duty_det > 0, settle_s >= 0 and "
                                 "lock_fraction + settle_s/period_s + duty_det <= 1")

    def slots(self):
        out = [(self.lock_slot, (), np.float64, "1", 0, 1.0), (self.meas_slot, (), np.float64, "1", 0, 0.0),
               (self.state_slot, (3,), np.float64, "1", 0, 0.0)]
        if self.phase_slot:
            out.append((self.phase_slot, (), np.float64, "rad", 0, 0.0))
        return out

    def init(self, layout, rng):
        super().init(layout, rng)
        if float(self.params["period_s"]) * layout.f_sim_hz / self.rate_div < 2.0:
            raise ValueError(f"block {self.name!r}: period_s must span >= 2 block steps")
        if self.toggles["jitter"]:
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self.toggles["jitter"]:
            self.draws.refill(state, n0, n_samples)

    def make_args(self, layout):
        if self.toggles["jitter"]:
            self.draws.bind(layout)
        a1 = self.det_mode == 1
        t = self.toggles
        return (self.o(layout, "t"), self.o(layout, self.lock_slot), self.o(layout, self.meas_slot),
                self.o(layout, self.phase_slot), self.o(layout, self.state_slot),
                self.p(layout, "period_s"), self.p(layout, "lock_fraction"), self.p(layout, "phase"),
                self.p(layout, "edge_s"), self.p(layout, "timing_jitter_s"), self.p(layout, "leak"),
                self.p(layout, "duty_det") if a1 else -1, self.p(layout, "settle_s") if a1 else -1,
                self.draws.buf_off, self.draws.idx_off, self.det_mode, t["gating"], t["jitter"], t["edges"], t["leak"])


@register_block(KIND)
def make_gate_schedule(params: Mapping[str, Any], *, name: str | None = None) -> GateSchedule:
    """Factory of the ``gate_schedule`` block (module docstring)."""
    return GateSchedule(params, name)
