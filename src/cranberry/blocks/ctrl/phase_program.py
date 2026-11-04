"""``phase_program``: per-port program of lock set points over time (Phase 6 WP6-P, D85-6; P6 Zhong et al. 2021,
notes/lit-zhong2021.md §3.1, §4.3 N6-1).

Time loop (CTRL), no inputs (reads ``t``). A table of program values p_k [rad] (rows = settings j, columns = ports k)
is played over time: setting j is active from the first sample with t ≥ times_s[j] (to 1e−9 of a step); before
times_s[0] every port is 0 (unprogrammed: the lock acquires at its nominal point first). Per step (h = rate_div/f_sim)::

    goal_k = table[j, k]                    (0 before the first switch)
    mode step:  c_k = goal_k
    mode ramp:  c_k += clip(goal_k − c_k, ±slew_rad_per_s·h)       (slew-limited transition, c_k(0) = 0)
    s_k = sign·c_k                          (the setting in the lock's coordinate)
    output rf_phase:     rf.<name>.port<k>   = q·round(s_k/q)      (toggle quantization, q = quantum_rad > 0: the DDS
                                                                    phase resolution; Moku:Delta 1e−6° = 1.745e−8 rad,
                                                                    notes/fact-p6-thermo-optic-moku.md §2) else s_k
    output dc_setpoint:  ctrl.<name>.vset<k> = v_fringe_v·sin(s_k)  (DC side-fringe set point, below; no quantisation)
    ctrl.<name>.target<k> = c_k             (the commanded program value: the set point that metric_ii setpoint_slot and
                                             lock_supervisor slip_center_slot follow)
    ctrl.<name>.index     = j               (−1 before the first switch)

Wiring (P6, lit-zhong2021 §3.1-2): ``rf_phase`` feeds a ``demodulator``'s ``in_program`` (φ = φ_d + rf). With the
heterodyne Liu probe, I = K cos(X + φ), the loop holds X = π/2 − φ: a program step moves the lock point by −Δrf, one to
one and independent of K (probe power, visibility, PD gain). ``sign = −1`` then makes the lock coordinate X − π/2 follow
+c_k (F_common: coupler_sign +1, the stretcher moves the probe and the squeezed light alike, so Θ_k follows c_k − c_0).
``dc_setpoint`` is the comparison (a bad configuration, D34): an ``adder`` subtracts vset from the DC fringe v = K cos X,
the loop holds K cos X = v_fringe_v·sin(s_k); with sign = −1 and v_fringe_v = K_cal that is X = π/2 + c_k only while
K = K_cal (power / visibility drift moves it: X = arccos(−(K_cal/K) sin c_k)) and the slope falls as cos c_k.

Table (structural): ``source: table`` with ``times_s`` (increasing, ≥ 0) and ``table_rad`` (len(times_s) rows ×
len(ports) columns); or ``source: random`` with ``n_settings``, ``t_first_s``, ``dwell_s`` (switching interval: times
t_first + j·dwell), ``step_rad`` and ``seed``: ``step_rad: null`` draws every value independently, uniform in
[−π/2, π/2] (Zhong's 30 random settings, Fig. 4(a), mapped to quadrature angles); a number 0 < step ≤ π/2 is a random
walk p_j = p_{j−1} ± step per port (sign drawn; reversed when it would leave [−π/2, π/2]), so every switch of every
port is a step of exactly ``step_rad`` (step-size axis). The random table is a deterministic function of ``seed``
(configuration, not noise: no draws at run time; ``table`` attribute).

Params (D34, all required): ``ports`` (list of port indices), ``source`` (+ its keys), ``mode`` (step | ramp),
``output`` (rf_phase | dc_setpoint), ``sign`` (±1); runtime ``slew_rad_per_s`` (ramp), ``quantum_rad`` (rf_phase),
``v_fringe_v`` (dc_setpoint). Toggles: ``enabled`` (off: every output 0, index −1), ``quantization`` (rf_phase only).
Slots: the outputs above, ``ctrl.<name>.x`` = [c_0 … c_{K−1}, j + 1] (internal).
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

__all__ = ["PhaseProgram", "make_phase_program", "random_table", "MODES", "OUTPUTS"]

KIND = "phase_program"
MODES: dict[str, int] = {"step": 0, "ramp": 1}
OUTPUTS: dict[str, int] = {"rf_phase": 0, "dc_setpoint": 1}
HALF_PI = 0.5 * math.pi


def random_table(n_settings: int, n_ports: int, step_rad: float | None, seed: int) -> np.ndarray:
    """(n_settings, n_ports) program table (module docstring: uniform for step None, ±step bounded walk otherwise)."""
    rng = np.random.default_rng(int(seed))
    if step_rad is None:
        return rng.uniform(-HALF_PI, HALF_PI, (n_settings, n_ports))
    out = np.zeros((n_settings, n_ports))
    prev = np.zeros(n_ports)
    for j in range(n_settings):
        s = np.where(rng.random(n_ports) < 0.5, -1.0, 1.0)
        nxt = prev + s * step_rad
        nxt = np.where(np.abs(nxt) > HALF_PI + 1e-12, prev - s * step_rad, nxt)
        out[j], prev = nxt, nxt
    return out


def _step(S, I, P, i_t, o_out, o_tgt, o_idx, st, times, table, sign, ps, pq, pv, h, mode, output,  # noqa: N803
          enabled, quant):
    n = o_out.shape[0]
    if not enabled:
        for k in range(n):
            S[o_out[k]] = 0.0
            S[o_tgt[k]] = 0.0
        S[o_idx] = -1.0
        return
    t = S[i_t]
    j = int(S[st + n]) - 1  # stored as j + 1 (the state slot starts at 0)
    while j + 1 < times.shape[0] and t >= times[j + 1] - 1e-9 * h:
        j += 1
    S[st + n] = j + 1
    S[o_idx] = j
    for k in range(n):
        goal = table[j, k] if j >= 0 else 0.0
        if mode == 0:
            c = goal
        else:
            c = S[st + k]
            lim = P[ps] * h
            d = goal - c
            c += lim if d > lim else (-lim if d < -lim else d)
        S[st + k] = c
        S[o_tgt[k]] = c
        s = sign * c
        if output == 0:
            if quant and P[pq] > 0.0:
                s = P[pq] * math.floor(s / P[pq] + 0.5)
            S[o_out[k]] = s
        else:
            S[o_out[k]] = P[pv] * math.sin(s)


class PhaseProgram(KitBlock):
    """The ``phase_program`` block (module docstring)."""

    KIND = KIND
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("enabled", "quantization")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "ports", "source", "mode", "output", "sign", block=name)
        for key, opts in (("mode", MODES), ("output", OUTPUTS)):
            if params[key] not in opts:
                raise ValueError(f"block {name!r}: {key} must be one of {sorted(opts)}")
        self.mode, self.output = MODES[params["mode"]], OUTPUTS[params["output"]]
        self.ports = [int(k) for k in params["ports"]]
        if not self.ports or len(set(self.ports)) != len(self.ports):
            raise ValueError(f"block {name!r}: ports must be a non-empty list of distinct port indices")
        if float(params["sign"]) not in (1.0, -1.0):
            raise ValueError(f"block {name!r}: sign must be +1 or -1")
        self.times, self.table = self._table(params, name)
        runtime = (("slew_rad_per_s",) if self.mode == 1 else ()) + (
            ("quantum_rad",) if self.output == 0 else ("v_fringe_v",))
        out = "rf.{n}.port{k}" if self.output == 0 else "ctrl.{n}.vset{k}"
        self.out_slots = [out.format(n=name, k=k) for k in self.ports]
        self.target_slots = [f"ctrl.{name}.target{k}" for k in self.ports]
        self.index_slot, self.state_slot = f"ctrl.{name}.index", f"ctrl.{name}.x"
        super().__init__(params, name, reads=("t",),
                         writes=(*self.out_slots, *self.target_slots, self.index_slot, self.state_slot), runtime=runtime)
        if self.mode == 1 and not float(params["slew_rad_per_s"]) > 0.0:
            raise ValueError(f"block {name!r}: slew_rad_per_s must be > 0")
        if self.output == 0 and not float(params["quantum_rad"]) >= 0.0:
            raise ValueError(f"block {name!r}: quantum_rad must be >= 0")
        unit = "rad" if self.output == 0 else "V"
        self.slot_units = {s: unit for s in self.out_slots} | {s: "rad" for s in self.target_slots}

    def _table(self, p: Mapping[str, Any], name: str) -> tuple[np.ndarray, np.ndarray]:
        if p["source"] == "table":
            require(p, "times_s", "table_rad", block=name)
            times = np.asarray(p["times_s"], dtype=float).ravel()
            table = np.asarray(p["table_rad"], dtype=float)
        elif p["source"] == "random":
            require(p, "n_settings", "t_first_s", "dwell_s", "step_rad", "seed", block=name)
            n = int(p["n_settings"])
            if n < 1 or not float(p["dwell_s"]) > 0.0:
                raise ValueError(f"block {name!r}: random needs n_settings >= 1 and dwell_s > 0")
            st = p["step_rad"]
            if st is not None and not 0.0 < float(st) <= HALF_PI:
                raise ValueError(f"block {name!r}: step_rad must be null or in (0, π/2]")
            times = float(p["t_first_s"]) + float(p["dwell_s"]) * np.arange(n)
            table = random_table(n, len(self.ports), None if st is None else float(st), int(p["seed"]))
        else:
            raise ValueError(f"block {name!r}: source must be 'table' or 'random'")
        if table.shape != (times.size, len(self.ports)) or times.size == 0:
            raise ValueError(f"block {name!r}: table_rad must have shape (len(times_s), len(ports)) = "
                             f"({times.size}, {len(self.ports)}), got {table.shape}")
        if np.any(times < 0.0) or np.any(np.diff(times) <= 0.0) or not np.all(np.isfinite(table)):
            raise ValueError(f"block {name!r}: times_s must be >= 0 and strictly increasing, table_rad finite")
        return times, table

    def slots(self):
        unit = "rad" if self.output == 0 else "V"
        return ([(s, (), np.float64, unit, 0, 0.0) for s in self.out_slots]
                + [(s, (), np.float64, "rad", 0, 0.0) for s in self.target_slots]
                + [(self.index_slot, (), np.float64, "1", 0, -1.0),
                   (self.state_slot, (len(self.ports) + 1,), np.float64, "1", 0, 0.0)])

    def make_args(self, layout):
        return (self.o(layout, "t"), np.array([self.o(layout, s) for s in self.out_slots], np.int64),
                np.array([self.o(layout, s) for s in self.target_slots], np.int64), self.o(layout, self.index_slot),
                self.o(layout, self.state_slot), self.times.copy(), np.ascontiguousarray(self.table), float(self.params["sign"]),
                self.p(layout, "slew_rad_per_s") if self.mode == 1 else -1,
                self.p(layout, "quantum_rad") if self.output == 0 else -1,
                self.p(layout, "v_fringe_v") if self.output == 1 else -1,
                float(self.h), self.mode, self.output, self.toggles["enabled"], self.toggles["quantization"])


@register_block(KIND)
def make_phase_program(params: Mapping[str, Any], *, name: str | None = None) -> PhaseProgram:
    """Factory of the ``phase_program`` block (module docstring)."""
    return PhaseProgram(params, name)
