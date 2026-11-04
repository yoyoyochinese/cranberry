"""``step_response``: per-switch step response of a programmed lock (Phase 6 WP6-P, D85-6; P6, notes/lit-zhong2021.md
§4.3 N6-3).

Measurement stage (hook ``measure.step_response``, after the time loop; no-op kernel). Reads the recorded series and
the switch table of a ``phase_program`` block (``program``). For every switch j (time t_j, segment [t_j, t_{j+1}) or
to the end of the record) and every programmed port k with goal g_jk = −sign·table_jk (where the program moves the
lock coordinate, the program's ``ctrl.<name>.target<k>``; g_{−1,k} = 0) and step Δ_jk = g_jk − g_{j−1,k}:

    y_k(t)     position_sign · positions[k] (the per-port lock coordinate, e.g. lock.L<k>.err = X_k − π/2, which a
               locked loop holds at the program value mod 2π), error e = wrap_2π(y − g_jk) (fringe slips: ``slips``)
    rise_s     time from 10 % to 90 % of the step ((u − g_{j−1})/Δ first ≥ 0.1, first ≥ 0.9); NaN for Δ = 0
    overshoot  max(0, max((u − g_jk)·sign Δ))/|Δ|; NaN for Δ = 0
    settle_s   from t_j to the start of the final run of samples with |e| ≤ tol; NaN if the segment ends outside
    bias_rad   mean of e over the last half of the segment (static programming bias); sigma_rad its std
    act_peak_frac  max |u|/range over the segment (``actuators`` / ``actuator_range``); act_step u_last − u_first
    slips      increase of the slip counter (``slips``, e.g. lock.sup<k>.n_slip) over the segment

Rise and overshoot read the PATH u(t): y unwrapped in time with the coordinate's period (2π, or π for Θ) and placed on
the branch nearest g_{j−1} at the segment start, so a step |Δ| ≥ period/4 (Θ: ≥ π/2, review P-2) is measured along the
way the lock actually went instead of folded around the target. A segment is ``ambiguous`` (rise / overshoot NaN) when
the record cannot resolve that path: a jump > period/4 between two recorded samples (decimation too coarse for the
lock, a slip) or a path ending > period/4 from g_jk (it went the other way round, or slipped: the quadrature may be
right mod period, the step was not taken as programmed). Settle and bias use the wrapped error e (unaffected).

Settling tolerance (review P-4): tol = max(``tol_rad``, ``tol_sigma``·σ_k), σ_k the std of the port's locked error before
the first switch (the second half of [0, t_first), the window of the program zero below). On a long port the lock
noise alone is of the order of a 10 mrad ``tol_rad``, so a fixed tolerance makes the settling time the segment length
or NaN, set by noise, not by the step; ``tol_sigma`` (e.g. 3) keeps the settling time a property of the transient.
``tol_sigma`` 0 = the fixed ``tol_rad`` (pre-review behaviour). The tolerances used are in ``tol_used_rad`` /
``theta_tol_used_rad``.

The same rise / overshoot / settle / bias, prefixed ``theta_``, for metric (ii) (``metric``: a ``metric_ii`` block,
or null): Θ_k(t) − Θ_set,k = m_k(t) + S[slot_k] − S[slot_0] (its ``setpoint_slot``s; the metric slot m_k holds the
deviation from the followed set point), error wrapped with period π (quadrature angle), against the programmed
difference g_jk − g_j0 (a port not in the program counts 0). This is the cross-port check of Θ_k − Θ_0 against the
program on the squeezed light. The program's zero (the static Θ_k − Θ_set,k without program, a per-port constant
calibrated once, lit-zhong2021 §3.1-4) is the circular mean (period π) over the second half of [0, t_first), restricted
to the metric's ``mask_slot`` > 0.5 when it has one (review P-6), and is subtracted first (``theta_zero_rad``).
Cross-port consistency on the lock coordinates: for every switch and every pair of programmed ports (k, l) with
the same non-zero step, ``cross_d_settle_s`` / ``cross_d_rise_s`` = (k) − (l); for every pair, ``pair_bias_rad`` =
mean over the last half of (y_k − y_l) − (g_jk − g_jl).

The resolution is the record's (``measure.record.decim``); the runner adds every slot this block reads to the record.
The record must be point-sampled: an ``antialias`` boxcar / fir filter averages the wrapped lock coordinates and metric
slots across their ±π (±π/2) folds and is refused (ValueError, review P-6).
Params (D34, all required): ``program``, ``positions`` ({k: slot}), ``position_sign`` (±1), ``metric`` (name or null),
``actuators`` ({k: slot}, may be {}), ``actuator_range`` ({k: number > 0} for the same keys), ``slips`` ({k: slot},
may be {}), ``tol_rad`` (> 0), ``tol_sigma`` (≥ 0).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block

__all__ = ["StepResponse", "make_step_response", "segment_response"]

KIND = "step_response"


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


def _keyed(spec: Any, what: str) -> dict[int, Any]:
    if not isinstance(spec, Mapping):
        raise ValueError(f"{what} must be a mapping {{k: ...}}")
    return {int(k): v for k, v in spec.items()}


def _wrap(x: np.ndarray, period: float) -> np.ndarray:
    """x folded to (−period/2, period/2]."""
    h = 0.5 * period
    return h - np.mod(h - x, period)


def unwrap_path(y: np.ndarray, g_prev: float, period: float) -> np.ndarray:
    """y unwrapped in time with ``period``, on the branch whose first sample is nearest ``g_prev`` (module docstring)."""
    if y.size == 0:
        return y.copy()
    w = 2.0 * math.pi / period
    u = np.unwrap(w * np.asarray(y, dtype=float)) / w
    return u + period * np.round((g_prev - u[0]) / period)


def segment_response(t: np.ndarray, y: np.ndarray, t0: float, g_prev: float, g: float, tol: float,
                     path: np.ndarray | None = None, period: float | None = None) -> dict[str, float]:
    """Rise / overshoot / settle / bias of one segment (module docstring); ``t``, ``y`` the segment's samples (y folded
    around g for settle / bias). ``path`` (default y): the unwrapped trajectory for rise / overshoot; with ``period`` a
    path that jumps or ends more than period/4 away marks the segment ambiguous (rise / overshoot NaN)."""
    nan = math.nan
    out = dict(rise_s=nan, overshoot=nan, settle_s=nan, bias_rad=nan, sigma_rad=nan, ambiguous=False)
    if t.size == 0:
        return out
    e = y - g
    tail = e[t.size // 2:]
    out["bias_rad"], out["sigma_rad"] = float(np.mean(tail)), float(np.std(tail))
    outside = np.abs(e) > tol
    if not outside[-1]:
        bad = np.flatnonzero(outside)
        out["settle_s"] = float(t[bad[-1] + 1] - t0) if bad.size else 0.0
    u = y if path is None else path
    if period is not None and u.size and (abs(u[-1] - g) > 0.25 * period
                                          or (u.size > 1 and np.max(np.abs(np.diff(u))) > 0.25 * period)):
        out["ambiguous"] = True
        return out
    d = g - g_prev
    if abs(d) > 1e-12:
        r = (u - g_prev) / d
        i10, i90 = np.flatnonzero(r >= 0.1), np.flatnonzero(r >= 0.9)
        if i10.size and i90.size:
            out["rise_s"] = float(t[i90[0]] - t[i10[0]])
        out["overshoot"] = float(max(0.0, np.max((u - g) * math.copysign(1.0, d))) / abs(d))
    return out


class StepResponse(BaseBlock):
    """``step_response`` (module docstring)."""

    KIND = KIND

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "program", "positions", "position_sign", "metric", "actuators", "actuator_range", "slips",
                "tol_rad", "tol_sigma", block=name)
        super().__init__(KIND, name, params)
        self.positions = _keyed(params["positions"], f"block {name!r}: positions")
        self.actuators = _keyed(params["actuators"], f"block {name!r}: actuators")
        self.act_range = {k: float(v) for k, v in _keyed(params["actuator_range"], f"block {name!r}: actuator_range").items()}
        self.slips = _keyed(params["slips"], f"block {name!r}: slips")
        if set(self.act_range) != set(self.actuators) or any(not r > 0.0 for r in self.act_range.values()):
            raise ValueError(f"block {name!r}: actuator_range needs a number > 0 for every actuators key")
        if float(params["position_sign"]) not in (1.0, -1.0):
            raise ValueError(f"block {name!r}: position_sign must be +1 or -1")
        if not float(params["tol_rad"]) > 0.0 or not float(params["tol_sigma"]) >= 0.0:
            raise ValueError(f"block {name!r}: tol_rad must be > 0 and tol_sigma >= 0")

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def _blocks(self, blocks: Mapping[str, Any]) -> tuple[Any, Any]:
        prog = blocks.get(self.params["program"])
        if prog is None or prog.kind != "phase_program":
            raise KeyError(f"block {self.name!r}: no phase_program block {self.params['program']!r}")
        met = None
        if self.params["metric"] is not None:
            met = blocks.get(self.params["metric"])
            if met is None or met.kind != "metric_ii":
                raise KeyError(f"block {self.name!r}: no metric_ii block {self.params['metric']!r}")
        unknown = sorted((set(self.positions) | set(self.actuators) | set(self.slips)) - set(prog.ports))
        if unknown:
            raise ValueError(f"block {self.name!r}: ports {unknown} are not programmed by {prog.name!r}")
        return prog, met

    def required_slots(self, blocks: Mapping[str, Any]) -> list[str]:
        """Slots this block reads from the record (the runner adds them to ``measure.record``)."""
        _, met = self._blocks(blocks)
        out = [*self.positions.values(), *self.actuators.values(), *self.slips.values()]
        if met is not None:
            out += [met.slot(p.port) for p in met.pairs] + list(met.setpoint_slots.values())
            out += [met.mask_slot] if met.mask_slot is not None else []
        return list(dict.fromkeys(out))

    def evaluate(self, series: Mapping[str, np.ndarray], t_s: np.ndarray, blocks: Mapping[str, Any],
                 record: Mapping[str, Any] | None = None) -> dict:
        """Per-switch results (module docstring) from the recorded ``series`` on the time axis ``t_s``; ``record`` the
        run's ``measure.record`` (an anti-aliased record is refused)."""
        rec = record or {}
        if int(rec.get("decim", 1) or 1) > 1 and str(rec.get("antialias", "none")) != "none":
            raise ValueError(f"block {self.name!r}: measure.record.antialias {rec.get('antialias')!r} filters the "
                             "wrapped lock coordinates / metric slots across their folds; step_response needs a "
                             "point-sampled record (antialias: none)")
        prog, met = self._blocks(blocks)
        t = np.asarray(t_s, dtype=float)

        def ser(slot: str) -> np.ndarray:
            if slot not in series:
                raise KeyError(f"block {self.name!r}: slot {slot!r} is not recorded")
            return np.asarray(series[slot], dtype=float).ravel()[: t.size]

        times = prog.times[prog.times <= t[-1]] if t.size else prog.times[:0]  # the switches the record reaches
        n = times.size
        ends = np.r_[times[1:], np.inf] if n else times
        col = {k: i for i, k in enumerate(prog.ports)}
        goals = -float(prog.params["sign"]) * prog.table[:n]  # where the lock coordinate goes (P-3)
        prev = np.vstack([np.zeros((1, goals.shape[1])), goals[:-1]]) if n else goals
        tol0, tol_sig = float(self.params["tol_rad"]), float(self.params["tol_sigma"])
        segs = [(t >= times[j]) & (t < ends[j]) for j in range(n)]
        pre = (t >= 0.5 * prog.times[0]) & (t < prog.times[0])  # unprogrammed, locked: noise level and program zero
        if met is not None and met.mask_slot is not None:
            pre &= ser(met.mask_slot) > 0.5

        def table_of(y: np.ndarray, g: np.ndarray, gp: np.ndarray, period: float) -> dict[str, np.ndarray]:
            sig = float(np.std(_wrap(y[pre], period))) if pre.sum() > 1 else 0.0
            tol = max(tol0, tol_sig * sig)
            rows = [segment_response(t[m], g[j] + _wrap(y[m] - g[j], period), times[j], gp[j], g[j], tol,
                                     path=unwrap_path(y[m], gp[j], period), period=period)
                    for j, m in enumerate(segs)]
            keys = ("rise_s", "overshoot", "settle_s", "bias_rad", "sigma_rad", "ambiguous")
            return {k: np.array([r[k] for r in rows]) for k in keys} | {"tol": tol}

        ports = sorted(self.positions)
        sgn = float(self.params["position_sign"])
        ys = {k: sgn * ser(self.positions[k]) for k in ports}
        per = {k: table_of(ys[k], goals[:, col[k]], prev[:, col[k]], 2.0 * math.pi) for k in ports}
        out: dict[str, Any] = {"ports": np.array(ports, dtype=np.int64), "t_switch_s": times.copy(),
                               "step_rad": np.array([[goals[j, col[k]] - prev[j, col[k]] for k in ports]
                                                     for j in range(n)]).reshape(n, len(ports)),
                               "tol_used_rad": np.array([per[k]["tol"] for k in ports])}
        for key in ("rise_s", "overshoot", "settle_s", "bias_rad", "sigma_rad", "ambiguous"):
            out[key] = np.array([per[k][key] for k in ports]).T.reshape(n, len(ports))
        nan = math.nan
        act = np.full((n, len(ports)), nan)
        act_step = np.full((n, len(ports)), nan)
        slips = np.full((n, len(ports)), nan)
        for i, k in enumerate(ports):
            u = ser(self.actuators[k]) if k in self.actuators else None
            c = ser(self.slips[k]) if k in self.slips else None
            for j, m in enumerate(segs):
                idx = np.flatnonzero(m)
                if not idx.size:
                    continue
                if u is not None:
                    act[j, i] = np.max(np.abs(u[idx])) / self.act_range[k]
                    act_step[j, i] = u[idx[-1]] - u[idx[0]]
                if c is not None:
                    slips[j, i] = c[idx[-1]] - (c[idx[0] - 1] if idx[0] > 0 else c[idx[0]])
        out.update(act_peak_frac=act, act_step=act_step, slips=slips)
        # cross-port consistency on the lock coordinates
        cs, cp, dset, drise, pairs, pbias = [], [], [], [], [], []
        for a in range(len(ports)):
            for b in range(a + 1, len(ports)):
                k, q = ports[a], ports[b]
                pairs.append((k, q))
                dy = ys[k] - ys[q]
                pbias.append([np.mean(_wrap(dy[m] - (goals[j, col[k]] - goals[j, col[q]]), 2.0 * math.pi)[m.sum() // 2:])
                              if m.any() else nan for j, m in enumerate(segs)])
                for j in range(n):
                    sk, sl = out["step_rad"][j, a], out["step_rad"][j, b]
                    if abs(sk) > 1e-12 and abs(sk - sl) < 1e-9:
                        cs.append(j)
                        cp.append((k, q))
                        dset.append(out["settle_s"][j, a] - out["settle_s"][j, b])
                        drise.append(out["rise_s"][j, a] - out["rise_s"][j, b])
        out.update(cross_switch=np.array(cs, dtype=np.int64), cross_ports=np.array(cp, dtype=np.int64).reshape(-1, 2),
                   cross_d_settle_s=np.array(dset), cross_d_rise_s=np.array(drise),
                   pair_ports=np.array(pairs, dtype=np.int64).reshape(-1, 2),
                   pair_bias_rad=np.array(pbias).reshape(len(pairs), n).T)
        # metric (ii): Θ_k − Θ_0 against the programmed difference
        if met is not None:
            tports = [p.port for p in met.pairs]
            sps = met.setpoint_slots
            th = {k: ser(met.slot(k)) + (ser(sps[k]) if k in sps else 0.0) - (ser(sps[0]) if 0 in sps else 0.0)
                  for k in tports}
            zero = {k: float(np.angle(np.mean(np.exp(2j * th[k][pre]))) / 2.0) if pre.any() else 0.0 for k in tports}
            th = {k: th[k] - zero[k] for k in tports}
            out["theta_zero_rad"] = np.array([zero[k] for k in tports])

            def g_of(k: int, g: np.ndarray) -> np.ndarray:
                return (g[:, col[k]] if k in col else 0.0) - (g[:, col[0]] if 0 in col else 0.0) + np.zeros(n)

            tper = {k: table_of(th[k], g_of(k, goals), g_of(k, prev), math.pi) for k in tports}
            out["theta_ports"] = np.array(tports, dtype=np.int64)
            out["theta_step_rad"] = np.array([g_of(k, goals) - g_of(k, prev) for k in tports]).T.reshape(n, len(tports))
            out["theta_tol_used_rad"] = np.array([tper[k]["tol"] for k in tports])
            for key in ("rise_s", "overshoot", "settle_s", "bias_rad", "sigma_rad", "ambiguous"):
                out[f"theta_{key}"] = np.array([tper[k][key] for k in tports]).T.reshape(n, len(tports))
        s = out["settle_s"]
        out["n_ambiguous"] = int(np.sum(out["ambiguous"])) + int(np.sum(out.get("theta_ambiguous", 0)))
        out["n_unsettled"] = int(np.sum(np.isnan(s)))
        out["settle_s_max"] = float(np.nanmax(s)) if np.any(np.isfinite(s)) else nan
        out["bias_rms_rad"] = float(np.sqrt(np.nanmean(out["bias_rad"] ** 2))) if n else nan
        return out


@register_block(KIND)
def make_step_response(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``step_response`` (module docstring)."""
    return StepResponse(params, name)
