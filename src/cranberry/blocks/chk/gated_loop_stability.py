"""``gated_loop_stability_check``: stability and residual error of a periodically gated (sampled-and-held) lock loop
(Phase 6 WP6-G, D85-2; notes/lit-madsen2022.md §4.2 N3; D75 rule: every point checked, none pruned; D80-2: an
unstable point emits ``loop_unstable`` so that ``RunResult.valid`` is False). Offline CHK block (G16-j, D34).

Why a separate check: ``loop_stability_check`` is the LTI Nyquist test of the continuous loop. A loop whose error
exists only in a lock window of D_lock·T per period T (``gate_schedule`` + ``gated_hold``) is linear periodically
time-varying; adding the hold as a delay (``lock_extra_delay_s``) is right only for f_u ≪ 1/T. Borealis' T = 100 µs
is in the kHz loop band, so the test here is exact in discrete time.

Model (lifted / Floquet analysis of the sampled loop at the time-loop rate, h = 1/f_sim). The open loop is
``loop_stability_check``'s (same params, :func:`~cranberry.blocks.chk.loop_stability.open_loop` for the ungated
reference):

    L(s) = k0 · C_1(s) · e^{−s·delay_s} · G_1(s) · [1 + Σ_k (Π_{j=2..k} C_j(s)) k_act,k A_k(s)]

discretised as the time loop runs it:

* C_1 (the gated controller) is the ``pid`` law of ``pid_discretisation`` (``tustin`` | ``dsc1`` | ``forward_euler``,
  ``blocks.ctrl.pid.first_order``; K_p, K_i, K_d, no ω_IS / ω_DS / K_i2), state [c[n−1], I, D];
* the gate acts on C_1's input c (``gated_hold`` between the error signal and the controller): in the active part of
  the lock window c = e; outside it c = off_gain·e (``hold_mode: zero``; off_gain = 0 for ``gated_hold`` zero /
  noise_only, √leak for a probe gated by ``gate_schedule`` without a hold block) or c = the last active e
  (``hold_mode: hold``). The controller's gate is the lock window delayed by ``gate_delay_s`` (the probe's flight to its
  detector when ``gated_hold`` reads a delayed ``gate.<name>.lock``) and its first ⌈settle_s·f_sim⌉ samples are blanked
  (``gated_hold.settle_s``);
* delay_s is the total loop delay of the time loop (D55-1, D61-2: tau_s + 1/f_sim + explicit delays; ≥ 1/f_sim), as
  integer samples + the linear fractional interpolation of ``delay_buffer``. No half-sample ZOH term (that belongs to
  the continuous approximation of ``loop_stability_check``);
* G_1, C_j (j ≥ 2: K_p + K_i/s, K_d = 0) and A_k are Tustin-discretised element by element (a mode prewarped at its own
  frequency, as ``pzt_fiber_stretcher``), so the stages run at the base rate (the TEC stage's own ``rate_div`` hold is
  neglected: it is ≪ the gate period's effect in the presets).

With z[n+1] = A_on z[n] + B_on d[n] (controller active) or A_off, B_off (otherwise), e[n] = C z[n] + d[n] (e: the lock
error, d: the free-running phase of the locked combination), the monodromy matrix over one period (N = round(T f_sim)
samples) is M = Π_n A_n. **Stable ⇔ ρ = max |eig M| < 1** (Floquet multipliers); ``growth_per_s`` = ln ρ / T;
``multipliers`` lists the four largest |eig M| (with a slow offload stage the largest is its slow mode, ρ ≈ 1 − 2πf_12 T,
and the fast loop's gated multiplier is the next one). A
window without an active sample gives ρ = 1 (the integrator is never corrected): reported as not stable.

Residual error (``disturbance`` given, stable points): the exact cyclostationary variance of e for a stationary d of
one-sided PSD S(f) [rad²/Hz] (f < f_sim/2): for d = e^{i2πfnh} the periodic steady state is solved from M
(x_0 = (e^{i2πfT} − M)⁻¹·forced response), and Var e[n] = ∫ S(f) |e_f[n]|² df at period position n. Reported:

* ``var_meas_rad2``   mean over the measurement window (what the GBS sees; the metric ``mask_slot`` window),
* ``var_lock_rad2``   mean over the lock window, ``var_period_rad2`` over the whole period,
* ``var_conversion_rad2``  the frequency-conversion (Dick-effect) part: writing e_f[n] = e^{i2πfnh} p_f[n] with p_f
  periodic, Σ_k |c_k|² = ⟨|p_f|²⟩ over the period (Parseval) and c_0 = ⟨p_f⟩ is the time-invariant part; ⟨|p_f|²⟩ −
  |⟨p_f⟩|² is the power the gating moves to f + k/T, k ≠ 0. Integrated over a stationary S this is the aliasing of
  the disturbance by the window (noise near k/T folded to low frequency and written onto the actuator, and the
  converse); it vanishes for an ungated loop,
* ``var_hold_drift_rad2``  the error accumulated during the hold, predicted from the free-running structure function
  D(τ) = ∫ S·4 sin²(πfτ) df averaged over the measurement window, τ from the last active sample (open loop during
  the hold; it neglects the residual at the window end, which ``var_meas_rad2`` contains),
* ``var_ungated_rad2`` = ∫ S |1/(1 + L)|² df of the continuous LTI loop (reference: the D_lock → 1 limit).

Backscatter round trip (optional ``round_trip``, notes/lit-nardelli2026.md §3.1-3): with a counter-propagating probe
the backscatter keeps reaching the GBS side for 2 n_g L/c after the probe is switched off; a detector window centred in
the probe-off time needs **f_gate ≤ (1 − D_pr − D_det) / [2 (2 n_g L/c + t_edge)]**; above it a CHK ``gate_round_trip``
warning (not invalidating: it is a stray-light, not a stability, problem; the size is ``window_stray_light``'s).

Params (all required, D34): ``loop_stability_check``'s ``k0``, ``kp``, ``ki_per_s``, ``kd_s``, ``delay_s``, ``poles_hz``,
``modes`` ({f_hz, q}), ``f_grid_hz`` (ungated Nyquist reference), optional ``stages``; ``pid_discretisation``,
``f_sim_hz``, ``period_s``, ``lock_fraction``, ``gate_delay_s``, ``settle_s``, ``hold_mode`` (``zero`` | ``hold``), ``off_gain``,
``meas_delay_s`` (measurement window start after the lock window end); toggles {enabled}. Optional (absent or null:
YAML overlays cannot carry a null, it deletes the key): ``meas_duty`` (fraction of the period; absent = to the end of
the period), ``disturbance`` ({f_hz: [...], psd_rad2_per_hz: [...]}; absent = stability only), ``round_trip`` ({n_g,
L_m, t_edge_s, duty_det}; absent = no check). Result in :attr:`result`.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.blocks.act.pzt import discretise, modes_ss, poles_ss
from cranberry.blocks.chk.loop_stability import (
    STAGE_KEYS,
    n_integrators,
    nyquist_stable,
    open_loop,
)
from cranberry.blocks.ctrl.pid import DISCRETISATIONS, first_order
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block

__all__ = ["GatedLoopStabilityCheck", "make_gated_loop_stability_check", "gated_loop_model", "floquet"]

KIND = "gated_loop_stability_check"
HOLD_MODES = ("zero", "hold")
C_LIGHT = 299792458.0
#: the residual (disturbance) needs two passes over one period per frequency: bounded work per check
MAX_RESPONSE_SAMPLES = 200_000

# --------------------------------------------------------------------------------------- discrete state-space algebra


def _ss(A=None, B=None, C=None, D=1.0):
    n = 0 if A is None else np.asarray(A).shape[0]
    return (np.zeros((0, 0)) if A is None else np.atleast_2d(np.asarray(A, float)),
            np.zeros((n, 1)) if B is None else np.asarray(B, float).reshape(n, 1),
            np.zeros((1, n)) if C is None else np.asarray(C, float).reshape(1, n),
            np.array([[float(D)]]))


def _series(s1, s2):
    """s1 then s2 (SISO)."""
    a1, b1, c1, d1 = s1
    a2, b2, c2, d2 = s2
    n1, n2 = a1.shape[0], a2.shape[0]
    a = np.block([[a1, np.zeros((n1, n2))], [b2 @ c1, a2]])
    return a, np.vstack([b1, b2 @ d1]), np.hstack([d2 @ c1, c2]), d2 @ d1


def _sum(s1, s2):
    """Same input, outputs added."""
    a1, b1, c1, d1 = s1
    a2, b2, c2, d2 = s2
    n1, n2 = a1.shape[0], a2.shape[0]
    a = np.block([[a1, np.zeros((n1, n2))], [np.zeros((n2, n1)), a2]])
    return a, np.vstack([b1, b2]), np.hstack([c1, c2]), d1 + d2


def _tustin(sys_c, h, prewarp_hz=None):
    ad, bd, cd, dd = discretise(sys_c, h, "tustin", prewarp_hz)
    if ad.size == 0:
        return _ss(D=dd)
    return _ss(ad, bd, cd, dd)


def _dynamics_d(q: Mapping[str, Any], h: float):
    """Discrete Π poles · Π modes (each element Tustin, a mode prewarped at its frequency)."""
    out = _ss()
    poles = [float(p) for p in q["poles_hz"] if math.isfinite(float(p))]
    if poles:
        out = _series(out, _tustin(poles_ss(poles), h))
    for m in q["modes"]:
        out = _series(out, _tustin(modes_ss([{"f_hz": m["f_hz"], "q": m["q"], "r": 1.0}]), h, float(m["f_hz"])))
    return out


def _stage_ctrl_d(st: Mapping[str, Any], h: float):
    """Discrete K_p + K_i/s of a stage controller (Tustin; dsc1 adds K_i·h_s/2 to K_p, as ``loop_stability``)."""
    if float(st["kd_s"]) != 0.0:
        raise ValueError("gated_loop_stability_check: stage controllers with kd_s != 0 are not modelled")
    kp, ki = float(st["kp"]), float(st["ki_per_s"])
    if st.get("discretisation", "continuous") == "dsc1":
        kp += ki * float(st["h_s"]) / 2.0
    if ki == 0.0:
        return _ss(D=kp)
    return _tustin((np.zeros((1, 1)), np.ones((1, 1)), np.array([[ki]]), np.array([[kp]])), h)


def plant_d(p: Mapping[str, Any], h: float):
    """Discrete G_1 · [1 + Σ_k chain_k · k_act,k · A_k] (module docstring; k0, C_1 and the delay excluded)."""
    inner, chain = _ss(), _ss()
    for st in p.get("stages") or ():
        chain = _series(chain, _stage_ctrl_d(st, h))
        path = _series(_series(chain, _ss(D=float(st["k_act"]))), _dynamics_d(st, h))
        inner = _sum(inner, path)
    return _series(inner, _dynamics_d(p, h))


# --------------------------------------------------------------------------------------- the gated loop


def gated_loop_model(p: Mapping[str, Any]) -> dict[str, Any]:
    """Discrete loop matrices {A_on, B_on, A_off, B_off, C, D, h} (module docstring)."""
    h = 1.0 / float(p["f_sim_hz"])
    ap, bp, cp, dp = plant_d(p, h)
    npl = ap.shape[0]
    dsm = float(p["delay_s"]) / h
    if dsm < 1.0 - 1e-9:
        raise ValueError(f"gated_loop_stability_check: delay_s must be >= 1/f_sim (the loop closure), got {p['delay_s']!r}")
    k = int(math.floor(dsm + 1e-9))
    mu = dsm - k
    if mu < 1e-9:
        mu = 0.0
    nh = k + (1 if mu > 0.0 else 0)
    w = np.zeros(nh)
    w[k - 1] = 1.0 - mu
    if mu > 0.0:
        w[k] = mu
    meth = DISCRETISATIONS[p["pid_discretisation"]]
    kp, ki, kd = float(p["kp"]), float(p["ki_per_s"]), float(p["kd_s"])
    b0, b1, a1 = first_order(meth, 0.0, ki, 0.0, h)
    k0 = float(p["k0"])
    hold = p["hold_mode"] == "hold"
    g_off = float(p["off_gain"])
    i_h, i_c = npl, npl + nh  # u history, controller [c1, I, D, eh]
    m = npl + nh + 4

    def step(z, d, on):
        xp, uh = z[:npl], z[i_h:i_c]
        c1, ii, dd, eh = z[i_c:i_c + 4]
        ud = w @ uh
        y = k0 * ((cp @ xp)[0] + dp[0, 0] * ud)
        e = d - y
        if hold:
            c = e if on else eh
            eh_new = c
        else:
            c = e if on else g_off * e
            eh_new = 0.0
        i_new = -a1 * ii + b0 * c + b1 * c1
        if kd == 0.0:
            d_new = 0.0
        elif meth == 1:
            d_new = -dd + kd * (2.0 / h) * (c - c1)
        else:
            d_new = kd * (c - c1) / h
        u = kp * c + i_new + d_new
        zn = np.empty(m)
        zn[:npl] = ap @ xp + bp[:, 0] * ud
        zn[i_h] = u
        zn[i_h + 1:i_c] = uh[:-1]
        zn[i_c:i_c + 4] = (c, i_new, d_new, eh_new)
        return zn, e

    eye = np.eye(m)
    out: dict[str, Any] = {"h": h, "m": m}
    for tag, on in (("on", True), ("off", False)):
        cols = [step(eye[:, j], 0.0, on) for j in range(m)]
        out[f"A_{tag}"] = np.column_stack([c[0] for c in cols])
        out[f"B_{tag}"] = step(np.zeros(m), 1.0, on)[0]
        out["C"] = np.array([c[1] for c in cols])
    out["D"] = 1.0
    return out


def schedule(p: Mapping[str, Any]) -> dict[str, Any]:
    """Sample schedule of one period: N, active[n] (controller fed), lock[n], meas[n] (module docstring)."""
    f = float(p["f_sim_hz"])
    T, dl = float(p["period_s"]), float(p["lock_fraction"])
    n = int(round(T * f))
    nl = min(n, int(math.ceil(dl * T * f - 1e-9)))
    ns = int(math.ceil(float(p["settle_s"]) * f - 1e-9))
    nd = int(math.ceil(float(p["gate_delay_s"]) * f - 1e-9))
    idx = np.arange(n)
    lock = idx < nl
    rel = (idx - nd) % n
    active = (rel < nl) & (rel >= ns)
    m0 = int(math.ceil((dl * T + float(p["meas_delay_s"])) * f - 1e-9))
    m1 = n if p.get("meas_duty") is None else int(math.ceil((dl * T + float(p["meas_delay_s"]) + float(p["meas_duty"]) * T)
                                                         * f - 1e-9))
    meas = (idx >= m0) & (idx < min(m1, n))
    return {"N": n, "N_lock": nl, "active": active, "lock": lock, "meas": meas}


def _runs(active: np.ndarray) -> list[tuple[bool, int]]:
    """Run-length encoding [(value, count), ...] of a boolean array."""
    a = np.asarray(active, bool)
    edges = np.flatnonzero(np.diff(a)) + 1
    starts = np.r_[0, edges]
    counts = np.diff(np.r_[starts, a.size])
    return [(bool(a[i]), int(c)) for i, c in zip(starts, counts)]


def floquet(model: Mapping[str, Any], active: np.ndarray) -> np.ndarray:
    """Monodromy matrix M = Π_n A_n over one period (sample 0 acts first)."""
    M = np.eye(model["m"])
    for on, cnt in _runs(active):
        M = np.linalg.matrix_power(model["A_on" if on else "A_off"], cnt) @ M
    return M


def periodic_response(model: Mapping[str, Any], active: np.ndarray, M: np.ndarray, f: np.ndarray) -> np.ndarray:
    """e_f[n] (shape (N, F)): periodic steady-state error for d[n] = e^{i2πf n h} (module docstring)."""
    h, n = model["h"], active.size
    w = 2.0 * np.pi * f * h
    A = {True: model["A_on"], False: model["A_off"]}
    B = {True: model["B_on"][:, None], False: model["B_off"][:, None]}
    x = np.zeros((model["m"], f.size), complex)
    for k in range(n):
        x = A[bool(active[k])] @ x + B[bool(active[k])] * np.exp(1j * w * k)[None, :]
    lhs = np.exp(1j * w * n)[:, None, None] * np.eye(model["m"])[None] - M[None]
    x = np.linalg.solve(lhs, x.T[:, :, None])[:, :, 0].T
    e = np.empty((n, f.size), complex)
    C = model["C"]
    for k in range(n):
        ph = np.exp(1j * w * k)
        e[k] = C @ x + model["D"] * ph
        x = A[bool(active[k])] @ x + B[bool(active[k])] * ph[None, :]
    return e


class GatedLoopStabilityCheck(BaseBlock):
    """The ``gated_loop_stability_check`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "k0", "kp", "ki_per_s", "kd_s", "delay_s", "poles_hz", "modes", "f_grid_hz",
                "pid_discretisation", "f_sim_hz", "period_s", "lock_fraction", "gate_delay_s", "settle_s", "hold_mode", "off_gain",
                "meas_delay_s", block=name)
        super().__init__(KIND, name, params)
        if params["pid_discretisation"] not in DISCRETISATIONS:
            raise ValueError(f"block {name!r}: pid_discretisation must be one of {sorted(DISCRETISATIONS)}")
        if params["hold_mode"] not in HOLD_MODES:
            raise ValueError(f"block {name!r}: hold_mode must be one of {HOLD_MODES}")
        if not 0.0 < float(params["lock_fraction"]) <= 1.0 or not float(params["period_s"]) > 0.0:
            raise ValueError(f"block {name!r}: need period_s > 0 and 0 < lock_fraction <= 1")
        for i, st in enumerate(params.get("stages") or ()):
            require(st, *STAGE_KEYS, block=f"{name}.stages[{i}]")
        dist = params.get("disturbance")
        if dist is not None:
            require(dist, "f_hz", "psd_rad2_per_hz", block=f"{name}.disturbance")
        if params.get("round_trip") is not None:
            require(params["round_trip"], "n_g", "L_m", "t_edge_s", "duty_det", block=f"{name}.round_trip")
        self.chk_log = ChkLog()
        self.result: dict[str, Any] | None = None
        if self._toggles["enabled"]:
            self.result = self._check()

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _check(self) -> dict[str, Any]:
        p = self.params
        model = gated_loop_model(p)
        sch = schedule(p)
        M = floquet(model, sch["active"])
        mags = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
        rho = float(mags[0])
        T = float(p["period_s"])
        lo, hi, n = p["f_grid_hz"]
        fg = np.geomspace(float(lo), float(hi), int(n))
        lg = open_loop(fg, p)
        u_stable, u_mx, u_f180 = nyquist_stable(fg, lg, n_integrators(p))
        res: dict[str, Any] = {
            "stable": bool(rho < 1.0), "rho": rho, "growth_per_s": math.log(rho) / T if rho > 0 else -math.inf,
            "multipliers": [float(x) for x in mags[:4]],
            "n_period": sch["N"], "n_lock": sch["N_lock"], "n_active": int(sch["active"].sum()),
            "n_meas": int(sch["meas"].sum()), "rho_ungated_per_sample": float(np.max(np.abs(np.linalg.eigvals(model["A_on"])))),
            "ungated_stable": bool(u_stable), "ungated_max_l_at_180": u_mx, "ungated_f_180_hz": u_f180,
        }
        if p.get("disturbance") is not None and res["stable"]:
            res.update(self._variances(model, sch, M))
        if not res["stable"]:
            self.chk_log.warn("loop_unstable", self.name,
                              f"gated loop unstable: Floquet multiplier ρ = {rho:.4g} per period (T = {T:.4g} s, "
                              f"lock fraction {float(p['lock_fraction']):.3g}, hold_mode {p['hold_mode']}); the time "
                              "loop will oscillate / run away from window to window", **res)
        rt = p.get("round_trip")
        if rt is not None:
            t_rt = 2.0 * float(rt["n_g"]) * float(rt["L_m"]) / C_LIGHT
            f_max = (1.0 - float(p["lock_fraction"]) - float(rt["duty_det"])) / (2.0 * (t_rt + float(rt["t_edge_s"])))
            res["round_trip_f_max_hz"] = f_max
            if 1.0 / T > f_max:
                self.chk_log.warn("gate_round_trip", self.name,
                                  f"gate rate {1.0 / T:.4g} Hz > (1 − D_pr − D_det)/[2(2n_gL/c + t_edge)] = {f_max:.4g} Hz: "
                                  "the probe's backscatter tail (2n_gL/c) reaches the detector window",
                                  gate_rate_hz=1.0 / T, f_max_hz=f_max, t_round_trip_s=t_rt)
        return res

    def _variances(self, model, sch, M) -> dict[str, float]:
        p = self.params
        d = p["disturbance"]
        f = np.asarray(d["f_hz"], float)
        s = np.asarray(d["psd_rad2_per_hz"], float)
        if f.ndim != 1 or f.shape != s.shape or f.size < 2 or np.any(np.diff(f) <= 0):
            raise ValueError(f"block {self.name!r}: disturbance needs increasing f_hz and psd_rad2_per_hz of equal length")
        if sch["N"] > MAX_RESPONSE_SAMPLES:
            raise ValueError(f"block {self.name!r}: disturbance response over {sch['N']} samples per period "
                             f"(> {MAX_RESPONSE_SAMPLES}); use it for periods ≤ {MAX_RESPONSE_SAMPLES} / f_sim")
        keep = f < 0.5 * float(p["f_sim_hz"])
        f, s = f[keep], s[keep]
        e = periodic_response(model, sch["active"], M, f)
        pw = np.abs(e) ** 2
        n = sch["N"]
        ph = np.exp(-2j * np.pi * np.outer(np.arange(n), f) * model["h"])
        c0 = np.mean(e * ph, axis=0)
        w_period = pw.mean(axis=0)
        act = np.flatnonzero(sch["active"])
        idx = np.arange(n)[sch["meas"]]
        # time since the last active sample before each measurement sample (cyclic)
        if act.size:
            pos = np.searchsorted(act, idx) - 1
            last = np.where(pos >= 0, act[np.maximum(pos, 0)], act[-1] - n)
        else:
            last = idx
        tau = (idx - last) * model["h"]
        w_drift = np.mean(4.0 * np.sin(np.pi * np.outer(tau, f)) ** 2, axis=0) if tau.size else np.zeros_like(f)
        lg = open_loop(f, p)

        def integ(wt):
            return float(np.trapezoid(s * wt, f))

        return {
            "var_meas_rad2": integ(pw[sch["meas"]].mean(axis=0)) if sch["meas"].any() else math.nan,
            "var_lock_rad2": integ(pw[sch["lock"]].mean(axis=0)) if sch["lock"].any() else math.nan,
            "var_period_rad2": integ(w_period),
            "var_conversion_rad2": integ(w_period - np.abs(c0) ** 2),
            "var_hold_drift_rad2": integ(w_drift),
            "var_ungated_rad2": integ(np.abs(1.0 / (1.0 + lg)) ** 2),
        }


@register_block(KIND)
def make_gated_loop_stability_check(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``gated_loop_stability_check`` block (module docstring)."""
    return GatedLoopStabilityCheck(params, name)
