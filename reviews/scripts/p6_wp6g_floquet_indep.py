"""Phase 6 review (WP6-G, D90): independent re-derivation of `gated_loop_stability_check`.

Independent of the block's matrix assembly (`gated_loop_model`, `plant_d`, `floquet`, `periodic_response`):

* the loop is re-written sample by sample from scipy.signal.bilinear sections (the PZT mode prewarped at f_res via
  fs_eff = ω0/(2 tan(ω0 h/2)), the TEC lag and the stage PI plain Tustin), a u-history delay line with linear
  fractional interpolation, the I-only Tustin / dsc1 PID and the gate (zero / hold, gate delay, settle blanking);
* the TEC path is modelled as the time loop wires it (from the undelayed fine output, not through the PZT mode or the
  loop delay) and optionally with its own rate_div hold (pid1s / tec1 at f_sim / tec_rate_div), in which case the
  monodromy is taken over lcm(T, T_tec) and compared as ρ^(lcm/T);
* the cyclostationary variance is computed in the time domain from the periodic Lyapunov recursion with an AR(1)
  disturbance state (exact, no frequency grid), against the block's frequency-domain var_meas / var_lock / var_period;
* the Dick (conversion) term is checked with a pure tone: a long time-domain run, the power left at the tone bin
  versus the total (the block's c_0 / ⟨|p_f|²⟩ at that frequency).

Run: uv run python reviews/scripts/p6_wp6g_floquet_indep.py
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import scipy.linalg
import scipy.signal

from cranberry.blocks.chk.gated_loop_stability import make_gated_loop_stability_check

FS = 1.0e6
H = 1.0 / FS

# F_liu20_nominal + gate/borealis_window + gate/borealis_si (printed by the preset, f_u = 1 kHz, 20 kHz Q 10)
BASE = {"k0": 0.4110960958218894, "kp": 0.0, "ki_per_s": 15283.981947378614, "kd_s": 0.0, "delay_s": 1e-06,
        "poles_hz": [], "modes": [{"f_hz": 20000.0, "q": 10.0}],
        "stages": [{"kp": 0.0, "ki_per_s": 5.489464870760138e-05, "kd_s": 0.0, "k_act": 3433.769295426536,
                    "poles_hz": [0.015915494309189534], "modes": []}],
        "f_grid_hz": [1e-6, 1e8, 40001], "pid_discretisation": "tustin", "f_sim_hz": FS, "period_s": 1e-4,
        "lock_fraction": 0.65, "gate_delay_s": 6.245738809079714e-06, "settle_s": 1e-06, "hold_mode": "zero",
        "off_gain": 0.0, "meas_delay_s": 1e-07, "meas_duty": None, "disturbance": None, "round_trip": None,
        "toggles": {"enabled": True}}


def block(p):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_gated_loop_stability_check(p, name="chk").result


class Sec:
    """One SISO DF-II-transposed section (b, a in z^-1, a[0] = 1)."""

    def __init__(self, b, a):
        b, a = np.atleast_1d(np.asarray(b, float)), np.atleast_1d(np.asarray(a, float))
        n = max(len(a), len(b))
        self.b = np.r_[b, np.zeros(n - len(b))] / a[0]
        self.a = np.r_[a, np.zeros(n - len(a))] / a[0]
        self.n = n - 1

    def step(self, z, x):
        y = self.b[0] * x + (z[0] if self.n else 0.0)
        zn = np.empty(self.n)
        for i in range(self.n):
            zn[i] = self.b[i + 1] * x - self.a[i + 1] * y + (z[i + 1] if i + 1 < self.n else 0.0)
        return y, zn


def mode_sec(f0, q, h):
    w0 = 2 * math.pi * f0
    fs_eff = w0 / (2 * math.tan(w0 * h / 2))  # Tustin prewarped at f0
    bz, az = scipy.signal.bilinear([w0 * w0], [1.0, w0 / q, w0 * w0], fs=fs_eff)
    return Sec(bz, az)


class Loop:
    """Sample-by-sample gated loop; state = flat vector. TEC stage optionally held at rate_div."""

    def __init__(self, p, tec_div=1):
        self.p, self.h = p, 1.0 / float(p["f_sim_hz"])
        h = self.h
        self.mode = [mode_sec(m["f_hz"], m["q"], h) for m in p["modes"]]
        st = (p.get("stages") or [None])[0]
        self.tec = st is not None
        self.tec_div = tec_div
        if self.tec:
            ht = h * tec_div
            bz, az = scipy.signal.bilinear([st["kp"], st["ki_per_s"]], [1.0, 0.0], fs=1 / ht)
            self.pi2 = Sec(bz, az)
            wp = 2 * math.pi * st["poles_hz"][0]
            bz, az = scipy.signal.bilinear([wp], [1.0, wp], fs=1 / h)
            self.lag = Sec(bz, az)
            self.kact = st["k_act"]
        dsm = float(p["delay_s"]) / h
        self.k = int(math.floor(dsm + 1e-9))
        self.mu = dsm - self.k if dsm - self.k > 1e-9 else 0.0
        self.nh = self.k + 1
        # layout
        self.sl = {}
        o = 0
        for name, n in [("uh", self.nh), ("mode", sum(s.n for s in self.mode)), ("I", 1), ("c1", 1), ("eh", 1),
                        ("pi2", self.pi2.n if self.tec else 0), ("vt", 1 if self.tec else 0),
                        ("lag", self.lag.n if self.tec else 0), ("cnt", 0)]:
            self.sl[name] = slice(o, o + n)
            o += n
        self.m = o

    def step(self, z, d, on, tec_sample):
        p, s = self.p, self.sl
        uh = z[s["uh"]]
        ud = (1 - self.mu) * uh[self.k - 1] + self.mu * uh[self.k] if self.mu else uh[self.k - 1]
        x, zm, o = ud, [], s["mode"].start
        par = bool(self.p["modes"]) and all("r" in m for m in self.p["modes"])
        acc = 0.0
        for sec, md in zip(self.mode, self.p["modes"]):
            y, zz = sec.step(z[o:o + sec.n], ud if par else x)
            if par:
                acc += float(md["r"]) * y  # parallel sum (pzt_fiber_stretcher modes_ss, D85-1)
            else:
                x = y  # series product
            zm.append(zz)
            o += sec.n
        psi = acc if par else x
        if self.tec:
            vt = z[s["vt"]][0]
            psi_t, zlag = self.lag.step(z[s["lag"]], vt)
            psi = psi + self.kact * psi_t
        e = d - float(p["k0"]) * psi
        if p["hold_mode"] == "hold":
            c = e if on else z[s["eh"]][0]
            eh = c
        else:
            c = e if on else float(p["off_gain"]) * e
            eh = 0.0
        ki, kp = float(p["ki_per_s"]), float(p["kp"])
        if p["pid_discretisation"] == "tustin":
            i_new = z[s["I"]][0] + ki * self.h / 2 * (c + z[s["c1"]][0])
        elif p["pid_discretisation"] == "dsc1":
            i_new = z[s["I"]][0] + ki * self.h * c
        else:
            i_new = z[s["I"]][0] + ki * self.h * z[s["c1"]][0]
        u = kp * c + i_new
        zn = np.zeros(self.m)
        zn[s["uh"]] = np.r_[u, uh[:-1]]
        zn[s["mode"]] = np.concatenate(zm) if zm else []
        zn[s["I"]], zn[s["c1"]], zn[s["eh"]] = i_new, c, eh
        if self.tec:
            zn[s["lag"]] = zlag
            if tec_sample:  # TEC controller reads the undelayed fine output u, holds its output
                v, zpi = self.pi2.step(z[s["pi2"]], u)
                zn[s["pi2"]], zn[s["vt"]] = zpi, v
            else:
                zn[s["pi2"]], zn[s["vt"]] = z[s["pi2"]], z[s["vt"]]
        return zn, e

    def matrices(self, on, tec_sample):
        eye = np.eye(self.m)
        A = np.column_stack([self.step(eye[:, j], 0.0, on, tec_sample)[0] for j in range(self.m)])
        C = np.array([self.step(eye[:, j], 0.0, on, tec_sample)[1] for j in range(self.m)])
        B, D = self.step(np.zeros(self.m), 1.0, on, tec_sample)
        return A, B, C, D


def schedule(p):
    n = int(round(p["period_s"] * FS))
    nl = int(math.ceil(p["lock_fraction"] * p["period_s"] * FS - 1e-9))
    nd = int(math.ceil(p["gate_delay_s"] * FS - 1e-9))
    ns = int(math.ceil(p["settle_s"] * FS - 1e-9))
    idx = np.arange(n)
    rel = (idx - nd) % n
    active = (rel < nl) & (rel >= ns)
    m0 = int(math.ceil((p["lock_fraction"] * p["period_s"] + p["meas_delay_s"]) * FS - 1e-9))
    return n, active, idx < nl, idx >= m0


def monodromy(lp, p, n_periods=1):
    n, active, _, _ = schedule(p)
    mats = {(on, ts): lp.matrices(on, ts)[0] for on in (True, False) for ts in (True, False)}
    M = np.eye(lp.m)
    for j in range(n * n_periods):
        ts = (j % lp.tec_div) == 0
        M = mats[(bool(active[j % n]), ts)] @ M
    return M


def rho_indep(p, tec_div=1):
    lp = Loop(p, tec_div)
    nper = 1 if tec_div == 1 else int(np.lcm(int(round(p["period_s"] * FS)), tec_div) // int(round(p["period_s"] * FS)))
    M = monodromy(lp, p, nper)
    ev = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    return ev[:4] ** (1.0 / nper)


def lyapunov_variances(p, a_ar):
    """Exact cyclostationary Var e[n] with d = AR(1) (d[n] = a d[n−1] + w, Var w = 1), periodic Lyapunov."""
    lp = Loop(p)
    n, active, lock, meas = schedule(p)
    m = lp.m + 1
    mats = {}
    for on in (True, False):
        A, B, C, D = lp.matrices(on, True)
        # augmented state [z, d_prev]; d[n] = a d_prev + w[n]
        Aa = np.zeros((m, m))
        Aa[:lp.m, :lp.m] = A
        Aa[:lp.m, lp.m] = B * a_ar
        Aa[lp.m, lp.m] = a_ar
        Ba = np.r_[B, 1.0]
        Ca = np.r_[C, D * a_ar]
        mats[on] = (Aa, Ba, Ca, D)
    # one-period accumulated noise and monodromy
    M, W = np.eye(m), np.zeros((m, m))
    for j in range(n):
        Aa, Ba, _, _ = mats[bool(active[j])]
        M = Aa @ M
        W = Aa @ W @ Aa.T + np.outer(Ba, Ba)
    P = scipy.linalg.solve_discrete_lyapunov(M, W)
    var = np.empty(n)
    for j in range(n):
        Aa, Ba, Ca, D = mats[bool(active[j])]
        var[j] = Ca @ P @ Ca + D * D  # e = C z + D (a d_prev + w)
        P = Aa @ P @ Aa.T + np.outer(Ba, Ba)
    return {"var_meas_rad2": var[meas].mean(), "var_lock_rad2": var[lock].mean(), "var_period_rad2": var.mean()}


def ar1_psd(a, f):
    return 2.0 / FS / np.abs(1 - a * np.exp(-2j * np.pi * f / FS)) ** 2


def tone_conversion(p, f_tone, n_per=400):
    """Fraction of the error power left at the tone bin (time domain) versus the block's |c0|²/⟨|p|²⟩."""
    lp = Loop(p)
    n, active, _, _ = schedule(p)
    z = np.zeros(lp.m)
    N = n * n_per
    e = np.empty(N)
    for j in range(N):
        z, e[j] = lp.step(z, math.cos(2 * math.pi * f_tone * j * H), bool(active[j % n]), True)
    e = e[N // 2:]  # steady state (whole periods)
    E = np.fft.rfft(e) / e.size
    k = int(round(f_tone * e.size * H))
    p_tone = 2 * abs(E[k]) ** 2
    p_tot = np.mean(e * e)
    return p_tone / p_tot


def main():
    print("== (a) Floquet multipliers: block vs independent re-implementation ==")
    cases = [
        ("F_liu20 borealis_si (zero, f_u 1k)", {}),
        ("f_u 5.5k, f_res 5k Q 100", {"ki_per_s": BASE["ki_per_s"] * 5.477, "modes": [{"f_hz": 5e3, "q": 100.0}]}),
        ("f_u 30k, f_res 60k Q 3", {"ki_per_s": BASE["ki_per_s"] * 30, "modes": [{"f_hz": 6e4, "q": 3.0}]}),
        ("T 10 µs, D 0.9, τ_pr 6.2 µs", {"period_s": 1e-5, "lock_fraction": 0.9}),
        ("T 1 ms, D 0.1", {"period_s": 1e-3, "lock_fraction": 0.1}),
        ("hold mode T 2 ms D 0.02 f_u 3k (unstable)", {"hold_mode": "hold", "period_s": 2e-3, "lock_fraction": 0.02,
                                                       "ki_per_s": BASE["ki_per_s"] * 3}),
        ("hold mode T 300 µs D 0.3 f_u 1k", {"hold_mode": "hold", "period_s": 3e-4, "lock_fraction": 0.3}),
        ("dsc1 law, delay 5.7 µs (frac), f_u 5.5k", {"pid_discretisation": "dsc1", "delay_s": 5.7e-6 + 0.3e-6,
                                                     "ki_per_s": BASE["ki_per_s"] * 5.477}),
        ("leak gate off_gain √1e-5, frac delay 2.4 µs", {"off_gain": math.sqrt(1e-5), "delay_s": 2.4e-6}),
    ]
    for name, upd in cases:
        p = {**BASE, **upd}
        r = block(p)
        mine = rho_indep(p)
        print(f"  {name:45s} block ρ {r['rho']:.9f} mult {np.round(r['multipliers'][:3], 7)} | indep "
              f"{np.round(mine[:3], 7)} | stable {r['stable']}")

    print("== (a') TEC stage rate_div 10000 hold (ignored by the check): ρ over lcm(T, 10 ms) ==")
    for name, upd in [("F_liu20 borealis_si", {}),
                      ("worst TEC: f_12 0.3 Hz, τ_th 1 s, T 1e-3",
                       {"period_s": 1e-3, "stages": [{**BASE["stages"][0], "ki_per_s": BASE["stages"][0]["ki_per_s"] * 10,
                                                      "poles_hz": [1 / (2 * math.pi * 1.0)]}]})]:
        p = {**BASE, **upd}
        r = block(p)
        a = rho_indep(p, 1)
        b = rho_indep(p, 10000)
        print(f"  {name:45s} block {r['multipliers'][:2]} | indep base-rate TEC {a[:2]} | TEC held 100 Hz {b[:2]}")

    print("== (b) cyclostationary variance: block (frequency domain) vs periodic Lyapunov (time domain) ==")
    f = np.r_[np.geomspace(0.01, 0.49999 * FS, 6000)]
    for name, upd, a in [("borealis_si, AR(1) 80 Hz", {}, 0.9995),
                         ("borealis_si, AR(1) 8 kHz", {}, 0.95),
                         ("T 1 ms D 0.1, AR(1) 800 Hz", {"period_s": 1e-3, "lock_fraction": 0.1}, 0.995),
                         ("hold T 300 µs D 0.3, AR(1) 80 Hz", {"hold_mode": "hold", "period_s": 3e-4, "lock_fraction": 0.3},
                          0.9995)]:
        p = {**BASE, **upd, "disturbance": {"f_hz": f.tolist(), "psd_rad2_per_hz": ar1_psd(a, f).tolist()}}
        r = block(p)
        ly = lyapunov_variances({**BASE, **upd}, a)
        print(f"  {name:36s} " + "  ".join(f"{k[4:-5]} {r[k]:.4e}/{ly[k]:.4e} ({r[k] / ly[k] - 1:+.2%})" for k in ly)
              + f"  conv {r['var_conversion_rad2']:.3e} drift {r['var_hold_drift_rad2']:.3e} ungated {r['var_ungated_rad2']:.3e}")

    print("== (c) Dick / conversion term: tone test (power left at the tone bin / total) ==")
    from cranberry.blocks.chk.gated_loop_stability import floquet, gated_loop_model, periodic_response
    from cranberry.blocks.chk.gated_loop_stability import schedule as bsched
    for ft in (2.0e3, 7.0e3, 1.0e4 + 37.0):
        p = {**BASE}
        n = int(round(p["period_s"] * FS))
        ft = round(ft * n * 200 * H) / (n * 200 * H)  # on a bin of the steady-state record
        model = gated_loop_model(p)
        sch = bsched(p)
        M = floquet(model, sch["active"])
        e = periodic_response(model, sch["active"], M, np.array([ft]))[:, 0]
        ph = np.exp(-2j * np.pi * np.arange(n) * ft * H)
        frac_block = abs(np.mean(e * ph)) ** 2 / np.mean(abs(e) ** 2)
        frac_time = tone_conversion(p, ft)
        print(f"  f = {ft:9.2f} Hz: block |c0|²/<|p|²> {frac_block:.5f}   time-domain tone fraction {frac_time:.5f}")


def multimode_with_gate():
    """lab/multimode915b_future + gate/borealis_window: the gated check keeps F_common's single mode {pzt_f_res_hz,
    pzt_q} (r = 1); the time-loop stretcher has the 4-mode parallel sum. ρ of the check as composed vs the loop."""
    from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop
    fit = [{"f_hz": 25200.0, "q": 25.0, "r": 0.67}, {"f_hz": 5400.0, "q": 15.0, "r": 0.11},
           {"f_hz": 10000.0, "q": 30.0, "r": 0.07}, {"f_hz": 40300.0, "q": 100.0, "r": 0.15}]
    print("== (d) lab/multimode915b_future + gate/borealis_window (fit point of D89-4) ==")
    ff = np.geomspace(1e-6, 1e8, 400001)
    for fu in (1e3, 2e3, 3e3, 5.48e3):
        p = {**BASE, "ki_per_s": BASE["ki_per_s"] * fu / 1e3}
        comp = block({**p, "modes": [{"f_hz": 25200.0, "q": 25.0}]})
        true = rho_indep({**p, "modes": fit})
        q = {**p, "modes": fit, "delay_s": p["delay_s"] + 0.5e-6}  # chk.lock_stability form (ungated, + ZOH)
        st, mx, f180 = nyquist_stable(ff, open_loop(ff, q), n_integrators(q))
        print(f"  f_u {fu:6.0f} Hz: chk.gated_stability as composed ρ {comp['rho']:.6f} mult "
              f"{np.round(comp['multipliers'][:3], 4)} stable {comp['stable']} | 4-mode loop mult {np.round(true[:3], 4)} | "
              f"chk.lock_stability (4 modes) stable {st}, max|L|@-180° {mx:.3f} @ {f180:.0f} Hz")
    # scan of the overlay axes: verdict of the check as composed vs the 4-mode loop
    opt = pes = agree = 0
    for fu in (1e3, 3e3, 5.48e3):
        for qm in (10.0, 25.0, 300.0):
            for r2 in (-0.15, -0.05, 0.05, 0.15):
                for r4 in (-0.3, 0.15, 0.3):
                    modes = [{"f_hz": 25200.0, "q": qm, "r": 1 - r2 - 0.07 - r4}, {"f_hz": 5400.0, "q": 15.0, "r": r2},
                             {"f_hz": 10000.0, "q": 30.0, "r": 0.07}, {"f_hz": 40300.0, "q": 100.0, "r": r4}]
                    p = {**BASE, "ki_per_s": BASE["ki_per_s"] * fu / 1e3}
                    a = block({**p, "modes": [{"f_hz": 25200.0, "q": qm}]})["stable"]
                    b = bool(rho_indep({**p, "modes": modes})[0] < 1.0)
                    agree += a == b
                    opt += a and not b
                    pes += b and not a
    print(f"  scan f_u × Q_main × mm2_r × mm4_r (108 points): agree {agree}, check says stable but loop unstable {opt}, "
          f"check says unstable but loop stable {pes}")


def dsc1_composed():
    """controller/dsc1 then gate/borealis_window (D90-4, loads silently): pid1f runs the dsc1 law, the gated check keeps
    tustin. Verdicts over f_u × f_res × Q at the dsc1 latency (delay 5.7 µs)."""
    import itertools
    print("== (e) controller/dsc1 + gate/borealis_window: check tustin (as composed) vs dsc1 (time loop) ==")
    diff = []
    for fu, fr, q in itertools.product([1e3, 5.477e3, 3e4], [5e3, 3.25e4, 6e4], [3.0, 17.32, 100.0]):
        p = {**BASE, "ki_per_s": BASE["ki_per_s"] * fu / 1e3, "modes": [{"f_hz": fr, "q": q}], "delay_s": 5.7e-6}
        a = block(p)
        b = block({**p, "pid_discretisation": "dsc1"})
        if a["stable"] != b["stable"]:
            diff.append((fu, fr, q, round(a["rho"], 4), round(b["rho"], 4)))
    print(f"  27 points: verdict differs at {len(diff)}: (f_u, f_res, Q, ρ_tustin, ρ_dsc1) {diff}")


if __name__ == "__main__":
    main()
    multimode_with_gate()
    dsc1_composed()
