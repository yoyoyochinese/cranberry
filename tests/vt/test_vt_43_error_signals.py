"""VT-43: error signals (S §4 VT-43 ← T14 §11 V1, V2, V3 (static part), V4–V7, V10–V12; T15 §7 V1, V2,
V3 (static part), V9, V11). WP-B3.

Tolerances (T14 §11 / T15 §7 / PLAN §6.3): shapes 1e−6 (V1, V4, V5), noise PSD 5 % (V2), radiometer
1/(BT) 5 %, κ_det linear 1.02 ± 0.1, log 1.32 ± 0.15 (V6), baseband vs RF-rate 10 % (V6) / 1e−3
deterministic, closed-loop sensor-limited residual 10 % (V7), Liu coefficients (V10), DC fringe 1e−12 /
1e−6 (V11), offsets 10 % (V12). Stochastic checks use fixed seeds and sample sizes whose 3σ error is
below the tolerance.

Ref B numbers are recomputed with x = √(P/P_th) = √(100/180) = 0.7454 (D51-10, T02), not T14's
x = 0.769 (G = 18.7): A−/A+ = (1 − x)/(1 + x) = 0.1459, so a lock-1 error δu = 0.1 rad moves the
lock-2 zero by 14.64 mrad (T14: 13.1 mrad).

Ref C quantum-noise lock (V7, D65-1): the radiometer noise is cyclostationary over the LO dither cycle;
the demodulated noise is 2⟨2 sin²θ S_v(Δ0 + m sinθ)⟩. For the square law this is 1.2257 × the value
T14 §6.3 uses (S_v at Δ0) at the squeezing point: 3.58 mrad in a continuous loop (3.62 mrad in the
test's 1 MHz discrete loop), not 3.24 mrad. For the log detector the weight is exactly 1 and the slope
falls instead: 4.11 mrad (κ_det = 1.32), not T14's 4.2. T14's formula is recovered where the weighting
is negligible (small m, or the anti-squeezing point: 78 mrad). V12: the LO RAM does not modulate the
dark noise, so the Ref C RAM shift is 2.54 mrad (D65-1, M1-A-03).
"""

import math

import numpy as np
import pytest
from scipy.optimize import brentq, minimize_scalar
from scipy.special import j0, j1

import cranberry.core.testing  # noqa: F401  (toy blocks for the loops)
from cranberry.blocks.err import fringe_lock_point, radiometer_mc
from cranberry.blocks.err.testing import (
    REF_A_DET,
    REF_A_SQ,
    REF_B,
    REF_C,
    X_REF_B,
    bhd,
    dcf,
    dem,
    liu,
    pbd,
    qnl,
    qnl_flat,
    record,
    rig,
)
from cranberry.core.fft import psd_onesided
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba
E = 1.602176634e-19
KV = 1e3 * 1.0 * 0.8 * 1e-3  # G ℜ η P_in of the pbd test parameters


def evaluator(objs, f_sim=16e6, external=None):
    kw = {} if external is None else {"external": external}
    root, s, p = rig(objs, f_sim, **kw)

    def ev(**phases):
        for k, v in phases.items():
            s[k] = v
        root(s, p)
        return s

    return ev, s, p


def band_mean(x, fs, lo, hi, nperseg=4096):
    f, s = psd_onesided(x, fs, nperseg=nperseg)
    sel = (f >= lo) & (f <= hi)
    return float(np.mean(s[sel]))


def dl(name, inp, out, unit="rad"):
    return {"block": "delay_buffer", "name": name, "params": {"in": inp, "out": out, "tau_s": 0.0,
                                                               "method": "linear", "unit": unit}}


def integ(name, inp, k, enabled=True):
    return toy("toy_integrator", name, in_=inp, out=f"ctrl.{name}.u", unit="rad", k=k, enabled=enabled)


# ------------------------------------------------------------------ V1: Ref A lock 1 shape at RF rate


def test_v1_refA_lock1_rf_rate_shape_and_no_fsh_component():
    phi1, phid = 0.7, 0.3
    n = 8 * 2000  # 2000 periods of the 2 MHz carrier at 16 MHz
    rec = record([pbd(rf=True), dem(rf=True, phi_d_rad=phid), dem("d1", "v.pbd", rf=True, f_ref_hz=1e6)],
                 16e6, n, ["dem.dem.i", "dem.d1.i", "dem.d1.q"], init={"beam.pr.phi": phi1 / 2})
    ap, am = REF_A_DET
    k = KV * 0.5 * (ap**2 - am**2)  # G ℜ η P_in sinh 2r_d
    assert np.mean(rec["dem.dem.i"]) == pytest.approx(k * math.cos(phi1 + phid), rel=1e-6)
    assert abs(np.mean(rec["dem.d1.i"])) < 1e-9 * k and abs(np.mean(rec["dem.d1.q"])) < 1e-9 * k


# ------------------------------------------------------------------ V2: noise through the demodulator


def test_v2_shot_noise_psd_2S_and_baseband_equals_rf_rate():
    ap, am = REF_A_DET
    s_v = 1e3**2 * 2 * E * 1.0 * 0.8e-3 * 0.5 * (ap**2 + am**2)  # G² 2eℜ P_DC
    tg = {"electronic_noise": False}
    n = 1 << 18
    out = {}
    for label, objs, fs in (
        ("rf", [pbd(rf=True, noise=True, toggles=tg), dem(rf=True)], 16e6),
        ("bb16", [pbd(noise=True, toggles=tg), dem()], 16e6),
        ("bb1", [pbd(noise=True, toggles=tg), dem()], 1e6),
    ):
        x = record(objs, fs, n, ["dem.dem.i"], seed=3)["dem.dem.i"]
        out[label] = band_mean(x, fs, 2e4, 2e5 if fs < 2e6 else 1e6, nperseg=2048 if fs < 2e6 else 4096)
    for label, v in out.items():
        assert v == pytest.approx(2 * s_v, rel=0.05), label


# ------------------------------------------------------------------ V3 (static): lock-2 reduction A−/A+


def zero_of(ev, var, f, lo, hi):
    return brentq(lambda x: f(ev(**{var: x})), lo, hi, xtol=1e-15)


def test_v3_lock2_zero_moves_by_A_minus_over_A_plus():
    du = 0.1
    # Ref A, heterodyne at f_sh: u′ = u + φ_d2 = δu′ (u = 0); zero of I vs Δ = φ_LO (φ_p = 0)
    ev, s, _ = evaluator([bhd(), dem(inp="v.bhd", f_ref_hz=1e6, phi_d_rad=du)])
    eps = zero_of(ev, "beam.lo.phi", lambda st: float(st["dem.dem.i"]), -0.05 + math.pi / 2, 0.05 + math.pi / 2) - math.pi / 2
    ratio_a = REF_A_SQ[1] / REF_A_SQ[0]
    assert 1 / ratio_a == pytest.approx(123.47, rel=1e-3)
    assert eps == pytest.approx(math.atan(math.tan(du) * ratio_a), rel=1e-3)
    assert eps == pytest.approx(0.813e-3, rel=1e-3)
    # Ref B, dither lock 2 with lock 1 at max de-amplification u0 = π/2 + δu (x = 0.7454, recomputed)
    ev, s, _ = evaluator([bhd(mode="dither")])
    ev(**{"beam.pr.phi": math.pi / 2 + du})
    eps_b = zero_of(ev, "beam.lo.phi", lambda st: float(st["v.bhd"][1]), math.pi / 2 - 0.05, math.pi / 2 + 0.05) - math.pi / 2
    ratio_b = (1 - X_REF_B) / (1 + X_REF_B)
    assert REF_B[1] / REF_B[0] == pytest.approx(ratio_b, rel=1e-12) and ratio_b == pytest.approx(0.1459, abs=1e-4)
    assert eps_b == pytest.approx(math.atan(ratio_b * math.tan(du)), rel=1e-3)
    assert eps_b == pytest.approx(14.64e-3, rel=1e-3)


# ------------------------------------------------------------------ V4: Ref B lock 1 (dither, DC extremum)


@pytest.mark.parametrize("beta", [0.1, 0.4, 1.0])  # β = 2m, m = 0.05 … 0.5 rad
def test_v4_refB_lock1_shape_bessel(beta):
    ev, _, _ = evaluator([pbd(mode="dither", dither_m_rad=beta)])
    ap, am = REF_B
    dg = 0.5 * (ap**2 - am**2)
    amp = 2 * KV * dg * j1(beta)
    for phi1 in np.linspace(-math.pi, math.pi, 13):
        a, b = ev(**{"beam.pr.phi": phi1 / 2})["v.pbd"]
        assert b == pytest.approx(-amp * math.sin(phi1), abs=1e-6 * amp)
        assert abs(a) < 1e-12 * amp
    h = 1e-6
    slope = lambda p0: (ev(**{"beam.pr.phi": (p0 + h) / 2})["v.pbd"][1] - ev(**{"beam.pr.phi": (p0 - h) / 2})["v.pbd"][1]) / (2 * h)  # noqa: E731
    s0, spi = slope(0.0), slope(math.pi)
    assert s0 == pytest.approx(-amp, rel=1e-6) and spi == pytest.approx(amp, rel=1e-6)  # one sign: one stable zero


# ------------------------------------------------------------------ V5: quantum noise lock error signal


def test_v5_qnl_square_law_shape_and_log_slope_ratio():
    m = REF_C["m"]
    q = qnl_flat()
    ev, _, _ = evaluator([qnl()], external={"beam.lo.phi": "rad", "beam.p.phi": "rad"})
    amp = 2 * q["p_delta"] * j1(2 * m)
    for d0 in np.linspace(-1.5, 1.5, 11):
        b = ev(**{"beam.lo.phi": d0 + math.pi / 2})["v.qnl"][1]
        assert b == pytest.approx(amp * math.sin(2 * d0), abs=1e-6 * amp)
    for d0 in (0.0, math.pi / 2):
        assert abs(ev(**{"beam.lo.phi": d0 + math.pi / 2})["v.qnl"][1]) < 1e-9 * amp
    # log detector: slopes at the two lock points differ by P+/P− (first order in m; m = 0.01 here)
    ev, _, _ = evaluator([qnl(law="log", dither_m_rad=0.01)], external={"beam.lo.phi": "rad", "beam.p.phi": "rad"})
    h = 1e-6

    def slope(d0):
        f = lambda d: ev(**{"beam.lo.phi": d + math.pi / 2})["v.qnl"][1]  # noqa: E731
        return (f(d0 + h) - f(d0 - h)) / (2 * h)

    p_minus, p_plus = q["p_mean"] - q["p_delta"], q["p_mean"] + q["p_delta"]
    assert abs(slope(0.0) / slope(math.pi / 2)) == pytest.approx(p_plus / p_minus, rel=2e-3)


# ------------------------------------------------------------------ V6: radiometer noise


@pytest.mark.parametrize("law,target,tol", [("square", 1.0, 0.05), ("linear", 1.02, 0.10), ("log", 1.32, 0.15)])
def test_v6_radiometer_mc_relative_variance(law, target, tol):
    """8–300 MHz flat band, f_s = 1 GHz, T = 0.5 µs (BT = 146), 12000 blocks (1σ of the variance ≈ 1.3 %)."""
    rng = np.random.default_rng(20261031)
    kappa = radiometer_mc(rng, f1_hz=8e6, f2_hz=300e6, fs_hz=1e9, t_avg_s=0.5e-6, n_blocks=12000, law=law)
    assert kappa == pytest.approx(target, abs=tol)


def cyclo_weight(d0, m, law="square", n=4096):
    """Independent numpy quadrature of 2⟨2 sin²θ S_v(Δ0 + m sinθ)⟩ / (2 S_v(Δ0)) for the flat Ref C band."""
    q = qnl_flat()
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    d = d0 + m * np.sin(th)
    p = q["p_mean"] - q["p_delta"] * np.cos(2 * d)
    qq = q["q0"] + q["q1"] * np.cos(2 * d) + q["q2"] * np.cos(2 * d) ** 2
    dv2 = {"square": 1.0, "linear": 0.25 / p, "log": (10 / np.log(10)) ** 2 / p**2}[law]
    sv = dv2 * 2 * qq
    p0 = q["p_mean"] - q["p_delta"] * math.cos(2 * d0)
    q0 = q["q0"] + q["q1"] * math.cos(2 * d0) + q["q2"] * math.cos(2 * d0) ** 2
    dv20 = {"square": 1.0, "linear": 0.25 / p0, "log": (10 / np.log(10)) ** 2 / p0**2}[law]
    return float(np.mean(2 * np.sin(th) ** 2 * sv) / (dv20 * 2 * q0)), dv20 * 2 * q0


def test_v6_qnl_baseband_equals_rf_rate_noise_psd():
    """Error-signal PSD of the equivalent Gaussian model: baseband (analytic dither) vs RF rate (time-series
    dither at 1 MHz, f_sim = 16 MHz) within 10 %, both equal to 2⟨2 sin²θ S_v⟩."""
    m, fs, n = REF_C["m"], 16e6, 1 << 18
    ext = {"beam.p.phi": "rad"}
    dith = toy("toy_sine", "dith", out="beam.lo.phi", unit="rad", amp=m, freq_hz=1e6)
    init = {"beam.p.phi": -math.pi}  # Δ = φ_LO
    rf = record([dith, qnl(rf=True, noise=True), dem(inp="v.qnl", rf=True, f_ref_hz=1e6, phi_d_rad=-math.pi / 2)],
                fs, n, ["dem.dem.i"], external=ext, init=init, seed=7)["dem.dem.i"]
    bb = record([qnl(noise=True), dem(inp="v.qnl", f_ref_hz=1e6, phi_d_rad=-math.pi / 2)], fs, n, ["dem.dem.i"],
                external={"beam.p.phi": "rad", "beam.lo.phi": "rad"}, init=init, seed=7)["dem.dem.i"]
    w, sv0 = cyclo_weight(0.0, m)
    s_rf, s_bb = band_mean(rf, fs, 1e4, 4e5), band_mean(bb, fs, 1e4, 4e5)
    assert w == pytest.approx(1.2265, rel=2e-3)
    assert s_bb == pytest.approx(2 * w * sv0, rel=0.05)
    assert s_rf == pytest.approx(s_bb, rel=0.10)


# ------------------------------------------------------------------ V7: sensor-limited lock residual


def qnl_loop(d0, m, n, seed, law="square", kappa=1.0):
    """QNL (baseband) → demod (I = b) → integrator → (delay 0) → φ_LO; Δ = φ_LO (φ_p = −π).
    First-order loop with noise bandwidth b_loop = π f_c/2 = 10 kHz at f_sim = 1 MHz.

    Returns (simulated rms, discrete-loop prediction, continuous-loop prediction, T14 stationary form).
    Square law: slope 4ΔP J1(2m), noise 2⟨2 sin²θ S_v⟩ (cyclostationary weight w). Log law (k = 1):
    slope (10/ln10)⟨2 sinθ ∂_Δ ln P(Δ0 + m sinθ)⟩, S_v = κ(10/ln10)²·2Q/P² — constant over a flat band, so
    the noise weight is exactly 1 (D65-1, M1-A-04)."""
    q = qnl_flat()
    th = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    sgn = 1 if d0 == 0.0 else -1
    w, sv0 = cyclo_weight(d0, m, law)
    if law == "square":
        kd = 4 * q["p_delta"] * j1(2 * m) * sgn  # dI/dΔ at the lock point
        kd_stat = kd
    else:
        dd = d0 + m * np.sin(th)
        pp = q["p_mean"] - q["p_delta"] * np.cos(2 * dd)
        c = 10 / np.log(10)
        kd = c * np.mean(2 * np.sin(th) * 2 * q["p_delta"] * np.sin(2 * dd) / pp)
        p0 = q["p_mean"] - q["p_delta"] * math.cos(2 * d0)
        kd_stat = c * 4 * q["p_delta"] * j1(2 * m) * sgn / p0  # T14 §6.2: square-law slope / P(Δ0)
    sv0 *= kappa
    f_c = 2 * 10e3 / math.pi
    objs = [qnl(law=law, noise=True, dither_m_rad=m, kappa_det=kappa),
            dem(inp="v.qnl", f_ref_hz=1e6, phi_d_rad=-math.pi / 2),
            integ("int", "dem.dem.i", -2 * math.pi * f_c / kd), dl("dl", "ctrl.int.u", "beam.lo.phi")]
    x = record(objs, 1e6, n, ["beam.lo.phi"], external={"beam.p.phi": "rad"},
               init={"beam.p.phi": -math.pi, "ctrl.int.u": d0}, seed=seed)["beam.lo.phi"]
    eps = 2 * math.pi * f_c / 1e6
    var_cont = 2 * (w * sv0) / kd**2 * 10e3  # continuous first-order loop, b_loop = 10 kHz
    var_pred = var_cont / (1 - eps / 2)  # exact AR(1) of the discrete loop (factor 1.0204)
    s_t14 = sv0 * 2 / kd_stat**2  # T14 §6.3: S_Δ,sens = 2 S_v(Δ0)/K_d² (rad²/Hz)
    return float(np.std(x[2000:] - d0)), math.sqrt(var_pred), math.sqrt(var_cont), math.sqrt(s_t14 * 10e3)


def test_v7_qnl_sensor_limited_residual():
    """Square law (D65-1): T14's 3.24 mrad → 3.58 mrad in the continuous loop (weight 1.2257), 3.62 mrad
    in this test's discrete loop at f_sim = 1 MHz (factor 1/(1 − ε/2) = 1.0204)."""
    rms, pred, cont, t14 = qnl_loop(0.0, REF_C["m"], 200_000, seed=11)
    assert t14 == pytest.approx(3.24e-3, rel=0.01)  # T14 §6.3 number (S_Δ,sens = 1.05e-9 rad²/Hz)
    assert cyclo_weight(0.0, REF_C["m"])[0] == pytest.approx(1.2257, rel=1e-3)  # square law only
    assert cont == pytest.approx(3.58e-3, rel=0.005)  # continuous loop, cyclostationary weighting
    assert pred == pytest.approx(3.62e-3, rel=0.005)  # discrete loop of this test
    assert rms == pytest.approx(pred, rel=0.10)
    rms, pred, cont, t14 = qnl_loop(math.pi / 2, REF_C["m"], 200_000, seed=12)  # anti-squeezing lock point
    assert t14 == pytest.approx(78e-3, rel=0.01) and rms == pytest.approx(t14, rel=0.10)
    rms, pred, cont, t14 = qnl_loop(0.0, 0.02, 200_000, seed=13)  # small dither: T14's formula itself holds
    assert rms == pytest.approx(t14, rel=0.10)


def test_v7_qnl_sensor_limited_residual_log_detector():
    """Log detector, κ_det = 1.32 (Ref C, D65-1 / M1-A-04): the cyclostationary noise weight is 1 and the
    slope is 0.905 × T14's stationary slope (1.222 in PSD): 4.11 mrad in the continuous loop, not T14's
    4.2 (rms × κ_det, an arithmetic slip) nor the stationary 3.24·√1.32 = 3.72 mrad."""
    rms, pred, cont, t14 = qnl_loop(0.0, REF_C["m"], 200_000, seed=14, law="log", kappa=1.32)
    assert cyclo_weight(0.0, REF_C["m"], "log")[0] == pytest.approx(1.0, rel=1e-9)
    assert t14 == pytest.approx(3.72e-3, rel=0.01)
    assert (t14 / cont) ** 2 == pytest.approx(1 / 1.222, rel=2e-3)
    assert cont == pytest.approx(4.11e-3, rel=0.005)
    assert rms == pytest.approx(pred, rel=0.10)


# ------------------------------------------------------------------ V10: Liu-type locks, integration schemes A/B


F_LIU, N_LIU, F_PSI = 1e5, 60_000, 10.0  # f_sim, samples (0.6 s), slow disturbance 10 Hz
F_C = 2e3  # loop crossover


def liu_rig(scheme, *, psi_amp=0.0, psi_diff=0.0, pump_amp=0.0, lock1=True):
    """Hand-assembled per-path phases (conv §9) for one port k with reference port 0 at φ_pr,0 = 0.

    ψ_f (forward, squeezed light), ψ_b = ψ_f + psi_diff (backward probe); a_k: port actuator (delay 0).
    Observable: θ_sq,k(GBS) − φ_pr,0 − π/2 = φ_p/2 + ψ_f + a_k (− φ_ref for scheme A).
    """
    ap, am = REF_A_DET
    k_het = KV * 0.5 * (ap**2 - am**2)
    objs = [
        toy("toy_sine", "pump", out="beam.p.phi", unit="rad", amp=pump_amp, freq_hz=F_PSI * 0.7),
        toy("toy_sine", "psif", out="psi.k.f", unit="rad", amp=psi_amp, freq_hz=F_PSI),
        toy("toy_sine", "psid", out="psi.k.d", unit="rad", amp=psi_diff, freq_hz=F_PSI * 1.3),
        toy("toy_sum", "psib", ins=["psi.k.f", "psi.k.d"], gains=[1.0, 1.0], out="psi.k.b", unit="rad"),
    ]
    if scheme == "B":  # backward probe returned to the detection OPA: Φ1 = 2(ψ_b + a_k) − φ_p
        objs += [
            toy("toy_sum", "prb", ins=["psi.k.b", "act.k.a"], gains=[1.0, 1.0], out="beam.pr.phi", unit="rad"),
            pbd(), dem(phi_d_rad=math.pi / 2),  # I = −K sin Φ1
            integ("ik", "dem.dem.i", math.pi * F_C / k_het), dl("dlk", "ctrl.ik.u", "act.k.a"),
            toy("toy_sum", "obs", ins=["beam.p.phi", "psi.k.f", "act.k.a"], gains=[0.5, 1.0, 1.0],
                out="metric.k.theta", unit="rad"),
        ]
    else:  # A: lock 1 = common reference vs pump (detection OPA); port k = Liu homodyne (probe vs reference)
        k_liu = 1e4 * 2 * math.sqrt(1e-5 * 1e-3)
        objs += [
            pbd("l1", in_pr="beam.ref.phi"), dem("d1", "v.l1", phi_d_rad=math.pi / 2),
            integ("i1", "dem.d1.i", math.pi * F_C / k_het, enabled=lock1), dl("dl1", "ctrl.i1.u", "beam.ref.phi"),
            toy("toy_sum", "prb", ins=["beam.ref.phi", "psi.k.b", "act.k.a"], gains=[1.0, 1.0, 1.0],
                out="beam.pr.phi", unit="rad"),
            liu(), integ("ik", "v.liu", 2 * math.pi * F_C / k_liu), dl("dlk", "ctrl.ik.u", "act.k.a"),
            toy("toy_sum", "obs", ins=["beam.p.phi", "psi.k.f", "act.k.a", "beam.ref.phi"],
                gains=[0.5, 1.0, 1.0, -1.0], out="metric.k.theta", unit="rad"),
        ]
    init = {"ctrl.ik.u": 0.0 if scheme == "B" else math.pi / 2}
    rec = record(objs, F_LIU, N_LIU, ["metric.k.theta", "beam.p.phi", "act.k.a", "psi.k.f", "psi.k.d"],
                 external=(), init=init)
    off = 0.0 if scheme == "B" else math.pi / 2  # lock constants C_k/2 or C_k (π/2 for the DC homodyne zero)
    sl = slice(10_000, None)
    return {k: v[sl] for k, v in rec.items()}, off


def coef(y, x):
    return float(np.dot(y - y.mean(), x - x.mean()) / np.dot(x - x.mean(), x - x.mean()))


@pytest.mark.parametrize("scheme", ["A", "B"])
def test_v10_port_path_noise_cancels(scheme):
    r, off = liu_rig(scheme, psi_amp=1.0)
    resid = r["metric.k.theta"] - off
    assert np.max(np.abs(resid - resid.mean())) < 1.5 * F_PSI / F_C  # suppressed by ~1/|L| (5e-3 rad)
    assert np.std(r["psi.k.f"]) > 0.6


def test_v10_scheme_B_pump_coefficient_one_actuators_follow_half():
    r, _ = liu_rig("B", pump_amp=0.1)
    assert coef(r["metric.k.theta"], r["beam.p.phi"]) == pytest.approx(1.0, abs=0.02)
    assert coef(r["act.k.a"], r["beam.p.phi"]) == pytest.approx(0.5, abs=0.02)


def test_v10_scheme_A_pump_cancels_with_lock1_half_without():
    r, _ = liu_rig("A", pump_amp=0.1)
    assert abs(coef(r["metric.k.theta"], r["beam.p.phi"])) < 0.02
    r, _ = liu_rig("A", pump_amp=0.1, lock1=False)
    assert coef(r["metric.k.theta"], r["beam.p.phi"]) == pytest.approx(0.5, abs=0.02)


@pytest.mark.parametrize("scheme", ["A", "B"])
def test_v10_forward_backward_difference_remains(scheme):
    r, off = liu_rig(scheme, psi_diff=0.2)
    assert coef(r["metric.k.theta"], r["psi.k.d"]) == pytest.approx(-1.0, abs=0.02)  # residual = ψ_f − ψ_b


# ------------------------------------------------------------------ V11: DC fringe (T14 V11, T15 V1/V2/V3/V9/V11)


DCF_EXT = {"beam.lo.phi": "rad", "beam.pr.phi": "rad", "v.dp": "1"}


def test_v11_dc_fringe_static_fringe_and_slope():
    ev, _, _ = evaluator([dcf(eps_r=0.0, toggles={"saturation": False})], external=DCF_EXT)
    phis = np.linspace(-8.0, 8.0, 41)
    v = np.array([float(ev(**{"beam.lo.phi": x})["v.dcf"]) for x in phis])
    np.testing.assert_allclose(v, 3.33 * np.cos(phis), rtol=1e-12, atol=1e-12)  # amplitude 3.33 V (T15 V1)
    ev, _, _ = evaluator([dcf(eps_r=0.03, v_off_v=0.01, mu_abs=0.9)], external=DCF_EXT)
    vv = lambda x: float(ev(**{"beam.lo.phi": x})["v.dcf"])  # noqa: E731
    np.testing.assert_allclose([vv(x) for x in phis[15:26]], 3.33 * (0.03 + 0.9 * np.cos(phis[15:26])) + 0.01, rtol=1e-12)
    h = 1e-6
    assert (vv(math.pi / 2 + h) - vv(math.pi / 2 - h)) / (2 * h) == pytest.approx(-3.33 * 0.9, rel=1e-6)


@pytest.mark.parametrize("kw,shift", [({"eps_r": 0.05}, 50.02e-3), ({"v_off_v": 0.015}, 4.505e-3)])
def test_v11_lock_point_shift_and_slope(kw, shift):
    ev, _, _ = evaluator([dcf(**kw)], external=DCF_EXT)
    vv = lambda x: float(ev(**{"beam.lo.phi": x})["v.dcf"])  # noqa: E731
    phi_lock = brentq(vv, 1.0, 2.2, xtol=1e-15)
    assert phi_lock - math.pi / 2 == pytest.approx(shift, abs=0.01e-3)
    want_phi, want_kd = fringe_lock_point(3.33, kw.get("eps_r", 0.0), 1.0, 0.0, kw.get("v_off_v", 0.0))
    assert phi_lock == pytest.approx(want_phi, rel=1e-12)
    h = 1e-6
    assert (vv(phi_lock + h) - vv(phi_lock - h)) / (2 * h) == pytest.approx(want_kd, rel=1e-6)


def test_v11_rin_does_not_move_lock_point_when_setpoint_is_electrical_zero():
    """T15 V3 (static part): v − S ∝ (1 + δp) at S = V_off, so δp leaves φ_lock unchanged."""
    for s_minus_voff, expect_zero in ((0.0, True), (0.015, False)):
        obj = dcf(eps_r=0.04, v_off_v=0.01, setpoint_v=0.01 + s_minus_voff, in_dp="v.dp", toggles={"rin": True})
        ev, _, _ = evaluator([obj], external=DCF_EXT)
        lock = lambda dp: brentq(lambda x: float(ev(**{"beam.lo.phi": x, "v.dp": dp})["v.dcf"]) - (0.01 + s_minus_voff),  # noqa: E731
                                 1.0, 2.2, xtol=1e-15)
        d = lock(1e-3) - lock(0.0)
        if expect_zero:
            assert abs(d) < 1e-12
        else:  # δφ = (S − V_off)/|K_d| · δp to first order
            assert d == pytest.approx(s_minus_voff / 3.33 * 1e-3, rel=0.05)


def test_v11_per_crosscoupling_lock_shift():
    """T15 V9: x = 0.01, Δχ over 0–2π: |μ| ∈ [0.98, 1], max lock shift arcsin(x/(1−x)) = 10.10 mrad (1e−6)."""
    x = 0.01

    ev, _, p = evaluator([dcf(per_x=x)], external=DCF_EXT)

    def lock(dchi):
        p["dcf.per_dchi_rad"] = dchi
        return brentq(lambda ph: float(ev(**{"beam.lo.phi": ph})["v.dcf"]), 1.0, 2.2, xtol=1e-15) - math.pi / 2

    dchi_max = math.acos(-x / (1 - x))
    res = minimize_scalar(lambda d: -abs(lock(d)), bounds=(dchi_max - 0.3, dchi_max + 0.3), method="bounded",
                          options={"xatol": 1e-9})
    assert abs(lock(res.x)) == pytest.approx(math.asin(x / (1 - x)), rel=1e-6)
    assert math.asin(x / (1 - x)) == pytest.approx(10.10e-3, abs=0.005e-3)
    ev, _, _ = evaluator([dcf(per_x=x, per_dchi_rad=math.pi)], external=DCF_EXT)
    assert float(ev(**{"beam.lo.phi": 0.0})["v.dcf"]) == pytest.approx(3.33 * (1 - 2 * x), rel=1e-12)  # |μ| = 0.98


@pytest.mark.parametrize("loop_sign,start,end", [(1.0, math.pi / 2 + 0.3, math.pi / 2),
                                                 (-1.0, math.pi / 2 + 0.01, -math.pi / 2)])
def test_v11_loop_sign_selects_the_lock_point(loop_sign, start, end):
    """T15 V11: e_v = loop_sign·(S − v) into an integrator on the arm-1 (915B, "pr") phase."""
    objs = [dcf(out_err="v.err", loop_sign=loop_sign), integ("i", "v.err", 2 * math.pi * 1e3 / 3.33),
            dl("d", "ctrl.i.u", "beam.pr.phi")]
    r = record(objs, 1e5, 5000, ["beam.pr.phi"], external={"beam.lo.phi": "rad"}, init={"beam.lo.phi": start})
    phid = start - r["beam.pr.phi"][-1]
    assert (phid - end + math.pi) % (2 * math.pi) - math.pi == pytest.approx(0.0, abs=1e-6)


# ------------------------------------------------------------------ V12: offsets


@pytest.mark.parametrize("m,exact_tol", [(REF_C["m"], 0.01), (0.02, 0.01)])
def test_v12_qnl_ram_offset(m, exact_tol):
    """a_RAM = 1e−2 (φ_a = 0) moves the QNL lock point by δΔ = a⟨2 sin²θ P_opt(m sinθ)⟩/(4ΔP J1(2m)).

    The LO RAM modulates only the LO-proportional (optical) part P_opt = P − p_dark of the band power,
    not the BHD dark / circuit noise (T14 §6.4; D65-1, M1-A-03). At Ref C's m = 0.0785 the shift is
    2.54 mrad (D65-1). T14 V12's 2.77 mrad (first-order form a R−/(2(R+−R−)J1(2m)) with R− incl. the dark
    noise) and D56-8's 3.06 mrad (cycle average, dark noise modulated too) are superseded; the first-order
    form with the optical R− holds for small m (checked at m = 0.02 within T14's 10 %)."""
    q = qnl_flat()
    ev, _, _ = evaluator([qnl(ram_a=1e-2, dither_m_rad=m, toggles={"ram_offset": True})],
                         external={"beam.lo.phi": "rad", "beam.p.phi": "rad"})
    d = brentq(lambda d0: float(ev(**{"beam.lo.phi": d0 + math.pi / 2})["v.qnl"][1]), -0.05, 0.05, xtol=1e-15)
    rp, rm = REF_C["r_plus"], REF_C["r_minus"]  # optical part (dark noise not modulated)
    t14 = 1e-2 * rm / (2 * (rp - rm) * j1(2 * m))
    th = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    p_opt = q["p_mean"] - q["p_dark"] - q["p_delta"] * np.cos(2 * m * np.sin(th))
    exact = 1e-2 * np.mean(2 * np.sin(th) ** 2 * p_opt) / (4 * q["p_delta"] * j1(2 * m))
    assert abs(d) == pytest.approx(exact, rel=exact_tol)
    if m == REF_C["m"]:
        assert exact == pytest.approx(2.54e-3, rel=0.01)
        # the superseded model (dark noise modulated too) gives D56-8's 3.06 mrad
        w_all = np.mean(2 * np.sin(th) ** 2 * (q["p_mean"] - q["p_delta"] * np.cos(2 * m * np.sin(th))))
        assert 1e-2 * w_all / (4 * q["p_delta"] * j1(2 * m)) == pytest.approx(3.06e-3, rel=0.01)
    else:
        assert abs(d) == pytest.approx(t14, rel=0.10)


@pytest.mark.parametrize("delta", [0.05, -0.2])
def test_v12_heterodyne_demod_phase_moves_lock_point_one_to_one(delta):
    ev, _, p = evaluator([pbd(), dem(phi_d_rad=math.pi / 2)])
    zero = lambda: brentq(lambda f1: float(ev(**{"beam.pr.phi": f1 / 2})["dem.dem.i"]), -0.5, 0.5, xtol=1e-15)  # noqa: E731
    z0 = zero()
    p["dem.phi_d_rad"] = math.pi / 2 + delta
    assert zero() - z0 == pytest.approx(-delta, abs=1e-12)


# ------------------------------------------------------------------ baseband == RF rate (deterministic, with LPF)


@pytest.mark.parametrize("case", ["pbd_het", "bhd_dither", "qnl_dither", "liu_het"])
def test_baseband_equals_rf_rate_after_lpf(case):
    fs, n = 16e6, 200_000
    slow = lambda out, amp: toy("toy_sine", "slow", out=out, unit="rad", amp=amp, freq_hz=200.0)  # noqa: E731
    fd = {"pbd_het": 2e6, "bhd_dither": 295e3, "qnl_dither": 1e6, "liu_het": 1e6}[case]
    m = {"bhd_dither": 0.1, "qnl_dither": REF_C["m"]}.get(case, 0.0)
    dithered = {"bhd_dither": "beam.pr.phi", "qnl_dither": "beam.lo.phi"}.get(case)
    gen = {"pbd_het": lambda rf: pbd(rf=rf), "bhd_dither": lambda rf: bhd(mode="dither", rf=rf),
           "qnl_dither": lambda rf: qnl(rf=rf), "liu_het": lambda rf: liu(mode="heterodyne", rf=rf)}[case]
    out = {"pbd_het": "v.pbd", "bhd_dither": "v.bhd", "qnl_dither": "v.qnl", "liu_het": "v.liu"}[case]
    target = dithered or "beam.pr.phi"
    ext = [k for k in ("beam.pr.phi", "beam.p.phi", "beam.lo.phi", "beam.ref.phi") if k != target]
    init = {"beam.p.phi": -math.pi} if case == "qnl_dither" else {}
    res = {}
    for rf in (False, True):
        objs = [slow("beam.slow.phi", 0.5)]
        if rf and dithered:
            objs += [toy("toy_sine", "dith", out="beam.dith.phi", unit="rad", amp=m, freq_hz=fd),
                     toy("toy_sum", "sum", ins=["beam.slow.phi", "beam.dith.phi"], gains=[1.0, 1.0], out=target, unit="rad")]
        else:
            objs += [toy("toy_sum", "sum", ins=["beam.slow.phi"], gains=[1.0], out=target, unit="rad")]
        objs += [gen(rf), dem(inp=out, rf=rf, order=3, f_ref_hz=fd, phi_d_rad=0.3, lpf_hz=3e3)]
        res[rf] = record(objs, fs, n, ["dem.dem.i", "dem.dem.q"], external=ext, init=init)
    scale = max(np.max(np.abs(res[False]["dem.dem.i"])), np.max(np.abs(res[False]["dem.dem.q"])))
    for key in ("dem.dem.i", "dem.dem.q"):
        diff = np.abs(res[True][key][20_000:] - res[False][key][20_000:])
        assert np.max(diff) < 2e-3 * scale, key
