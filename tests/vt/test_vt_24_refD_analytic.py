"""VT-24 (WP-I2): Ref D pre-blind analytic checks on the refD preset (D04; T15 §7; notes/00_summary.md
§4 VT-24). No measured Ref D data are used (`blind/` is never read).

Every check runs the preset (``refD_nominal`` + the overrides named in each test) through
``runner.run`` and compares with the closed forms of T15:

1.  static fringe (DSC1 output swept ±4 V, T_opt = 1): v = V_0 cos φ_D (eq. 3) to 1e−12 V, amplitude
    ±3.33 V, span k_act·g·8 V (8.26 fringes);
2.  lock point / slope: ε = 0.05 → φ_lock − π/2 = 50.02 mrad; V_off = 15 mV → 4.505 mrad (eq. 4);
    K_d = −V_0𝒱√(1 − c²) (eq. 5) against a numerical derivative of the block to 1e−6;
3.  RIN: S = V_off ⇒ no RIN in the true phase (< 1e−6 of the S − V_off = 15 mV case); S − V_off = 15 mV
    ⇒ S_φ = |T|²(15 mV/|K_d|)² S_RIN below f_c within the Welch error (eq. 6, 11);
4.  P-only instability, rigid 915B, T15 §2.9 model (914 one pole 28 kHz, 120 / 100 kHz poles,
    τ = 5.67 µs, K_0 = 21.6): K_p,u = 0.093 and 40.9 kHz, within 2 % (0.98 K_p,u decays, 1.02 K_p,u
    oscillates at 40.9 kHz);
5.  integrator loop residual σ² = πS₀/(2 f_c) for S_env = S₀/f² within 5 % (T15 test 7);
6.  DAC limit cycle (rigid 915B; D65-3, which replaces D56-7): the quantiser is a relay on the code
    threshold, so the limit cycle runs at the loop's phase crossover f₁₈₀ (≈ 17 kHz). Fixed operating
    point: band ≤ 1 LSB = 12.7 mrad; rms ≤ 1.2·(√2/π)·k·q·|P(i2πf₁₈₀)| (≈ 5.7 mrad; P = DAC output →
    optical phase, DC 1, from the preset's discrete chain); slow ramp: band ≤ 2 LSB, same rms bound.
    ±LSB/2 bounds neither the instantaneous error nor (with actuator dynamics) the rms. With the nominal
    915B resonance the limit cycle is amplified (recorded only, T15 test 10);
7.  PER: x = 0.01, maximum lock-point shift arcsin(x/(1 − x)) = 10.10 mrad to 1e−6, |μ| ∈ [0.98, 1];
    the per_fringe → dc_fringe slot path gives the shift −arg μ(x, Δχ);
8.  loop sign flipped ⇒ +π/2 unstable, lock at −π/2 (T15 test 11);
9.  laser × ΔL: in-loop variance + ledger out-of-band variance = 2πΔντ = 0.0277 rad² (Δν = 100 kHz,
    ΔL = 9 m) within 5 % (T15 test 6);
10. the preset's DSC1 chain (P only) reaches −180° at 60 ± 3 kHz (DSC1 p.2: > 58 kHz, typ 60 kHz).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parents[1] / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

FS = 1.0e6
T_OPT1 = {"coupler_excess_db": 0.0}  # T_opt = 1: V_0 = G ℜ P_in = 3.33 V (T15 §6.1)
LINEAR = {**rs.NOISE_OFF, "pdb450c.per_crosscoupling": False, "pdb450c.visibility_hf_averaging": False}
REC = lambda slots, decim=1: {"record": {"slots": slots, "decim": decim}}  # noqa: E731


def _k_act_g(spec):
    leaves, params = rs.built_leaves(spec)
    return leaves["stretcher915b"].gain_rad_per_v(params), leaves["dsc1.dac"].gain()


# --------------------------------------------------------------------------------------------- 1


def test_1_static_fringe():
    dur = 0.05
    spec = rs.load(params=T_OPT1, run={"duration_s": dur}, toggles={**LINEAR, "stretcher915b.resonance": False},
                   measure=REC(["lock.lockD.err", "v.pdb", "v.dsc1.cv"]))
    ramp = {"block": "ramp", "name": "dsc1.ramp", "params": {"out": "v.dsc1.cv", "waveform": "sawtooth",
                                                             "amplitude_v": 4.0, "offset_v": 0.0,
                                                             "period_s": dur * 1.0001, "toggles": {"enabled": True}}}
    spec = rs.replace_block(spec, "dsc1.pid", [ramp])
    res = run(spec)
    phid = res.series["lock.lockD.err"] + math.pi / 2
    v = res.series["v.pdb"]
    v0 = rs.fringe(spec)["v0"]
    assert v0 == pytest.approx(3.33, rel=1e-12)
    assert np.max(np.abs(v - v0 * np.cos(phid))) < 1e-12
    assert v.max() == pytest.approx(3.33, rel=1e-5) and v.min() == pytest.approx(-3.33, rel=1e-5)
    k, g = _k_act_g(spec)
    cv = res.series["v.dsc1.cv"]
    mid = slice(int(0.1 * cv.size), int(0.9 * cv.size))  # the ramp after the start step, before the end
    slope = np.polyfit(cv[mid], phid[mid], 1)[0]  # rad/V, independent of the (constant) chain lag
    assert slope == pytest.approx(-k * g, rel=1e-4)  # φ_D = φ_0 − k_act·g·u (conv §5: arm 1 = pr)
    assert k * g * 8.0 == pytest.approx(51.91, abs=0.01)  # swept span ±4 V (g = DAC load divider)
    assert k * 8.0 / (2 * math.pi) == pytest.approx(8.26, abs=0.005)  # T15 test 1 (g = 1)


# --------------------------------------------------------------------------------------------- 2


def _locked_err(spec, t_avg=0.005):
    res = run(spec)
    e = res.series["lock.lockD.err"]
    n = int(round(t_avg * FS / spec.measure["record"]["decim"]))
    return float(np.mean(e[-n:])), float(np.ptp(e[-n:])), res


@pytest.mark.parametrize("case,expected", [("eps", 50.02e-3), ("voff", 4.505e-3)])
def test_2_lock_point_shift(case, expected):
    params = {**T_OPT1, "eps_r": 0.05} if case == "eps" else {**T_OPT1, "v_off0_v": 0.015}
    tog = dict(LINEAR)
    if case == "voff":
        tog["pdb.offset.offset"] = True
    spec = rs.load(params=params, run={"duration_s": 0.04}, toggles=tog, measure=REC(["lock.lockD.err"]))
    mean, ptp, _ = _locked_err(spec)
    exact = math.asin(0.05) if case == "eps" else math.asin(0.015 / 3.33)
    assert mean == pytest.approx(exact, abs=1e-7)
    assert mean == pytest.approx(expected, abs=5e-6)  # T15 §7 test 2 values (rounded)
    assert ptp < 1e-9


def test_2_slope_matches_eq5():
    from cranberry.blocks.err.dc_fringe import fringe_lock_point

    spec = rs.load(params={**T_OPT1, "eps_r": 0.05, "v_off0_v": 0.015})
    p = {**rs.block_params(spec, "pdb450c"), "in_lo": "v.lo", "in_pr": "v.pr", "in_dp": None, "in_voff": None,
         "in_per_x": None, "in_per_dchi": None, "v_off_v": 0.015,
         "toggles": {**rs.block_params(spec, "pdb450c")["toggles"], "shot_noise": False, "electronic_noise": False,
                     "rin": False, "offset_drift": False, "per_crosscoupling": False,
                     "visibility_hf_averaging": False}}
    lay = SlotLayout(FS)
    lay.add("v.lo", (), np.float64, "rad")
    lay.add("v.pr", (), np.float64, "rad")
    blk = build(spec_from_obj([{"block": "dc_fringe", "name": "pdb", "params": p}]), lay, 0, external=("v.lo", "v.pr"))
    st, pr = State(lay), Params(lay)
    v0 = rs.fringe(spec)["v0"]
    phi_lock, k_d = fringe_lock_point(v0, 0.05, 1.0, 0.0, 0.015)

    def v_at(phi):
        st["v.lo"], st["v.pr"] = phi, 0.0
        blk(st, pr)
        return float(st["v.pdb"])

    assert abs(v_at(phi_lock)) < 1e-15  # the lock point of eq. (4)
    dphi = 1e-6
    k_num = (v_at(phi_lock + dphi) - v_at(phi_lock - dphi)) / (2 * dphi)
    c = -0.015 / v0 - 0.05
    assert k_d == pytest.approx(-v0 * math.sqrt(1 - c * c), rel=1e-14)
    assert k_num == pytest.approx(k_d, rel=1e-6)


# --------------------------------------------------------------------------------------------- 3


def test_3_rin_coupling_vanishes_at_s_equal_voff():
    s_rin = 1.0e-8  # 1/Hz, raised from −140 dB/Hz so the 15 mV case is far above rounding
    dur, nperseg = 2.0, 32768
    out = {}
    for case, (sp_v, voff) in {"zero": (0.015, 0.015), "fifteen": (0.015, 0.0)}.items():
        tog = {**LINEAR, "dx1.rin": True, "pdb.offset.offset": True}
        spec = rs.load(params={**T_OPT1, "setpoint_v": sp_v, "v_off0_v": voff}, run={"duration_s": dur},
                       blocks=[{"name": "dx1", "params": {"rin_psd": [[2.0e7, s_rin], [1.0e9, s_rin / 10]]}}],
                       toggles=tog, measure=REC(["lock.lockD.err"]))
        res = run(spec)
        x = res.series["lock.lockD.err"][int(0.05 * FS):]
        f, p = signal.welch(x, fs=FS, window="hann", nperseg=nperseg, detrend="constant")
        out[case] = (f, p, spec, (x.size - nperseg) // (nperseg // 2) + 1)
    f, p15, spec, n_avg = out["fifteen"]
    band = (f >= 30.0) & (f <= 300.0)
    v0 = rs.fringe(spec)["v0"]
    sv = 0.015 / math.sqrt(v0**2 - 0.015**2)  # (S − V_off)/|K_d| (eq. 5 with c = 0.015/V_0)
    pid = rs.block_params(spec, "dsc1.pid")
    k, g = _k_act_g(spec)
    f_c = v0 * k * g * abs(pid["ki_per_s"]) / (2 * math.pi)
    t2 = 1.0 / (1.0 + (f[band] / f_c) ** 2)  # |L/(1+L)|² of the integrator-dominated loop (VT-13)
    ratio = p15[band].mean() / (sv**2 * s_rin * t2.mean())
    tol = 3.0 / math.sqrt(n_avg * band.sum() / 1.5)  # Hann 50 % overlap: ~1.5 bins per independent estimate
    assert abs(ratio - 1) < max(tol, 0.03), (ratio, tol)
    assert out["zero"][1][band].mean() < 1e-6 * p15[band].mean()


# --------------------------------------------------------------------------------------------- 4


def _t15_p_only(k0):
    """K_p,u and f_180 of T15 §2.9 (P only, rigid 915B, 914 = 28 kHz pole, τ = 5.67 µs)."""
    def phase(f):
        return (math.atan(f / 1.2e5) + math.atan(f / 1.0e5) + math.atan(f / 2.8e4)
                + 2 * math.pi * f * 5.67e-6)
    from scipy.optimize import brentq

    f180 = brentq(lambda f: phase(f) - math.pi, 1e4, 1e5)
    mag = 1 / math.sqrt((1 + (f180 / 1.2e5) ** 2) * (1 + (f180 / 1.0e5) ** 2) * (1 + (f180 / 2.8e4) ** 2))
    return 1.0 / (k0 * mag), f180


def test_4_p_only_instability_rigid():
    fs = 1.0e7
    params = {**T_OPT1, "theta_coupler_rad": 0.0, "latency_s": 5.67e-6, "dsc1_explicit_delay_s": 0.0,
              "dsc1_i": 0.0, "dsc1_d": 0.0}
    blocks = [{"name": "amp914", "params": {"poles_hz": [2.8e4]}},
              {"name": "dsc1.adc", "params": {"avg_samples": 1}},
              {"name": "dsc1.pid", "rate_div": 1}, {"name": "dsc1.sup", "rate_div": 1},
              {"name": "dsc1.dac", "rate_div": 1}]
    base = rs.load(params=params, run={"duration_s": 0.02, "f_sim_hz": fs}, blocks=blocks,
                   toggles={**LINEAR, "stretcher915b.resonance": False}, measure=REC(["lock.lockD.err"]))
    k, g = _k_act_g(base)
    k0 = rs.fringe(base)["v0"] * k * g
    kpu, f180 = _t15_p_only(k0)
    kpu_t15, _ = _t15_p_only(3.33 * 6.491)
    assert kpu_t15 == pytest.approx(0.093, rel=0.01) and f180 == pytest.approx(40.9e3, rel=0.005)
    offset = -(math.pi / 2) / (k * g)  # P-only equilibrium exactly at +π/2 (θ = 0: φ_D = 0 at rest)
    rms, fpk = {}, {}
    for fac in (0.98, 1.02):
        spec = rs.load(params={**params, "dsc1_p": fac * kpu, "offset_v": offset},
                       run={"duration_s": 0.02, "f_sim_hz": fs}, blocks=blocks,
                       toggles={**LINEAR, "stretcher915b.resonance": False}, measure=REC(["lock.lockD.err"]))
        e = run(spec).series["lock.lockD.err"]
        tail = e[-int(0.005 * fs):]
        rms[fac] = float(np.std(tail))
        spec_f = np.abs(np.fft.rfft((tail - tail.mean()) * np.hanning(tail.size)))
        fpk[fac] = float(np.fft.rfftfreq(tail.size, 1 / fs)[np.argmax(spec_f)])
    assert rms[0.98] < 1e-3, rms  # decays (from the π/2 acquisition step)
    assert rms[1.02] > 0.05, rms  # self-sustained oscillation (limit cycle set by the fringe, ~0.4 rad)
    assert fpk[1.02] == pytest.approx(f180, rel=0.02), (fpk, f180)


# --------------------------------------------------------------------------------------------- 5


def test_5_integrator_residual_closed_form():
    f_c = 2.0e3
    base_params = {"latency_s": 1.0e-6, "dsc1_explicit_delay_s": 0.0, "in_out_poles": False,
                   "dsc1_p": 0.0, "dsc1_d": 0.0, "mech_f_c_hz": 2.0e5}
    blocks = [{"name": "dsc1.adc", "params": {"avg_samples": 1}}, {"name": "dsc1.pid", "rate_div": 1},
              {"name": "dsc1.sup", "rate_div": 1}, {"name": "dsc1.dac", "rate_div": 1}]
    tog = {**LINEAR, "env.mech.mech": True, "stretcher915b.resonance": False, "amp914.poles": False}
    probe = rs.load(params=base_params, blocks=blocks, toggles=tog)
    k, g = _k_act_g(probe)
    k0 = rs.fringe(probe)["v0"] * k * g
    ki = 2 * math.pi * f_c / k0  # L = K_0 K_i / s → unity gain at f_c
    dur = 2.0
    spec = rs.load(params={**base_params, "dsc1_i": ki * probe.offline["params"]["dsc1_h_s"]},
                   run={"duration_s": dur + 0.02}, blocks=blocks, toggles=tog, measure=REC(["lock.lockD.err"], 10))
    res = run(spec)
    e = res.series["lock.lockD.err"][int(0.02 * FS / 10):]
    var_sim = float(np.var(e))
    mech = rs.block_params(spec, "env.mech")
    a1, a2 = rs.segment(spec, "a1"), rs.segment(spec, "a2")
    s0 = mech["h_ref_rad2_hz"] * (a1["L_m"] + a2["L_m"]) / mech["L_ref_m"]  # γ = 1
    closed = math.pi * s0 / (2 * f_c)
    # exact prediction: |S_z|² of the discrete loop × the two-pole S_env of T08 eq. 2.7
    f = np.geomspace(1e-3, FS / 2, 200_000)
    s_env = s0 / mech["f_lo_hz"] ** 2 / ((1 + (f / mech["f_lo_hz"]) ** 2) * (1 + (f / mech["f_c_hz"]) ** 2))
    k_d = -rs.fringe(spec)["v0"]
    s2 = np.abs(rs.sensitivity_z(spec, f, k_d=k_d)) ** 2
    var_pred = float(np.trapezoid(s2 * s_env, f))
    assert closed == pytest.approx(var_pred, rel=0.05)  # T15 test 7 closed form incl. finite corrections
    assert var_sim == pytest.approx(var_pred, rel=0.05), (var_sim, var_pred, closed)


# --------------------------------------------------------------------------------------------- 6


def _limit_cycle_rms_bound(spec, f_lc):
    """D65-3 at the limit-cycle frequency ``f_lc``: (|arg L + π| [deg], |P|, 1.2·(√2/π)·k·q·|P|), from the
    preset's discrete loop (refD_support). P = DAC output → optical phase (DC 1)."""
    v0 = rs.fringe(spec)["v0"]
    a, b, c, g = rs.loop_parts(spec, np.array([f_lc]), k_d=-v0)  # K_d = −V_0 at +π/2 (T15 eq. 5, θ = 0)
    lg = complex(-g * c[0] * b[0] * a[0])  # closed loop 1/(1 + L) (refD_support: 1/(1 − gCBA))
    leaves, params = rs.built_leaves(spec)
    k_act = leaves["stretcher915b"].gain_rad_per_v(params)
    p_abs = float(abs(a[0]) / (k_act * leaves["dsc1.pid"].rate_div))  # A = hold (DC r)·LPF_out·delays·914·915B·k
    lsb = k_act * g * 8.0 / 2**12  # k·q: 12-bit DAC on ±4 V (DSC1 p.2) as phase
    off180 = abs(math.degrees(math.remainder(np.angle(lg) + math.pi, 2 * math.pi)))
    return off180, p_abs, 1.2 * math.sqrt(2) / math.pi * lsb * p_abs


@pytest.mark.parametrize("ramp", [False, True])
def test_6_dac_limit_cycle_bound(ramp, record_property):
    """Rigid 915B (D65-3). The nominal 20 kHz / Q 10 resonance amplifies the 2-code limit cycle: that
    case is recorded (no pass/fail, T15 test 10), see test_6_dac_limit_cycle_with_resonance_recorded."""
    tog = {**LINEAR, "dsc1.adc.quantization": True, "dsc1.dac.quantization": True,
           "stretcher915b.resonance": False}
    params = dict(T_OPT1)
    if ramp:
        tog["env.thermal.drift"] = True
        params.update(r_T_k_per_s=1.0e-3, tau_th_s=0.5)  # φ ramp ≈ 0.4 rad/s: ~30 LSB crossed in 1 s
    dur = 1.2 if ramp else 0.1
    spec = rs.load(params=params, run={"duration_s": dur}, toggles=tog, measure=REC(["lock.lockD.err"]))
    res = run(spec)
    e = res.series["lock.lockD.err"][int(0.05 * FS):]
    k, g = _k_act_g(spec)
    lsb = k * g * 8.0 / 2**12  # 12-bit DAC on ±4 V (DSC1 p.2): 12.7 mrad
    adc = 8.0 / 2**16 / rs.fringe(spec)["v0"]  # ADC LSB as phase (37 µrad scale)
    f, psd = signal.welch(e - e.mean(), fs=FS, nperseg=e.size)
    f_lc = float(f[np.argmax(psd)])  # limit-cycle frequency (measured 21 kHz fixed point, 18 kHz ramp)
    off180, p_abs, rms_max = _limit_cycle_rms_bound(spec, f_lc)
    rms = float(np.std(e))
    record_property("band_over_lsb", float(np.ptp(e) / lsb))
    record_property("max_abs_over_lsb", float(np.max(np.abs(e)) / lsb))
    record_property("rms_rad", rms)
    record_property("rms_bound_rad", rms_max)
    record_property("f_limit_cycle_hz", f_lc)
    assert lsb == pytest.approx(12.67e-3, rel=0.002)
    # relay describing function: the cycle sits where the loop phase is near −180°. The discrete loop
    # comes closest to −180° (−171° at 22 kHz) below the PID's I/D notch at √(K_i/K_d)/2π = 28 kHz; the
    # sampled relay adds up to ~h/2 of lag, so allow 20°
    assert 15e3 < f_lc < 25e3 and off180 < 20.0, (f_lc, off180)
    assert 5.0e-3 < rms_max < 6.0e-3, rms_max  # ≈ 5.7 mrad (review: |P| 0.83 at 17.5 kHz)
    assert rms <= rms_max, (rms, rms_max)  # measured 3.1 mrad (fixed point), 2.8 mrad (ramp)
    if not ramp:
        # fixed operating point: the error stays within one 1-LSB band (measured 0.74 LSB, max |e| 0.53 LSB;
        # the max |e| itself lies anywhere between LSB/2 and 1 LSB depending on the operating point)
        assert np.ptp(e) <= lsb + 2 * adc
    else:
        # the operating point crosses the code thresholds: band ≤ 2 LSB (measured 1.24 LSB)
        assert np.ptp(e) <= 2 * lsb


def test_6_dac_limit_cycle_with_resonance_recorded(record_property):
    tog = {**LINEAR, "dsc1.adc.quantization": True, "dsc1.dac.quantization": True}
    spec = rs.load(params=T_OPT1, run={"duration_s": 0.1}, toggles=tog, measure=REC(["lock.lockD.err"]))
    e = run(spec).series["lock.lockD.err"][int(0.05 * FS):]
    k, g = _k_act_g(spec)
    lsb = k * g * 8.0 / 2**12
    record_property("resonant_band_over_lsb", float(np.ptp(e) / lsb))
    record_property("resonant_rms_rad", float(np.std(e)))
    assert np.isfinite(e).all() and abs(float(np.mean(e))) < lsb  # still locked; amplitude recorded only


# --------------------------------------------------------------------------------------------- 7


def test_7_per_lock_point_shift():
    x = 0.01
    dchi_star = math.acos(-x / (1 - x))
    shifts = {}
    for dchi in (0.0, dchi_star, math.pi):
        blocks = [{"name": "pdb450c", "params": {"in_per_x": None, "in_per_dchi": None, "per_x": x,
                                                 "per_dchi_rad": dchi}}]
        spec = rs.load(params=T_OPT1, run={"duration_s": 0.04}, blocks=blocks,
                       toggles={**LINEAR, "pdb450c.per_crosscoupling": True}, measure=REC(["lock.lockD.err"]))
        shifts[dchi], ptp, _ = _locked_err(spec)
        mu = (1 - x) + x * np.exp(1j * dchi)
        assert shifts[dchi] == pytest.approx(-np.angle(mu), abs=1e-8)
        assert 0.98 - 1e-15 <= abs(mu) <= 1.0
    assert abs(shifts[dchi_star]) == pytest.approx(math.asin(x / (1 - x)), abs=1e-6)
    assert abs(shifts[dchi_star]) == pytest.approx(10.10e-3, abs=5e-6)


def test_7_per_fringe_slot_path():
    tog = {**LINEAR, "pdb450c.per_crosscoupling": True}
    spec = rs.load(run={"duration_s": 0.04}, toggles=tog, measure=REC(["lock.lockD.err", "pol.per.x", "pol.per.dchi"]))
    mean, ptp, res = _locked_err(spec)
    x, dchi = float(res.series["pol.per.x"][-1]), float(res.series["pol.per.dchi"][-1])
    assert 1e-4 < x < 0.02  # coupler-port PER 25 dB at both ends of both arms (nominal)
    assert mean == pytest.approx(-np.angle((1 - x) + x * np.exp(1j * dchi)), abs=1e-7)


# --------------------------------------------------------------------------------------------- 8


def test_8_loop_sign_flip_locks_at_minus_half_pi():
    for sign, target in ((1.0, 0.0), (-1.0, -math.pi)):
        spec = rs.load(params={**T_OPT1, "loop_sign": sign}, run={"duration_s": 0.04}, toggles=LINEAR,
                       measure=REC(["lock.lockD.err"]))
        mean, ptp, res = _locked_err(spec)
        wrapped = (mean - target + math.pi) % (2 * math.pi) - math.pi
        assert abs(wrapped) < 1e-7, (sign, mean)
        if sign < 0:  # +π/2 is unstable: the loop runs away from φ_D(0) ≈ 0 to −π/2
            assert abs(mean) > 3.0


# --------------------------------------------------------------------------------------------- 9


def test_9_laser_delta_l_total_variance():
    from cranberry.blocks.noise.laser import laser_oob_variance

    h0 = 1.0e5 / math.pi  # Δν = 100 kHz
    tog = {k: v for k, v in rs.NOISE_OFF.items() if k != "dx1.fm_white"}
    spec = rs.load(params={"h0_hz2_per_hz": h0, "dL_m": 9.0}, run={"duration_s": 0.5},
                   toggles={**tog, "dsc1.pid.feedback": False}, measure=REC(["lock.lockD.err"]))
    res = run(spec)
    var_in = float(np.var(res.series["lock.lockD.err"]))
    tau = spec.offline["params"]["n_g"] * 9.0 / rs.C_M_PER_S
    decl = [d for d in res.offline.band_ledger.decls if d.source == "dx1" and d.quantity == "freq"]
    assert len(decl) == 1 and decl[0].targets == ("lockD",)
    f_n = decl[0].f_hi_hz
    var_out = laser_oob_variance(lambda f: h0 * np.ones_like(np.asarray(f, float)), tau, f_n, 1e10) + 2 * h0 / 1e10
    total = 2 * math.pi * math.pi * h0 * tau  # 2πΔντ
    assert total == pytest.approx(0.0277, rel=0.01)  # T15 test 6 (n_g = 1.468)
    assert var_in == pytest.approx((2 * math.pi * tau) ** 2 * h0 * f_n, rel=0.05)
    assert var_in + var_out == pytest.approx(total, rel=0.05)


# -------------------------------------------------------------------------------------------- 10


def test_10_dsc1_chain_minus_180_deg_frequency():
    spec = rs.load()
    f = np.linspace(4e4, 8e4, 4001)
    a, b, c, g = rs.loop_parts(spec, f, k_d=1.0)
    st = rs.built_leaves(spec)[0]["stretcher915b"]
    leaves, params = rs.built_leaves(spec)
    k_act = leaves["stretcher915b"].gain_rad_per_v(params) * float(params[st.param_name("polarity")])
    # DSC1 input → output (P only): sensor part (LPF_in, average, decimation ½) × hold × LPF_out × closure
    # × delay buffer, i.e. A·B without the 914, 915B, transit and plant sign
    amp = leaves["amp914"].freq_response(f)
    h_str = leaves["stretcher915b"].freq_response(f)
    a1 = rs.segment(spec, "a1")
    from cranberry.path.segments import frac_delay

    k1, mu1 = frac_delay(0.5 * a1["n_g"] * a1["L_m"] / rs.C_M_PER_S, 1.0 / FS)
    zi = np.exp(-2j * np.pi * f / FS)
    # the 1 MHz sequence drives the zoh-discretised 914: its half-sample ZOH interface delay is part
    # of the DSC1 output timing (physical 500 kHz ZOH = discrete hold (1 + z⁻¹) + this half sample)
    zoh = np.exp(-1j * np.pi * f / FS)
    h_dsc1 = zoh * a * b / (amp * k_act * h_str * (-1.0) * zi**k1 * ((1 - mu1) + mu1 * zi))
    ph = np.unwrap(np.angle(h_dsc1))
    f180 = float(np.interp(-np.pi, ph[::-1], f[::-1]))
    assert 57e3 <= f180 <= 63e3, f180
