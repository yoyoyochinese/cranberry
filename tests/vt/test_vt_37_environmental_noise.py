"""VT-37 (WP-B2 part; path-model parts V2, V6, V11, V12 (driver handling) are WP-F6's
``test_vt_37_f6_counterprop.py``): environmental noise sources through ``segment_phase`` (T08 §9).

* V1 thermal coefficient and sign: k_T = (2π/λ)(dn/dT + nα_L) = 40.60 / 43.84 rad/(m·K) (dn/dT = 9.2e−6 /
  1.0e−5, rel 1e−4 on the 4-digit values); a warming segment gives ψ = +k_T L ΔT > 0 (conv §4).
* V3 mechanical PSD: h = 1e−3 rad²·Hz, f_lo = 0.1 Hz, f_c = 1 kHz: Welch estimate within ±10 % of the
  model (T08 eq. 2.7 with the f_lo flattening of ``segment_mech_noise``) in every 1/3-octave band
  1–300 Hz; variance = ∫S df ±5 % (separate run at 20 Hz: the exact discretisation keeps the marginal
  variance at any step).
* V4 length scaling: L × 10 gives PSD × 10 for γ = 1 with independent subsegments (κ = 0) and × 100
  for γ = 2 lumped (±10 %).
* V5 correlation: ρ = 0.9 (constant) on two equal segments: difference PSD 0.2 S, coherence 0.81;
  ρ = 1: difference ≡ 0; correlation off: 2 S (±10 %). Frequency-dependent ρ(f) = ρ0/(1 + (f/f_ρ)²):
  coherence follows ρ(f)².
* V7 range exhaustion: lumped 63 m, τ_th = 10 s, r_T = 1e−4 K/s, k_T = 43.84 rad/(m·K): from steady drift,
  64.9 rad are used up in 235 s (±1 %).
* V9 thermodynamic: Dong Table I, L = 1308 m: S_TM(1 Hz) = 4.36e−10, S_TC(1 Hz) = 2.86e−12 rad²/Hz (±1 %);
  327 → 134 / 68 / 34 m: −3.9 / −6.8 / −9.8 dB; the time series reproduces S_TM + S_TC (±10 %, 1–30 Hz).
* V10 first-order lag: ramp response ΔT = r(t − τ(1 − e^{−t/τ})) exact to 1e−9 (a step response is the
  derivative of the same lag; the block's ambient model has ramps, not steps).
* V11 toggles: every env source off ⇒ ψ ≡ 0 and τ_g = nominal (exact).
* V12 axis coefficients: the same ΔT on a pm_diagonal segment gives ψ_slow − ψ_fast = K_T L ΔT exactly.
"""

import math

import numpy as np
import pytest
import scipy.signal

import cranberry.blocks.noise.testing as nt
import cranberry.path.testing as pt
from cranberry.blocks.noise.mech import s_tc, s_tm, two_pole_psd
from cranberry.blocks.noise.thermal import k_T_rad_per_m_K
from cranberry.core.fft import psd_onesided

pytestmark = pytest.mark.numba
KT = k_T_rad_per_m_K(1548.8e-9, 1.0e-5, 1.468, 5.5e-7)


def sp(topo, **tg):
    t = {"phase": True, "group_delay": True, "pm_axis_drift": True}
    t.update(tg)
    return {"block": "segment_phase", "name": "sp", "params": {"topology": topo, "segments": "all", "toggles": t}}


def run(blocks, slots, n, fs, seed=21):
    rec, _, _ = pt.run_and_record([*blocks, pt.record(slots, n)], fs, n, compiled=True, seed=seed)
    return rec


def band_means(f, s, model, centres):
    out = []
    for fc in centres:
        m = (f >= fc * 2 ** (-1 / 6)) & (f < fc * 2 ** (1 / 6))
        out.append(np.mean(s[m]) / np.mean(model(f[m])))
    return np.array(out)


THIRD_OCT = 10 ** (np.arange(0, 25) / 10)  # 1 … 251 Hz


def test_v1_thermal_coefficient_and_sign():
    assert k_T_rad_per_m_K(1548.8e-9, 9.2e-6, 1.468, 5.5e-7) == pytest.approx(40.60, rel=1e-4)
    assert KT == pytest.approx(43.84, rel=1e-4)
    topo = nt.seg_topo({"a": (1.0, 1, "uniform")}, coef={"dT": KT})
    rec = run([nt.thermal(topo, ["a"], {"a": 1.0}, on=("drift",)), sp(topo)], ["drv.a.dT", "psi.a.b_a"], 2000, 10.0)
    assert rec["drv.a.dT"][-1] > 0 and rec["psi.a.b_a"][-1] > 0
    # ψ reads the driver at the subsegment transit time (segment_phase interpolation, δ = 4.9 ns ≪ 0.1 s)
    np.testing.assert_allclose(rec["psi.a.b_a"], KT * rec["drv.a.dT"], rtol=1e-6)


def _mech_psd(on_corr=False, fs=4096.0, n=(1 << 17) * 60, **over):
    L = 1000.0
    topo = nt.seg_topo({"a": (L, 1, "uniform")}, coef={"strain": nt.K_EPS})
    on = ("mech", "correlation") if on_corr else ("mech",)
    blk = nt.mech(topo, ["a"], on=on, h_ref_rad2_hz=1e-3, L_ref_m=L, **over)
    return run([blk, sp(topo)], ["psi.a.b_a"], n, fs)["psi.a.b_a"]


def test_v3_mech_psd_third_octaves_and_variance():
    fs = 4096.0
    psi = _mech_psd(fs=fs)
    f, s = psd_onesided(psi, fs, nperseg=1 << 17)
    r = band_means(f, s, lambda x: two_pole_psd(x, 1e-3, 0.1, 1e3), THIRD_OCT)
    np.testing.assert_allclose(r, 1.0, atol=0.10)
    psi20 = _mech_psd(fs=20.0, n=2_000_000)
    var = math.pi / 2 * 1e-3 * 1e3 / (0.1 * (0.1 + 1e3))
    assert np.var(psi20) == pytest.approx(var, rel=0.05)


def _scaling_pair(gamma, place, n_sub, kappa=0.0):
    topo = nt.seg_topo({"a": (6.3, n_sub, place), "b": (63.0, n_sub, place)}, coef={"strain": nt.K_EPS})
    sc = {"a": kappa, "b": kappa} if place == "per_sub" else {}
    blk = nt.mech(topo, ["a", "b"], on=("mech",), gamma=gamma, spatial_coherence=sc)
    rec = run([blk, sp(topo)], ["psi.a.b_a", "psi.b.b_b"], 1 << 20, 1024.0)
    f, sa = psd_onesided(rec["psi.a.b_a"], 1024.0, nperseg=8192)
    _, sb = psd_onesided(rec["psi.b.b_b"], 1024.0, nperseg=8192)
    m = (f >= 1.0) & (f <= 100.0)
    return np.mean(sb[m]) / np.mean(sa[m])


def test_v4_length_scaling():
    assert _scaling_pair(1.0, "per_sub", 8) == pytest.approx(10.0, rel=0.10)
    assert _scaling_pair(2.0, "uniform", 1) == pytest.approx(100.0, rel=0.10)


def _pair(rho, f_rho=math.inf, corr=True, n=1 << 20, fs=1024.0):
    L = 9.0
    topo = nt.seg_topo({"a": (L, 1, "uniform"), "b": (L, 1, "uniform")}, coef={"strain": nt.K_EPS})
    on = ("mech", "correlation") if corr else ("mech",)
    blk = nt.mech(topo, ["a", "b"], on=on, rho_0=rho, f_rho_hz=f_rho)
    rec = run([blk, sp(topo)], ["psi.a.b_a", "psi.b.b_b"], n, fs)
    return rec["psi.a.b_a"], rec["psi.b.b_b"]


def test_v5_correlation():
    fs = 1024.0
    a, b = _pair(0.9)
    f, sa = psd_onesided(a, fs, nperseg=8192)
    _, sd = psd_onesided(a - b, fs, nperseg=8192)
    m = (f >= 1.0) & (f <= 100.0)
    assert np.mean(sd[m]) / np.mean(sa[m]) == pytest.approx(0.2, rel=0.10)
    fc, coh = scipy.signal.coherence(a, b, fs, nperseg=8192)
    assert np.mean(coh[m]) == pytest.approx(0.81, abs=0.02)
    a1, b1 = _pair(1.0, n=1 << 14)
    np.testing.assert_array_equal(a1, b1)
    a0, b0 = _pair(0.9, corr=False)
    _, s0 = psd_onesided(a0, fs, nperseg=8192)
    _, sd0 = psd_onesided(a0 - b0, fs, nperseg=8192)
    assert np.mean(sd0[m]) / np.mean(s0[m]) == pytest.approx(2.0, rel=0.10)
    # frequency-dependent ρ(f) = ρ0/(1 + (f/f_ρ)²)
    a2, b2 = _pair(0.95, f_rho=10.0)
    fc, coh2 = scipy.signal.coherence(a2, b2, fs, nperseg=8192)
    rho = 0.95 / (1 + (fc / 10.0) ** 2)
    for lo, hi in ((1.0, 3.0), (8.0, 12.0), (40.0, 80.0)):
        mm = (fc >= lo) & (fc <= hi)
        assert np.mean(coh2[mm]) == pytest.approx(np.mean(rho[mm] ** 2), abs=0.03)


def test_v7_range_exhaustion():
    L, fs = 63.0, 10.0
    topo = nt.seg_topo({"a": (L, 1, "uniform")}, coef={"dT": KT})
    n = 5000
    rec = run([nt.thermal(topo, ["a"], {"a": 10.0}, on=("drift",)), sp(topo)], ["psi.a.b_a"], n, fs)
    t = np.arange(n) / fs
    psi = rec["psi.a.b_a"]
    i0 = int(100 * fs)  # 10 τ_th: steady drift
    i1 = np.argmax(psi - psi[i0] >= 64.9)
    t_sat = t[i1] - t[i0]
    assert t_sat == pytest.approx(64.9 / (KT * L * 1e-4), rel=0.01)
    assert 64.9 / (KT * L * 1e-4) == pytest.approx(235.0, rel=0.01)


def test_v9_thermodynamic():
    assert float(s_tm(1.0, 1308.0, **nt.DONG)) == pytest.approx(4.36e-10, rel=0.01)
    assert float(s_tc(1.0, 1308.0, **nt.DONG)) == pytest.approx(2.86e-12, rel=0.01)
    for L, db in ((134.0, -3.9), (68.0, -6.8), (34.0, -9.8)):
        assert 10 * math.log10(float(s_tm(1.0, 4 * L, **nt.DONG) / s_tm(1.0, 4 * 327.0, **nt.DONG))) == pytest.approx(db, abs=0.05)
    L, fs = 1308.0, 2048.0  # f_N ≫ the conduction corner (190 Hz): no aliasing of its OU
    topo = nt.seg_topo({"a": (L, 1, "uniform")}, coef={"strain": nt.K_EPS})
    rec = run([nt.thermo(topo, ["a"], accumulate=False), sp(topo)], ["psi.a.b_a"], 1 << 22, fs)
    f, s = psd_onesided(rec["psi.a.b_a"], fs, nperseg=1 << 15)
    r = band_means(f, s, lambda x: s_tm(x, L, **nt.DONG) + s_tc(x, L, **nt.DONG), THIRD_OCT[:16])
    np.testing.assert_allclose(r, 1.0, atol=0.10)


def test_v10_first_order_lag():
    topo = nt.seg_topo({"a": (1.0, 1, "uniform")})
    fs, n, tau = 50.0, 3000, 7.0
    rec = run([nt.thermal(topo, ["a"], {"a": tau}, on=("drift",), r_T_k_per_s=1.0)], ["drv.a.dT"], n, fs)
    t = np.arange(n) / fs
    np.testing.assert_allclose(rec["drv.a.dT"], t - tau * (1 - np.exp(-t / tau)), rtol=1e-9, atol=1e-12)


def test_v11_all_env_off_is_identity():
    topo = nt.seg_topo({"a": (63.0, 1, "uniform"), "b": (9.0, 4, "per_sub")}, coef={"dT": KT, "strain": nt.K_EPS})
    blocks = [nt.thermal(topo, ["a", "b"], {"a": 3.0, "b": 3.0}, on=()),
              nt.mech(topo, ["a", "b"], on=(), spatial_coherence={"b": 0.5}),
              nt.thermo(topo, ["a", "b"], on=()), sp(topo)]
    rec = run(blocks, ["psi.a.b_a", "psi.b.b_b", "tau_g.a.b_a"], 3000, 1e3)
    assert np.all(rec["psi.a.b_a"] == 0.0) and np.all(rec["psi.b.b_b"] == 0.0)
    assert np.all(rec["tau_g.a.b_a"] == rec["tau_g.a.b_a"][0])


def test_v12_axis_coefficients():
    k_t, L = -2.05, 9.0
    pol = {"K_T_rad_per_m_K": k_t, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.3e-12, "axis_drivers": ["dT"]}
    seg = pt.seg_spec("p", L, drivers={"dT": "uniform"}, pol_model="pm_diagonal", pol_params=pol)
    beams = [pt.beam_spec("s", [pt.step_spec("p", 1, axis="slow")], "o"),
             pt.beam_spec("f", [pt.step_spec("p", 1, axis="fast")], "o")]
    topo = pt.topo_spec([seg], beams, ["o"], pt.coef_spec("p", "dT", "*", KT * L))
    rec = run([nt.thermal(topo, ["p"], {"p": 2.0}, on=("drift", "random")), sp(topo)],
              ["drv.p.dT", "psi.p.s", "psi.p.f"], 3000, 100.0)
    # ψ_slow = c·ΔT̃ (ΔT̃: the driver sample at the transit time, segment_phase interpolation) → ΔT̃ = ψ_slow/c
    np.testing.assert_allclose(rec["psi.p.s"] - rec["psi.p.f"], k_t * L * rec["psi.p.s"] / (KT * L), rtol=1e-10,
                               atol=1e-16)
    assert np.std(rec["drv.p.dT"]) > 0
