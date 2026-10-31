"""VT-32 (WP-B11): pulsed OPA supermodes and temporal-mode overlaps (00_summary §4 VT-32; T03 §7 V1–V11).

* V1 Gaussian overlap: numerical κ (FFT delay + GDD) vs the closed form, Δt ∈ [−3T0, 3T0],
  χ ∈ [−5, 5]: |κ_num − κ_ana| < 1e−6 (sign of Δt and of arg κ included).
* V2 symplectic: C C† − S S† = 𝟙, C Sᵀ = S Cᵀ to 1e−8 at N = 512.
* V3 low gain: S/g → JSA F (relative 1e−3), K of the supermodes = K of F (1e−3); Gaussian JSA
  K = ½(σ_−/σ_+ + σ_+/σ_−) (1e−3).
* V4 multimode BHD: u_LO = ψ_0 → V_− = e^{−2r_0}/2; u_LO ⟂ every ψ_n → V = 1/2 (1e−9).
* V6 quasi-CW: τ_p = 10, 30, 100 τ_GVM (β2 = 0): r_0/√(αP_pk) increases towards 1, within 1e−2 at 100.
* V7 single-mode mismatch = loss: V = η_t V_sq + (1 − η_t)/2 (1e−12, prompt §7-3 precision).
* V8 jitter: Monte Carlo (10⁶) ⟨η_t⟩ vs the closed form within 3σ; Gauss–Hermite table average 1e−6.
* V9 identical paths ⇒ η_t = 1, θ_off = 0 (1e−12).
* V10 thermal timing with an ideal phase lock: Δτ_th(1 − R_τ,st/R_τ,th); TEC → 0 (1e−6).
* V11 table interpolation: (V_−, V_+, η_t) at half-grid points vs direct evaluation, relative 1e−3.
"""

import numpy as np
import pytest

from cranberry.blocks.fib.mode_overlap import (gaussian_mode, gaussian_overlap, gaussian_source, jitter_averaged_eta,
                                               make_mode_overlap_table, net_dispersion, timing_after_phase_lock)
from cranberry.blocks.source import temporal_modes as tm
from cranberry.blocks.source.opa_supermodes import make_opa_supermodes, schmidt_number_jsa
from cranberry.blocks.source.pulsed_testing import BOREALIS

TOG = {"multimode": True, "gvm": True, "gvd": True}
OPA = dict(length_m=1e-3, alpha_per_w=1.0, peak_pump_w=0.25, pump_fwhm_s=1e-12, gvm_s_per_m=3e-10,
           beta2_signal_s2_per_m=1e-25, beta2_pump_s2_per_m=4e-25, n_points=128, time_window_s=8e-12, n_z=20,
           phi_opa0_rad=0.0, segment="xtal", psi_p_xtal_rad=0.0, toggles=TOG)


def opa(**kw):
    tg = dict(TOG, **kw.pop("toggles", {}))
    return make_opa_supermodes({**OPA, **kw, "toggles": tg})


def test_v1_gaussian_overlap_numeric_vs_closed_form():
    t0 = 1.0
    g = (np.arange(4096) - 2048) * (160 * t0 / 4096)
    u = gaussian_mode(g, t0)
    for dt in np.linspace(-3, 3, 7) * t0:
        for chi in (-5.0, -1.0, 0.0, 0.5, 5.0):
            k = tm.overlap(u, tm.transform_mode(u, g, dt_s=dt, phi2_s2=2 * chi * t0**2), g)
            assert abs(k - complex(gaussian_overlap(dt, 2 * chi * t0**2, t0))) < 1e-6
    k = tm.overlap(u, tm.transform_mode(u, g, dt_s=0.5, phi2_s2=2.0), g)  # χ = 1: arg κ = ½ arctan χ − χΔt²/(4T0²(1+χ²))
    assert np.angle(k) == pytest.approx(0.5 * np.arctan(1.0) - 0.25 / 8, abs=1e-9)
    us = tm.transform_mode(u, g, dt_s=0.5)  # Δt > 0: the squeezed pulse arrives later
    assert np.sum(g * np.abs(us) ** 2) * (g[1] - g[0]) == pytest.approx(0.5, abs=1e-9)


def test_v2_symplectic_n512():
    c, s = opa(n_points=512, time_window_s=16e-12, n_z=3).propagator
    assert np.max(np.abs(c @ c.conj().T - s @ s.conj().T - np.eye(512))) < 1e-8
    assert np.max(np.abs(c @ s.T - s @ c.T)) < 1e-8


def test_v3_low_gain_schmidt():
    b = opa(peak_pump_w=1e-8, n_z=200)
    c, s = b.propagator
    f = b.jsa()
    assert np.max(np.abs(s / b.g - f)) / np.max(np.abs(f)) < 1e-3
    assert b.schmidt_number == pytest.approx(schmidt_number_jsa(f), rel=1e-3)
    x = np.linspace(-12, 12, 401)
    w1, w2 = np.meshgrid(x, x, indexing="ij")
    for sp, sm in ((1.0, 1.0), (0.7, 2.0), (0.5, 3.0)):
        jsa = np.exp(-((w1 + w2) / np.sqrt(2)) ** 2 / (2 * sp**2) - ((w1 - w2) / np.sqrt(2)) ** 2 / (2 * sm**2))
        assert schmidt_number_jsa(jsa) == pytest.approx(0.5 * (sm / sp + sp / sm), rel=1e-3)


def test_v4_multimode_bhd():
    b = opa()
    vm, _, th = b.source_output.project(b.leading_mode)
    assert float(vm) == pytest.approx(0.5 * np.exp(-2 * b.supermode_r[0]), rel=1e-9)
    assert float(th) == pytest.approx(0.0, abs=1e-9)
    so = gaussian_source(1.0, 1.0, 512, 40.0, 0.0)
    h1 = so.grid * so.modes[0]  # first Hermite–Gauss: orthogonal to the only supermode
    vm, vp, _ = so.project(h1)
    assert float(vm) == pytest.approx(0.5, rel=1e-9) and float(vp) == pytest.approx(0.5, rel=1e-9)
    vm, _, _ = so.project(so.modes[0])
    assert float(vm) == pytest.approx(0.5 * np.exp(-2.0), rel=1e-9)


def test_v6_quasi_cw():
    tg, ratios = 1e-10 * 0.01, []
    for m in (10, 30, 100):
        tau = m * tg
        t0 = tau / (2 * np.sqrt(np.log(2)))
        b = opa(length_m=0.01, gvm_s_per_m=1e-10, pump_fwhm_s=tau, n_points=64, time_window_s=10 * t0, n_z=40,
                toggles={"gvd": False})
        ratios.append(b.supermode_r[0] / np.sqrt(1.0 * 0.25))  # √(αP_pk)
    assert ratios[0] < ratios[1] < ratios[2] <= 1.0 + 1e-12
    assert abs(1 - ratios[2]) < 1e-2


def _overlap_block(**kw):
    p = dict(source={"gaussian": dict(T0_s=0.6e-12, r=1.0, n_points=512, window_s=40 * 0.6e-12, phi_opa0_rad=0.0)},
             dt_grid_s=[0.0], dL_grid_m=[0.0], beta2_s2_per_m=-21.7e-27, beta3_s3_per_m=0.0, pad=2,
             toggles={"dispersion": True, "jitter_fast": False})
    tg = dict(p["toggles"], **kw.pop("toggles", {}))
    return make_mode_overlap_table({**p, **kw, "toggles": tg})


def test_v7_single_mode_mismatch_is_loss():
    b = _overlap_block()
    for dt, dl in ((0.3e-12, 0.0), (0.5e-12, 30.0), (-0.4e-12, -50.0)):
        vm, vp, th, eta = b.evaluate(dt, dl)
        k = gaussian_overlap(dt, -21.7e-27 * dl, 0.6e-12)
        assert eta == pytest.approx(abs(k) ** 2, abs=1e-12)
        assert vm == pytest.approx(eta * 0.5 * np.exp(-2.0) + 0.5 * (1 - eta), abs=1e-12)
        assert vp == pytest.approx(eta * 0.5 * np.exp(2.0) + 0.5 * (1 - eta), abs=1e-12)
        assert th == pytest.approx(np.angle(k), abs=1e-12)  # GDD turns timing into θ_off = arg κ


def test_v8_jitter_average():
    t0, mu, sg = 1.0, 0.3, 0.4
    rng = np.random.default_rng(8)
    eta = np.abs(gaussian_overlap(rng.normal(mu, sg, 1_000_000), 0.0, t0)) ** 2
    se = eta.std() / np.sqrt(eta.size)
    assert abs(eta.mean() - jitter_averaged_eta(mu, sg, t0)) < 3 * se
    t0s = 0.6e-12
    b = _overlap_block(jitter_fast_s=0.4 * t0s, n_jitter=40, toggles={"jitter_fast": True, "dispersion": False})
    assert b.evaluate(0.3 * t0s, 0.0)[3] == pytest.approx(float(jitter_averaged_eta(0.3 * t0s, 0.4 * t0s, t0s)), abs=1e-6)


def test_v9_identical_paths():
    path = [{"L_m": 63.0, "beta2_s2_per_m": -21.7e-27, "beta3_s3_per_m": 0.127e-39},
            {"L_m": 2.0, "beta2_s2_per_m": 5e-27, "beta3_s3_per_m": 0.0}]
    p2, p3 = net_dispersion(path, list(path))
    assert p2 == 0.0 and p3 == 0.0
    vm, vp, th, eta = _overlap_block().evaluate(0.0, 0.0)
    assert abs(1 - eta) < 1e-12 and abs(th) < 1e-12


def test_v10_thermal_timing_cocompensation():
    r_th, r_st = 0.886e-15, 1.069e-15  # T03 §5.5 [AE]: fs/rad
    assert timing_after_phase_lock(2.5e-12, r_st, r_th, tec=False) == pytest.approx(2.5e-12 * (1 - r_st / r_th), rel=1e-6)
    assert timing_after_phase_lock(2.5e-12, r_st, r_th, tec=True) == 0.0


def test_v11_table_interpolation():
    from cranberry.blocks.source.pulsed_testing import TOGGLES_ON

    kc = np.pi * BOREALIS["fwhm_hz"]
    src = {"block": "pulsed_opo_source", "name": "opo", "params": dict(BOREALIS, x_ss=2.4, toggles=TOGGLES_ON)}
    dts = np.linspace(-0.3, 0.3, 61) / kc
    b = make_mode_overlap_table(dict(source=src, dt_grid_s=dts.tolist(), dL_grid_m=[0.0, 1000.0],
                                     beta2_s2_per_m=-21.7e-27, beta3_s3_per_m=0.0, pad=2,
                                     toggles={"dispersion": True, "jitter_fast": False}))
    tab = b.table
    assert tab.axes == ("dt_s", "dL_m") and tab.value_shape == (4,)
    for x in 0.5 * (dts[:-1] + dts[1:])[::6]:
        direct = np.array(b.evaluate(x, 500.0))
        got = tab.interp(dt_s=x, dL_m=500.0)
        for k in (0, 1, 3):
            assert got[k] == pytest.approx(direct[k], rel=1e-3)
