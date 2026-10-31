"""VT-35 (PLAN §6.2 "GNLSE extras", T06 §9 V6–V9, V11).

* V6 dispersion: linear only, Gaussian, L = 3 L_D. Output/input spectral phase = β2Ω²L/2 + β3Ω³L/6
  (≤ 1e−9 rad where |Ã| > 1e−3 max; float64 rounding of the phase ≈ 1e−12·|phase|); β2 < 0 ⇒ blue
  (Ω > 0) half arrives first, at T ≈ β2 L ⟨Ω⟩_blue (5 %); rms width × √(1 + (L/L_D)²) (1e−6, β3 = 0).
* V7 loss: E(L)/E(0) = e^{−αL} (1e−12) with α from dB/km (0.2 dB/km → 4.605e−5 /m).
* V8 pure SPM (β = 0, f_R = 0, shock off, with loss): |A|e^{αL/2} = |A0|, phase γ|A0|²L_eff (rel 1e−8,
  abs 1e−11 rad for FFT rounding in the wings).
* V9 SSFS (Blow–Wood, f_R = 0.18, shock off, β3 = 0, N = 1, T0 = 0.5 ps, L = L_D): centroid shift
  vs Gordon −8T_R|β2|L/(15T0⁴), T_R = f_R ∫t h_R dt = 1.46 fs, within 3 % (prototype 1.004).
  Not applied to Lin–Agrawal (T06 §2.3: h_b has zero first moment).
* V11 grid edge: energy in the outer 5 % of the window > 1e−6 → CHK ``grid_edge`` warning, none otherwise.
"""

import warnings

import numpy as np
import pytest

from cranberry.blocks.fib.gnlse import (
    GnlseModel,
    db_per_km_to_alpha_per_m,
    l_eff,
    pulse_field,
    raman_first_moment_s,
    solve_pulsed,
    time_axis,
)
from cranberry.core.checks import ChkLog, ChkWarning

B2, B3, GAMMA = -21.7e-27, 0.126e-39, 1.3e-3


def _model(n, dt, **kw):
    d = dict(lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1.0, c_raman=0.0,
             raman_model=None, self_steepening=False)
    d.update(kw)
    return GnlseModel(n, dt, **d)


def test_v6_dispersion_phase_sign_and_broadening():
    t0 = 100e-15
    ld = t0**2 / abs(B2)
    t, dt = time_axis(4096, 200 * t0)
    a0 = pulse_field("gaussian", 1.0, t0, t)
    m = _model(4096, dt, betas=(B2, B3), gamma_per_w_m=0.0)
    r = solve_pulsed(m, a0, t, 3 * ld, n_steps=10, n_save=1, control="fixed", rtol=None)
    s0, s1 = r.spectrum(0), r.spectrum(-1)
    w = m.omega
    ok = np.abs(s0) > 1e-3 * np.abs(s0).max()
    dphi = np.angle(s1[ok] / s0[ok]) - np.angle(np.exp(1j * (B2 * w[ok] ** 2 * 3 * ld / 2 + B3 * w[ok] ** 3 * 3 * ld / 6)))
    assert np.max(np.abs(np.angle(np.exp(1j * dphi)))) <= 1e-9

    # sign: blue half first (β2 < 0), at T ≈ β2 L ⟨Ω⟩_blue
    m2 = _model(4096, dt, gamma_per_w_m=0.0)
    r2 = solve_pulsed(m2, a0, t, 3 * ld, n_steps=10, n_save=1, control="fixed", rtol=None)
    s = r2.spectrum(-1)
    blue = np.where(w > 0, s, 0.0)
    from cranberry.core.fft import ifft_agrawal

    ab = ifft_agrawal(blue, dt, t0=float(t[0]))
    t_blue = np.sum(t * np.abs(ab) ** 2) / np.sum(np.abs(ab) ** 2)
    w_blue = np.sum(w * np.abs(blue) ** 2) / np.sum(np.abs(blue) ** 2)
    assert t_blue < 0
    assert t_blue == pytest.approx(B2 * 3 * ld * w_blue, rel=0.05)

    rms = lambda a: np.sqrt(np.sum(t**2 * np.abs(a) ** 2) / np.sum(np.abs(a) ** 2))  # noqa: E731
    assert rms(r2.a_out) / rms(a0) == pytest.approx(np.sqrt(1 + 9), rel=1e-6)


def test_v7_loss():
    alpha = db_per_km_to_alpha_per_m(0.2)
    assert alpha == pytest.approx(4.60517e-5, rel=1e-5)
    t, dt = time_axis(512, 20e-12)
    m = _model(512, dt, alpha_per_m=alpha, betas=(), gamma_per_w_m=0.0)
    r = solve_pulsed(m, pulse_field("sech", 1.0, 1e-12, t), t, 5000.0, n_steps=7, n_save=1, control="fixed", rtol=None)
    assert r.energy_j[1] / r.energy_j[0] == pytest.approx(np.exp(-alpha * 5000.0), rel=1e-12)


def test_v8_pure_spm_with_loss():
    alpha, length, p0 = 1e-3, 500.0, 2.0
    t, dt = time_axis(1024, 40e-12)
    a0 = pulse_field("gaussian", p0, 2e-12, t)
    m = _model(1024, dt, alpha_per_m=alpha, betas=())
    r = solve_pulsed(m, a0, t, length, n_steps=1000, n_save=1, control="fixed", rtol=None)
    big = np.abs(a0) > 1e-3 * np.sqrt(p0)
    np.testing.assert_allclose(np.abs(r.a_out[big]) * np.exp(alpha * length / 2), np.abs(a0[big]), rtol=1e-8)
    phase = np.angle(r.a_out[big] / a0[big])
    np.testing.assert_allclose(phase, GAMMA * np.abs(a0[big]) ** 2 * l_eff(alpha, length), rtol=1e-8, atol=1e-11)


def test_v9_ssfs_vs_gordon_blow_wood():
    t0, f_r = 0.5e-12, 0.18
    ld = t0**2 / abs(B2)
    p1 = abs(B2) / (GAMMA * t0**2)
    t, dt = time_axis(4096, 80 * t0)
    m = _model(4096, dt, c_kerr=1 - f_r, c_raman=f_r, raman_model="blow_wood")
    r = solve_pulsed(m, pulse_field("sech", p1, t0, t), t, ld, n_steps=200, n_save=1, control="fixed", rtol=None)
    shift = m.centroid(r.spectrum(-1)) - m.centroid(r.spectrum(0))
    t_r = f_r * raman_first_moment_s("blow_wood")
    assert t_r == pytest.approx(1.46e-15, rel=0.01)
    gordon = -8 * t_r * abs(B2) * ld / (15 * t0**4)
    assert shift / gordon == pytest.approx(1.0, abs=0.03), shift / gordon


def test_v11_grid_edge_monitor():
    t0 = 1e-12
    t, dt = time_axis(256, 40 * t0)
    m = _model(256, dt, gamma_per_w_m=0.0)
    log = ChkLog()
    with warnings.catch_warnings():
        warnings.simplefilter("error", ChkWarning)
        r = solve_pulsed(m, pulse_field("sech", 1.0, t0, t), t, 1.0, n_steps=2, n_save=1, control="fixed", rtol=None,
                         log=log)
    assert not log.events and max(r.edge_fractions) <= 1e-6
    t, dt = time_axis(256, 8 * t0)  # window too narrow for a sech of T0
    m = _model(256, dt, gamma_per_w_m=0.0)
    with pytest.warns(ChkWarning, match="grid_edge"):
        r = solve_pulsed(m, pulse_field("sech", 1.0, t0, t), t, 1.0, n_steps=2, n_save=1, control="fixed", rtol=None,
                         log=log)
    assert log.events and log.events[-1].code == "grid_edge" and r.edge_fractions[0] > 1e-6
