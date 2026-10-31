"""VT-36 fiber quantum extras (PLAN §6.2; T07 §8 V2–V5, V8–V11; T10 §7 test 6; T17 §6).

* V2 symplectic (no loss/noise; Kerr + β2, soliton order 1, L = L_D): ‖MΩMᵀ − Ω‖_max < 1e−9
  (400 RK4IP steps; the RK4IP map is symplectic only to O(h⁴)), also with self-steepening in the
  photon-number quadratures; tangent columns vs central finite
  differences of the B9 map < 1e−6; step halving reduces the error of M by ~16 (12–22).
* V3 physicality with loss and every noise source: V + iΩ/2 ≥ −1e−12.
* V4 linearized Kerr (no dispersion): V_min = ½(√(1+φ²) − φ)² → −4.18/−7.66/−12.54 dB at φ = 0.5/1/2,
  anti-squeezing the sign flip, 0.01 dB (CW sideband path exact; CW on the grid with RK4IP).
* V5 pure loss: −10 dB input, V_out = ηV_in + (1−η)/2 to 1e−12 (also with dispersion: η × lossless + vacuum).
* V8 SBS (g_B = 5e−11 m/W, A_eff = 80 µm², 0.2 dB/km): P_th = 0.534 W / 54.1 mW / 34.4 mW at
  63 m / 630 m / 1 km (1 %); G_B > 5 → CHK warning (15.4 mW CW LO at 630 m), none at 63 m.
* V9 TW vs linear (optional block): φ_NL = 0.5, N_traj = 1e4, V at the linear θ_min within 3σ
  (σ = V√(2/N)) for 1e6 photons per bin; at 3 photons per bin TW departs by > 5σ (breakdown recorded).
* V10 Brillouin backscatter: b_B = 4.3e−9 /m and Φ_B(1 mW, 63 m) = 2.1e9 /s (1 %); N_B ∝ P·D, L, τ_eff.
* V11 phase-insensitive added noise: the backscatter covariance projected on a LO mode is N·𝟙 (every θ
  gains N, ΔV/V_SN = 2N), N = ∫ n|ũ|² df = the stray-table photons per mode (1e−5).
* Stray light: CW 1 mW, b_R = 6.18e−8 /m, 63 m, same polarization: 302.8 photons per 10 ns (0.5 %, T10
  test 6); Borealis-type mode (τ_eff = 14.5 ns): 435 photons per mode (2 %, T17 §6).
* All toggles off ⇒ identity (bright pulse, squeezed input).
"""

import numpy as np
import pytest

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.gnlse import model_from_params, pulse_field, solve_pulsed
from cranberry.blocks.fib.quantum import TW_TOGGLES, fiber_quantum_linear, fiber_quantum_tw, project_cov
from cranberry.blocks.fib.quantum_testing import BRILLOUIN, LAMBDA0, PROBE, linear_params, toggles
from cranberry.core.block import BlockSpec
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import make_block
from cranberry.core.rng import stream
from cranberry.offline import gaussian

KERR = toggles(beta2=True, kerr=True, kerr_fluct=True)
HNU = sc.H_PLANCK * 299792458.0 / LAMBDA0


def _db(v):
    return 10 * np.log10(v / 0.5)


def test_v2_symplectic_tangent_and_step_halving():
    r = fiber_quantum_linear(linear_params(steps=dict(n_steps=400, n_save=4, control="fixed")), KERR)
    om = gaussian.symplectic_form(64)
    assert np.max(np.abs(r.transfer @ om @ r.transfer.T - om)) < 1e-9
    # self-steepening: symplectic in the photon-number quadratures â(Ω) ∝ Ã/√(1 + Ω/ω0)
    ss = fiber_quantum_linear(linear_params(steps=dict(n_steps=400, n_save=4, control="fixed")),
                              dict(KERR, self_steepening=True)).transfer
    assert np.max(np.abs(ss @ om @ ss.T - om)) < 1e-9
    # finite differences of the classical B9 map (centred grid, δA = e_n → q column, i e_n → p column)
    p = linear_params(steps=dict(n_steps=400, n_save=4, control="fixed"))
    model = model_from_params(p, KERR, 64, r.dt)
    base = pulse_field("sech", p["pulse"]["peak_power_w"], 1e-12, r.t_s).astype(complex)
    eps = 1e-4 * np.abs(base).max()
    for j in (20, 32, 64 + 32, 64 + 40):
        d = np.zeros(64, dtype=complex)
        d[j % 64] = 1.0 if j < 64 else 1j

        def out(s):
            return solve_pulsed(model, base + s * eps * d, r.t_s, p["L_m"], n_steps=400, n_save=1, control="fixed",
                                rtol=None, monitor_edges=False).a_out

        fd = (out(1) - out(-1)) / (2 * eps)
        col = np.concatenate([fd.real, fd.imag])
        assert np.linalg.norm(col - r.transfer[:, j]) < 1e-6 * np.linalg.norm(col)

    def m(n):
        return fiber_quantum_linear(linear_params(steps=dict(n_steps=n, n_save=1, control="fixed")), KERR).transfer

    ref = m(1600)
    e1, e2 = np.abs(m(50) - ref).max(), np.abs(m(100) - ref).max()
    assert 12 < e1 / e2 < 22


def test_v3_physical_with_loss_and_all_noise():
    tg = toggles(True, cw_shortcut=False)
    r = fiber_quantum_linear(linear_params(alpha_per_m=0.02, xpm_b=2.0, xpm_power_w=0.01), tg)
    assert gaussian.is_physical(r.cov, tol=1e-12)
    assert all(np.max(np.abs(v)) > 0 for v in r.noise_cov.values())


@pytest.mark.parametrize("phi, db", [(0.5, -4.18), (1.0, -7.66), (2.0, -12.54)])
def test_v4_linearized_kerr(phi, db):
    power = phi / (1.3e-3 * 100.0)
    exact = _db(0.5 * (np.sqrt(1 + phi**2) - phi) ** 2)
    assert exact == pytest.approx(db, abs=0.005)
    cw = fiber_quantum_linear(linear_params(laser_mode="cw", power_w=power, L_m=100.0, f_s_hz=[1e6, 1e9],
                                            steps=dict(n_steps=1, n_save=1, control="fixed")),
                              toggles(kerr=True, kerr_fluct=True, cw_shortcut=True))
    vm, vp, _ = gaussian.noise_levels_from_cov(cw.cov)
    np.testing.assert_allclose(_db(vm), exact, atol=0.01)
    np.testing.assert_allclose(_db(vp), -exact, atol=0.01)
    grid = fiber_quantum_linear(linear_params(laser_mode="cw", power_w=power, L_m=100.0,
                                              grid=dict(n_t=8, t_window_s=2e-12),
                                              steps=dict(n_steps=100, n_save=1, control="fixed")),
                                toggles(kerr=True, kerr_fluct=True))
    vm2, vp2, _ = gaussian.noise_levels_from_cov(project_cov(grid.cov, np.ones(8), grid.dt))
    assert _db(vm2) == pytest.approx(exact, abs=0.01) and _db(vp2) == pytest.approx(-exact, abs=0.01)


def test_v5_pure_loss_exact():
    r_db = 10 / (20 / np.log(10))  # −10 dB
    sq = dict(kind="squeezed", r=r_db, theta_sq_rad=np.pi / 2, mode=dict(shape="sech", t0_s=1e-12))
    dark = dict(shape="sech", peak_power_w=0.0, t0_s=1e-12)
    alpha, L = 0.02, linear_params()["L_m"]
    eta = np.exp(-alpha * L)
    v_in = fiber_quantum_linear(linear_params(input_state=sq, pulse=dark), toggles()).cov
    out = fiber_quantum_linear(linear_params(input_state=sq, pulse=dark, alpha_per_m=alpha),
                               toggles(loss=True, loss_noise=True)).cov
    np.testing.assert_allclose(out, gaussian.pure_loss(v_in, eta), rtol=0, atol=1e-12)
    disp = fiber_quantum_linear(linear_params(input_state=sq, pulse=dark), toggles(beta2=True)).cov
    both = fiber_quantum_linear(linear_params(input_state=sq, pulse=dark, alpha_per_m=alpha),
                                toggles(beta2=True, loss=True, loss_noise=True)).cov
    np.testing.assert_allclose(both, gaussian.pure_loss(disp, eta), rtol=0, atol=1e-12)
    cw = fiber_quantum_linear(linear_params(laser_mode="cw", power_w=0.0, f_s_hz=[1e6], input_state=sq,
                                            alpha_per_m=alpha, steps=dict(n_steps=5, n_save=1, control="fixed")),
                              toggles(loss=True, loss_noise=True, cw_shortcut=True)).cov[0]
    np.testing.assert_allclose(cw, gaussian.pure_loss(gaussian.squeezed_vacuum(r_db, np.pi / 2), eta), atol=1e-12)
    vm = gaussian.noise_levels_from_cov(cw)[0]
    assert _db(vm) == pytest.approx(_db(eta * 0.05 + (1 - eta) / 2), abs=1e-9)


def test_v8_sbs_threshold_and_warning():
    a = 0.2 * np.log(10) / 1e4
    for L, pth in ((63.0, 0.534), (630.0, 54.1e-3), (1000.0, 34.4e-3)):
        assert sc.sbs_threshold_w(g_b_m_per_w=5e-11, a_eff_m2=80e-12, alpha_per_m=a, length_m=L,
                                  laser_linewidth_hz=0.0, dnu_b_hz=30e6) == pytest.approx(pth, rel=0.01)
    p = dict(brillouin=dict(BRILLOUIN), alpha_per_m=a, L_m=630.0, power_w=15.4e-3, laser_linewidth_hz=0.0,
             gain_warn=5.0)
    with pytest.warns(ChkWarning, match="sbs_gain"):
        make_block(BlockSpec("sbs_check", "sbs", p)).check()
    quiet = make_block(BlockSpec("sbs_check", "sbs", dict(p, L_m=63.0)))
    quiet.check()
    assert quiet.chk_log.events == []


def test_v9_truncated_wigner_vs_linear():
    dt = 0.25e-12
    sig = {}
    for nbar in (1e6, 3.0):
        power = nbar * HNU / dt
        p = linear_params(laser_mode="cw", power_w=power, L_m=1.0, gamma_per_w_m=0.5 / power,
                          grid=dict(n_t=4, t_window_s=4 * dt), steps=dict(n_steps=50, n_save=1, control="fixed"))
        lin = project_cov(fiber_quantum_linear(p, toggles(kerr=True, kerr_fluct=True)).cov, np.eye(4)[0], dt)
        th = gaussian.noise_levels_from_cov(lin)[2]
        tw = fiber_quantum_tw(p, dict(dict.fromkeys(TW_TOGGLES, False), kerr=True), 10_000, stream(7, "tw"))
        v_lin = gaussian.homodyne_variance(lin, th)
        v_tw = gaussian.homodyne_variance(project_cov(tw, np.eye(4)[0], dt), th)
        sig[nbar] = abs(v_tw - v_lin) / (v_lin * np.sqrt(2 / 10_000))
    assert v_lin == pytest.approx(0.5 * (np.sqrt(1.25) - 0.5) ** 2, rel=1e-7)  # RK4IP, 50 steps
    assert sig[1e6] < 3 and sig[3.0] > 5


def test_v10_brillouin_backscatter():
    b, s, a = sc.brillouin_backscatter_coeff(**BRILLOUIN, temperature_k=295.0, lambda0_m=LAMBDA0)
    assert b == pytest.approx(4.3e-9, rel=0.01) and s + a == pytest.approx(1.0)
    p = linear_params(L_m=63.0)
    tg = dict(rayleigh=False, brillouin=True, raman=False)
    lines, _ = sc.backscatter_sources(p, tg, 1e-3)
    assert sum(ln[2] for ln in lines) == pytest.approx(2.1e9, rel=0.01)
    assert {round(ln[1] / 1e9, 6) for ln in lines} == {-10.8, 10.8}

    def nb(pw, L, tau):
        q = linear_params(L_m=L, probe=dict(PROBE), gbs_filter=dict(shape="flat", fwhm_hz=1e12), t_gate_s=1e-9,
                          reflections=[], wdm_isolation_lin=1.0, mode=dict(shape="gaussian", t_fwhm_s=tau))
        return sc.stray_photons(q, dict(stray_rayleigh=False, stray_brillouin=True, stray_raman=False,
                                        stray_reflections=False), pw)["brillouin_mode"]

    n0 = nb(1e-3, 63.0, 100e-15)
    assert nb(0.5e-3, 63.0, 100e-15) == pytest.approx(0.5 * n0, rel=1e-12)  # P̄ = P·D
    assert nb(1e-3, 126.0, 100e-15) == pytest.approx(2 * n0, rel=1e-12)
    assert nb(1e-3, 63.0, 200e-15) == pytest.approx(2 * n0, rel=0.01)  # ∝ τ_eff while ν_B ≪ mode bandwidth


def test_v11_phase_insensitive_and_stray_consistency():
    p = linear_params(pulse=dict(shape="sech", peak_power_w=0.0, t0_s=1e-12), grid=dict(n_t=128, t_window_s=20e-12),
                      L_m=63.0)
    r = fiber_quantum_linear(p, toggles(rayleigh_backscatter=True))
    u = np.exp(-2 * np.log(2) * (r.t_s / 1e-12) ** 2)  # intensity FWHM 1 ps
    n2 = project_cov(r.noise_cov["rayleigh"], u, r.dt)
    n = n2[0, 0]
    np.testing.assert_allclose(n2, n * np.eye(2), atol=1e-12 * n)
    out = gaussian.add_noise(gaussian.squeezed_vacuum(1.0, 0.3), n)
    th = np.linspace(0, np.pi, 7)
    np.testing.assert_allclose(gaussian.homodyne_variance(out, th) - gaussian.homodyne_variance(
        gaussian.squeezed_vacuum(1.0, 0.3), th), n, rtol=1e-12)
    q = dict(p, gbs_filter=dict(shape="flat", fwhm_hz=1e12), t_gate_s=1e-9, reflections=[], wdm_isolation_lin=1.0,
             mode=dict(shape="gaussian", t_fwhm_s=1e-12))
    st = sc.stray_photons(q, dict(stray_rayleigh=True, stray_brillouin=False, stray_raman=False,
                                  stray_reflections=False), 1e-3)
    assert n == pytest.approx(st["rayleigh_mode"], rel=1e-5)


def test_stray_light_references():
    q = linear_params(L_m=63.0, rayleigh=dict(b_r_per_m=6.18e-8), gbs_filter=dict(shape="flat", fwhm_hz=1e12),
                      t_gate_s=10e-9, reflections=[], wdm_isolation_lin=1.0,
                      mode=dict(shape="exponential", tau_eff_s=14.5e-9))
    tg = dict(stray_rayleigh=True, stray_brillouin=False, stray_raman=False, stray_reflections=False)
    st = sc.stray_photons(q, tg, 1e-3)
    assert st["clicks_per_gate"] == pytest.approx(302.8, rel=0.005)
    st2 = sc.stray_photons(dict(q, rayleigh=dict(b_r_per_m=6.2e-8)), tg, 1e-3)
    assert st2["photons_per_mode"] == pytest.approx(435, rel=0.02)


def test_all_toggles_off_identity_bright_pulse():
    sq = dict(kind="squeezed", r=0.7, theta_sq_rad=0.4, mode=dict(shape="sech", t0_s=1e-12))
    r = fiber_quantum_linear(linear_params(input_state=sq, alpha_per_m=0.05), toggles())
    r0 = fiber_quantum_linear(linear_params(input_state=sq, L_m=0.0), toggles(True, cw_shortcut=False))
    np.testing.assert_allclose(r.cov, r0.cov, atol=1e-13)
