"""offline.source_output: SourceOutput (G07), projection, §7-1 single-mode reference (WP-F3)."""

import numpy as np
import pytest

from cranberry.core.units import noise_db
from cranberry.offline import gaussian as g
from cranberry.offline.source_output import SourceOutput, wrap_half_pi

REF_A = dict(loss=0.08, alpha_per_w=9.06, pump_w=0.64)  # prompt §7-1: L = 8 %, α = 906 %/W, P = 640 mW


def test_prompt_7_1_numbers_through_jitter_mapping():
    """§7-1 (T04 V1): L = 8 %, α = 906 %/W, P = 0.64 W, θ = 9 mrad → R′− = −10.148 dB, R′+ = +20.556 dB
    (±1e−3 dB); small-angle vs exact < 1e−3 dB. Prompt: ≈ −10.1 / +20.5 dB."""
    so = SourceOutput.single_mode_reference(**REF_A, phi_opa0_rad=0.0)
    vm, vp, th = so.project(None)
    assert vm.shape == () and float(th) == pytest.approx(0.0, abs=1e-12)
    # R± = L + (1 − L) e^{±2√(αP)}
    r = np.sqrt(9.06 * 0.64)
    assert float(vm) / 0.5 == pytest.approx(0.08 + 0.92 * np.exp(-2 * r), rel=1e-12)
    assert float(vp) / 0.5 == pytest.approx(0.08 + 0.92 * np.exp(2 * r), rel=1e-12)
    chi = g.chi_gaussian(9e-3)
    sq_db = float(noise_db(g.jitter_average(vm, vp, chi, delta0_rad=0.0)))
    asq_db = float(noise_db(g.jitter_average(vm, vp, chi, delta0_rad=np.pi / 2)))
    assert sq_db == pytest.approx(-10.148, abs=1e-3)
    assert asq_db == pytest.approx(20.556, abs=1e-3)
    assert sq_db == pytest.approx(-10.15, abs=0.10) and asq_db == pytest.approx(20.56, abs=0.05)  # D48
    sm = float(noise_db(g.jitter_average_small_angle(vm, vp, 9e-3)))
    sm_anti = float(noise_db(g.jitter_average_small_angle(vp, vm, 9e-3)))
    assert abs(sm - sq_db) < 1e-3 and abs(sm_anti - asq_db) < 1e-3
    assert so.provenance["model"] == "single_mode_reference"


def test_ref_a_fig5_optimum_point():
    """T04 V2 / §7-2 sanity: (8 %, 9 mrad) optimum pump ≈ 612 mW, −10.152 dB."""
    p = np.linspace(0.55, 0.68, 1301)
    vals = []
    for pw in p:
        vm, vp, _ = SourceOutput.single_mode_reference(loss=0.08, alpha_per_w=9.06, pump_w=pw, phi_opa0_rad=0.0).project(None)
        vals.append(float(noise_db(g.jitter_average(vm, vp, g.chi_gaussian(9e-3), delta0_rad=0.0))))
    k = int(np.argmin(vals))
    assert p[k] == pytest.approx(0.612, abs=0.002)
    assert vals[k] == pytest.approx(-10.152, abs=2e-3)


def test_reference_cov_frame_and_theta_off():
    """cov is in the lab frame at φ_p = 0: minimum at θ_sq,0 = φ_OPA,0/2 + π/2 (conv §4); θ_off = 0."""
    for phi0 in (0.0, 0.8, -2.5):
        so = SourceOutput.single_mode_reference(**REF_A, phi_opa0_rad=phi0)
        _, _, th_min = g.noise_levels_from_cov(so.cov)
        assert float(wrap_half_pi(th_min - (phi0 / 2 + np.pi / 2))) == pytest.approx(0.0, abs=1e-12)
        assert so.theta_sq0_rad == pytest.approx(phi0 / 2 + np.pi / 2)
        assert float(so.project(None)[2]) == pytest.approx(0.0, abs=1e-12)
    assert not so.cov.flags.writeable


def test_wrap_half_pi():
    a = np.array([0.0, np.pi / 2, -np.pi / 2, np.pi, 3.0, -3.0])
    w = wrap_half_pi(a)
    assert np.all((w > -np.pi / 2 - 1e-15) & (w <= np.pi / 2 + 1e-15))
    np.testing.assert_allclose(np.cos(2 * (w - a)), 1.0, atol=1e-12)  # equal mod π
    assert wrap_half_pi(-np.pi / 2) == pytest.approx(np.pi / 2)


def test_cw_sideband_projection_arrays():
    f = np.array([1e5, 1e6, 1e7])
    vm = np.array([0.05, 0.1, 0.3])
    vp = np.array([20.0, 10.0, 1.0])
    rot = np.array([0.0, -0.05, -0.2])  # e.g. OPO detuning rotation vs f_s
    phi0 = 0.4
    cov = g.cov_from_noise_levels(vm, vp, phi0 / 2 + np.pi / 2 + rot)
    so = SourceOutput(basis="cw_sideband", grid=f, cov=cov, phi_opa0_rad=phi0)
    a, b, th = so.project(None)
    assert a.shape == b.shape == th.shape == (3,)
    np.testing.assert_allclose(a, vm, rtol=1e-12)
    np.testing.assert_allclose(b, vp, rtol=1e-12)
    np.testing.assert_allclose(th, rot, atol=1e-12)
    with pytest.raises(ValueError):
        so.project(np.ones(3))
    with pytest.raises(ValueError):
        SourceOutput(basis="cw_sideband", grid=f, cov=cov[:2], phi_opa0_rad=0.0)
    with pytest.raises(ValueError):
        SourceOutput(basis="cw_sideband", grid=f[::-1], cov=cov, phi_opa0_rad=0.0)


def _two_mode_temporal():
    t = np.linspace(-10, 10, 2001)
    dt = t[1] - t[0]
    f0 = np.exp(-(t**2) / 2)
    f1 = t * np.exp(-(t**2) / 2)
    f0 = f0 / np.sqrt(np.sum(f0**2) * dt)
    f1 = f1 / np.sqrt(np.sum(f1**2) * dt)
    modes = np.stack([f0, f1])
    cov = g.vacuum(2)
    cov = g.squeeze(cov, 1.0, np.pi, mode=0)  # θ_sq = π/2 in mode 0 (φ_OPA,0 = 0 nominal)
    return t, modes, cov


def test_temporal_projection():
    t, modes, cov = _two_mode_temporal()
    so = SourceOutput(basis="temporal", grid=t, cov=cov, phi_opa0_rad=0.0, modes=modes)
    assert so.n_modes == 2
    vm, vp, th = so.project(modes[0])
    assert float(vm) == pytest.approx(0.5 * np.exp(-2), rel=1e-9)
    assert float(vp) == pytest.approx(0.5 * np.exp(2), rel=1e-9)
    assert float(th) == pytest.approx(0.0, abs=1e-9)
    # vacuum mode
    vm1, vp1, _ = so.project(modes[1])
    assert float(vm1) == pytest.approx(0.5, rel=1e-9) and float(vp1) == pytest.approx(0.5, rel=1e-9)
    # LO overlapping half with the squeezed mode → pure loss 1/2 (T12 §1.3); unnormalised input is fine
    vmh, vph, _ = so.project(3.0 * (modes[0] + modes[1]))
    assert float(vmh) == pytest.approx(0.5 * 0.5 * np.exp(-2) + 0.25, rel=1e-9)
    assert float(vph) == pytest.approx(0.5 * 0.5 * np.exp(2) + 0.25, rel=1e-9)
    # LO with a phase e^{iφ}: c = e^{−iφ} → the ellipse appears rotated by −φ
    _, _, thp = so.project(np.exp(0.3j) * modes[0])
    assert float(thp) == pytest.approx(-0.3, abs=1e-9)
    # LO outside span{f_n} sees vacuum
    u_out = t**2 * np.exp(-(t**2) / 2)
    u_out = u_out - (np.sum(u_out * modes[0]) * (t[1] - t[0])) * modes[0]
    vmo, vpo, _ = so.project(u_out)
    assert float(vmo) == pytest.approx(0.5, rel=1e-9) and float(vpo) == pytest.approx(0.5, rel=1e-9)
    with pytest.raises(ValueError):
        so.project(None)
    with pytest.raises(ValueError):
        so.project(np.zeros_like(t))


def test_temporal_validation():
    t, modes, cov = _two_mode_temporal()
    with pytest.raises(ValueError):
        SourceOutput(basis="temporal", grid=t, cov=cov, phi_opa0_rad=0.0)  # no modes
    with pytest.raises(ValueError):
        SourceOutput(basis="temporal", grid=t, cov=cov, phi_opa0_rad=0.0, modes=np.stack([modes[0], modes[0]]))
    with pytest.raises(ValueError):
        SourceOutput(basis="temporal", grid=t, cov=g.vacuum(3), phi_opa0_rad=0.0, modes=modes)


def test_single_mode_validation():
    with pytest.raises(ValueError):
        SourceOutput(basis="nope", grid=None, cov=g.vacuum(1), phi_opa0_rad=0.0)
    with pytest.raises(ValueError):
        SourceOutput(basis="single_mode", grid=np.arange(2.0), cov=g.vacuum(1), phi_opa0_rad=0.0)
    with pytest.raises(ValueError):
        SourceOutput(basis="single_mode", grid=None, cov=np.array([[0.5, 0.1], [0.0, 0.5]]), phi_opa0_rad=0.0)
    so = SourceOutput(basis="single_mode", grid=None, cov=g.vacuum(1), phi_opa0_rad=0.0)
    with pytest.raises(ValueError):
        so.project(np.ones(3))
    with pytest.raises(ValueError):
        SourceOutput.single_mode_reference(loss=1.5, alpha_per_w=9.06, pump_w=0.64, phi_opa0_rad=0.0)


def test_zero_pump_is_vacuum():
    """Prompt §7-4 for the reference: P = 0 → vacuum regardless of loss."""
    for loss in (0.0, 0.08, 1.0):
        vm, vp, _ = SourceOutput.single_mode_reference(loss=loss, alpha_per_w=9.06, pump_w=0.0, phi_opa0_rad=0.0).project(None)
        assert float(vm) == pytest.approx(0.5, rel=1e-12) and float(vp) == pytest.approx(0.5, rel=1e-12)
