"""VT-31 (WP-B6): CW OPO (00_summary §4 VT-31; T02 V-B1, V-B2, V-B7, V-B8, V-C1).

* V-B1: Langevin 2×2 matrix (``cov``) = closed form R± = 1 ± ρ4x/((1∓x)² + (Ω/κ)²) at Δ = 0, Ref B
  values, f_s = 0, 1, 5, 20 MHz, relative 1e−12.
* V-B2: lossless (ρ = 1) R+R− = 1 for every Ω (1e−12); x = 0 → vacuum; Ω → ∞ → vacuum; R± − 1 linear
  in the downstream efficiency (VT-05 / VT-06 OPO parts, run here).
* V-B7: detuning, Ref B x = 0.745, Δ/κ = 1e−3: δθ_sq = −(Δ/κ)/(1 + x²) (numerical ellipse angle,
  relative 1e−3); probe output phase at the lock-1 extremum −Δ/(2κ) (singular vector of the
  transmission, both extrema); residual (Δ/κ)(1 − x²)/(2(1 + x²)). Signs: δl > 0 → Δ < 0 → θ_sq
  increases. Δ/κ = 0.1, f_s = 1 kHz, no detection loss (ρ = 0.970): minimum −12.88 dB, rotation
  −3.69° (±0.01).
* V-B8: L_ic as a callable ≡ 0.0038 gives the same covariance as the constant (bit for bit).
* V-C1: Ref C, FWHM 1.75 GHz, η_tot = 0.79, θ̃ = 0.1 rad, x from 10.0 dB anti-squeezing at 30 MHz:
  anti-squeezing 8.7 / 5.0 / 2.6 dB at 0.25 / 0.69 / 1.21 GHz (±0.05 dB, T02 computation); Ref C
  Fig. 4 values 8.5 / 5.2 / 2.8 dB within ±0.5 dB.
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from cranberry.blocks.source.opo_cw import closed_form, langevin_spectrum, probe_phase_coef, residual_coef, squeezing_angle_coef
from cranberry.blocks.source.opo_detuning import make_opo_detuning_noise
from cranberry.blocks.source.opo_testing import REF_B, REF_C, REF_C_ETA_D, db, measured_db, opo
from cranberry.offline import gaussian as g
from cranberry.offline.source_output import wrap_half_pi

X = 0.745


def r_pm(block):
    vm, vp, _ = block.source_output.project(None)
    return 2 * vm, 2 * vp


def test_v_b1_langevin_equals_closed_form():
    fs = np.array([0.0, 1e6, 5e6, 20e6])
    b = opo(REF_B, f_s_hz=fs)
    x, kappa = b.x, b.kappa(b.x)
    rho = b.source_output.provenance["escape_efficiency"]
    rm, rp = closed_form(x, rho, 2 * np.pi * fs / kappa)
    cov = b.source_output.cov
    np.testing.assert_allclose(2 * cov[:, 0, 0], rp, rtol=1e-12)  # φ_OPA,0 = 0: q amplified, p squeezed
    np.testing.assert_allclose(2 * cov[:, 1, 1], rm, rtol=1e-12)
    np.testing.assert_allclose(cov[:, 0, 1], 0.0, atol=1e-15)


def test_v_b2_purity_and_limits():
    fs = np.geomspace(1e3, 1e9, 25)
    lossless = opo(REF_B, loss_ic=0.0, f_s_hz=fs)
    rm, rp = r_pm(lossless)
    np.testing.assert_allclose(rm * rp, 1.0, rtol=1e-12)
    assert np.all(rm < 1) and rp[-1] - 1 < 1e-3  # Ω ≫ κ → vacuum
    np.testing.assert_allclose(opo(REF_B, pump_w=0.0, f_s_hz=fs).source_output.cov, 0.5 * np.broadcast_to(np.eye(2), (25, 2, 2)),
                               atol=1e-15)
    b = opo(REF_B, f_s_hz=fs)
    cov = b.source_output.cov
    for eta in (1.0, 0.5, 0.1):  # R± − 1 linear in the downstream efficiency (pure loss, conv §2)
        vm, vp, _ = g.noise_levels_from_cov(g.pure_loss(cov, eta))
        np.testing.assert_allclose(2 * vp - 1, eta * (r_pm(b)[1] - 1), rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(2 * vm - 1, eta * (r_pm(b)[0] - 1), rtol=1e-12, atol=1e-15)


def _ellipse_angle(delta_over_kappa, f_s=0.0, **kw):
    b = opo(REF_B, loss_ic=0.0038, pump_mode="x", x=X, f_s_hz=[f_s], toggles={"detuning": True}, **kw)
    k = b.kappa(X)
    b = opo(REF_B, loss_ic=0.0038, pump_mode="x", x=X, f_s_hz=[f_s], toggles={"detuning": True},
            detuning_rad_s=delta_over_kappa * k, **kw)
    vm, _, th = b.source_output.project(None)
    return float(th[0]), float(vm[0]), k


def _probe_angle(delta_over_kappa):
    """Output direction of the transmitted probe at the lock-1 extrema (singular vectors of M(0)⁻¹)."""
    b = opo(REF_B, loss_ic=0.0038, pump_mode="x", x=X)
    k = b.kappa(X)
    d, e = delta_over_kappa * k, X * k
    t = np.linalg.inv(np.array([[k - e, -d], [d, k + e]]))
    u, _, _ = np.linalg.svd(t)
    amp = float(wrap_half_pi(np.arctan2(u[1, 0], u[0, 0])))           # maximal amplification
    deamp = float(wrap_half_pi(np.arctan2(u[1, 1], u[0, 1]) - np.pi / 2))  # maximal de-amplification
    return amp, deamp, k


def test_v_b7_detuning_rotations_and_signs():
    d = 1e-3
    th, _, k = _ellipse_angle(d)
    assert th == pytest.approx(-d / (1 + X**2), rel=1e-3)
    assert th == pytest.approx(squeezing_angle_coef(X, k) * d * k, rel=1e-3)
    amp, deamp, _ = _probe_angle(d)
    for phi in (amp, deamp):
        assert phi == pytest.approx(-d / 2, rel=1e-3)
        assert phi == pytest.approx(probe_phase_coef(k) * d * k, rel=1e-3)
    assert amp - th == pytest.approx(d * (1 - X**2) / (2 * (1 + X**2)), rel=1e-3)
    assert residual_coef(X, k) * d * k == pytest.approx(d * (1 - X**2) / (2 * (1 + X**2)), rel=1e-12)
    # signs: δl > 0 lowers the resonance (Δ < 0, gain −ω0/l < 0) and θ_sq increases
    noise = make_opo_detuning_noise({"out": "opo.dn.d", "segment": "opo", "nu0_hz": 3.49e14, "round_trip_m": 0.5,
                                     "toggles": {"length_noise": False, "peaks": False, "thermal": False}})
    assert noise.gain == pytest.approx(-2 * np.pi * 3.49e14 / 0.5, rel=1e-15)
    assert _ellipse_angle(-d)[0] > 0
    # Δ/κ = 0.1 at 1 kHz, no detection loss
    th, vm, _ = _ellipse_angle(0.1, f_s=1e3)
    assert db(2 * vm) == pytest.approx(-12.88, abs=0.01)
    assert np.degrees(th) == pytest.approx(-3.69, abs=0.01)


def test_v_b7_frequency_dependent_rotation():
    """T02 §2.4 [AE]: Δ/κ = 0.01 rotates by −0.369° / −0.360° / −0.252° at 1 kHz / 1 MHz / 5 MHz.

    Computed here −0.3685 / −0.3621 / −0.2560° (ρ = 0.970); the note's values sit between this and the
    lossless case (−0.3617 / −0.2512°), so the check is the frequency dependence within 0.005°.
    """
    got = [np.degrees(_ellipse_angle(0.01, f_s=f)[0]) for f in (1e3, 1e6, 5e6)]
    np.testing.assert_allclose(got, [-0.369, -0.360, -0.252], atol=5e-3)
    assert got[0] < got[1] < got[2] < 0  # |rotation| falls as Ω approaches κ


def test_v_b8_callable_loss_equals_constant():
    const = opo(REF_B, loss_ic=0.0038)
    call = opo(REF_B, loss_ic=lambda x: 0.0038)
    np.testing.assert_array_equal(const.source_output.cov, call.source_output.cov)
    rm, rp = measured_db(call, 0.99 * 0.988**2 * 0.998, np.radians(1.5))
    assert (float(rm[0]), float(rp[0])) == pytest.approx((-9.44, 14.95), abs=0.01)


def test_v_c1_ref_c_linewidth_fwhm():
    theta = 0.1

    def anti(x):
        return measured_db(opo(REF_C, x=x), REF_C_ETA_D, theta)[1]

    x = brentq(lambda x: anti(x)[0] - 10.0, 0.1, 0.95, xtol=1e-12)
    assert x == pytest.approx(0.56, abs=0.005)
    a = anti(x)[1:]
    np.testing.assert_allclose(a, [8.7, 5.0, 2.6], atol=0.05)
    np.testing.assert_allclose(a, [8.5, 5.2, 2.8], atol=0.5)  # Ref C Fig. 4 (read values)
    # HWHM reading would give 9.6 / 7.8 / 5.6 dB: discriminated
    hw = measured_db(opo(REF_C, x=x, fwhm_hz=2 * 1.75e9), REF_C_ETA_D, theta)[1][1:]
    assert np.all(hw - a > 0.8)
    assert opo(REF_C).kappa(0.0) == pytest.approx(np.pi * 1.75e9, rel=1e-15)


def test_langevin_detuned_spectrum_hermitian_physical():
    s = langevin_spectrum(0.8, 3e7, 1e6, 5e6, 2 * np.pi * np.array([0.0, 1e6, 1e7]))
    np.testing.assert_allclose(s, np.conj(np.swapaxes(s, -1, -2)), atol=1e-12)
    b = opo(REF_B, f_s_hz=[0.0, 1e6, 1e7], detuning_rad_s=5e6, toggles={"detuning": True})
    assert np.all(g.is_physical(b.source_output.cov))
