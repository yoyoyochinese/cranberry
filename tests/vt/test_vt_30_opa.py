"""VT-30 (WP-B1): waveguide OPA.

Sources: S §4 VT-30; T01 V-T01-3..10 (§9), §2.2–§2.7, §7.3; T14 §3.4 (lock-2 reduction).
Tolerances (T01 §9): uncertainty R+R− = 1 lossless (1e−9), > 1 lossy; sign θ_sq = φ_p/2 + π/2 to
1e−4 rad; optimum pump numeric = closed form 1e−6 (P_opt = 612.3 mW); mismatch Δkℓ/2 = 0.01 at
r = 2.408 → δθ_mm = −2.04 ± 0.05 mrad (sign fixed here: negative for Δk > 0); coupled-wave energy
conservation 1e−9, weak-seed gain e^r 1e−4, SHG P = αP² and sinc² half width 1.392; +1 % pump →
R′+ +0.104 ± 0.005 dB, R′− +0.002 dB; lock 2: ε = arctan((A−/A+) tan δu′) = 8.1e−4 rad for
δu′ = 0.1 rad (T14); crystal phase enters only as ½ψ_p,xtal via the pump-path segment.
Not here: the PD-OPA 2Ω_m beat of V-T01-10 (error-signal block, WP-B3, from this gain table).
"""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from cranberry.blocks.source.opa_waveguide import make_opa_waveguide, propagate
from cranberry.blocks.source.single_mode import make_single_mode_reference, r_prime
from cranberry.core.units import lin_to_db
from cranberry.offline import gaussian as g

R_A = np.sqrt(9.06 * 0.64)  # 2.408
OFF = dict(pump_loss=False, phase_mismatch=False, pump_depletion=False, pump_rin=False)
BASE = dict(alpha_per_w=9.06, pump_w=0.64, pump_ref="output", length_m=0.045, phi_opa0_rad=0.0,
            segment="xtal", psi_p_xtal_rad=0.0, gain_table_pump_w=None, toggles=OFF)
LUMPED = dict(BASE, model="lumped", waveguide_loss=0.0, loss_source="test")
DIST = dict(BASE, model="distributed", loss_signal_db_per_cm=0.1)


def opa(params, **kw):
    toggles = dict(params["toggles"], **kw.pop("toggles", {}))
    return make_opa_waveguide(dict(params, **kw, toggles=toggles))


def lossless(**kw):
    defaults = dict(alpha_per_w=9.06, pump_in_w=0.64, length_m=0.045, a_s_per_m=0.0, a_p_per_m=0.0,
                    delta_k_per_m=0.0, phi_opa0_rad=0.0)
    return propagate(**(defaults | kw))


def test_uncertainty_v_t01_4():
    for p_w in (0.05, 0.64, 1.2):
        vm, vp, _ = g.noise_levels_from_cov(opa(LUMPED, pump_w=p_w).source_output.cov)
        assert 4 * vm * vp == pytest.approx(1.0, rel=1e-9)
        for a_s in (0.03, 0.1, 0.5):
            vm, vp, _ = g.noise_levels_from_cov(opa(DIST, pump_w=p_w, loss_signal_db_per_cm=a_s).source_output.cov)
            assert 4 * vm * vp > 1.0


@pytest.mark.parametrize("phi_p", [0.7, -1.1, 2.9])
def test_sign_v_t01_5(phi_p):
    """θ_sq = φ_p/2 + π/2: φ_p and φ_OPA,0 enter G together, so the crystal is run with φ_OPA,0 = φ_p."""
    so = opa(LUMPED, phi_opa0_rad=phi_p).source_output
    _, _, th = g.noise_levels_from_cov(so.cov)
    assert np.mod(th - (phi_p / 2 + np.pi / 2) + np.pi / 2, np.pi) - np.pi / 2 == pytest.approx(0.0, abs=1e-4)
    # same state as conv §4 H_int: ξ = r e^{i(φ_p + π)}  ⇒ squeezed vacuum at φ_p/2 + π/2
    np.testing.assert_allclose(so.cov, g.squeeze(g.vacuum(1), R_A, phi_p + np.pi), rtol=1e-9, atol=1e-12)
    # seed on θ = φ_p/2 is amplified by e^r and keeps its direction
    t = lossless(phi_opa0_rad=phi_p)["T"]
    e = np.array([np.cos(phi_p / 2), np.sin(phi_p / 2)])
    np.testing.assert_allclose(t @ e, np.exp(R_A) * e, rtol=1e-9)
    # LO on θ_sq reads R−
    assert 2 * g.homodyne_variance(so.cov, phi_p / 2 + np.pi / 2) == pytest.approx(np.exp(-2 * R_A), rel=1e-9)


def test_phi_opa0_shift_and_literature_offset():
    """φ_OPA,0 moves θ_sq by φ_OPA,0/2; the raw coupled-wave form (arg A_p = 0) is φ_OPA,0 = π/2 (T01 §2.2)."""
    for phi0 in (0.0, 0.4):
        assert opa(LUMPED, phi_opa0_rad=phi0).source_output.project(None)[2] == pytest.approx(0.0, abs=1e-9)
    _, _, th = g.noise_levels_from_cov(lossless(phi_opa0_rad=np.pi / 2)["cov"])
    assert th == pytest.approx(np.pi / 4 + np.pi / 2, abs=1e-9)  # θ_amp = (arg A_p + π/2)/2 = π/4


def test_optimum_pump_v_t01_6():
    theta, loss = 0.009, 0.08

    def rm(p_w):
        b = make_single_mode_reference(dict(loss=loss, alpha_per_w=9.06, pump_w=p_w, theta_rms_rad=theta,
                                            phi_opa0_rad=0.0, average="small_angle"))
        return b.r_prime()[0]

    res = minimize_scalar(rm, bounds=(0.3, 1.0), method="bounded", options={"xatol": 1e-10})
    r_opt = 0.5 * np.log(1 / np.tan(theta))
    assert res.x == pytest.approx(r_opt**2 / 9.06, rel=1e-6)
    assert res.x * 1e3 == pytest.approx(612.3, abs=0.05)
    assert res.fun == pytest.approx(loss + (1 - loss) * np.sin(2 * theta), rel=1e-6)


def r_eff_formula(gg, dd, ell):
    """T01 §2.5: cosh(2r_eff) = cosh²(sℓ) + sinh²(sℓ)(g² + δ²)/(g² − δ²), s = √(g² − δ²)."""
    s = np.sqrt(complex(gg**2 - dd**2))
    c = np.cosh(s * ell) ** 2 + np.sinh(s * ell) ** 2 * (gg**2 + dd**2) / (gg**2 - dd**2)
    return 0.5 * np.arccosh(c.real)


def test_phase_mismatch_v_t01_7():
    ell = 0.045
    gg = R_A / ell
    # Δkℓ/2 = 0.01 → δθ_mm = −2.04 mrad (sign fixed: negative for Δk > 0)
    so = opa(LUMPED, delta_k_per_m=0.02 / ell, toggles=dict(phase_mismatch=True)).source_output
    assert so.project(None)[2] * 1e3 == pytest.approx(-2.04, abs=0.05)
    assert so.provenance["delta_theta_mm_rad"] == pytest.approx(float(so.project(None)[2]))
    for half in (0.01, 0.1, 0.3, 1.392, 1.2 * R_A, 3 * R_A):
        dd = half / ell
        vm, vp, _ = g.noise_levels_from_cov(lossless(delta_k_per_m=2 * dd)["cov"])
        assert 0.25 * np.log(vp / vm) == pytest.approx(r_eff_formula(gg, dd, ell), rel=1e-7)
        if dd > gg:  # no amplification: cosh(2r_eff) ≤ C = (g² + δ²)/(δ² − g²) for every ℓ, so R+ ≤ C + √(C² − 1)
            c = (gg**2 + dd**2) / (dd**2 - gg**2)
            assert 2 * vp <= c + np.sqrt(c * c - 1) + 1e-9
    # δ = 0 is §2.1
    np.testing.assert_array_equal(opa(LUMPED, delta_k_per_m=0.0, toggles=dict(phase_mismatch=True)).source_output.cov,
                                  opa(LUMPED).source_output.cov)


def shg(p_fh, dk_half, ell=0.045):
    out = propagate(alpha_per_w=9.06, pump_in_w=0.0, length_m=ell, a_s_per_m=0.0, a_p_per_m=0.0,
                    delta_k_per_m=2 * dk_half / ell, phi_opa0_rad=0.0, seed=np.sqrt(p_fh) + 0j, depletion=True)
    return abs(out["pump_out"]) ** 2


def test_coupled_wave_v_t01_8():
    # (a) weak seed: gain e^r, r = √(αP)
    s0 = np.sqrt(1e-12)
    out = lossless(seed=s0 * np.exp(0j), depletion=True)
    assert abs(out["seed_out"]) / s0 == pytest.approx(np.exp(R_A), rel=1e-4)
    # (b) energy conservation with strong seed, depletion and mismatch (lossless)
    for seed_w, dk in ((0.05, 0.0), (0.3, 30.0)):
        out = lossless(seed=np.sqrt(seed_w) * np.exp(0.3j), depletion=True, delta_k_per_m=dk)
        assert abs(out["seed_out"]) ** 2 + abs(out["pump_out"]) ** 2 == pytest.approx(0.64 + seed_w, rel=1e-9)
    # (c) SHG at low conversion: P_SH = α P_FH², sinc² tuning with half maximum at Δkℓ/2 = 1.392
    p_fh = 1e-7  # conversion αP ≈ 1e−6
    assert shg(p_fh, 0.0) == pytest.approx(9.06 * p_fh**2, rel=1e-4)
    assert shg(p_fh, 1.392) / shg(p_fh, 0.0) == pytest.approx(0.5, abs=1e-3)


def test_pump_sensitivity_v_t01_9():
    def db(p_w):
        return make_single_mode_reference(dict(loss=0.08, alpha_per_w=9.06, pump_w=p_w, theta_rms_rad=0.009,
                                               phi_opa0_rad=0.0, average="exact")).r_prime_db()

    m0, p0 = db(0.64)
    m1, p1 = db(0.64 * 1.01)
    assert p1 - p0 == pytest.approx(0.104, abs=0.005)
    assert p1 - p0 == pytest.approx(16.3 * 0.0064, abs=0.005)  # analytic slope 8.686 r/(2P) = 16.3 dB/W
    assert m1 - m0 == pytest.approx(0.002, abs=0.001)
    # the OPA's RIN hook: δr = (r/2)·δP/P
    so0, so1 = opa(LUMPED).source_output, opa(LUMPED, pump_w=0.64 * 1.01).source_output
    assert so1.provenance["r_eff"] - so0.provenance["r_eff"] == pytest.approx(so0.provenance["dr_d_rin"] * 0.01, rel=3e-3)


def test_locks_v_t01_10_partial():
    """Lock 2 (T14 §3.4): tan ε = (A−/A+) tan δu′ with A± from the gain table (A−/A+ = e^{−2r} = 1/123)."""
    t = opa(LUMPED, gain_table_pump_w=[0.1, 0.64, 1.0]).source_output.gain_table
    a_p, a_m = t.interp(pump_w=0.64)
    assert np.arctan(a_m / a_p * np.tan(0.1)) == pytest.approx(8.1e-4, abs=0.05e-4)
    # crystal temperature: ψ_p,xtal is a pump-path segment phase (½ψ into θ_sq via the time loop), not in cov
    a, b = opa(LUMPED).source_output, opa(LUMPED, psi_p_xtal_rad=0.5).source_output
    np.testing.assert_array_equal(a.cov, b.cov)
    assert b.pump_path_phase_rad["xtal"] - a.pump_path_phase_rad["xtal"] == 0.5


def test_distributed_vs_lumped_t01_7_3():
    """Full budget with downstream lumped losses (1 %, 1 %, 0.2 %, 4 %) and θ = 9 mrad:
    distributed waveguide at α = 9.06 /W gives R′− ≈ −10.17, R′+ ≈ 20.29 dB (0.27 dB below lumped, D48)."""
    eta = 0.99 * 0.99 * 0.998 * 0.96

    def budget(cov):
        vm, vp, _ = g.noise_levels_from_cov(g.pure_loss(cov, eta))
        return [float(lin_to_db(x)) for x in r_prime(2 * vm, 2 * vp, 0.009, "exact")]

    dist = budget(opa(DIST).source_output.cov)
    lump = budget(g.pure_loss(opa(LUMPED).source_output.cov, 0.98))
    assert dist[0] == pytest.approx(-10.17, abs=0.01) and dist[1] == pytest.approx(20.29, abs=0.01)
    assert lump[1] - dist[1] == pytest.approx(0.27, abs=0.01)
