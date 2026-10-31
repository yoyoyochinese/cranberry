"""offline.gaussian: conv §2–§3 / §5–§6 operations, symplectic and physicality checks (WP-F3)."""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from cranberry.core.rng import stream
from cranberry.core.units import noise_db
from cranberry.offline import gaussian as g

TOL = 1e-12


def _random_gate_sequence(n_modes, rng, n_gates=12):
    s = np.eye(2 * n_modes)
    for _ in range(n_gates):
        kind = rng.integers(3 if n_modes > 1 else 2)
        if kind == 0:
            m = g.rotation(rng.uniform(-np.pi, np.pi), n_modes, int(rng.integers(n_modes)))
        elif kind == 1:
            m = g.squeezer(rng.uniform(0, 1.5), rng.uniform(-np.pi, np.pi), n_modes, int(rng.integers(n_modes)))
        else:
            i, j = rng.choice(n_modes, 2, replace=False)
            m = g.beam_splitter_matrix(rng.uniform(0, 1), n_modes, int(i), int(j))
        s = m @ s
    return s


def test_vacuum_and_symplectic_form():
    np.testing.assert_array_equal(g.vacuum(2), 0.5 * np.eye(4))
    om = g.symplectic_form(2)
    np.testing.assert_array_equal(om, [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    with pytest.raises(ValueError):
        g.vacuum(0)
    with pytest.raises(ValueError):
        g.n_modes_of(np.eye(3))


def test_squeezer_matches_conv3_variance():
    """Var(x_θ) = ½[e^{−2r} cos²(θ − φ_ξ/2) + e^{2r} sin²(θ − φ_ξ/2)], θ_sq = φ_ξ/2 (conv §3)."""
    theta = np.linspace(-np.pi, np.pi, 37)
    for r in (0.0, 0.3, 1.2, 2.408):
        for phi in (-2.0, 0.0, 0.7, np.pi):
            cov = g.squeeze(g.vacuum(1), r, phi)
            v = g.homodyne_variance(cov, theta)
            expect = 0.5 * (np.exp(-2 * r) * np.cos(theta - phi / 2) ** 2 + np.exp(2 * r) * np.sin(theta - phi / 2) ** 2)
            np.testing.assert_allclose(v, expect, rtol=1e-12, atol=1e-15)
            vm, vp, th = g.noise_levels_from_cov(cov)
            assert vm == pytest.approx(0.5 * np.exp(-2 * r), rel=1e-12)
            assert vp == pytest.approx(0.5 * np.exp(2 * r), rel=1e-12)
            if r > 0:
                assert np.isclose(np.mod(th - phi / 2 + np.pi / 2, np.pi) - np.pi / 2, 0.0, atol=1e-12)


def test_squeezer_heisenberg_form_phi0():
    """φ_ξ = 0: S†âS = â cosh r − â† sinh r ⇒ q → e^{−r} q, p → e^{r} p (conv §3)."""
    np.testing.assert_allclose(g.squeezer(0.4, 0.0), np.diag([np.exp(-0.4), np.exp(0.4)]), atol=1e-15)
    with pytest.raises(ValueError):
        g.squeezer(-0.1, 0.0)


def test_rotation_adds_to_squeezing_angle():
    """A → A e^{iψ} (conv §4) rotates the squeezing angle by +ψ (T04 §4)."""
    for theta_sq in (0.0, 0.4, 1.3):
        for psi in (-0.9, 0.25, 2.0):
            cov = g.rotate(g.squeezed_vacuum(1.0, theta_sq), psi)
            _, _, th = g.noise_levels_from_cov(cov)
            assert np.isclose(np.mod(th - (theta_sq + psi), np.pi), 0.0, atol=1e-12) or np.isclose(
                np.mod(th - (theta_sq + psi), np.pi), np.pi, atol=1e-12
            )
    # per-mode angles and single-mode embedding agree
    cov = g.vacuum(2)
    cov = g.squeeze(cov, 0.5, 0.0, mode=0)
    cov = g.squeeze(cov, 0.8, 0.0, mode=1)
    a = g.rotate(cov, np.array([0.3, -0.7]))
    b = g.rotate(g.rotate(cov, 0.3, mode=0), -0.7, mode=1)
    np.testing.assert_allclose(a, b, atol=1e-14)


def test_gates_are_symplectic_to_1e12():
    rng = stream(11, "test.gaussian.symplectic")
    for n in (1, 2, 3, 5):
        for _ in range(5):
            s = _random_gate_sequence(n, rng)
            assert g.is_symplectic(s, tol=TOL)
    assert not g.is_symplectic(np.diag([2.0, 2.0]), tol=TOL)
    assert g.is_symplectic(g.beam_splitter_matrix(0.5, 2, 0, 1), tol=TOL)


def test_symplectic_invariance_of_symplectic_eigenvalues_and_physicality():
    """Symplectic eigenvalues (and so purity / physicality) are invariant under gates (1e-12)."""
    rng = stream(12, "test.gaussian.invariance")
    for n in (1, 2, 4):
        cov = g.vacuum(n)
        cov = g.add_noise(g.pure_loss(cov, rng.uniform(0.2, 1.0, n)), rng.uniform(0, 0.3, n))
        nu0 = g.symplectic_eigenvalues(cov)
        s = _random_gate_sequence(n, rng)
        cov2 = g.apply_symplectic(cov, s)
        np.testing.assert_allclose(g.symplectic_eigenvalues(cov2), nu0, rtol=1e-12)
        assert g.is_physical(cov2)
        # pure states stay pure: all ν = 1/2
        pure = g.apply_symplectic(g.vacuum(n), s)
        np.testing.assert_allclose(g.symplectic_eigenvalues(pure), 0.5, rtol=1e-12)


def test_is_physical():
    assert g.is_physical(g.vacuum(3))
    assert g.is_physical(g.squeezed_vacuum(2.0, 0.3))
    assert not g.is_physical(np.diag([0.1, 0.1]))  # V_q V_p < 1/4
    assert not g.is_physical(np.diag([0.4, 0.6]))  # 0.24 < 1/4
    assert not g.is_physical(np.array([[0.5, 0.1], [0.0, 0.5]]))  # not symmetric
    # lossy squeezed state is physical; batch input
    batch = np.stack([g.pure_loss(g.squeezed_vacuum(r, 0.0), 0.7) for r in (0.1, 1.0, 3.0)])
    assert g.is_physical(batch)


def test_noise_levels_round_trip_and_batch():
    vm = np.array([0.02, 0.3, 0.5])
    vp = np.array([12.0, 0.9, 0.5])
    th = np.array([0.1, 2.9, 0.0])
    cov = g.cov_from_noise_levels(vm, vp, th)
    assert cov.shape == (3, 2, 2)
    a, b, c = g.noise_levels_from_cov(cov)
    np.testing.assert_allclose(a, vm, rtol=1e-12)
    np.testing.assert_allclose(b, vp, rtol=1e-12)
    np.testing.assert_allclose(c, np.mod(th, np.pi), atol=1e-12)
    # minimum variance is at θ_sq
    np.testing.assert_allclose(g.homodyne_variance(cov, th), vm, rtol=1e-12)
    np.testing.assert_allclose(g.homodyne_variance(cov, th + np.pi / 2), vp, rtol=1e-12)
    with pytest.raises(ValueError):
        g.cov_from_noise_levels(1.0, 0.5, 0.0)
    assert g.noise_levels_from_cov(g.vacuum(1)) == (0.5, 0.5, 0.0)


def test_homodyne_equals_quadrature_variance_multimode():
    rng = stream(13, "test.gaussian.homodyne")
    cov = g.apply_symplectic(g.vacuum(3), _random_gate_sequence(3, rng))
    for k in range(3):
        for th in (0.0, 0.4, 2.0):
            w = np.zeros(6)
            w[k], w[k + 3] = np.cos(th), np.sin(th)
            assert g.homodyne_variance(cov, th, mode=k) == pytest.approx(g.quadrature_variance(cov, w), rel=1e-12)


def test_beam_splitter_conv5_and_loss_equivalence():
    """η = 1/2: b1 = (a_s + a_LO)/√2, b2 = (a_s − a_LO)/√2 (conv §5); BS with vacuum = pure loss."""
    m = g.beam_splitter_matrix(0.5, 2, 0, 1)
    h = 1 / np.sqrt(2)
    np.testing.assert_allclose(m[:2, :2], [[h, h], [h, -h]], atol=1e-15)
    np.testing.assert_allclose(m[2:, 2:], [[h, h], [h, -h]], atol=1e-15)
    sq = g.squeezed_vacuum(1.3, 0.2)
    for eta in (0.0, 0.3, 0.92, 1.0):
        two = np.zeros((4, 4))
        two[np.ix_([0, 2], [0, 2])] = sq
        two[np.ix_([1, 3], [1, 3])] = g.vacuum(1)
        out = g.beam_splitter(two, eta, 0, 1)
        reduced = out[np.ix_([0, 2], [0, 2])]
        np.testing.assert_allclose(reduced, g.pure_loss(sq, eta), rtol=1e-12, atol=1e-15)
    with pytest.raises(ValueError):
        g.beam_splitter_matrix(1.2, 2, 0, 1)
    with pytest.raises(ValueError):
        g.beam_splitter_matrix(0.5, 2, 1, 1)


def test_pure_loss_per_mode_and_batch():
    cov = g.vacuum(2)
    cov = g.squeeze(cov, 1.0, 0.0, mode=0)
    cov = g.squeeze(cov, 0.5, 1.0, mode=1)
    cov = g.beam_splitter(cov, 0.3, 0, 1)  # correlations between the modes
    out = g.pure_loss(cov, np.array([0.8, 0.5]))
    # equals embedding of each loss as a BS with its own vacuum ancilla
    big = np.zeros((8, 8))
    idx_sys = [0, 1, 4, 5]
    big[np.ix_(idx_sys, idx_sys)] = cov
    big[np.ix_([2, 3, 6, 7], [2, 3, 6, 7])] = g.vacuum(2)
    big = g.beam_splitter(big, 0.8, 0, 2)
    big = g.beam_splitter(big, 0.5, 1, 3)
    np.testing.assert_allclose(out, big[np.ix_(idx_sys, idx_sys)], atol=1e-14)
    # scalar eta on a (n_f, 2, 2) batch
    batch = g.cov_from_noise_levels(np.array([0.1, 0.2]), np.array([3.0, 4.0]), 0.0)
    np.testing.assert_allclose(g.pure_loss(batch, 0.5), 0.5 * batch + 0.25 * np.eye(2))
    with pytest.raises(ValueError):
        g.pure_loss(cov, 1.1)
    with pytest.raises(ValueError):
        g.pure_loss(cov, np.array([0.5, 0.5, 0.5]))


def test_add_noise_phase_insensitive():
    """V → V + N for every θ, ΔV/V_SN = 2N (conv §2, T07 V11)."""
    sq = g.squeezed_vacuum(1.0, 0.3)
    out = g.add_noise(sq, 0.05)
    th = np.linspace(0, np.pi, 13)
    np.testing.assert_allclose(g.homodyne_variance(out, th) - g.homodyne_variance(sq, th), 0.05, rtol=1e-12)
    with pytest.raises(ValueError):
        g.add_noise(sq, -0.1)


# ------------------------------------------------------------------------------------------------
# jitter mapping (conv §6, T04)
# ------------------------------------------------------------------------------------------------


def test_jitter_average_gaussian_vs_monte_carlo():
    """T04 V4: exact formula equals the average of V(Δ) over Gaussian δ (statistical tolerance)."""
    vm, vp, sigma = 0.0437, 56.83, 0.1
    rng = stream(14, "test.gaussian.mc")
    delta = rng.normal(0.0, sigma, 1_000_000)
    v_mc = np.mean(vm * np.cos(delta) ** 2 + vp * np.sin(delta) ** 2)
    v_ex = g.jitter_average(vm, vp, g.chi_gaussian(sigma), delta0_rad=0.0)
    # std of the MC mean: vp·std(sin²δ)/√N ≈ vp·√2 σ²/√N
    assert abs(v_mc - v_ex) < 4 * vp * np.sqrt(2) * sigma**2 / 1e3
    # small-angle form differs by O(σ⁴): 0.027 dB at σ = 0.1 with the Ref A source (T04 §3.1)
    d_db = noise_db(g.jitter_average_small_angle(vm, vp, sigma)) - noise_db(v_ex)
    assert 0.02 < abs(d_db) < 0.035


def test_jitter_average_delta0_and_timeseries_chi():
    vm, vp = 0.1, 5.0
    # Δ0 = π/2 → anti-squeezing, χ = 1 → no jitter
    assert g.jitter_average(vm, vp, 1.0, delta0_rad=np.pi / 2) == pytest.approx(vp)
    assert g.jitter_average(vm, vp, 1.0, delta0_rad=0.0) == pytest.approx(vm)
    # static offset only: V(Δ0) = V− cos²Δ0 + V+ sin²Δ0
    assert g.jitter_average(vm, vp, 1.0, delta0_rad=0.3) == pytest.approx(vm * np.cos(0.3) ** 2 + vp * np.sin(0.3) ** 2)
    # time-series χ equals Gaussian χ for a deterministic Gaussian-quantile series
    rng = stream(15, "test.gaussian.chi")
    d = rng.normal(0.0, 0.05, 400_000)
    assert g.chi_timeseries(d) == pytest.approx(g.chi_gaussian(0.05), abs=3e-5)


def test_dither_ref_c():
    """T04 V5: ±4.5° dither → χ = J0(2m) = 0.99384; Ref C source (L = 21 %, R+ = 10 dB):
    −5.63 dB (no jitter) → −5.18 dB (dither only) → −4.32 dB (0.1 rad total)."""
    m = np.deg2rad(4.5)
    assert m / np.sqrt(2) == pytest.approx(0.0555, abs=5e-5)
    assert g.chi_sinusoid(m) == pytest.approx(0.99384, abs=1e-5)
    # sinusoidal time average
    t = np.linspace(0, 1, 200_000, endpoint=False)
    assert g.chi_timeseries(m * np.sin(2 * np.pi * 50 * t)) == pytest.approx(g.chi_sinusoid(m), abs=1e-10)
    loss, rp = 0.21, 10.0
    g2r = (rp - loss) / (1 - loss)
    rm = loss + (1 - loss) / g2r
    assert noise_db(rm, v_sn=1.0) == pytest.approx(-5.63, abs=0.005)
    v_dither = g.jitter_average(rm, rp, g.chi_sinusoid(m), delta0_rad=0.0)
    assert noise_db(v_dither, v_sn=1.0) == pytest.approx(-5.18, abs=0.005)
    v_total = g.jitter_average(rm, rp, g.chi_gaussian(0.1), delta0_rad=0.0)
    assert noise_db(v_total, v_sn=1.0) == pytest.approx(-4.32, abs=0.005)
    # independent components multiply
    both = g.jitter_average(rm, rp, g.chi_sinusoid(m) * g.chi_gaussian(0.083), delta0_rad=0.0)
    assert noise_db(both, v_sn=1.0) == pytest.approx(noise_db(v_total, v_sn=1.0), abs=0.01)


def test_ref_b_lossless_bound():
    """T04 V3: loss 0, optimal pump: 3.9° → −8.67 dB (small angle) / −8.68 dB (exact); 1.5° → −12.81 dB."""

    def best(sigma, exact):
        def f(r):
            vm, vp = np.exp(-2 * r), np.exp(2 * r)
            if exact:
                return g.jitter_average(vm, vp, g.chi_gaussian(sigma), delta0_rad=0.0)
            return g.jitter_average_small_angle(vm, vp, sigma)

        return noise_db(minimize_scalar(f, bounds=(0, 6), method="bounded", options={"xatol": 1e-10}).fun, v_sn=1.0)

    assert best(np.deg2rad(3.9), exact=False) == pytest.approx(-8.67, abs=0.006)
    assert best(np.deg2rad(3.9), exact=True) == pytest.approx(-8.68, abs=0.006)
    assert best(np.deg2rad(1.5), exact=False) == pytest.approx(-12.81, abs=0.006)
