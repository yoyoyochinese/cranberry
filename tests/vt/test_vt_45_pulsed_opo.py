"""VT-45 (WP-B11): pulse-pumped OPO, ``pulsed_opo_source`` (00_summary §4 VT-45; T17 §7 V-P1–V-P15).

Also the pulsed-OPO parts of VT-05 (pure loss, V-P13) and VT-06 (vacuum and the initial-cavity
vacuum term, V-P2 / V-P3). Tolerances from T17 §7 unless noted:

* V-P1 CW limit (ε const, window [0, 40/κ_c], flat mode of 20/κ_c in the middle, x = 0.5, ρ = 0.95):
  R_p = 0.183704, R_q = 7.84003; (R − 1) relative error < 1e−4 at δt = 1/(80κ_c); halving δt divides
  the step-to-step change by ≈ 4 (δt²); kernel element C(t, t′) = δt·2κ_T(κ_c/λ − 1)e^{−λ|τ|}.
* V-P2 vacuum: x_ss = 0 → every eigenvalue 1 within 2(κ_c δt)². V-P3: without the initial-vacuum
  columns the minimum eigenvalue at ρ = 1 is ≈ 0 (regression); with them it is 1.
* V-P4 purity: ρ = 1 with detuning and a resonant detuned pump: λ_i λ_{2N+1−i} = 1 within 1.3e−3
  at δt = 1/(40κ_c), decreasing ≈ δt².
* V-P5 Borealis K = 1.10 ± 0.01; V-P6 Liu 1/K = 0.975 ± 0.005.
* V-P7 long pump → CW (non-resonant, x = 0.5): r_0 increases with T_pκ_c = 10, 30 towards 2 artanh x;
  at 30: r_0 within 3 %, detuning rotation within 1 % of −1/(1 + x²) = −0.8 (slow).
* V-P8 phase convention: squeezed quadrature p at φ = 0; φ_OPA,0 = 0.4 → θ_sq = π/2 + 0.2 (1e−4 rad).
* V-P9 detuning (Borealis, ρ = 0.95, x_ss = 1.3): −0.638 / −0.425 / −1.063 / −0.543 per Δ_s/κ_c
  (±0.005), same at d = 0.05 (±0.002); d > 0 lowers θ_sq.
* V-P10 timing: overlap 0.99 / 0.999 at |Δt| = 0.090/κ_c / 0.028/κ_c (±3 %), both presets.
* V-P11 tail (Borealis, x_ss = 2.4): ε_tail(10 ns) = 0.138, ε_tail(20 ns) = 2.3e−3, ≈ 8.8e^{−2κ_c T} (±10 %).
* V-P12 short pump T_pκ_c = 0.1, r_0 ≈ 0.5: non-resonant K = 1.001, doubly resonant 1.086 (±0.005).
* V-P13 pure loss: optimal-mode R_− ≈ ρe^{−2r_0} + 1 − ρ (±0.05 dB); T17 §2.6 x_ss ↔ r_0 table (1 %).
* V-P14 noise in the mode: N = Φ|ũ(f₀)|² (narrow band), N = n₀ (flat): |ũ(0)|² = τ_eff, ∫|ũ|²df = 1
  (1e−6), for the optimal mode and a 1 ps Gaussian (τ_eff = 2√π T0 = 2.13 τ_FWHM).
* V-P15 FSR comb: T_f = 0 → n_FSR = 0, T_f = 1 → Σ 2 sinh² r_m (exact).
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from cranberry.blocks.fib.mode_overlap import gaussian_mode
from cranberry.blocks.source import temporal_modes as tm
from cranberry.blocks.source.pulsed_opo import transfer_matrix
from cranberry.blocks.source.pulsed_testing import BOREALIS, LIU, pulsed_opo
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline import gaussian


# ---------------------------------------------------------------------------------------------- V-P1


def _cw_r(nper):
    kappa, x, rho = 1.0, 0.5, 0.95
    dt = 1.0 / (nper * kappa)
    n = int(round(40.0 / dt))
    t = transfer_matrix(np.full(n, x * kappa + 0j), dt, kappa, rho, 0.0, True)
    tc = (np.arange(n) + 0.5) * dt
    w = ((tc > 10) & (tc < 30)).astype(float)
    w /= np.linalg.norm(w)
    return np.sum((w @ t[n:]) ** 2), np.sum((w @ t[:n]) ** 2), t, dt, n


def _flat_analytic(lam, kt=0.95, kappa=1.0, length=20.0):
    return 1 + 2 * kt * (kappa / lam - 1) * (2 / lam) * (1 - (1 - np.exp(-lam * length)) / (lam * length))


def test_v_p1_cw_limit():
    assert _flat_analytic(1.5) == pytest.approx(0.183704, abs=1e-6)
    assert _flat_analytic(0.5) == pytest.approx(7.84003, abs=1e-5)
    res = {n: _cw_r(n) for n in (20, 40, 80)}
    rp, rq, t, dt, n = res[80]
    assert abs(rp - _flat_analytic(1.5)) / abs(_flat_analytic(1.5) - 1) < 1e-4
    assert abs(rq - _flat_analytic(0.5)) / abs(_flat_analytic(0.5) - 1) < 1e-4
    for k in (0, 1):  # Richardson: the change from δt → δt/2 shrinks ≈ 4× (δt²)
        d1, d2 = res[20][k] - res[40][k], res[40][k] - res[80][k]
        assert 3.0 < d1 / d2 < 5.0
    # kernel element C(t, t') for t ≠ t' in the stationary part: p quadrature, λ = κ(1 + x)
    c = t @ t.T
    i, j = int(15 / dt), int(16 / dt)
    lam = 1.5
    ana = dt * 2 * 0.95 * (1 / lam - 1) * np.exp(-lam * (j - i) * dt)
    assert c[n + i, n + j] == pytest.approx(ana, rel=1e-3)


# ---------------------------------------------------------------------------------------- V-P2, V-P3


def test_v_p2_vacuum_and_v_p3_initial_vacuum_term():
    b = pulsed_opo(BOREALIS, x_ss=0.0, escape_efficiency=1.0)
    lam = np.linalg.eigvalsh(b.cov_kernel)
    assert np.max(np.abs(lam - 1)) <= 2 * BOREALIS["dt_kappa"] ** 2
    bad = pulsed_opo(BOREALIS, x_ss=0.0, escape_efficiency=1.0, toggles={"include_initial_vacuum": False})
    assert np.linalg.eigvalsh(bad.cov_kernel)[0] < 1e-3  # spurious "perfect squeezing" (regression)
    rho = pulsed_opo(BOREALIS, x_ss=0.0, escape_efficiency=0.9, toggles={"include_initial_vacuum": False})
    assert np.linalg.eigvalsh(rho.cov_kernel)[0] == pytest.approx(0.1, abs=2e-3)  # 1 − ρ without the term


# ---------------------------------------------------------------------------------------------- V-P4


def _purity_err(dt_kappa):
    kc = np.pi * BOREALIS["fwhm_hz"]
    b = pulsed_opo(BOREALIS, x_ss=2.4, escape_efficiency=1.0, dt_kappa=dt_kappa, delta_s_rad_s=0.05 * kc,
                   delta_p_rad_s=0.1 * kc, toggles={"detuning": True})
    lam = np.linalg.eigvalsh(b.cov_kernel)
    return np.max(np.abs(lam * lam[::-1] - 1))


def test_v_p4_purity():
    e40, e80 = _purity_err(1 / 40), _purity_err(1 / 80)
    assert e40 <= 1.3e-3
    assert 3.0 < e40 / e80 < 5.0


# ---------------------------------------------------------------------------------------- V-P5, V-P6


def test_v_p5_borealis_schmidt():
    b = pulsed_opo(BOREALIS, escape_efficiency=1.0)
    assert b.schmidt_number == pytest.approx(1.10, abs=0.01)
    assert b.supermode_r[0] == pytest.approx(0.533, abs=0.01)  # T17 §2.6: x_ss = 1.297 ↔ r_0 = 0.533


def test_v_p6_liu_purity():
    b = pulsed_opo(LIU, escape_efficiency=1.0)
    assert 1.0 / b.schmidt_number == pytest.approx(0.975, abs=0.005)
    assert b.supermode_r[0] == pytest.approx(1.8, rel=0.02)


# ---------------------------------------------------------------------------------------------- V-P7


@pytest.mark.slow
def test_v_p7_long_pump_to_cw():
    kc = np.pi * BOREALIS["fwhm_hz"]
    r0, rot = [], []
    for tpk in (10.0, 30.0):
        b = pulsed_opo(BOREALIS, x_ss=0.5, escape_efficiency=1.0, pump_duration_s=tpk / kc, window_kappa=tpk + 10.0,
                       toggles={"pump_resonant": False})
        r0.append(b.supermode_r[0])
        rot.append(b.angle_coefficients(0.5)[0] * kc)
    cw_r, cw_rot = 2 * np.arctanh(0.5), -1 / (1 + 0.25)
    assert r0[0] < r0[1] < cw_r and abs(r0[1] - cw_r) / cw_r < 0.03
    assert abs(rot[1] - cw_rot) < abs(rot[0] - cw_rot) and abs(rot[1] - cw_rot) / abs(cw_rot) < 0.01
    assert r0 == pytest.approx([0.996, 1.078], abs=0.01) and rot == pytest.approx([-0.775, -0.794], abs=0.005)


# ---------------------------------------------------------------------------------------------- V-P8


def test_v_p8_phase_convention():
    b0 = pulsed_opo(BOREALIS)
    u, _ = b0.optimal_mode
    assert np.max(np.abs(u.imag)) < 1e-12  # real mode at φ = 0, Δ = 0
    th0 = gaussian.noise_levels_from_cov(b0.source_output.lo_mode_cov(u))[2]
    assert th0 == pytest.approx(np.pi / 2, abs=1e-4)
    b1 = pulsed_opo(BOREALIS, phi_opa0_rad=0.4)
    th1 = gaussian.noise_levels_from_cov(b1.source_output.lo_mode_cov(u))[2]
    assert th1 == pytest.approx(np.pi / 2 + 0.2, abs=1e-4)
    assert float(b1.source_output.project(b1.optimal_mode[0])[2]) == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------------------------- V-P9


@pytest.mark.parametrize("d, tol", [(0.01, 0.005), (0.05, 0.005)])
def test_v_p9_detuning_coefficients(d, tol):
    b = pulsed_opo(BOREALIS, x_ss=1.3, angle_probe_d=d)
    cs, cp, cb = (c * b.kappa_c for c in b.angle_coefficients(1.3))
    both = cs + 2 * cp  # Δ_p = 2Δ_s for a cavity-length change
    locked = both - 0.5 * 2 * cb  # transmitted-pump lock removes ½ arg∫b dt
    assert cs == pytest.approx(-0.638, abs=tol)
    assert 2 * cp == pytest.approx(-0.425, abs=tol)
    assert both == pytest.approx(-1.063, abs=tol)
    assert locked == pytest.approx(-0.543, abs=tol)
    assert cs < 0 and cp < 0  # d > 0 (shorter cavity) lowers θ_sq


def test_v_p9_linear_in_d():
    c1 = np.array(pulsed_opo(BOREALIS, x_ss=1.3, angle_probe_d=0.01).angle_coefficients(1.3))
    c5 = np.array(pulsed_opo(BOREALIS, x_ss=1.3, angle_probe_d=0.05).angle_coefficients(1.3))
    kc = np.pi * BOREALIS["fwhm_hz"]
    assert np.max(np.abs(c1 - c5)) * kc < 0.002


# --------------------------------------------------------------------------------------------- V-P10


@pytest.mark.parametrize("base", [BOREALIS, LIU], ids=["borealis", "liu"])
def test_v_p10_timing_overlap(base):
    b = pulsed_opo(base, x_ss=2.4)
    u, g, kc = b.optimal_mode[0], b.grid, b.kappa_c

    def eta(s):
        return abs(tm.overlap(u, tm.transform_mode(u, g, dt_s=s / kc), g)) ** 2

    assert brentq(lambda s: eta(s) - 0.99, 1e-3, 1.0) == pytest.approx(0.090, rel=0.03)
    assert brentq(lambda s: eta(s) - 0.999, 1e-3, 1.0) == pytest.approx(0.028, rel=0.03)


# --------------------------------------------------------------------------------------------- V-P11


def test_v_p11_tail():
    b = pulsed_opo(BOREALIS, x_ss=2.4)
    assert b.tail(10e-9) == pytest.approx(0.138, rel=0.10)
    assert b.tail(20e-9) == pytest.approx(2.3e-3, rel=0.10)
    for kt in (2.0, 3.0, 4.0, 5.0):
        assert b.tail(kt / b.kappa_c) == pytest.approx(8.8 * np.exp(-2 * kt), rel=0.10)
    assert b.source_output.tail_into_next_slot == pytest.approx(b.tail(1 / 6e6), abs=1e-15)
    assert pulsed_opo(BOREALIS, x_ss=2.4, rep_rate_hz=100e6).source_output.tail_into_next_slot == pytest.approx(0.138, rel=0.1)


# --------------------------------------------------------------------------------------------- V-P12


def test_v_p12_short_pump():
    kc = np.pi * BOREALIS["fwhm_hz"]
    common = dict(escape_efficiency=1.0, pump_duration_s=0.1 / kc, dt_kappa=1 / 80, window_kappa=10.0)
    nonres = pulsed_opo(BOREALIS, x_ss=5.2, toggles={"pump_resonant": False}, **common)
    res = pulsed_opo(BOREALIS, x_ss=7.3, **common)
    assert nonres.supermode_r[0] == pytest.approx(0.5, abs=0.02) and res.supermode_r[0] == pytest.approx(0.5, abs=0.02)
    assert nonres.schmidt_number == pytest.approx(1.001, abs=0.005)
    assert res.schmidt_number == pytest.approx(1.086, abs=0.005) and res.schmidt_number >= 1.06


# --------------------------------------------------------------------------------------------- V-P13


@pytest.mark.parametrize("rho", [0.95, 0.8, 0.5])
def test_v_p13_pure_loss(rho):
    b = pulsed_opo(BOREALIS, x_ss=2.5, escape_efficiency=rho)
    r0 = b.supermode_r[0]
    single = 10 * np.log10(rho * np.exp(-2 * r0) + 1 - rho)
    assert 10 * np.log10(b.optimal_mode[1]) == pytest.approx(single, abs=0.05)
    if rho == 0.95:  # T17 §5.2 quotes −7.66 dB with r_0 = 1.029; its own §2.6 table gives r_0 ≈ 1.04 at 2.5
        assert 10 * np.log10(b.optimal_mode[1]) == pytest.approx(-7.66, abs=0.1)


@pytest.mark.parametrize("base, table", [
    (BOREALIS, [(1.297, 0.533), (2.136, 0.886), (2.401, 1.0), (3.531, 1.5), (4.180, 1.8)]),
    (LIU, [(1.279, 0.533), (2.107, 0.886), (2.369, 1.0), (3.486, 1.5), (4.129, 1.8)]),
], ids=["borealis", "liu"])
def test_r0_vs_x_ss_table(base, table):
    """T17 §2.6 [sim] x_ss ↔ r_0 (ρ = 1) within 1 % (this discretisation is −0.5 % / +0.8 % systematic)."""
    for x, r0 in table:
        assert pulsed_opo(base, x_ss=x, escape_efficiency=1.0).supermode_r[0] == pytest.approx(r0, rel=0.01)


# --------------------------------------------------------------------------------------------- V-P14


def _check_mode_noise(u, grid):
    f, s = tm.mode_spectrum(u, grid, pad=8)
    df = f[1] - f[0]
    assert np.sum(s) * df == pytest.approx(1.0, abs=1e-6)  # flat n₀: N = n₀ ∫|ũ|² df = n₀
    assert s[np.argmin(np.abs(f))] == pytest.approx(tm.tau_eff_s(u, grid), rel=1e-6)  # Φ at f₀ = 0: N = Φ τ_eff
    return f, s


def test_v_p14_narrow_and_broadband_noise():
    b = pulsed_opo(BOREALIS, x_ss=2.4)
    _check_mode_noise(b.optimal_mode[0], b.grid)
    assert b.tau_eff_s * b.kappa_c == pytest.approx(3.00, abs=0.06)  # T17 §5.2: 14.5 ns
    tau = 1e-12
    t0 = tau / (2 * np.sqrt(np.log(2)))
    g = (np.arange(4096) - 2048) * (40 * t0 / 4096)
    u = gaussian_mode(g, t0)
    u /= np.sqrt(np.sum(u**2) * (g[1] - g[0]))
    _check_mode_noise(u, g)
    assert tm.tau_eff_s(u, g) == pytest.approx(2 * np.sqrt(np.pi) * t0, rel=1e-6)
    assert tm.tau_eff_s(u, g) / tau == pytest.approx(2.13, abs=0.01)


# --------------------------------------------------------------------------------------------- V-P15


def _fsr(stages, fsr_modes=True):
    return make_block(BlockSpec("fsr_comb_filter", "fsr", dict(
        fsr_hz=3.7e9, stages=stages, r_m=lambda m: 0.5 + 0.0 * m, m_max=20, insertion_eta=1.0, segment="opo.out",
        applies_to=["gbs"], source_tag="[T17 §2.7]", toggles={"fsr_modes": fsr_modes})))


def test_v_p15_fsr_comb_limits():
    assert _fsr([]).n_fsr == 20 * 2 * np.sinh(0.5) ** 2  # T_f = 1
    block_all = [{"type": "flat_top", "center_hz": 0.0, "bandwidth_hz": 1e9, "order": 50, "floor": 0.0}]
    assert _fsr(block_all).n_fsr == 0.0  # T_f(m·FSR) = 0 for every m ≥ 1, passband kept
    assert _fsr(block_all).eta_passband == 1.0
    assert _fsr([], fsr_modes=False).n_fsr == 0.0  # toggle off = ideal comb removal
