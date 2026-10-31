"""VT-40: time → space demultiplexing (T11 §9 tests 1–10; 00_summary §4 VT-40; WP-B8).

1. Pure loss: tree, t = 0.99, ε = 0.005, 5 stages: η_sw = (0.99·0.995)⁵ = 0.92745; V_in = 0.05 →
   V_out = 0.082647 (−7.817 dB) through the LossRegistry and ``pure_loss`` (1e−12 relative to η V + (1 − η)/2).
2. Crosstalk ⟂ target (K = 4, 32; ER 10/20/30 dB): covariance of all (port, bin) modes after the switch
   network: (a) target bins of different ports uncorrelated (exactly 0), (b) target variance = pure loss
   with η_sw,k (< 1e−12), (c) photons in non-target bins = Σ_j (1 − (1 − ε)^S) sinh²r (t = 1; < 1e−12),
   and the leak output ≈ T11 §7.6 n̄_leak = 0.035 (ε = 1/200, r = 1; first order, 3 %).
3. Dilution (B5 ``expected_r_meas``; the time-domain check is VT-42): K = 32, η = 1, −10 dB, integrating
   −0.124 dB, pulse-resolved −10.00 dB.
4. Squeezing-angle sign: a switch phase ψ_sw = 0.3 rad (``demux_switch`` → segment driver → F6 path)
   moves the minimum-variance LO phase to θ_sq + 0.3 rad (conv §3–4), 1e−6.
5. Drive noise: 5 stages, κ_pm = 1, δV/V_π = 1.4235e−3 white per pulse → σ_sw = 10.0 mrad;
   R− = 0.1, R+ = 100: R′ = (R− + R+)/2 + (R− − R+)/2 e^{−2σ²} = 0.10999 (−9.587 dB). ``timeseries``:
   measured from the block's samples within 3σ_stat; ``declared``: from ``bands()`` (exact formula, 1e−12).
6. Bias drift 0.1 rad: stage transmission cos²(0.05) = 0.997502, leak 2.4979e−3 (−26.02 dB), 1e−6.
7. Delay lines: 100 MHz, n_g = 1.468, N_b = 1, L_0 = 0: L_1 = 2.0422 m, L_31 = 63.31 m; N_b = 50:
   L_31 = 3165 m; all target pulses reach the GBS together (1e−18 s); δT/T = 1e−6 → port 31 vs 0: 310 fs.
8. Overlap placeholder: δt = 0.0849 τ_p → η_mm = 0.990 (±5e−4).
9. Rise time: erf edge, t_r = 10 ns: leak 0.5 at the edge centre, < 1e−6 at ≥ 2 t_r; PC with
   f_tog,max = 1 MHz at f_rep = 100 MHz, N_b = 1 → ChkLog "demux_rise_time" and "demux_toggle_rate".
10. Block off ⇒ identity: design T = 𝟙 and no enabled loss, switch writes ψ = 0 (beam phase unchanged),
    delay line writes δt = 0.
"""

import math

import numpy as np
import pytest

import cranberry.blocks.dmx.testing as dt
import cranberry.path.testing as pt
from cranberry.blocks.dmx.demux import edge_leak, gaussian_overlap
from cranberry.blocks.err.testing import record
from cranberry.blocks.meas.pulse_sample import expected_r_meas
from cranberry.core.checks import ChkWarning
from cranberry.core.compile import run_reference
from cranberry.offline.gaussian import (chi_gaussian, cov_from_noise_levels, homodyne_variance, jitter_average,
                                        pure_loss, rotate)
from cranberry.offline.loss_registry import BeamPath, LossRegistry


def design_block(**kw):
    return dt.make(dt.design(**kw), 1e6)[0]


def test_1_pure_loss_tree():
    b = design_block()
    reg = LossRegistry(paths={"sq3": BeamPath("sq3", ("UP", "D3", "F3"))})
    reg.add_all(b.losses())
    eta = reg.eta("sq3", "gbs")
    assert eta == pytest.approx((0.99 * 0.995) ** 5, rel=1e-12) and eta == pytest.approx(0.92745, abs=5e-6)
    v = pure_loss(cov_from_noise_levels(0.05, 20.0, 0.0), eta)
    assert v[0, 0] == pytest.approx(eta * 0.05 + (1 - eta) / 2, rel=1e-12)
    assert v[0, 0] == pytest.approx(0.082647, abs=5e-7)
    assert 10 * math.log10(v[0, 0] / 0.5) == pytest.approx(-7.817, abs=1e-3)


def network_cov(T, v_minus, v_plus):
    """Output covariance (xxpp) of the K² (port k′, slot j) modes: Σ = M(V − 𝟙/2)Mᵀ + 𝟙/2, M[(k′, j), j] =
    √T[k′, target(j)], target(j) = K − 1 − j (T11 §3.1 numbering; leak phases do not matter)."""
    K = T.shape[0]
    M = np.zeros((K * K, K))
    for kp in range(K):
        for j in range(K):
            M[kp * K + j, j] = math.sqrt(T[kp, K - 1 - j])
    d = M @ M.T
    return np.block([[(v_minus - 0.5) * d, np.zeros_like(d)], [np.zeros_like(d), (v_plus - 0.5) * d]]) + 0.5 * np.eye(2 * K * K)


@pytest.mark.parametrize("K", [4, 32])
@pytest.mark.parametrize("er_db", [10.0, 20.0, 30.0])
def test_2_crosstalk_orthogonal_to_target(K, er_db):
    r = 1.0
    vm, vp = 0.5 * math.exp(-2 * r), 0.5 * math.exp(2 * r)
    b = design_block(K=K, er_db=er_db, toggles={"insertion_loss": False})
    T = b.matrix()
    cov = network_cov(T, vm, vp)
    n = K * K
    tgt = [k * K + (K - 1 - k) for k in range(K)]  # port k's target bin holds slot K − 1 − k
    sub = cov[np.ix_(tgt, tgt)]
    assert np.all(sub[~np.eye(K, dtype=bool)] == 0.0)  # (a)
    eta = b.design()["eta_sw"]
    np.testing.assert_allclose(np.diag(sub), eta * vm + (1 - eta) / 2, rtol=1e-12)  # (b)
    other = np.setdiff1d(np.arange(n), tgt)
    photons = 0.5 * (cov[other, other] + cov[other + n, other + n]) - 0.5
    S, eps = math.log2(K), 10 ** (-er_db / 10)
    assert photons.sum() == pytest.approx(K * (1 - (1 - eps) ** S) * math.sinh(r) ** 2, rel=1e-12)  # (c)
    np.testing.assert_allclose(b.design()["leak_in"], 1 - (1 - eps) ** S, rtol=1e-12)


def test_2_leak_photons_vs_first_order():
    b = design_block(er_db=10 * math.log10(200), toggles={"insertion_loss": False})
    n_leak = b.design()["leak_in"][7] * math.sinh(1.0) ** 2
    assert n_leak == pytest.approx(0.035, rel=0.03)
    assert n_leak == pytest.approx(5 / 200 * math.sinh(1.0) ** 2, rel=0.02)  # Σε sinh²r up to O(ε²)


def test_3_dilution():
    r = np.ones(32)
    r[0] = 0.1
    db = lambda x: 10 * math.log10(x)  # noqa: E731
    assert db(expected_r_meas(r, "integrating", [0], [1.0], 0.0)[0]) == pytest.approx(-0.124, abs=5e-4)
    assert db(expected_r_meas(r, "pulse_resolved", [0], [1.0], 0.0)[0]) == pytest.approx(-10.00, abs=1e-9)


def switch_path(sw_obj):
    seg = pt.seg_spec("dmx", 0.0, drivers={"eom": "per_sub"})
    topo = pt.topo_spec([seg], [pt.beam_spec("sq", [pt.step_spec("dmx", 1)], "gbs")], ["gbs"],
                        pt.coef_spec("dmx", "eom", "*", 1.0))
    return [sw_obj, pt.pathtest("pathtest_table", "src", "laser.las.phi", (), "rad", n_rows=4),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "act.sw.dpsi"], 4)]


def test_4_squeezing_angle_sign():
    theta_sq = 0.2
    sw = dt.switch(n_stages=1, toggles={"drive_phase_noise": False, "bias_drift": True}, bias_offset_rad=0.3)
    rec, _, _ = pt.run_and_record(switch_path(sw), 1e6, 4, {"laser.las.phi_table": [theta_sq] * 4})
    psi = rec["beam.sq.phi"][-1] - theta_sq
    assert psi == pytest.approx(0.3, abs=1e-15)
    cov = rotate(cov_from_noise_levels(0.05, 20.0, theta_sq), psi)  # â → â e^{iψ} (conv §4)
    th = np.linspace(-0.5, 1.0, 150001)
    assert th[np.argmin(homodyne_variance(cov, th))] == pytest.approx(theta_sq + 0.3, abs=1e-5)
    off = dict(sw, params={**sw["params"], "toggles": {**sw["params"]["toggles"], "enabled": False}})
    rec, _, _ = pt.run_and_record(switch_path(off), 1e6, 4, {"laser.las.phi_table": [theta_sq] * 4})
    assert rec["beam.sq.phi"][-1] == theta_sq  # test 10: off ⇒ identity


R_M, R_P = 0.1, 100.0
R_PRIME = (R_M + R_P) / 2 + (R_M - R_P) / 2 * math.exp(-2 * 0.0100**2)


@pytest.mark.numba
def test_5_drive_noise_timeseries():
    assert R_PRIME == pytest.approx(0.10999, abs=5e-6)
    assert 10 * math.log10(R_PRIME) == pytest.approx(-9.587, abs=1e-3)
    n = 400_000
    psi = record([dt.switch()], 1e6, n, ["act.sw.dpsi"], external=())["act.sw.dpsi"]
    assert psi.std() == pytest.approx(0.0100, rel=0.01)
    r = (R_M + R_P) / 2 + (R_M - R_P) / 2 * np.cos(2 * psi)
    assert abs(r.mean() - R_PRIME) < 3 * r.std() / math.sqrt(n)


def test_5_drive_noise_declared():
    b = dt.make(dt.switch(mode="declared"), 1e6)[0]
    sigma = math.sqrt(b.bands()[0].variance_outside)
    assert sigma == pytest.approx(0.0100, abs=1e-6)
    r = jitter_average(R_M, R_P, chi_gaussian(sigma), delta0_rad=0.0)
    assert r == pytest.approx((R_M + R_P) / 2 + (R_M - R_P) / 2 * math.exp(-2 * sigma**2), rel=1e-12)
    assert r == pytest.approx(0.10999, abs=5e-6)


def test_6_bias_drift():
    b, st, pr = dt.make(dt.switch(n_stages=1, toggles={"drive_phase_noise": False, "bias_drift": True},
                                  bias_offset_rad=0.1), 1e6)
    run_reference(b, st, pr, 2)
    assert st["pow.sw.trans"] == pytest.approx(0.997502, abs=1e-6)
    assert st["leak.sw.frac"] == pytest.approx(2.4979e-3, rel=1e-4)
    assert 10 * math.log10(st["leak.sw.frac"]) == pytest.approx(-26.02, abs=0.01)


def test_7_delay_lines():
    d = design_block().design()
    assert d["L_m"][1] == pytest.approx(2.0422, abs=5e-5) and d["L_m"][31] == pytest.approx(63.31, abs=5e-3)
    assert design_block(N_b=50).design()["L_m"][31] == pytest.approx(3165, abs=0.5)
    t_rep = 1e-8
    arrive = (31 - np.arange(32)) * t_rep + d["tau_s"]  # slot σ(k) = K − 1 − k reaches the switch at σ(k)·T_rep
    assert np.ptp(arrive) < 1e-18
    shift = {}
    for k in (0, 31):
        b, st, pr = dt.make(dt.delay(port=k, rep_rel_drift=1e-6, toggles={"static_mismatch": False}), 1e6)
        run_reference(b, st, pr, 1)
        shift[k] = st["timing.dl.dt"]
    assert shift[31] - shift[0] == pytest.approx(310e-15, rel=1e-9)


def test_8_overlap_placeholder():
    assert gaussian_overlap(0.0849, 1.0) == pytest.approx(0.990, abs=5e-4)
    assert gaussian_overlap(0.0, 1.0) == 1.0


def test_9_rise_time():
    assert edge_leak(0.0, 10e-9) == pytest.approx(0.5, rel=1e-12)
    assert np.all(edge_leak(np.array([20e-9, 30e-9, -25e-9]), 10e-9) < 1e-6)
    b = design_block()
    with pytest.warns(ChkWarning):
        b.check()
    assert {"demux_rise_time", "demux_toggle_rate"} <= {e.code for e in b.chk_log.events}


def test_10_off_is_identity():
    b = design_block(toggles={"enabled": False})
    np.testing.assert_array_equal(b.matrix(), np.eye(32))
    assert not any(e.enabled for e in b.losses())
    assert np.all(b.design()["eta_sw"] == 1.0) and np.all(b.design()["leak_in"] == 0.0)
    b, st, pr = dt.make(dt.delay(dtau_s=1e-12, rep_rel_drift=1e-6, toggles={
        "static_mismatch": False, "delay_thermal": False, "rep_rate_drift": False, "fiber_loss": False}), 1e6)
    run_reference(b, st, pr, 1)
    assert st["timing.dl.dt"] == 0.0 and not b.losses()[0].enabled
