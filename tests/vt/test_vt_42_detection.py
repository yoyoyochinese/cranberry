"""VT-42 detection (WP-B5): T13 §6 tests 1 (Monte Carlo part), 3–8.

Sources: notes/bhd-detection.md (T13) §1, §4.1, §6; notes/00_summary.md VT-42; PDB450C values from
devices/pdb450c.yaml (PDB manual v2.4, printed pages cited there). Tolerances: analytic 1e−12 relative,
quoted rounded values to their last digit, Monte Carlo to a few standard errors (PLAN §6.3).
Tests 5–6 drive B3's ``dc_fringe`` with the device preset (the Ref D balanced detector, G16-d).
"""

import math

import numpy as np
import pytest

import cranberry.blocks.det.testing as T
from cranberry.blocks.det import bhd
from cranberry.blocks.err.dc_fringe import fringe_lock_point
from cranberry.blocks.err.testing import dcf, rig
from cranberry.core.block import BlockSpec
from cranberry.core.compile import run_reference
from cranberry.core.registry import make_block
from cranberry.core.rng import stream
from cranberry.core.units import lin_to_db
from cranberry.hardware.devices import device_values, load_device

pytestmark = pytest.mark.vt


def mk(obj):
    return make_block(BlockSpec.from_obj(obj))


def test_1_subtract_dark_monte_carlo():
    """R_corr = (S − S_dark)/(S_SN − S_dark) → R; Monte Carlo within a few √(2/N) (T13 test 1)."""
    c, r, n = 10.0, 0.1, 1_000_000
    fe = mk(T.frontend(clearance_lin=c, toggles={"imbalance_leak": False, "freq_response": False, "dip": False}))
    assert fe.psd(0.0, r)["r_meas"] == pytest.approx((c * r + 1) / (c + 1), rel=1e-12)
    fe_sub = mk(T.frontend(clearance_lin=c, toggles={"imbalance_leak": False, "freq_response": False, "dip": False,
                                                    "subtract_dark": True}))
    assert fe_sub.psd(0.0, r)["r_meas"] == pytest.approx(r, rel=1e-12)
    rec = fe.samples(r, n, stream(42, "vt42.mc"))
    v = {k: x.var() for k, x in rec.items()}
    r_corr = (v["signal"] - v["dark"]) / (v["shot"] - v["dark"])
    tol = 5 * math.sqrt(2 / n) * (1 + 1 / (c * r)) * math.sqrt(3)
    assert abs(r_corr / r - 1) < tol


def test_3_imbalance_leak():
    """δ = 0: R_meas = R for any RIN; |δ|² = 2.5e−4, N_ex = 60: −10 dB → −9.45 dB; vacuum stays 1 (T13 test 3)."""
    base = {"circuit_noise": False, "freq_response": False, "dip": False}
    for n_ex in (0.0, 60.0, 1e6):
        fe = mk(T.frontend(imbalance_abs2=0.0, n_ex=n_ex, toggles=base))
        assert fe.psd(0.0, 0.1)["r_meas"] == pytest.approx(0.1, rel=1e-12)
    fe = mk(T.frontend(toggles=base))
    r_meas = fe.psd(0.0, 0.1)["r_meas"]
    assert r_meas == pytest.approx(0.1135, abs=5e-5) and lin_to_db(r_meas) == pytest.approx(-9.45, abs=0.005)
    assert fe.psd(0.0, 1.0)["r_meas"] == pytest.approx(1.0, rel=1e-12)
    assert r_meas == pytest.approx(bhd.r_meas_imbalance(0.1, 2.5e-4, 60.0), rel=1e-12)
    # N_ex = 60 is the LO RIN −150 dB/Hz at 15.4 mW, η_q = 1 (T13 §4.2); CMRR 30 dB ↔ |δ|² = 2.5e−4
    assert bhd.n_ex_from_rin(0.0154, 1e-15, 1.0, 1.55e-6) == pytest.approx(60.0, abs=0.1)
    assert bhd.imbalance_from_cmrr(1000.0) == pytest.approx(2.5e-4)
    det = mk(T.det_eff(s_rin_per_hz=1e-15, n_ex=None, wavelength_m=1.55e-6, qe=1.0,
                       toggles={"circuit_noise": False, "freq_response": False, "dip": False, "qe": False}))
    assert det.breakdown()["eta_bal"] == pytest.approx(1 / (1 + 2.5e-4 * 61), rel=2e-3)


def test_4_frequency_response_and_dip():
    """Dip at 0.96 GHz, r = 0.9: |H| min at f_dip with depth (1 − r)/(1 + r); without circuit noise R(f)
    is unchanged, with it the equivalent loss grows only near the dip (T13 test 4)."""
    f = np.linspace(0.5e9, 1.4e9, 90001)
    fe = mk(T.frontend(toggles={"freq_response": False, "imbalance_leak": False, "circuit_noise": False}))
    h = np.abs(fe.h(f))
    assert f[np.argmin(h)] == pytest.approx(0.96e9, abs=1e4)
    assert h.min() == pytest.approx(0.1 / 1.9, rel=1e-6)
    np.testing.assert_allclose(fe.psd(f, 0.3)["r_meas"], 0.3, rtol=1e-12)
    fe_el = mk(T.frontend(clearance_lin=100.0, toggles={"freq_response": False, "imbalance_leak": False}))
    loss = fe_el.psd(f, 0.0)["r_meas"]  # R = 0 → R_meas = 1/(1 + C) = 1 − η_e
    ref = 1 - bhd.eta_electronic(100.0)
    assert loss[np.argmin(h)] == pytest.approx(1 / (1 + 100 * (0.1 / 1.9) ** 2), rel=1e-9)
    np.testing.assert_allclose(loss, 1 / (1 + 100 * h**2), rtol=1e-12)
    assert f[np.argmax(loss)] == f[np.argmin(h)] and loss[np.argmin(h)] > 30 * ref
    assert fe_el.psd(0.0, 0.0)["r_meas"] == pytest.approx(ref, rel=1e-12)  # no extra loss away from the dip (DC)
    # the same C(f) in the registered loss of detector_efficiency
    det = mk(T.det_eff(clearance_lin=100.0, f_s_hz=0.96e9, toggles={"freq_response": False, "qe": False,
                                                                     "imbalance_leak": False}))
    assert 1 - det.breakdown()["eta_det"] == pytest.approx(loss[np.argmin(h)], rel=1e-6)


def _pdb(**over):
    """dc_fringe with the PDB450C device values (Ref D, high-Z into the DSC1)."""
    dev = device_values(load_device("pdb450c"))
    vals = {k: dev["rf"][k] for k in ("gain_v_per_a", "v_sat_v", "v_off_v", "i_n_a_rthz")}
    vals["responsivity_a_per_w"] = dev["photodiode"]["responsivity_a_per_w"]
    return dcf(**{**vals, **over})


def _fringe(obj, phi):
    root, s, p = rig([obj], 1e6)
    out = []
    for x in phi:
        s["beam.lo.phi"] = x
        run_reference(root, s, p, 1)
        out.append(float(s["v.dcf"]))
    return np.array(out)


def test_5_ref_d_fringe_and_clipping():
    """P_in 3.33 mW, ℜ 1.0, G 10³: ±3.33 V, slope −3.33 V/rad at π/2 (s = +1); P_in 5 mW clips at
    ±4.6 V (PDB) and ±4.0 V (DSC1 input); 50 Ω load halves amplitude and slope (T13 test 5)."""
    phi = np.linspace(0, 2 * np.pi, 721)
    v = _fringe(_pdb(), phi)
    assert v.max() == pytest.approx(3.33, rel=1e-12) and v.min() == pytest.approx(-3.33, rel=1e-12)
    d = 1e-6
    slope = (_fringe(_pdb(), [np.pi / 2 + d])[0] - _fringe(_pdb(), [np.pi / 2 - d])[0]) / (2 * d)
    assert slope == pytest.approx(-3.33, rel=1e-6)
    v5 = _fringe(_pdb(p_in_w=5e-3), phi)
    assert v5.max() == 4.6 and v5.min() == -4.6
    lo, hi = device_values(load_device("dsc1"))["adc"]["range_v"]
    v_dsc1 = np.clip(v5, lo, hi)
    assert v_dsc1.max() == 4.0 and v_dsc1.min() == -4.0
    assert np.clip(v, lo, hi).max() == pytest.approx(3.33)  # 14–17 % margin at 3.33 mW (T13 §4.1)
    div = device_values(load_device("pdb450c"))["rf"]["load_divisor_50ohm"]
    v50 = _fringe(_pdb(gain_v_per_a=1e3 / div, v_sat_v=4.6 / div), phi)
    np.testing.assert_allclose(v50, v / 2, rtol=1e-12)
    fe = mk(T.frontend(device="pdb450c", load="50ohm"))  # the same through bhd_frontend
    assert fe.dc_voltage(3.33e-3, 0.0) == pytest.approx(3.33 / 2)


def test_6_null_shift_and_rin():
    """ε = 0.03, V = 1: null at π/2 + arcsin(0.03); V_off = 0: 1 % P_in step does not move it; V_off =
    15 mV (PDB offset bound): shift 4.5 mrad, moved 45 µrad by the step (T13 test 6)."""
    phi0, _ = fringe_lock_point(3.33, 0.03, 1.0, 0.0, 0.0)
    assert phi0 == pytest.approx(np.pi / 2 + math.asin(0.03), abs=1e-12)
    assert phi0 - np.pi / 2 == pytest.approx(0.030, abs=1e-4)
    assert _fringe(_pdb(eps_r=0.03), [phi0])[0] == pytest.approx(0.0, abs=1e-12)
    phi1, _ = fringe_lock_point(3.33 * 1.01, 0.03, 1.0, 0.0, 0.0)
    assert phi1 == pytest.approx(phi0, abs=1e-15)
    v_off = device_values(load_device("pdb450c"))["rf"]["v_off_v"]
    v_off_max = load_device("pdb450c")["rf"]["v_off_v"]["sweep"][1]
    assert v_off == 0.0 and v_off_max == 0.015
    a, _ = fringe_lock_point(3.33, 0.0, 1.0, 0.0, v_off_max)
    b, _ = fringe_lock_point(3.33 * 1.01, 0.0, 1.0, 0.0, v_off_max)
    assert a - np.pi / 2 == pytest.approx(4.5e-3, abs=0.05e-3)
    assert abs(b - a) == pytest.approx(45e-6, abs=1e-6)


def test_7_inter_pulse_crosstalk():
    """First-order h(t), short window: adjacent correlation of an independent vacuum train = Σc_m c_{m+1}/Σc_m²
    (= c_1); → 0 for B → ∞; per-pulse normalized variance = R whatever c_m for a uniform state (T13 test 7)."""
    pk = mk(T.pulsed_kernel(toggles={"circuit_noise": False}))
    assert pk.rho_adj == pytest.approx(math.exp(-math.pi), rel=1e-9)
    ps = mk(T.pulse_sample(target_slots=list(range(32)), n_frames=8000, toggles={"electronic_noise": False}))
    rng = stream(7, "vt42.ps")
    out = ps.evaluate(np.ones(32), c=pk.c, clearance_lin=math.inf, rng=rng)
    n = 32 * 8000
    assert out["rho_adj_sample"] == pytest.approx(pk.rho_adj, abs=4 / math.sqrt(n))
    fast = mk(T.pulsed_kernel(response={"model": "first_order", "bandwidth_hz": 1e10}, toggles={"circuit_noise": False}))
    assert fast.rho_adj < 1e-100
    out = ps.evaluate(np.ones(32), c=fast.c, clearance_lin=math.inf, rng=rng)
    assert abs(out["rho_adj_sample"]) < 4 / math.sqrt(n)
    for c in (pk.c, mk(T.pulsed_kernel(response={"model": "second_order", "f0_hz": 1e8, "q": 3.0})).c):
        out = ps.evaluate(np.full(32, 0.3), c=c, clearance_lin=math.inf, rng=rng)
        np.testing.assert_allclose(out["expected"], 0.3, rtol=1e-12)
        assert out["sample"].mean() == pytest.approx(0.3, rel=4 * math.sqrt(4 / n))


def test_8_shot_noise_scaling():
    """Vacuum input: dark-subtracted output PSD ∝ P_LO; with circuit noise C ∝ P_LO (T13 test 8)."""
    p_lo = np.array([1e-3, 2e-3, 5e-3, 15.4e-3])
    sn, cl = [], []
    for p in p_lo:
        fe = mk(T.frontend(p_lo_w=p, clearance_lin=None, i_n_a_rthz=69e-12,
                           toggles={"imbalance_leak": False, "freq_response": False, "dip": False}))
        s = fe.psd(1e6, 1.0)
        sn.append(float(s["shot"] - s["dark"]))
        cl.append(fe.clearance(1e6))
    np.testing.assert_allclose(np.array(sn) / p_lo, sn[0] / p_lo[0], rtol=1e-12)
    np.testing.assert_allclose(np.array(cl) / p_lo, cl[0] / p_lo[0], rtol=1e-12)
    assert sn[0] == pytest.approx(1e3**2 * 2 * 1.602176634e-19 * 1e-3, rel=1e-12)  # G² 2eℜP_LO


@pytest.mark.parametrize("mode", ["integrating", "pulse_resolved", "gated"])
def test_t11_integrating_dilution(mode):
    """K = 32, η = 1, one slot at −10 dB, LO over the whole train: integrating BHD −0.124 dB; per-pulse
    (resolved / gated, no crosstalk) −10.00 dB (T11 §3.7, §7.7 test 3; D45 comparison configuration)."""
    r = np.ones(32)
    r[0] = 0.1
    ps = mk(T.pulse_sample(mode=mode, n_frames=20000, toggles={"crosstalk": False, "electronic_noise": False}))
    out = ps.evaluate(r, c=[1.0], clearance_lin=math.inf, rng=stream(8, "vt42.t11"))
    want = -0.124 if mode == "integrating" else -10.00
    assert out["expected_db"][0] == pytest.approx(want, abs=1e-3)
    assert out["sample_db"][0] == pytest.approx(want, abs=0.15 if mode == "integrating" else 0.25)
