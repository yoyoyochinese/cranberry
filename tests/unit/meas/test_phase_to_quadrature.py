"""Unit tests of ``blocks.meas.phase_to_quadrature`` (conv §6, T04 §3.1, §13, G05)."""

import math

import numpy as np
import pytest

from cranberry.blocks.meas.phase_to_quadrature import PhaseStats, phase_to_quadrature, quadrature_for_metric
from cranberry.core.units import V_VAC
from cranberry.offline import gaussian
from cranberry.offline.band_ledger import BandDecl, BandLedger
from cranberry.offline.loss_registry import LossElement, LossRegistry
from cranberry.offline.source_output import SourceOutput

VM, VP = 0.0875 * V_VAC, 113.7 * V_VAC


def gauss(sigma, n=200_001, mean=0.0):
    from scipy.special import ndtri

    x = ndtri((np.arange(n) + 0.5) / n)
    return mean + sigma * x / np.sqrt(np.mean(x * x))  # exact sample rms σ


def test_phase_stats_series_acc_and_shift():
    d = gauss(0.05, 1001, mean=0.01)
    st = PhaseStats.from_series(d)
    acc = [d.size, d.sum(), (d * d).sum(), np.cos(2 * d).sum(), np.sin(2 * d).sum()]
    st2 = PhaseStats.from_acc(acc)
    for f in ("mean", "mean_sq", "mean_cos2", "mean_sin2"):
        assert getattr(st, f) == pytest.approx(getattr(st2, f), rel=1e-12)
    sh, ref = st.shifted(0.3), PhaseStats.from_series(d - 0.3)
    for f in ("mean", "mean_sq", "mean_cos2", "mean_sin2"):
        assert getattr(sh, f) == pytest.approx(getattr(ref, f), rel=1e-12, abs=1e-14)
    assert st.std == pytest.approx(np.std(d), rel=1e-12)
    with pytest.raises(ValueError):
        PhaseStats.from_acc([0, 0, 0, 0, 0])


def test_exact_matches_jitter_average_and_gaussian_chi():
    sigma = 0.03
    out = phase_to_quadrature(gauss(sigma), VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0)
    chi = out["chi_loop"]
    assert chi == pytest.approx(gaussian.chi_gaussian(sigma), rel=1e-6)
    assert out["r_minus"] == pytest.approx(gaussian.jitter_average(VM, VP, chi, delta0_rad=0.0), rel=1e-12)
    assert out["r_plus"] == pytest.approx(gaussian.jitter_average(VM, VP, chi, delta0_rad=math.pi / 2), rel=1e-12)
    assert out["v_prime"] == out["r_minus"]
    assert out["sigma_loop_rad"] == pytest.approx(sigma, rel=1e-12)
    assert np.mean(out["v_prime_t"]) == pytest.approx(out["v_prime"], rel=1e-12)  # V′(t) averages to V′


def test_setpoint_offset_and_theta_off_rotate_the_same_way():
    d = gauss(0.02, 2001)
    a = phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=math.pi / 2)
    assert a["v_prime"] == pytest.approx(a["r_plus"], rel=1e-12)
    b = phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0, theta_off_rad=0.1)
    c = phase_to_quadrature(d - 0.1, VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0)
    assert b["r_minus"] == pytest.approx(c["r_minus"], rel=1e-12)
    assert b["mean_loop_rad"] == pytest.approx(-0.1, abs=1e-12)
    # static offset Δ0 enters as cos 2Δ0 (conv §6)
    z = phase_to_quadrature(np.zeros(4), VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.2)
    assert z["v_prime"] == pytest.approx(gaussian.jitter_average(VM, VP, 1.0, delta0_rad=0.2), rel=1e-12)


def test_hf_and_dither_multiply_chi():
    d = gauss(0.01, 2001)
    base = phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0)
    out = phase_to_quadrature(d, VM, VP, hf_variance_rad2=4e-4, average="exact", setpoint_rad=0.0, dither_m_rad=[0.05])
    assert out["chi_hf"] == pytest.approx(math.exp(-8e-4))
    assert out["chi_dither"] == pytest.approx(float(gaussian.chi_sinusoid(0.05)))
    assert out["chi"] == pytest.approx(base["chi_loop"] * out["chi_hf"] * out["chi_dither"], rel=1e-12)
    assert out["sigma_total_rad"] == pytest.approx(math.sqrt(1e-4 + 4e-4 + 0.05**2 / 2), rel=1e-9)
    assert out["r_minus_db"] > base["r_minus_db"]


def test_small_angle_and_errors():
    d = gauss(0.009, 2001)
    sa = phase_to_quadrature(d, VM, VP, hf_variance_rad2=1e-6, average="small_angle", setpoint_rad=math.pi / 2)
    s = math.sqrt(0.009**2 + 1e-6)
    assert sa["r_minus"] == pytest.approx(gaussian.jitter_average_small_angle(VM, VP, s), rel=1e-9)
    assert sa["v_prime"] == sa["r_plus"]
    with pytest.raises(ValueError, match="small_angle"):
        phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="small_angle", setpoint_rad=0.3)
    with pytest.raises(ValueError):
        phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="median", setpoint_rad=0.0)
    with pytest.raises(ValueError):
        phase_to_quadrature(d, VM, VP, hf_variance_rad2=-1.0, average="exact", setpoint_rad=0.0)


def test_cw_sideband_arrays_broadcast():
    grid = np.array([1e6, 2e6, 3e6])
    vm, vp = np.array([0.1, 0.2, 0.3]) * V_VAC, np.array([10.0, 8.0, 6.0]) * V_VAC
    out = phase_to_quadrature(PhaseStats.from_series(gauss(0.02, 2001)), vm, vp, hf_variance_rad2=0.0,
                              average="exact", setpoint_rad=0.0, theta_off_rad=np.array([0.0, 0.01, 0.02]))
    assert np.shape(out["r_minus_db"]) == grid.shape and "v_prime_t" not in out
    assert out["r_minus"][1] < out["r_minus"][2]


def test_quadrature_for_metric_uses_registry_eta_and_ledger_once():
    so = SourceOutput.single_mode_reference(loss=0.0, alpha_per_w=9.06, pump_w=0.64, phi_opa0_rad=0.0)
    reg = LossRegistry(paths={"sqz": {"segments": ["s"], "observers": ["bhd"]}})
    reg.add(LossElement(id="pd", eta=0.9, kind="detector", applies_to={"bhd"}, source="test", enabled=True,
                        observer="bhd"))
    led = BandLedger(beams={"sqz": ["s"], "lo": ["l"]})
    led.define_metric("i.port0", f_max_hz=1e6, weight=None)
    led.declare(BandDecl(source="n", quantity="phase", targets=("l",), f_lo_hz=1.0, f_hi_hz=1e3, unit="rad",
                         variance_outside=2e-4))
    d = gauss(0.01, 2001)
    out = quadrature_for_metric(d, metric="i.port0", source=so, lo_mode=None, loss_registry=reg, path="sqz",
                                target="bhd", band_ledger=led, beams={"lo": 1, "sqz": -1}, average="exact",
                                setpoint_rad=0.0, include_hf_phase_noise=True)
    vm, vp, _ = so.project(None)
    assert out["eta"] == pytest.approx(0.9)
    assert out["v_minus"] == pytest.approx(0.9 * vm + 0.1 * V_VAC)
    assert out["hf_variance_rad2"] == pytest.approx(2e-4)
    gbs = quadrature_for_metric(d, metric="i.port0", source=so, lo_mode=None, loss_registry=reg, path="sqz",
                                target="gbs", band_ledger=led, beams={"lo": 1, "sqz": -1}, average="exact",
                                setpoint_rad=0.0, include_hf_phase_noise=False)
    assert gbs["eta"] == 1.0 and gbs["hf_variance_rad2"] == 0.0  # BHD-only element, HF switched off
    with pytest.raises(ValueError, match="BandLedger"):
        quadrature_for_metric(d, metric="i.port0", source=so, lo_mode=None, loss_registry=None, path=None,
                              target="bhd", band_ledger=None, beams={}, average="exact", setpoint_rad=0.0,
                              include_hf_phase_noise=True)
