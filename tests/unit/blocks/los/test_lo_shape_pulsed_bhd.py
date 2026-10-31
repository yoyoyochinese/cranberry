"""lo_temporal_shape and pulsed_bhd_kernel (WP-B11): T17 §6.3 LO overlaps, registered losses, toggles;
T13 §1.8 crosstalk, T11 §7.7 integrating dilution, D46 gate check."""

import warnings

import numpy as np
import pytest

from cranberry.blocks.source.pulsed_testing import BOREALIS, LIU, TOGGLES_ON
from cranberry.core.block import BlockSpec
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import make_block
from cranberry.offline.loss_registry import LossRegistry

KC = np.pi * BOREALIS["fwhm_hz"]


def lo(shape, base=BOREALIS, fit=True, on=True, **kw):
    src = {"block": "pulsed_opo_source", "name": "opo", "params": dict(base, x_ss=2.4, toggles=TOGGLES_ON)}
    return make_block(BlockSpec("lo_temporal_shape", "lo", dict(
        source=src, shape=shape, fit=fit, observer="bhd0", source_tag="[T17 §6.3]", toggles={"shape_mismatch": on}, **kw)))


@pytest.mark.parametrize("base", [BOREALIS, LIU], ids=["borealis", "liu"])
def test_t17_lo_overlaps(base):
    kc = np.pi * (base.get("fwhm_hz") or base["fsr_hz"] / base["finesse"])
    # T17 §6.3 [sim]; the continuous optimum of the Gaussian is 0.955 (T17's coarser search: 0.951)
    assert lo("gauss", base).eta_t == pytest.approx(0.951, abs=0.005)
    rect = lo("rect", base)
    assert rect.eta_t == pytest.approx(0.864, abs=0.005)
    assert abs(rect.shape_params[1]) * kc == pytest.approx(2.1, abs=0.15)  # best width ≈ 2/κ_c
    assert lo("one_sided_exp", base, decay_rate_rad_s=kc).eta_t == pytest.approx(0.843, abs=0.005)


def test_given_params_and_cw_digital_limit():
    g = lo("gauss", fit=False, center_s=1.4 / KC, t0_s=0.79 / KC)
    assert g.shape_params == (1.4 / KC, 0.79 / KC) and g.eta_t == pytest.approx(0.955, abs=0.002)
    assert lo("cw_digital", bhd_bandwidth_hz=1e13).eta_t == pytest.approx(1.0, abs=1e-6)  # B → ∞: weight = u_opt
    slow = lo("cw_digital", bhd_bandwidth_hz=30e6).eta_t
    assert 0.5 < slow < lo("cw_digital", bhd_bandwidth_hz=300e6).eta_t < 1.0


def test_projection_matches_registry_for_matched_shape_and_toggle_off():
    b = lo("gauss")
    (el,) = b.losses()
    assert el.kind == "lo_shape" and el.applies_to == {"bhd"} and el.observer == "bhd0" and el.eta == pytest.approx(b.eta_t)
    off = lo("gauss", on=False)
    assert off.eta_t == pytest.approx(1.0, abs=1e-12) and off.projection[2] == pytest.approx(0.0, abs=1e-12)
    reg = LossRegistry(off.losses(), paths={"sqz": {"segments": [], "observers": ["bhd0"]}})
    assert reg.eta("sqz", "bhd") == 1.0  # disabled element (identity)
    # multimode source (K ≈ 1.1): projection ≈ single-mode registry reading
    vm, _, _ = b.projection
    vopt = 0.5 * lo("optimal").source.project(lo("optimal").u_opt)[0]
    assert 2 * vm == pytest.approx(b.eta_t * 4 * vopt + (1 - b.eta_t), abs=0.01)


def bhd(mode="pulse_resolved", resp="first_order", on=(True, True), **kw):
    p = dict(response=resp, bandwidth_hz=200e6, gate_s=25e-9, rep_rate_hz=20e6, m_max=5, n_lo_photons=1e8,
             sigma_el_electrons=1e3, quantum_efficiency=0.99, mode=mode, observer="bhd0", source_tag="[T13 §1.8]",
             toggles={"crosstalk": on[0], "electronic_noise": on[1]})
    p.update(kw)
    return make_block(BlockSpec("pulsed_bhd_kernel", "pbhd", p))


def test_crosstalk_first_order_and_adjacent_correlation():
    b = bhd()
    np.testing.assert_allclose(b.c_m, np.exp(-2 * np.pi * 200e6 * np.arange(6) / 20e6), rtol=1e-12)
    rng = np.random.default_rng(1)
    q = rng.standard_normal(400_000)
    b2 = bhd(bandwidth_hz=5e6)
    x = np.convolve(q, b2.c_m)[: q.size]
    rho = np.corrcoef(x[1:], x[:-1])[0, 1]
    assert rho == pytest.approx(b2.adjacent_correlation, abs=4 / np.sqrt(q.size))
    assert bhd(on=(False, True)).c_m.tolist() == [1, 0, 0, 0, 0, 0]
    assert bhd(bandwidth_hz=0.74 * 20e6).c_m[1] < 1e-2  # B/f_rep > 0.73 ⇒ c_1 < 1e-2 (T13 §1.8)


def test_second_order_reduces_to_first_order_shape_and_rings():
    over = bhd(resp="second_order", q=0.3, bandwidth_hz=10e6)
    assert np.all(over.c_m > 0) and np.all(np.diff(over.c_m) < 0)
    ring = bhd(resp="second_order", q=5.0, bandwidth_hz=10e6)
    assert np.any(ring.c_m[1:] < 0)  # ringing: sign changes (T13 §1.8)


def test_integrating_dilution_and_electronic_loss():
    b = bhd(mode="integrating", n_slots_integrated=32)
    r = 10 ** (-10 / 10)
    assert 10 * np.log10(b.measured_r(r)) == pytest.approx(-0.124, abs=1e-3)  # T11 §7.7
    assert bhd().measured_r(r) == r
    (el,) = bhd().losses()
    c = 0.99 * 1e8 / 1e6
    assert el.kind == "electronic" and el.eta == pytest.approx(c / (1 + c)) and el.enabled
    assert not bhd(on=(True, False)).losses()[0].enabled


def test_gate_shorter_than_mode_warns():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        short = bhd(mode="gated", gate_s=5e-9, mode_length_s=25e-9)
        ok = bhd(mode="gated", gate_s=30e-9, mode_length_s=25e-9)
    assert [e.code for e in short.chk_log.events] == ["gate_shorter_than_mode"] and not ok.chk_log.events
    assert any(issubclass(x.category, ChkWarning) for x in w)
