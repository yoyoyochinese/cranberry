"""lo_temporal_shape (WP-B11; pulsed_bhd_kernel is owned by WP-B5, blocks/det/pulsed.py): T17 §6.3 LO overlaps, registered losses, toggles;
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


def lo(shape, base=BOREALIS, fit=True, on=True, carrier="registry", src=None, **kw):
    src = src or {"block": "pulsed_opo_source", "name": "opo", "params": dict(base, x_ss=2.4, toggles=TOGGLES_ON)}
    return make_block(BlockSpec("lo_temporal_shape", "lo", dict(
        source=src, shape=shape, fit=fit, observer="bhd0", source_tag="[T17 §6.3]", eta_t_carrier=carrier,
        toggles={"shape_mismatch": on}, **kw)))


@pytest.mark.parametrize("base", [BOREALIS, LIU], ids=["borealis", "liu"])
def test_t17_lo_overlaps(base):
    kc = np.pi * (base.get("fwhm_hz") or base["fsr_hz"] / base["finesse"])
    # T17 §6.3 [sim] as corrected by D60-5 / D70-2: Gaussian 0.955, rectangular 0.866 (T17: 0.951 / 0.864)
    assert lo("gauss", base).eta_t == pytest.approx(0.955, abs=0.002)
    rect = lo("rect", base)
    assert rect.eta_t == pytest.approx(0.866, abs=0.002)
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
    vm, _, _ = lo("gauss", carrier="projection").projection
    vopt = 0.5 * lo("optimal").source.project(lo("optimal").u_opt)[0]
    assert 2 * vm == pytest.approx(b.eta_t * 4 * vopt + (1 - b.eta_t), abs=0.01)


def test_eta_t_single_source():
    """M5-09: η_t is carried either by the projection (element disabled) or by the registry element
    (projection η_t-free), never both; for a single-mode source the two carriers give the same V."""
    g = {"gaussian": {"T0_s": 1.0 / KC, "r": 1.0, "n_points": 256, "window_s": 16.0 / KC, "phi_opa0_rad": 0.0}}
    kw = dict(fit=False, center_s=0.3 / KC, t0_s=1.4 / KC, src=g)
    proj, regi = lo("gauss", carrier="projection", **kw), lo("gauss", carrier="registry", **kw)
    eta = proj.eta_t
    assert 0.5 < eta < 0.99 and regi.eta_t == eta
    paths = {"sqz": {"segments": [], "observers": ["bhd0"]}}
    assert LossRegistry(proj.losses(), paths=paths).eta("sqz", "bhd") == 1.0  # disabled
    assert LossRegistry(regi.losses(), paths=paths).eta("sqz", "bhd") == pytest.approx(eta, rel=1e-12)
    vm_r, vp_r, _ = regi.projection  # η_t-free: the u_opt projection
    assert [vm_r, vp_r] == pytest.approx([float(v) for v in proj.source.project(proj.u_opt)[:2]], rel=1e-12)
    vm_p, vp_p, _ = proj.projection
    assert vm_p == pytest.approx(eta * vm_r + 0.5 * (1 - eta), rel=1e-9)
    assert vp_p == pytest.approx(eta * vp_r + 0.5 * (1 - eta), rel=1e-9)
    with pytest.raises(ValueError):
        lo("gauss", carrier="both")
