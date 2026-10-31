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
