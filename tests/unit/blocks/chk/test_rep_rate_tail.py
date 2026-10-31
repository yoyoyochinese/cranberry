"""rep_rate_tail_check (WP-B11, T17 §5.3, D50): repetition-rate limits, warning, source model, toggle."""

import warnings

import numpy as np
import pytest

from cranberry.blocks.source.pulsed_testing import BOREALIS, TOGGLES_ON
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block


def chk(f_rep, model="analytic", on=True, **kw):
    p = dict(rep_rate_hz=f_rep, threshold=1e-3, tail_model=model, toggles={"enabled": on})
    p.update(kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_block(BlockSpec("rep_rate_tail_check", "chk", p))


def test_limits_borealis_and_liu():
    b = chk(6e6, fwhm_hz=66e6)
    assert b.result["ok"] and not b.chk_log.events
    assert b.result["f_rep_max_hz"] == pytest.approx(45e6, rel=0.02)  # T17 §5.3: ≲ 45 MHz
    assert chk(20e6, fwhm_hz=7.4e9 / 59).result["f_rep_max_hz"] == pytest.approx(87e6, rel=0.02)  # ≲ 87 MHz


def test_bad_bor_100mhz_warns_but_continues():
    b = chk(100e6, fwhm_hz=66e6)
    assert not b.result["ok"] and b.result["tail"] == pytest.approx(0.138, rel=0.05)  # 13.8 % (T17 §5.3)
    (ev,) = b.chk_log.events
    assert ev.code == "rep_rate_tail" and ev.context["tail"] == b.result["tail"]
    liu = chk(100e6, fwhm_hz=7.4e9 / 59)
    assert liu.result["tail"] == pytest.approx(3.4e-3, rel=0.05)  # 0.34 %


def test_source_model_and_toggle_off():
    src = {"block": "pulsed_opo_source", "name": "opo", "params": dict(BOREALIS, x_ss=2.4, toggles=TOGGLES_ON)}
    b = chk(100e6, model="source", source=src)
    assert b.result["tail"] == pytest.approx(0.138, rel=0.10) and b.chk_log.events
    off = chk(100e6, on=False, fwhm_hz=66e6)
    assert off.result is None and not off.chk_log.events
    assert np.isfinite(chk(1e9, fwhm_hz=66e6).result["tail"])
