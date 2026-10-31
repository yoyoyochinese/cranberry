"""fsr_comb_filter (WP-B11): FBG pair / MZI cascade transmission, n_FSR, registered passband loss, toggle."""

import numpy as np
import pytest

from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline.loss_registry import LossRegistry

FBG = [{"type": "flat_top", "center_hz": c, "bandwidth_hz": 5e9, "order": 4, "floor": 1e-3} for c in (-0.5e9, 0.5e9)]


def fsr(stages, on=True, eta=0.8, r_m=0.886):
    return make_block(BlockSpec("fsr_comb_filter", "fsr", dict(
        fsr_hz=3.7e9, stages=stages, r_m=r_m, m_max=10, insertion_eta=eta, segment="opo.out", applies_to=["gbs", "bhd"],
        source_tag="[M-SI p.19]", toggles={"fsr_modes": on})))


def test_fbg_pair_borealis_like():
    b = fsr(FBG)
    t = b.transmission(np.array([-2, -1, 0, 1, 2]) * 3.7e9)
    assert t[2] > 0.999
    assert np.all(10 * np.log10(t[[0, 1, 3, 4]]) < -30)  # ≥ 30 dB on ±1, ±2 [M-SI p.19]
    m = np.arange(1, 11)
    expect = 0.8 * np.sum((b.transmission(m * 3.7e9) + b.transmission(-m * 3.7e9)) * np.sinh(0.886) ** 2)
    assert b.n_fsr == pytest.approx(expect, rel=1e-12) and b.n_fsr < 1e-4


def test_mzi_stage_formula():
    st = {"type": "mzi", "period_hz": 8 * 7.4e9, "offset_hz": 0.0, "extinction": 1e-4, "eta": 0.998}
    b = fsr([st])
    assert b.transmission(0.0) == pytest.approx(0.998, rel=1e-12)
    assert b.transmission(4 * 7.4e9) == pytest.approx(0.998 * 1e-4, rel=1e-9)  # half period: extinction floor


def test_passband_loss_registered_and_toggle():
    b = fsr(FBG)
    reg = LossRegistry(b.losses(), paths={"sqz": {"segments": ["opo.out"]}})
    assert reg.eta("sqz", "gbs") == pytest.approx(0.8 * float(b.transmission(0.0)), rel=1e-12)
    off = fsr(FBG, on=False)
    assert off.n_fsr == 0.0 and off.eta_passband == b.eta_passband  # ideal removal, loss kept


def test_r_m_callable_and_bad_stage():
    b = fsr([], eta=1.0, r_m=lambda m: 0.5 / m)
    assert b.n_fsr == pytest.approx(np.sum(2 * np.sinh(0.5 / np.arange(1, 11)) ** 2), rel=1e-12)
    with pytest.raises(ValueError):
        fsr([{"type": "etalon"}])
