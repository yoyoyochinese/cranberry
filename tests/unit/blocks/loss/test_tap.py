"""tap block (WP-B1): registered loss, P_tap output, toggle off ⇒ identity, compiled = reference."""

import pytest

from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.offline.loss_registry import LossRegistry


def spec(enabled=True, r_t=3e-3):
    return [{"block": "tap", "name": "tap1", "params": {
        "tap_ratio": r_t, "p_probe_w": 10e-6, "out": "v.tap1.p", "segment": "sqz.tap", "applies_to": ["gbs", "bhd"],
        "source": "[T12 §6.2]", "toggles": {"enabled": enabled}}}]


def make(enabled=True):
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(spec(enabled)), lay, 3)
    return blk, State(lay), Params(lay)


def registry(blk):
    (leaf,) = iter_leaves(blk)
    return LossRegistry(leaf.losses(), paths={"sqz": {"segments": ["sqz.tap"]}})


def test_loss_and_output():
    blk, s, p = make()
    run_reference(blk, s, p, 3)
    assert float(s["v.tap1.p"]) == pytest.approx(3e-3 * 10e-6)
    assert registry(blk).eta("sqz", "bhd") == pytest.approx(1 - 3e-3)
    p["tap1.p_probe_w"] = 20e-6
    run_reference(blk, s, p, 1)
    assert float(s["v.tap1.p"]) == pytest.approx(3e-3 * 20e-6)


def test_toggle_off_is_identity():
    blk, s, p = make(enabled=False)
    run_reference(blk, s, p, 3)
    assert float(s["v.tap1.p"]) == 0.0
    reg = registry(blk)
    assert reg.eta("sqz", "gbs") == 1.0 and len(reg.breakdown("sqz", "gbs")) == 1  # kept in the record


@pytest.mark.numba
def test_compiled_equals_reference():
    (b1, s1, p1), (b2, s2, p2) = make(), make()
    run_reference(b1, s1, p1, 100)
    run_compiled(b2, s2, p2, 100)
    assert float(s1["v.tap1.p"]) == float(s2["v.tap1.p"]) and s1.n == s2.n == 100


def test_bad_ratio():
    with pytest.raises(ValueError):
        build(spec_from_obj(spec(r_t=1.0)), SlotLayout(1e6), 0)
