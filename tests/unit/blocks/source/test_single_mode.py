"""single_mode_reference block contract (WP-B1). Numbers: tests/vt/test_vt_01_*.py, test_vt_04_*.py."""

import numpy as np
import pytest

from cranberry.blocks.source.single_mode import make_single_mode_reference
from cranberry.core.block import Block, BlockSpec, MissingParameterError
from cranberry.core.registry import make_block

P = dict(loss=0.08, alpha_per_w=9.06, pump_w=0.64, theta_rms_rad=0.009, phi_opa0_rad=0.0, average="exact")


def test_factory_registry_and_contract():
    b = make_block(BlockSpec("single_mode_reference", "sqz", dict(P)))
    assert isinstance(b, Block) and b.kind == "single_mode_reference" and b.name == "sqz"
    assert b.reads == () and b.writes == () and b.toggles == {}
    assert b.losses() == [] and b.bands() == []  # L is inside R± (§7-1), never registered twice
    assert b(state := object(), None) is state  # offline block: identity in the time loop


def test_missing_and_bad_params():
    for k in P:
        with pytest.raises(MissingParameterError):
            make_single_mode_reference({x: v for x, v in P.items() if x != k})
    with pytest.raises(ValueError):
        make_single_mode_reference(dict(P, average="mean"))


def test_frame_rule_4():
    for phi0 in (0.0, 0.7, -2.0):
        so = make_single_mode_reference(dict(P, phi_opa0_rad=phi0)).source_output
        vm, vp, th_off = so.project(None)
        assert th_off == pytest.approx(0.0, abs=1e-12)
        assert so.theta_sq0_rad == pytest.approx(phi0 / 2 + np.pi / 2)
        r = np.sqrt(9.06 * 0.64)
        assert 2 * vm == pytest.approx(0.08 + 0.92 * np.exp(-2 * r), rel=1e-12)
        assert 2 * vp == pytest.approx(0.08 + 0.92 * np.exp(2 * r), rel=1e-12)
