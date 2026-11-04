"""Phase 5 WS-01: the OPO-mode overlays (``stray/opo_mode_*``, ``ref/opo_mode_*``) point winstray1 / frefstray1 at the
preset's ``pulsed_opo_source`` mode, the same spec as stray1. winstray1's ungated ("on") value then equals stray1's
on the OPO presets, also with a frequency-offset probe (``probe/freq_80MHz``, where the Gaussian stand-in is ~10³× low).
Specs are loaded at a fixed point and the blocks built directly (no time loop)."""

import pytest
import yaml

from cranberry import config
from cranberry.blocks.fib.scattering import stray_photons
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.config.overlay import REPO_ROOT

pytestmark = [pytest.mark.integration]

POINT = dict(win_period_s=1.0e-4, win_duty_lock=0.65, win_edge_s=1.0e-7, win_eps_off=1.0e-5, win_limit_statistic="mean",
             win_refl_r=0.0, fref_p_peak_w=1.0e-6, fref_t_pulse_s=1.0e-9, fref_eps_mod=1.0e-5, fref_eps_gate=1.0,
             fref_edge_tau_s=1.0e-9, det_p_afterpulse=1.0e-3, det_tau_afterpulse_s=1.0e-7, det_recovery_s=4.3e-8)
SHORT = {"F_bor6": "bor6", "F_liu20": "liu20"}


def load(tmp_path, base, overlays):
    pt = tmp_path / "point.yaml"
    pt.write_text(yaml.safe_dump({"offline": {"params": POINT}}))
    spec = config.load(REPO_ROOT / f"presets/future/{base}_nominal.yaml", overlays=[*overlays, str(pt)])
    assert not config.sweep_markers(spec)
    return spec


def blocks(spec):
    return {b["name"]: b for b in spec.offline["blocks"]}


def build(b):
    return make_block(BlockSpec(b["block"], b["name"], b["params"]))


@pytest.mark.parametrize("base", ["F_bor6", "F_liu20"])
@pytest.mark.parametrize("probe", [(), ("probe/freq_80MHz",)])
def test_winstray_on_value_equals_stray1(tmp_path, base, probe):
    spec = load(tmp_path, base, [*probe, "stray/window_tail", f"stray/opo_mode_{SHORT[base]}"])
    bl = blocks(spec)
    assert bl["winstray1"]["params"]["mode"] == bl["stray1"]["params"]["mode"]
    win, s1 = build(bl["winstray1"]), build(bl["stray1"])
    on = win.check()["photons_per_mode_on"]
    ref = stray_photons(s1.params, s1._toggles, float(spec.offline["params"]["p_pr_avg_w"]))["photons_per_mode"]
    assert on == pytest.approx(ref, rel=1e-9)


@pytest.mark.parametrize("base", ["F_bor6", "F_liu20"])
@pytest.mark.parametrize("ref_overlay", ["ref/aurora_calframe", "ref/aurora_guardslot"])
def test_frefstray_mode_equals_stray1(tmp_path, base, ref_overlay):
    bl = blocks(load(tmp_path, base, [ref_overlay, f"ref/opo_mode_{SHORT[base]}"]))
    assert bl["frefstray1"]["params"]["mode"] == bl["stray1"]["params"]["mode"]
