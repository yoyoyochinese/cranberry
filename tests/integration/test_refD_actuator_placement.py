"""WP-I2 integration check: B4 actuator drivers × F6 placement weights give exactly the intended δψ.

B4 (``blocks.ctrl._kit.DrvTarget``) writes ``drv.<seg>.<driver>`` with shape (n_sub,) as an intensive
density δψ·N_sub/(hi − lo) on ``sub = [lo, hi)`` (D61). F6
(``path.segments``) reads a driver with its placement weights: scalar placements (``uniform`` 1/N on
every subsegment, ``actuator`` / ``{point: k}`` weight 1) read a 1-element slot; ``per_sub`` reads
drv[m]/N_sub (an intensive per-subsegment value such as strain or temperature).

Finding (WP-I2): with the earlier "sum = δψ" convention, a split actuator on a ``per_sub`` placement
gave δψ/N_sub (a factor-N error). Fixed at the root in ``DrvTarget`` (D61); both the scalar placement
used by the Ref D preset (``n_sub: 1, sub: [0, 1]``) and ``per_sub`` placements now give exact δψ.
"""

import pytest

from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.path.testing import beam_spec, coef_spec, seg_spec, step_spec, topo_spec

pytestmark = pytest.mark.numba

FS = 1.0e6
K = 2.0  # rad/V (direct gain)
U = 0.75  # V (static command)
N_SUB = 4


def _pipeline(placement, act_n_sub, act_sub, coef=1.0, actuator_sub=None):
    seg = seg_spec("a1", 8.4, n_sub=N_SUB, drivers={"stretcher": placement},
                   spatial="distributed" if placement == "per_sub" else "lumped", actuator_sub=actuator_sub)
    topo = topo_spec([seg], [beam_spec("pr", [step_spec("a1")], "det")], ["det"],
                     coef_spec("a1", "stretcher", "pr", coef, 0.0))
    return [
        {"block": "ramp", "name": "cmd", "params": {"out": "v.u", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                    "period_s": 1.0, "offset_v": U, "toggles": {"enabled": True}}},
        {"block": "pzt_fiber_stretcher", "name": "str", "params": {
            "in": "v.u", "gain_mode": "direct", "k_rad_per_v": K, "polarity": 1.0, "modes": [],
            "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -5.0, "drive_max_v": 5.0,
            "hysteresis_frac": 0.0, "segment": "a1", "driver": "stretcher", "n_sub": act_n_sub, "sub": act_sub,
            "ring_len": 0, "insertion_loss_lin": None,
            "toggles": {"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": False}}},
        {"block": "segment_phase", "name": "sp", "params": {
            "topology": topo, "segments": "all", "toggles": {"phase": True, "group_delay": True,
                                                            "pm_axis_drift": False}}},
    ]


def _psi(obj, compiled):
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(obj), lay, 1)
    st, pr = State(lay), Params(lay)
    (run_compiled if compiled else run_reference)(blk, st, pr, 50)  # the segment transit (≈ 41 ns) is long settled
    return float(st["psi.a1.pr"]), float(st["act.str.dpsi"])


@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("placement,actuator_sub", [("uniform", None), ("actuator", 2), ({"point": 1}, None)])
def test_scalar_placement_gives_exact_dpsi(placement, actuator_sub, compiled):
    psi, dpsi = _psi(_pipeline(placement, 1, [0, 1], actuator_sub=actuator_sub), compiled)
    assert dpsi == pytest.approx(K * U, rel=1e-15)
    assert psi == pytest.approx(K * U, rel=1e-12)  # δψ on the beam = actuator δψ, no factor N


@pytest.mark.parametrize("compiled", [False, True])
def test_per_sub_placement_gives_exact_dpsi(compiled):
    # D61: actuators write intensive densities δψ·N/(hi−lo); F6 per_sub reads drv[m]/N → exact δψ
    psi, dpsi = _psi(_pipeline("per_sub", N_SUB, [0, N_SUB]), compiled)
    assert dpsi == pytest.approx(K * U, rel=1e-15)
    assert psi == pytest.approx(K * U, rel=1e-12)  # was δψ/N_sub before D61 (factor-N bug, fixed in _kit)
    psi_part, _ = _psi(_pipeline("per_sub", N_SUB, [1, 3]), compiled)
    assert psi_part == pytest.approx(K * U, rel=1e-12)  # a partial actuator range also gives δψ


def test_refD_preset_uses_the_scalar_convention():
    from cranberry import config

    spec = config.load("presets/reference/refD_nominal.yaml")
    seg = next(s for s in spec.topology["segments"] if s["id"] == "a1")
    assert seg["drivers"]["stretcher"] in ("uniform", "actuator")
    act = next(b for b in spec.pipeline.children if getattr(b, "name", "") == "stretcher915b")
    assert act.params["n_sub"] == 1 and list(act.params["sub"]) == [0, 1]
    coef = [c for c in spec.topology["coef"] if c["segment"] == "a1" and c["driver"] == "stretcher"
            and c["effect"] == "phase"]
    assert len(coef) == 1 and coef[0]["value"] == 1.0
