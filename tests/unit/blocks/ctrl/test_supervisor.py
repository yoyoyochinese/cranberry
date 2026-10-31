"""``lock_supervisor``: ramp re-acquisition on a fringe, reset to centre, error-threshold loss, fringe-slip
counting, toggle off = identity (T16 §3.6, V16)."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.block import BaseBlock
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

FS = 1.0e5


class Fringe(BaseBlock):
    """Test plant (Python only): v = sin(φ0 + φ_a) with φ0 from a slot, δ = φ0 + φ_a (true phase)."""

    def __init__(self, params):
        super().__init__("test_fringe", "fringe", params, reads=("v.phi0", "act.a.dpsi"), writes=("v.e", "v.delta"))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add("v.e", (), np.float64, "V")
        layout.add("v.delta", (), np.float64, "rad")

    def __call__(self, state, params):
        d = float(state["v.phi0"]) + float(state["act.a.dpsi"])
        state["v.delta"], state["v.e"] = d, math.sin(d)
        return state


def const(name, out, value):
    return {"block": "ramp", "name": name, "params": {"out": out, "waveform": "sawtooth", "amplitude_v": 0.0,
                                                      "period_s": 1.0, "offset_v": value, "toggles": {"enabled": True}}}


def sup(relock, on=True, true_in="v.delta", **kw):
    p = {"u_in": "ctrl.pid.u", "e_in": "v.e", "true_in": true_in, "out": "ctrl.sup.u", "relock": relock,
         "u_min_v": -4.0, "u_max_v": 4.0, "t_sat_s": 0.0, "e_th": 0.3, "t_e_s": 2e-4, "toggles": {"relock": on}}
    p.update({"reset_center": {"t_hop_s": 1e-3}, "fringe_hop": {"t_hop_s": 1e-3, "u_2pi_v": 2 * math.pi},
              "ramp_search": {"ramp_period_s": 0.02, "e_cap": 0.05, "slope_sign": 1.0}}.get(relock, {}))
    p.update(kw)
    return {"block": "lock_supervisor", "name": "sup", "params": p}


def loop(phi0, supervisor, ki=-2e3):
    """φ0 → fringe v = sin(φ0 + K u) → pid (I) → supervisor → delay → static actuator K = 1 rad/V.
    Stable lock points: cos(δ) > 0 (ki < 0: u' = ki·v)."""
    pid = {"in": "v.e", "out": "ctrl.pid.u", "device": None, "discretisation": "dsc1", "anti_windup": "none",
           "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
           "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1e9, "u_max_v": 1e9,
           "toggles": {"feedback": True, "p_term": False, "i_term": True, "d_term": False}}
    act = {"in": "delay.dl.out", "gain_mode": "direct", "k_rad_per_v": 1.0, "polarity": 1.0, "modes": [],
           "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -4.0, "drive_max_v": 4.0,
           "hysteresis_frac": 0.0, "segment": None, "insertion_loss_db": None,
           "toggles": {"resonance": False, "hysteresis": False, "saturation": True, "insertion_loss": False}}
    obj = [phi0, {"block": "pid", "name": "pid", "params": pid}, supervisor,
           {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.sup.u", "out": "delay.dl.out",
                                                              "tau_s": 0.0, "method": "linear"}},
           {"block": "pzt_fiber_stretcher", "name": "a", "params": act}]
    lay = SlotLayout(FS)
    blk = build(_with_fringe(obj), lay, 9)
    return blk, State(lay), Params(lay)


def _with_fringe(obj):
    from cranberry.core.pipeline import Seq, from_spec

    blocks = [from_spec(spec_from_obj(o)) for o in obj]
    return Seq([blocks[0], Fringe({}), *blocks[1:]])


def test_ramp_search_reacquires_on_the_correct_slope():
    """Start on the wrong side of the fringe (|e| > e_th): loss → ramp over the range → engage at a zero
    crossing with positive slope (stable) → locked with |e| → 0 and cos δ > 0."""
    blk, st, pr = loop(const("p0", "v.phi0", 2.5), sup("ramp_search"))
    run_reference(blk, st, pr, int(0.05 * FS))
    assert int(st["lock.sup.n_loss"]) == 1 and int(st["lock.sup.n_relock"]) == 1
    assert int(st["lock.sup.state"]) == 1
    assert abs(float(st["v.e"])) < 1e-6 and math.cos(float(st["v.delta"])) > 0.99
    assert -4.0 <= float(st["ctrl.sup.u"]) <= 4.0


def test_reset_center_and_error_threshold():
    """Loss by |e| > e_th for t_e: the output jumps to the range centre and the state holds for t_hop;
    with the threshold then disabled (runtime parameter) the loop re-locks by itself."""
    blk, st, pr = loop(const("p0", "v.phi0", 2.5), sup("reset_center"))
    run_reference(blk, st, pr, 20)  # |e| = sin 2.5 > e_th from n = 0: lost when 20 samples · h ≥ t_e
    assert int(st["lock.sup.n_loss"]) == 1 and float(st["lock.sup.t_loss"]) == pytest.approx(19 / FS, rel=1e-12)
    assert float(st["ctrl.sup.u"]) == 0.0 and int(st["lock.sup.state"]) == 2 and int(st["lock.sup.n_relock"]) == 1
    pr["sup.e_th"] = float("inf")
    run_reference(blk, st, pr, int(0.02 * FS))
    assert int(st["lock.sup.state"]) == 1 and abs(float(st["v.e"])) < 1e-6 and int(st["lock.sup.n_loss"]) == 1


def test_fringe_slip_count_and_toggle_off_identity():
    """Without relock the output is the controller output (identity); the unlocked loop slides over
    fringes (δ = φ0 + K u keeps growing while the actuator is at its range limit) and slips are counted
    from the true phase."""
    drift = {"block": "ramp", "name": "p0", "params": {"out": "v.phi0", "waveform": "sawtooth",
                                                       "amplitude_v": 50.0, "period_s": 2.0, "offset_v": 50.0,
                                                       "toggles": {"enabled": True}}}
    blk, st, pr = loop(drift, sup("ramp_search", on=False, e_th=float("inf")))
    us = []
    for _ in range(int(0.3 * FS) // 100):
        run_reference(blk, st, pr, 100)
        us.append((float(st["ctrl.pid.u"]), float(st["ctrl.sup.u"])))
    assert all(a == b for a, b in us)  # relock off → out = u_in
    assert int(st["lock.sup.n_loss"]) == 1  # saturation at the ±4 V range
    k_end = math.floor(float(st["v.delta"]) / (2 * math.pi) + 0.5)
    assert int(st["lock.sup.n_slip"]) == abs(k_end) >= 1


def test_requires_mode_params():
    p = sup("fringe_hop")
    del p["params"]["u_2pi_v"]
    with pytest.raises(KeyError):
        loop(const("p0", "v.phi0", 0.0), p)
    with pytest.raises(ValueError):
        loop(const("p0", "v.phi0", 0.0), sup("magic"))
