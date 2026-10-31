"""VT-21 (block part, WP-B4): lock duration set by the actuator range (prompt §7-10 second half;
notes/00_summary.md §4 VT-21; T16 V12, T15 §7 test 8).

(a) Linear drift r = 0.1 rad/s, range R = 25.96 rad (DSC1 ±4 V × 6.49 rad/V), start at the centre: the
    ``lock_supervisor`` reports the loss at t = R/r = 259.6 s ± 1 sample. The controller is PI + double
    integrator (type-2 loop: zero ramp error, so the actuator command equals the drift exactly).
(b) Random walk (Var = 2 D t): mean exit time of (−R, R) from the centre is R²/(2D) within 5 %. Each trial
    is a fresh State on the same build (compiled once; the noise stream continues, so trials are
    independent). Exit-time spread is √(2/3)·mean (Brownian exit from a symmetric interval), so
    N trials give a 3σ error of 3·0.816/√N; N = 3000 → 4.5 %. R/σ_step = 200 keeps the discrete
    overshoot bias ≈ 2·0.58/200 ≈ 0.6 %.

The thermal-ramp case (63 m, 1e−4 K/s, 64.9 rad → 235 s) needs the segment coefficients (WP-F6 /
WP-I5) and is not part of this file.
"""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.compile import compile_loop, run_compiled
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba


def chain(err_src, k_act, u_max, gains, relock="none"):
    """err_src → adder e = d − φ_a → pid → lock_supervisor → dac (no quantisation) → delay → actuator."""
    pid = {"in": "v.e", "out": "ctrl.pid.u", "device": None, "discretisation": "dsc1", "anti_windup": "none",
           "kp": 0.0, "ki_per_s": 0.0, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
           "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1e9, "u_max_v": 1e9,
           "toggles": {"feedback": True, "p_term": True, "i_term": True, "d_term": False}}
    pid.update(gains)
    return [
        err_src,
        {"block": "adder", "name": "err", "params": {"ins": ["v.d", "act.a.dpsi"], "gains": [1.0, -1.0],
                                                     "out": "v.e", "unit": "rad"}},
        {"block": "pid", "name": "pid", "params": pid},
        {"block": "lock_supervisor", "name": "sup", "params": {
            "u_in": "ctrl.pid.u", "e_in": "v.e", "true_in": None, "out": "ctrl.sup.u", "relock": relock,
            "u_min_v": -u_max, "u_max_v": u_max, "t_sat_s": 0.0, "e_th": float("inf"), "t_e_s": 0.0,
            "toggles": {"relock": False}}},
        {"block": "dac", "name": "dac", "params": {
            "in": "ctrl.sup.u", "out": "v.dac", "device": None, "bits": 12, "range_v": [-u_max, u_max],
            "noise_vrms": 0.0, "r_out_ohm": 220.0, "r_load_ohm": float("inf"),
            "toggles": {"quantization": False, "noise": False, "saturation": True}}},
        {"block": "delay_buffer", "name": "dl", "params": {"in": "v.dac", "out": "delay.dl.out", "tau_s": 0.0,
                                                           "method": "linear"}},
        {"block": "pzt_fiber_stretcher", "name": "a", "params": {
            "in": "delay.dl.out", "gain_mode": "direct", "k_rad_per_v": k_act, "polarity": 1.0, "modes": [],
            "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -u_max, "drive_max_v": u_max,
            "hysteresis_frac": 0.0, "segment": None, "insertion_loss_db": None,
            "toggles": {"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": False}}},
    ]


def test_vt21a_linear_drift_lock_time():
    fs, r, k_act, u_max = 1.0e3, 0.1, 6.49, 4.0  # R = 25.96 rad
    drift = {"block": "ramp", "name": "drift", "params": {"out": "v.d", "waveform": "sawtooth",
                                                          "amplitude_v": r * 500.0, "period_s": 1000.0,
                                                          "offset_v": r * 500.0, "toggles": {"enabled": True}}}
    obj = chain(drift, k_act, u_max, {"kp": 0.05, "ki_per_s": 2.0, "ki2_per_s2": 4.0})
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, 21)
    st, pr = State(lay), Params(lay)
    run_compiled(blk, st, pr, int(262 * fs))
    assert int(st["lock.sup.n_loss"]) == 1 and int(st["lock.sup.state"]) == 0
    assert float(st["lock.sup.t_loss"]) == pytest.approx(25.96 / r, abs=1.0 / fs)
    assert int(st["ctrl.dac.sat"]) == 1  # the actuator sits at its range after the loss


def test_vt21b_random_walk_mean_exit_time():
    fs, sigma, u_max = 1.0e4, 0.005, 1.0  # K = 1 rad/V → R = 1 rad = 200 σ
    walk = toy("toy_noise_ar1", "walk", out="v.d", unit="rad", a=1.0, b=sigma)
    obj = chain(walk, 1.0, u_max, {"kp": 0.0, "ki_per_s": 0.2 * fs})  # g = K·K_i/f_sim = 0.2
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, 2121)
    loop = compile_loop(blk)
    n_trials = 3000
    exits = np.empty(n_trials)
    mean_samples = (u_max / sigma) ** 2  # R²/(2D) in samples, 2 D h = σ²
    for i in range(n_trials):
        st, pr = State(lay), Params(lay)
        while int(st["lock.sup.n_loss"]) == 0:
            run_compiled(blk, st, pr, int(mean_samples // 2), loop)
        exits[i] = float(st["lock.sup.t_loss"])
    d = sigma**2 * fs / 2
    expected = u_max**2 / (2 * d)
    rel_3sigma = 3 * math.sqrt(2 / 3) / math.sqrt(n_trials)
    assert rel_3sigma < 0.05
    assert exits.mean() == pytest.approx(expected, rel=0.05)
