"""WP6-P (D85-6): ``phase_program`` (step / ramp / quantisation / dc_setpoint / random tables, reference = kernel),
the backward-compatible slots ``demodulator.in_program``, ``metric_ii.setpoint_slot``, ``lock_supervisor.slip_center_slot``,
and the ``step_response`` segment analysis."""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.ctrl.phase_program import make_phase_program, random_table
from cranberry.blocks.meas.step_response import segment_response

pytestmark = pytest.mark.numba
FS = 1.0e6
TIMES = [5e-6, 20e-6, 40e-6]
TABLE = [[0.3, -0.2], [-0.5, 0.4], [1.0, 1.0]]


def prog(mode="step", output="rf_phase", enabled=True, quant=False, **kw):
    p = {"ports": [0, 1], "source": "table", "times_s": TIMES, "table_rad": TABLE, "mode": mode, "output": output,
         "sign": -1.0, "toggles": {"enabled": enabled, "quantization": quant}}
    if mode == "ramp":
        p["slew_rad_per_s"] = 5e4
    if output == "rf_phase":
        p["quantum_rad"] = 0.07
    else:
        p["v_fringe_v"] = 2.0
    p.update(kw)
    return {"block": "phase_program", "name": "pg", "params": p}


def run(compiled=False, n=60, **kw):
    out = "rf.pg.port" if kw.get("output", "rf_phase") == "rf_phase" else "ctrl.pg.vset"
    slots = [f"{out}0", f"{out}1", "ctrl.pg.target0", "ctrl.pg.target1", "ctrl.pg.index"]
    rec, _, _ = pt.run_and_record([prog(**kw), pt.record(slots, n)], FS, n, compiled=compiled)
    return [rec[s] for s in slots]


def goal(t):
    j = np.searchsorted(TIMES, t + 1e-15, side="right") - 1
    return np.where(j[:, None] >= 0, np.array(TABLE)[np.maximum(j, 0)], 0.0), j


def test_step_table_sign_and_index():
    rf0, rf1, c0, c1, idx = run()
    t = np.arange(60) / FS
    g, j = goal(t)
    np.testing.assert_array_equal(np.c_[c0, c1], g)
    np.testing.assert_array_equal(np.c_[rf0, rf1], -g)  # sign −1
    np.testing.assert_array_equal(idx, j)
    assert idx[4] == -1 and idx[5] == 0  # the switch lands on the sample t = 5 µs


def test_ramp_slew_limit_and_quantisation():
    _, _, c0, c1, _ = run(mode="ramp")
    lim = 5e4 / FS
    for c in (c0, c1):
        assert np.max(np.abs(np.diff(c))) <= lim + 1e-15
    assert c0[19] == pytest.approx(0.3) and c0[20] == pytest.approx(0.3 - lim)  # reached, then slews
    rf0, _, c0, _, _ = run(mode="ramp", quant=True)
    np.testing.assert_allclose(rf0, 0.07 * np.floor(-c0 / 0.07 + 0.5), rtol=0, atol=1e-15)
    assert np.max(np.abs(rf0 + c0)) <= 0.035 + 1e-12


def test_dc_setpoint_and_disabled():
    v0, v1, c0, c1, _ = run(output="dc_setpoint")
    np.testing.assert_allclose(v0, 2.0 * np.sin(-c0), rtol=0, atol=1e-15)
    rf0, rf1, c0, c1, idx = run(enabled=False)
    assert not np.any(np.r_[rf0, rf1, c0, c1]) and np.all(idx == -1)


@pytest.mark.parametrize("kw", [{}, {"mode": "ramp", "quant": True}, {"output": "dc_setpoint", "mode": "ramp"}])
def test_reference_equals_kernel(kw):
    a, b = run(compiled=False, **kw), run(compiled=True, **kw)
    for x, y in zip(a, b):
        np.testing.assert_allclose(x, y, rtol=0, atol=1e-12)


def test_random_tables():
    walk = random_table(30, 3, math.pi / 4, 7)
    steps = np.abs(np.diff(np.vstack([np.zeros(3), walk]), axis=0))
    np.testing.assert_allclose(steps, math.pi / 4, rtol=0, atol=1e-12)  # every switch is exactly step_rad
    assert np.all(np.abs(walk) <= math.pi / 2 + 1e-12)
    uni = random_table(30, 3, None, 7)
    assert np.all(np.abs(uni) <= math.pi / 2) and uni.std() > 0.6
    np.testing.assert_array_equal(random_table(30, 3, None, 7), uni)  # deterministic in the seed
    blk = make_phase_program({**prog()["params"], "source": "random", "n_settings": 30, "t_first_s": 1e-3,
                              "dwell_s": 2e-3, "step_rad": None, "seed": 7}, name="r")
    assert blk.table.shape == (30, 2) and blk.times[-1] == pytest.approx(1e-3 + 29 * 2e-3)


def test_validation():
    base = prog()["params"]
    for bad, msg in (({"mode": "slow"}, "mode"), ({"sign": 2.0}, "sign"), ({"table_rad": [[0.0, 0.0]]}, "shape"),
                     ({"times_s": [2e-6, 1e-6, 3e-6]}, "increasing"), ({"source": "random", "n_settings": 3,
                      "t_first_s": 0.0, "dwell_s": 1e-3, "step_rad": 2.0, "seed": 0}, "step_rad")):
        with pytest.raises(ValueError, match=msg):
            make_phase_program({**base, **bad})
    with pytest.raises(KeyError):
        make_phase_program({**base, "mode": "ramp"})  # slew_rad_per_s required (D34)


# ---- backward-compatible slots ------------------------------------------------------------------------------------------

def test_demodulator_in_program():
    """I = a cos φ − b sin φ with φ = φ_d + rf: (a, b) = (cos X, sin X) gives I = cos(X + φ_d + rf)."""
    x = 0.4
    gen = pt.pathtest("pathtest_table", "ab", "v.ab", shape=[2], unit="V", n_rows=60)
    dem = {"block": "demodulator", "name": "dm",
           "params": {"in": "v.ab", "in_phase": None, "in_program": "rf.pg.port1", "out_i": "v.i", "out_q": "v.q",
                      "lpf_order": 0, "f_ref_hz": 1e5, "phi_d_rad": 0.1, "v_off_v": 0.0, "lpf_hz": 1e5,
                      "toggles": {"rf_rate_mode": False, "demod_offset": False, "rf_reference_phase_noise": False}}}
    tab = np.tile([math.cos(x), math.sin(x)], (60, 1))
    for compiled in (False, True):
        rec, _, _ = pt.run_and_record([prog(), gen, dem, pt.record(["v.i", "rf.pg.port1"], 60)], FS, 60,
                                      {"v.ab_table": tab}, compiled=compiled)
        np.testing.assert_allclose(rec["v.i"], np.cos(x + 0.1 + rec["rf.pg.port1"]), rtol=0, atol=1e-14)
    old = {**dem, "params": {k: v for k, v in dem["params"].items() if k != "in_program"}}
    rec, _, _ = pt.run_and_record([gen, old, pt.record(["v.i"], 10)], FS, 10, {"v.ab_table": tab})
    np.testing.assert_allclose(rec["v.i"], math.cos(x + 0.1), rtol=0, atol=1e-14)  # absent: unchanged


def test_supervisor_slip_center_slot():
    """The true error follows the program; with the centre following it too, no slip is counted, while a fixed centre
    counts the excursion across δ_c + π as a slip."""
    from cranberry.core.testing import toy

    def pipe(slot):
        p = {"u_in": "v.zero", "e_in": "v.zero", "true_in": "v.true", "out": "ctrl.sup.u", "relock": "none",
             "u_min_v": -1.0, "u_max_v": 1.0, "t_sat_s": 0.0, "e_th": 1e9, "t_e_s": 0.0, "slip_center_rad": 0.0,
             "toggles": {"relock": False}}
        if slot:
            p["slip_center_slot"] = "ctrl.pg.target1"
        big = {**prog()["params"], "table_rad": [[0.0, 2.5], [0.0, 3.5], [0.0, 0.5]]}
        return [{"block": "phase_program", "name": "pg", "params": big},
                toy("toy_gain", "z", in_="ctrl.pg.target0", out="v.zero", unit="V", gain=0.0),
                toy("toy_gain", "tr", in_="ctrl.pg.target1", out="v.true", unit="rad", gain=1.0),
                {"block": "lock_supervisor", "name": "sup", "params": p}, pt.record(["v.true"], 60)]

    for compiled in (False, True):
        _, fixed, _ = pt.run_and_record(pipe(False), FS, 60, compiled=compiled)
        _, moving, _ = pt.run_and_record(pipe(True), FS, 60, compiled=compiled)
        assert int(fixed["lock.sup.n_slip"]) == 2 and int(moving["lock.sup.n_slip"]) == 0


def test_segment_response_first_order():
    """First-order step y = Δ(1 − e^{−t/τ}): rise 10–90 % = τ ln 9, settle to tol = τ ln(Δ/tol), no overshoot."""
    tau, d, tol = 1e-4, 0.5, 5e-3
    t = np.arange(0, 2e-3, 1e-7)
    r = segment_response(t, d * (1 - np.exp(-t / tau)), 0.0, 0.0, d, tol)
    assert r["rise_s"] == pytest.approx(tau * math.log(9), rel=1e-3)
    assert r["settle_s"] == pytest.approx(tau * math.log(d / tol), rel=1e-3)
    assert r["overshoot"] == 0.0 and abs(r["bias_rad"]) < 1e-5  # tail e^{−10}·Δ
    r = segment_response(t, d * (1 - np.exp(-t / tau) * np.cos(2e4 * t)), 0.0, 0.0, d, tol)
    assert r["overshoot"] > 0.05
    assert math.isnan(segment_response(t, np.full(t.size, 0.1), 0.0, 0.0, d, tol)["settle_s"])
