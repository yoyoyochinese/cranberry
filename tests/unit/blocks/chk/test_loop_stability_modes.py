"""``loop_stability_check`` mode lists (D85-1): with ``r`` on every mode the actuator response is the
parallel sum of ``pzt_fiber_stretcher`` (``modes_ss``), at the top level and per stage; without ``r`` it
stays the series product (presets before D85)."""

import numpy as np

import cranberry.core.testing  # noqa: F401
from cranberry.blocks.act.pzt import modes_ss
from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj

# main mode + a low mode with an anti-resonance just below it (r < 0) + a high mode; Σr = 1
MODES = [{"f_hz": 25.0e3, "q": 10.0, "r": 1.15}, {"f_hz": 5.3e3, "q": 18.0, "r": -0.05},
         {"f_hz": 40.0e3, "q": 20.0, "r": -0.10}]
F = np.geomspace(100.0, 2.0e5, 3000)


def _check(modes, **kw):
    p = {"k0": 1.0, "kp": 1.0, "ki_per_s": 0.0, "kd_s": 0.0, "delay_s": 0.0, "poles_hz": [], "modes": modes}
    p.update(kw)
    return p


def _ss_response(modes, f):
    A, B, C, _ = modes_ss(modes)
    return np.array([(C @ np.linalg.solve(2j * np.pi * x * np.eye(A.shape[0]) - A, B))[0, 0] for x in f])


def test_parallel_sum_matches_pzt_block():
    want = _ss_response(MODES, F)
    np.testing.assert_allclose(open_loop(F, _check(MODES)), want, rtol=1e-10, atol=1e-12)
    # r < 0 puts an anti-resonance (a zero pair) just below the 5.3 kHz mode, absent from the product
    mag = np.abs(want)
    dip, peak = mag[(F > 4.5e3) & (F < 5.3e3)].min(), mag[(F > 5.3e3) & (F < 6.0e3)].max()
    assert dip < 0.75 < 1.5 < peak  # DC gain 1 (Σr = 1); the product of the same sections has no dip
    # and the time-loop block itself (ZOH at 50 MHz; half-sample lag ≤ 2.5 mrad at 40 kHz)
    lay = SlotLayout(50.0e6)
    obj = {"block": "pzt_fiber_stretcher", "name": "a", "params": {
        "in": "v.u", "gain_mode": "direct", "k_rad_per_v": 1.0, "polarity": 1.0, "modes": MODES,
        "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -5.0, "drive_max_v": 5.0,
        "hysteresis_frac": 0.0, "segment": None, "insertion_loss_db": None,
        "toggles": {"resonance": True, "hysteresis": False, "saturation": False, "insertion_loss": False}}}
    ramp = {"block": "ramp", "name": "u", "params": {"out": "v.u", "waveform": "triangle", "amplitude_v": 1.0,
                                                     "period_s": 1e-3, "offset_v": 0.0, "toggles": {"enabled": True}}}
    act = build(spec_from_obj([ramp, obj]), lay, 1).children[1]
    fs = F[::50]
    np.testing.assert_allclose(act.freq_response(fs), open_loop(fs, _check(MODES)), rtol=5e-3, atol=1e-3)


def test_stage_parallel_sum():
    """A stage with an r-mode list: Σ = 1 + C_2·k_act·M_2 with M_2 the parallel sum."""
    st = {"kp": 2.0, "ki_per_s": 0.0, "kd_s": 0.0, "k_act": 0.5, "poles_hz": [], "modes": MODES}
    got = open_loop(F, _check([], stages=[st]))
    np.testing.assert_allclose(got, 1.0 + 2.0 * 0.5 * _ss_response(MODES, F), rtol=1e-10, atol=1e-12)


def test_without_r_product_unchanged():
    modes = [{"f_hz": m["f_hz"], "q": m["q"]} for m in MODES]
    s = 2j * np.pi * F
    want = np.ones_like(s)
    for m in modes:
        w = 2 * np.pi * m["f_hz"]
        want = want * w * w / (s * s + s * w / m["q"] + w * w)
    np.testing.assert_allclose(open_loop(F, _check(modes)), want, rtol=1e-12)
    # mixed list (one r missing) is also the product; one mode with r = 1 is the same either way
    mixed = [dict(MODES[0]), *modes[1:]]
    np.testing.assert_allclose(open_loop(F, _check(mixed)), want, rtol=1e-12)
    one = [{"f_hz": 25e3, "q": 10.0, "r": 1.0}]
    np.testing.assert_allclose(open_loop(F, _check(one)), open_loop(F, _check([{"f_hz": 25e3, "q": 10.0}])),
                               rtol=1e-12)
    np.testing.assert_allclose(open_loop(F, _check([])), 1.0)
