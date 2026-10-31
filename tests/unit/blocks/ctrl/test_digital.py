"""``adc`` / ``dac``: quantiser, range, averaging, noise, divider, device values, toggles; DAC
quantisation inside the 915B loop (12.7 mrad/LSB bound, limit cycle, white-noise model; T16 V11, T15 §7-10)."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy
from cranberry.runner.run import run

FS = 1.0e6
ALL = {"quantization": True, "noise": False, "saturation": True}


def conv(kind, toggles=None, inp="v.x", **kw):
    p = {"in": inp, "out": f"v.{kind}", "device": None, "bits": 12, "range_v": [-4.0, 4.0], "noise_vrms": 0.0,
         "toggles": dict(toggles or ALL)}
    p.update({"avg_samples": 1} if kind == "adc" else {"r_out_ohm": 220.0, "r_load_ohm": float("inf")})
    p.update(kw)
    return {"block": kind, "name": kind, "params": p}


def make(obj, seed=3, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def run_trace(obj, n, slots, fs=FS, compiled=False):
    blk, st, pr = make(obj, fs=fs)
    out = {s: np.empty(n) for s in slots}
    step = run_compiled if compiled else run_reference
    for i in range(n):
        step(blk, st, pr, 1)
        for s in slots:
            out[s][i] = float(st[s])
    return out


def sweep_src():
    return {"block": "ramp", "name": "x", "params": {"out": "v.x", "waveform": "sawtooth", "amplitude_v": 5.0,
                                                     "period_s": 1e-3, "offset_v": 0.0, "toggles": {"enabled": True}}}


@pytest.mark.parametrize("kind", ["adc", "dac"])
def test_quantiser_and_range(kind):
    tr = run_trace([sweep_src(), conv(kind)], 1000, ["v.x", f"v.{kind}"])
    x, y = tr["v.x"], tr[f"v.{kind}"]
    q = 8.0 / 2**12
    assert q == pytest.approx(1.953e-3, rel=1e-3)
    inside = (x > -4.0) & (x < 4.0 - q)
    assert np.all(np.abs(y[inside] - x[inside]) <= q / 2 + 1e-15)
    np.testing.assert_allclose(y / q, np.round(y / q), atol=1e-9)  # on the code grid
    assert y.max() == pytest.approx(4.0 - q) and y.min() == pytest.approx(-4.0)  # two's complement codes


@pytest.mark.parametrize("kind", ["adc", "dac"])
def test_toggles_off_identity(kind):
    off = {"quantization": False, "noise": False, "saturation": False}
    tr = run_trace([sweep_src(), conv(kind, off)], 500, ["v.x", f"v.{kind}"])
    np.testing.assert_array_equal(tr[f"v.{kind}"], tr["v.x"])
    sat_only = {"quantization": False, "noise": False, "saturation": True}
    tr = run_trace([sweep_src(), conv(kind, sat_only)], 500, ["v.x", f"v.{kind}"])
    np.testing.assert_array_equal(tr[f"v.{kind}"], np.clip(tr["v.x"], -4, 4))


def test_saturation_flag():
    blk, st, pr = make([sweep_src(), conv("dac")])
    flags = []
    for _ in range(1000):
        run_reference(blk, st, pr, 1)
        flags.append((float(st["v.x"]), int(st["ctrl.dac.sat"])))
    assert all(f == (1 if x > 4.0 - 8 / 2**13 else (-1 if x < -4.0 - 8 / 2**13 else 0)) or abs(abs(x) - 4) < 3e-3
               for x, f in flags)


def test_adc_average_and_noise():
    tr = run_trace([sweep_src(), conv("adc", {"quantization": False, "noise": False, "saturation": False},
                                      avg_samples=2)], 200, ["v.x", "v.adc"])
    np.testing.assert_allclose(tr["v.adc"][1:], 0.5 * (tr["v.x"][1:] + tr["v.x"][:-1]), rtol=1e-13, atol=1e-15)
    const = {"block": "ramp", "name": "x", "params": {"out": "v.x", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                      "period_s": 1.0, "offset_v": 0.1, "toggles": {"enabled": True}}}
    tr = run_trace([const, conv("adc", {"quantization": False, "noise": True, "saturation": False},
                                noise_vrms=3e-4)], 40_000, ["v.adc"], compiled=True)
    assert np.mean(tr["v.adc"]) == pytest.approx(0.1, abs=5 * 3e-4 / 200)
    assert np.std(tr["v.adc"]) == pytest.approx(3e-4, rel=3 * math.sqrt(0.5 / 40_000) + 0.005)


def test_dac_divider_and_device_values():
    tr = run_trace([sweep_src(), conv("dac", {"quantization": False, "noise": False, "saturation": False},
                                      r_load_ohm=1e4)], 50, ["v.x", "v.dac"])
    np.testing.assert_allclose(tr["v.dac"], tr["v.x"] * 1e4 / (1e4 + 220.0), rtol=1e-15)
    obj = {"block": "dac", "name": "dac", "params": {"in": "v.x", "out": "v.dac", "device": "dsc1",
                                                     "noise_vrms": 0.0, "r_load_ohm": float("inf"),
                                                     "toggles": ALL}}
    blk, _, _ = make([sweep_src(), obj])
    leaf = blk.children[1]
    assert leaf.values["bits"] == 12 and leaf.values["r_out_ohm"] == 220.0 and leaf.lsb_v == pytest.approx(8 / 4096)
    obj = {"block": "adc", "name": "adc", "params": {"in": "v.x", "out": "v.adc", "device": "dsc1",
                                                     "noise_vrms": 0.0, "bits": 14, "toggles": ALL}}
    blk, _, _ = make([sweep_src(), obj])
    assert blk.children[1].values["bits"] == 14 and blk.children[1].values["avg_samples"] == 2  # spec overrides


# ----------------------------------------------------------------- DAC quantisation in the 915B loop

K_S = 2 * math.pi * 16e-6 / (1548.8e-9 * 10.0)  # 6.49 rad/V (8669-1, δOPL reading)


def loop(dist, ki, quant=True, chain914=True):
    """dist → e = d − φ_a → pid (I, DSC1 law, 500 kHz) → dac (DSC1, 12 bit) → delay → [914] → 915B."""
    acts = []
    inp = "delay.dl.out"
    if chain914:
        acts.append({"block": "hv_amp_914", "name": "amp", "params": {
            "in": inp, "out": "v.amp", "poles_hz": [29.5e3, 112e3], "discretisation": "zoh", "prewarp_hz": None,
            "input_min_v": -5.0, "input_max_v": 5.0, "noise_v_per_rthz": 0.0, "slew_v_per_s": 1e9,
            "toggles": {"poles": True, "slew": False, "saturation": True, "noise": False}}})
        inp = "v.amp"
    return [dist,
            {"block": "adder", "name": "err", "params": {"ins": ["v.d", "act.s.dpsi"], "gains": [1.0, -1.0],
                                                         "out": "v.e", "unit": "rad"}},
            {"block": "pid", "name": "pid", "rate_div": 2 if chain914 else 1, "params": {
                "in": "v.e", "out": "ctrl.pid.u", "device": "dsc1" if chain914 else None, "discretisation": "dsc1",
                "anti_windup": "none", "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "omega_is_rad_per_s": 0.0,
                "omega_ds_rad_per_s": 0.0, "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0,
                "u_min_v": -4.0, "u_max_v": 4.0,
                "toggles": {"feedback": True, "p_term": False, "i_term": True, "d_term": False}}},
            {"block": "dac", "name": "dac", "rate_div": 2 if chain914 else 1, "params": {
                "in": "ctrl.pid.u", "out": "v.dac", "device": "dsc1", "noise_vrms": 0.0, "r_load_ohm": float("inf"),
                "toggles": {"quantization": quant, "noise": False, "saturation": True}}},
            {"block": "delay_buffer", "name": "dl", "params": {"in": "v.dac", "out": "delay.dl.out",
                                                               "tau_s": 4.7e-6 if chain914 else 0.0,
                                                               "method": "linear"}},
            *acts,
            {"block": "pzt_fiber_stretcher", "name": "s", "params": {
                "in": inp, "gain_mode": "dopl", "dopl_full_m": 16e-6, "v_full_v": 10.0, "lambda0_m": 1548.8e-9,
                "polarity": 1.0, "modes": [], "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -5.0,
                "drive_max_v": 5.0, "hysteresis_frac": 0.0, "segment": None, "insertion_loss_db": None,
                "toggles": {"resonance": False, "hysteresis": False, "saturation": True, "insertion_loss": False}}}]


def record(obj, n, slot, fs=FS):
    spec = RunSpec(run=RunSettings(duration_s=n / fs, f_sim_hz=fs, chunk_samples=1 << 16, seed=5),
                   pipeline=spec_from_obj(obj), measure={"record": {"slots": [slot], "decim": 1}})
    return run(spec, keep_state=False).series[slot]


def dist(r, offset):
    """d(t) = offset + r·t (sawtooth of period 10 s, first period only)."""
    return {"block": "ramp", "name": "d", "params": {"out": "v.d", "waveform": "sawtooth", "amplitude_v": 5.0 * r,
                                                     "period_s": 10.0, "offset_v": 5.0 * r + offset,
                                                     "toggles": {"enabled": True}}}


LSB_RAD = K_S * 8.0 / 2**12


@pytest.mark.numba
@pytest.mark.parametrize("offset", [0.0031, 0.007, 0.0105, -0.0222])
def test_dac_lsb_bound_and_limit_cycle(offset):
    """No noise, constant disturbance (not on the code grid), I loop f_c ≈ 2 kHz with the 914 and the
    DSC1 delay: the true phase error stays within K_s·q = 12.7 mrad (V11a) and the limit cycle between
    neighbouring codes stays inside a 1-LSB band at a fixed operating point (half peak-to-peak ≤ K_s·q/2;
    D65-3, which replaces D56-7 / T15 §7-10: ±LSB/2 about the set point is not a bound, only the band)."""
    assert LSB_RAD == pytest.approx(12.7e-3, abs=0.05e-3)
    e = record(loop(dist(0.0, offset), 2 * math.pi * 2e3 / K_S), 60_000, "v.e")[20_000:]
    assert np.max(np.abs(e)) <= LSB_RAD
    assert 0.5 * (e.max() - e.min()) <= LSB_RAD / 2


@pytest.mark.numba
def test_dac_lsb_bound_slow_ramp():
    """Slow phase ramp (1 rad/s ≈ 8 LSB in 0.1 s): |error| ≤ K_s·q throughout. (The limit-cycle centre
    moves with the fractional code position; D65-3: over a ramp the band may approach 2 LSB, while
    |error| < 1 LSB still holds here.)"""
    e = record(loop(dist(1.0, 0.0031), 2 * math.pi * 2e3 / K_S), 100_000, "v.e")[5_000:]
    assert np.max(np.abs(e)) <= LSB_RAD


@pytest.mark.numba
def test_dac_quantisation_white_model_in_loop():
    """V11b: with a disturbance ≫ q the DAC quantisation adds the variance K²(q²/12)·‖z⁻¹S‖² =
    K² q²/12 · 2/(2 − G) (I loop, G = K·K_i·h, static actuator, all at f_sim) within 10 %."""
    g = 0.05
    ki = g * FS / K_S
    white = toy("toy_noise_ar1", "d", out="v.d", unit="rad", a=0.0, b=1.0)
    n = 400_000
    e_on = record(loop(white, ki, True, False), n, "v.e")
    e_off = record(loop(white, ki, False, False), n, "v.e")
    q = 8.0 / 2**12
    pred = K_S**2 * q**2 / 12 * 2 / (2 - g)
    assert np.var(e_on - e_off) == pytest.approx(pred, rel=0.10)
