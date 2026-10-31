"""``iir`` / ``lpf`` / ``hpf`` / ``band_split`` against scipy, complementarity, toggles; ``adder``,
``ramp``, ``sine_injection``."""

import math

import numpy as np
import pytest
from scipy import signal

import cranberry.core.testing  # noqa: F401
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

FS = 1.0e6


def make(obj, seed=4, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def trace(obj, n, slots, fs=FS):
    blk, st, pr = make(obj, fs=fs)
    out = {s: np.empty(n) for s in slots}
    for i in range(n):
        run_reference(blk, st, pr, 1)
        for s in slots:
            out[s][i] = float(st[s])
    return blk, out


SRC = toy("toy_noise_ar1", "x", out="v.x", unit="V", a=0.5, b=1.0)


def filt(kind, enabled=True, **kw):
    p = {"in": "v.x", "unit": "V", "toggles": {"enabled": enabled}, **kw}
    if kind == "band_split":
        p.update(out_low="v.lo", out_high="v.hi")
    else:
        p["out"] = "v.y"
    return {"block": kind, "name": "f", "params": p}


@pytest.mark.parametrize("kind, method", [("lpf", "bilinear"), ("hpf", "bilinear"), ("lpf", "zoh")])
def test_designed_filter_equals_scipy(kind, method):
    btype = "lowpass" if kind == "lpf" else "highpass"
    _, tr = trace([SRC, filt(kind, order=3, f_c_hz=40e3, method=method)], 2000, ["v.x", "v.y"])
    if method == "bilinear":
        sos = signal.butter(3, 40e3, btype=btype, fs=FS, output="sos")
        assert abs(signal.sosfreqz(sos, [40e3], fs=FS)[1][0]) == pytest.approx(1 / math.sqrt(2), rel=1e-9)
    else:
        b, a = signal.butter(3, 2 * np.pi * 40e3, analog=True)
        bd, ad, _ = signal.cont2discrete((b, a), 1 / FS, method="zoh")
        np.testing.assert_allclose(tr["v.y"], signal.lfilter(np.squeeze(bd), ad, tr["v.x"]), rtol=1e-8, atol=1e-11)
        return
    np.testing.assert_allclose(tr["v.y"], signal.sosfilt(sos, tr["v.x"]), rtol=1e-9, atol=1e-12)


def test_first_order_zoh_is_exact_held_input_response():
    """1st-order ZOH LPF: y[n] = a y[n−1] + (1 − a) x[n−1], a = e^{−2π f_c h}."""
    _, tr = trace([SRC, filt("lpf", order=1, f_c_hz=10e3, method="zoh")], 300, ["v.x", "v.y"])
    a = math.exp(-2 * math.pi * 10e3 / FS)
    y = signal.lfilter([0.0, 1 - a], [1.0, -a], tr["v.x"])
    np.testing.assert_allclose(tr["v.y"], y, rtol=1e-9, atol=1e-13)


def test_iir_given_sos_and_toggle_identity():
    sos = signal.cheby1(4, 0.5, 0.1, output="sos")
    blk, tr = trace([SRC, filt("iir", sos=sos.tolist())], 1000, ["v.x", "v.y"])
    np.testing.assert_allclose(tr["v.y"], signal.sosfilt(sos, tr["v.x"]), rtol=1e-9, atol=1e-12)
    f = np.array([1e3, 5e4])
    np.testing.assert_allclose(blk.children[1].freq_response(f), signal.sosfreqz(sos, f, fs=FS)[1])
    _, tr = trace([SRC, filt("iir", enabled=False, sos=sos.tolist())], 100, ["v.x", "v.y"])
    np.testing.assert_array_equal(tr["v.y"], tr["v.x"])


def test_band_split_is_complementary():
    _, tr = trace([SRC, filt("band_split", order=2, f_c_hz=5e3, method="bilinear")], 1000, ["v.x", "v.lo", "v.hi"])
    np.testing.assert_allclose(tr["v.lo"] + tr["v.hi"], tr["v.x"], rtol=0, atol=1e-14)
    _, tr = trace([SRC, filt("band_split", enabled=False, order=2, f_c_hz=5e3, method="bilinear")], 100,
                  ["v.x", "v.lo", "v.hi"])
    np.testing.assert_array_equal(tr["v.lo"], tr["v.x"])
    np.testing.assert_array_equal(tr["v.hi"], 0.0)


def test_filter_validation():
    with pytest.raises(ValueError, match="f_c_hz"):
        make([SRC, filt("lpf", order=1, f_c_hz=600e3, method="bilinear")])
    with pytest.raises(ValueError, match="method"):
        make([SRC, filt("lpf", order=1, f_c_hz=1e3, method="matched")])


def test_adder_ramp_sine():
    obj = [SRC,
           {"block": "ramp", "name": "saw", "params": {"out": "v.saw", "waveform": "sawtooth", "amplitude_v": 2.0,
                                                       "period_s": 1e-4, "offset_v": 0.5, "toggles": {"enabled": True}}},
           {"block": "ramp", "name": "tri", "params": {"out": "v.tri", "waveform": "triangle", "amplitude_v": 1.0,
                                                       "period_s": 1e-4, "offset_v": 0.0, "toggles": {"enabled": True}}},
           {"block": "sine_injection", "name": "inj", "params": {"in": "v.x", "out": "v.inj", "unit": "V",
                                                                 "freq_hz": 3e3, "amplitude": 0.2,
                                                                 "toggles": {"enabled": True}}},
           {"block": "adder", "name": "sum", "params": {"ins": ["v.saw", "v.inj"], "gains": [2.0, -1.0],
                                                        "out": "v.sum", "unit": "V"}}]
    _, tr = trace(obj, 300, ["t", "v.x", "v.saw", "v.tri", "v.inj", "v.sum"])
    t = np.arange(300) / FS  # value written at sample n sees t = n / f_sim
    ph = t / 1e-4 - np.floor(t / 1e-4)
    np.testing.assert_allclose(tr["v.saw"], 0.5 + 2.0 * (2 * ph - 1), atol=1e-12)
    np.testing.assert_allclose(tr["v.tri"], 1.0 - 4.0 * np.abs(ph - 0.5), atol=1e-12)
    np.testing.assert_allclose(tr["v.inj"], tr["v.x"] + 0.2 * np.sin(2 * np.pi * 3e3 * t), atol=1e-14)
    np.testing.assert_allclose(tr["v.sum"], 2 * tr["v.saw"] - tr["v.inj"], atol=1e-14)


def test_signal_toggles_off():
    obj = [SRC,
           {"block": "ramp", "name": "r", "params": {"out": "v.r", "waveform": "sawtooth", "amplitude_v": 2.0,
                                                     "period_s": 1e-4, "offset_v": 0.5, "toggles": {"enabled": False}}},
           {"block": "sine_injection", "name": "inj", "params": {"in": "v.x", "out": "v.inj", "unit": "V",
                                                                 "freq_hz": 3e3, "amplitude": 0.2,
                                                                 "toggles": {"enabled": False}}}]
    _, tr = trace(obj, 100, ["v.x", "v.r", "v.inj"])
    np.testing.assert_array_equal(tr["v.r"], 0.0)
    np.testing.assert_array_equal(tr["v.inj"], tr["v.x"])
