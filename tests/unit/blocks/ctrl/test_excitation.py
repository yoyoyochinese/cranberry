"""``excitation`` (D85-7): reference = compiled (1e−12), waveform properties, toggle, runtime amplitude, validation."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.blocks.ctrl.excitation import excitation_schedule, prbs_bits
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

FS = 1.0e5
SRC = toy("toy_noise_ar1", "src", out="v.x", unit="V", a=0.5, b=1.0)

WAVES = {
    "stepped_sine": dict(freqs_hz=[1.0e3, 2.3e3, 4.1e3], dwell_s=2e-3, settle_s=5e-4, amp_f_corner_hz=2.0e3),
    "multisine": dict(f_start_hz=1.0e3, f_stop_hz=2.0e4, n_freqs=8, spacing="log", period_s=2e-3, settle_s=0.0,
                      amp_f_corner_hz=math.inf),
    "chirp": dict(f_start_hz=500.0, f_stop_hz=2.0e4, spacing="log", period_s=3e-3, settle_s=0.0),
    "prbs": dict(prbs_order=5, bit_s=2e-5, settle_s=0.0),
    "white": dict(settle_s=0.0),
}


def exc(waveform, enabled=True, rate_div=1, src="v.x", amplitude=0.3, t_start_s=1e-3, **over):
    p = {"in": src, "out": "v.y", "unit": "V", "waveform": waveform, "amplitude": amplitude, "t_start_s": t_start_s,
         "toggles": {"enabled": enabled}, **WAVES[waveform], **over}
    return {"block": "excitation", "name": "exc", "rate_div": rate_div, "params": p}


def make(obj, seed=3, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def trace(obj, n, slots=("v.x", "v.y", "ctrl.exc.x", "ctrl.exc.f_hz", "ctrl.exc.step"), fs=FS):
    blk, st, pr = make(obj, fs=fs)
    out = {s: np.empty(n) for s in slots}
    for i in range(n):
        run_reference(blk, st, pr, 1)
        for s in slots:
            out[s][i] = float(st[s])
    return blk, out


@pytest.mark.numba
@pytest.mark.parametrize("waveform", sorted(WAVES))
@pytest.mark.parametrize("rate_div", [1, 2])
def test_reference_equals_compiled(waveform, rate_div):
    obj = [SRC, exc(waveform, rate_div=rate_div)]
    (b1, s1, p1), (b2, s2, p2) = make(obj), make(obj)
    run_reference(b1, s1, p1, 1500)
    run_compiled(b2, s2, p2, 1500)
    for k in ("v.y", "ctrl.exc.x", "ctrl.exc.f_hz"):
        np.testing.assert_allclose(s2[k], s1[k], rtol=1e-12, atol=1e-12)
    assert int(s1["ctrl.exc.step"]) == int(s2["ctrl.exc.step"])


def test_stepped_sine_schedule_and_waveform():
    t0, a = 1e-3, 0.3
    _, tr = trace([SRC, exc("stepped_sine", t_start_s=t0, amplitude=a)], 800)
    sc = excitation_schedule({**WAVES["stepped_sine"], "waveform": "stepped_sine"})
    tm = sc["dwell_s"] - sc["settle_s"]
    np.testing.assert_allclose(sc["freqs_hz"] * tm, np.round(sc["freqs_hz"] * tm))   # whole periods per window
    np.testing.assert_allclose(sc["gains"], np.maximum(1.0, sc["freqs_hz"] / 2.0e3))
    n = np.arange(800)
    rel = n - round(t0 * FS)
    k = np.where(rel >= 0, rel // round(sc["dwell_s"] * FS), -1)
    k = np.where(k >= len(sc["freqs_hz"]), len(sc["freqs_hz"]), k)
    np.testing.assert_array_equal(tr["ctrl.exc.step"], k)
    on = (k >= 0) & (k < len(sc["freqs_hz"]))
    kk = np.clip(k, 0, len(sc["freqs_hz"]) - 1)
    tau = (rel - kk * round(sc["dwell_s"] * FS)) / FS
    x = np.where(on, a * sc["gains"][kk] * np.sin(2 * np.pi * sc["freqs_hz"][kk] * tau), 0.0)
    np.testing.assert_allclose(tr["ctrl.exc.x"], x, atol=1e-13)
    np.testing.assert_allclose(tr["ctrl.exc.f_hz"], np.where(on, sc["freqs_hz"][kk], 0.0))
    np.testing.assert_allclose(tr["v.y"], tr["v.x"] + tr["ctrl.exc.x"], atol=1e-14)


def test_multisine_tones_rms_and_crest():
    a = 0.3
    blk, tr = trace([exc("multisine", src=None, t_start_s=0.0, amplitude=a)], 400, slots=("ctrl.exc.x",))
    sc = excitation_schedule({**WAVES["multisine"], "waveform": "multisine"})
    per = tr["ctrl.exc.x"][:200]                                       # one period (2 ms at 100 kHz)
    np.testing.assert_allclose(tr["ctrl.exc.x"][200:], per, atol=1e-12)  # periodic
    assert np.sqrt(np.mean(per**2)) == pytest.approx(a / math.sqrt(2), rel=1e-9)   # amplitude = √2 RMS
    spec = np.fft.rfft(per) * 2 / per.size
    m = np.round(sc["freqs_hz"] * sc["period_s"]).astype(int)
    np.testing.assert_allclose(np.abs(spec[m]), a / math.sqrt(m.size), rtol=1e-9)
    others = np.delete(np.abs(spec), m)
    assert others.max() < 1e-9
    # Schroeder phases give a lower crest factor than equal phases
    zero = np.cos(2 * np.pi * np.outer(np.arange(200), m) / 200).sum(axis=1)
    crest = lambda y: np.abs(y).max() / np.sqrt(np.mean(y**2))  # noqa: E731
    assert crest(per) < 0.7 * crest(zero)


def test_chirp_instantaneous_frequency():
    _, tr = trace([exc("chirp", src=None, t_start_s=0.0, spacing="lin")], 300, slots=("ctrl.exc.x", "ctrl.exc.f_hz"))
    tau = np.arange(300) / FS
    f0, f1, T = 500.0, 2.0e4, 3e-3
    np.testing.assert_allclose(tr["ctrl.exc.f_hz"], f0 + (f1 - f0) * tau / T, rtol=1e-12)
    np.testing.assert_allclose(tr["ctrl.exc.x"], 0.3 * np.sin(2 * np.pi * (f0 * tau + 0.5 * (f1 - f0) * tau**2 / T)),
                               atol=1e-12)
    _, tr = trace([exc("chirp", src=None, t_start_s=0.0)], 300, slots=("ctrl.exc.f_hz",))   # log: geometric
    np.testing.assert_allclose(tr["ctrl.exc.f_hz"], f0 * (f1 / f0) ** (tau / T), rtol=1e-12)


def test_prbs_maximal_length_and_levels():
    for order in (3, 5, 7, 10):
        b = prbs_bits(order)
        n = 2**order - 1
        assert b.size == n and np.sum(b > 0) == 2 ** (order - 1)
        ac = np.array([np.dot(b, np.roll(b, k)) for k in range(1, n)])
        np.testing.assert_array_equal(ac, -1.0)                       # two-valued autocorrelation: m-sequence
    _, tr = trace([exc("prbs", src=None, t_start_s=0.0)], 124, slots=("ctrl.exc.x", "ctrl.exc.step"))
    np.testing.assert_allclose(np.abs(tr["ctrl.exc.x"]), 0.3 / math.sqrt(2))
    np.testing.assert_allclose(tr["ctrl.exc.x"][:62:2], 0.3 / math.sqrt(2) * prbs_bits(5))   # bit = 2 samples
    np.testing.assert_array_equal(tr["ctrl.exc.step"][[0, 61, 62]], [0, 0, 1])


def test_white_std_and_draw_stream():
    _, tr = trace([exc("white", src=None, t_start_s=0.0, amplitude=2.0)], 20000, slots=("ctrl.exc.x",))
    assert np.std(tr["ctrl.exc.x"]) == pytest.approx(2.0 / math.sqrt(2), rel=0.03)


def test_toggle_off_is_identity_and_amplitude_is_runtime():
    _, tr = trace([SRC, exc("stepped_sine", enabled=False)], 300)
    np.testing.assert_array_equal(tr["v.y"], tr["v.x"])
    np.testing.assert_array_equal(tr["ctrl.exc.x"], 0.0)
    np.testing.assert_array_equal(tr["ctrl.exc.step"], -1)
    blk, st, pr = make([exc("stepped_sine", src=None, t_start_s=0.0, amplitude=1.0)])
    run_reference(blk, st, pr, 17)
    x1 = float(st["ctrl.exc.x"])
    blk, st, pr = make([exc("stepped_sine", src=None, t_start_s=0.0, amplitude=1.0)])
    pr.set({"exc.amplitude": 0.25})
    run_reference(blk, st, pr, 17)
    assert float(st["ctrl.exc.x"]) == pytest.approx(0.25 * x1)


def test_validation():
    no_dwell = exc("stepped_sine", src=None)
    del no_dwell["params"]["dwell_s"]
    with pytest.raises(MissingParameterError):
        make([no_dwell])
    with pytest.raises(ValueError, match="dwell_s must exceed"):
        make([exc("stepped_sine", dwell_s=4e-4)])
    with pytest.raises(ValueError, match="positive integer"):
        make([exc("multisine", src=None, period_s=2.000005e-3)])
    with pytest.raises(ValueError, match="Nyquist"):
        make([exc("stepped_sine", src=None, freqs_hz=[6.0e4])])
    with pytest.raises(ValueError, match="waveform"):
        bad = exc("stepped_sine", src=None)
        bad["params"]["waveform"] = "square"
        make([bad])
