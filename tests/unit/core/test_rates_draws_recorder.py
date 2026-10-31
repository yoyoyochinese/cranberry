"""core.rates step counting, core.rng.DrawBuffer (G11 sequential draws), runner Recorder."""

import numpy as np
import pytest

from cranberry.core.block import BaseBlock
from cranberry.core.compile import compile_loop, run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import Seq, build
from cranberry.core.rates import PIECE_SAMPLES, max_steps_in, n_steps_in
from cranberry.core.rng import DrawBuffer, stream
from cranberry.core.state import Params, State
from cranberry.runner.run import Recorder


def test_n_steps_in_brute_force():
    for r in (1, 2, 3, 7):
        for n0 in range(0, 15):
            for m in range(0, 20):
                expect = sum(1 for n in range(n0, n0 + m) if n % r == 0)
                assert n_steps_in(n0, m, r) == expect
                assert expect <= max_steps_in(m, r)
    with pytest.raises(ValueError):
        n_steps_in(0, 1, 0)


@pytest.mark.parametrize("rate_div", [1, 3])
def test_draw_buffer_sequential_for_any_pieces(rate_div):
    lay = SlotLayout(1.0)
    db = DrawBuffer("n1", 2, "standard_normal")
    db.init(lay, stream(5, "n1"), rate_div)
    lay.freeze()
    db.bind(lay)
    st = State(lay)
    assert db.capacity == max_steps_in(PIECE_SAMPLES, rate_div)
    got, n0 = [], 0
    for m in (5, 1, PIECE_SAMPLES, 13, 0, 400):
        db.refill(st, n0, m)
        for _ in range(n_steps_in(n0, m, rate_div) * 2):
            got.append(db.take(st))
        n0 += m
    expect = stream(5, "n1").standard_normal(len(got))
    np.testing.assert_array_equal(got, expect)
    with pytest.raises(ValueError):
        db.refill(st, 0, PIECE_SAMPLES + 1)


def test_draw_buffer_validation():
    with pytest.raises(ValueError):
        DrawBuffer("x", 1, "gumbel")
    with pytest.raises(ValueError):
        DrawBuffer("x", 0, "random")
    db = DrawBuffer("x", 1, "random", label="phase")
    assert db.slot_names == ("rng.x.phase.buf", "rng.x.phase.idx")


def test_draw_buffer_refill_checks_consumption():
    """Phase 5 F1-5: refill raises when the previous piece did not consume exactly what was drawn for it
    (unread numbers would be discarded at the piece boundary, making results chunk-length dependent)."""
    from cranberry.core.rng import DrawCountError

    lay = SlotLayout(1.0)
    db = DrawBuffer("n1", 1, "standard_normal")
    db.init(lay, stream(5, "n1"), 1)
    lay.freeze()
    db.bind(lay)
    st = State(lay)
    db.refill(st, 0, 3)
    for _ in range(3):
        db.take(st)
    db.refill(st, 3, 3)  # exact consumption: fine
    db.take(st)
    with pytest.raises(DrawCountError, match="consumed 1 of the 3"):
        db.refill(st, 6, 3)
    st2 = State(lay)  # the check is per state: a fresh state starts clean
    db.refill(st2, 0, 3)


class Writer(BaseBlock):
    """Writes a float vector, an int and a complex 2×2 slot each step (Python and kernel)."""

    def __init__(self):
        super().__init__("writer", "w", {}, writes=("v.w", "lock.w.state", "J.w"))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add("v.w", (2,))
        layout.add("lock.w.state", (), np.int64)
        layout.add("J.w", (2, 2), np.complex128)

    def bind(self, layout):
        import numba

        o_v, o_i, o_z, o_n = (layout.slot(s).offset for s in ("v.w", "lock.w.state", "J.w", "n"))

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            n = I[o_n]
            S[o_v] = n
            S[o_v + 1] = -n
            I[o_i] = 2 * n
            Z[o_z + 3] = n + 1j * n

        self.kernel = k
        self._o = (o_v, o_i, o_z, o_n)

    def __call__(self, state, params):
        o_v, o_i, o_z, o_n = self._o
        n = state.I[o_n]
        state.S[o_v], state.S[o_v + 1], state.I[o_i], state.Z[o_z + 3] = n, -n, 2 * n, n + 1j * n
        return state


@pytest.mark.parametrize("compiled", [False, True])
def test_recorder_float_int_complex(compiled):
    rec = Recorder(["J.w", "v.w", "lock.w.state"], decim=3, chunk_samples=10)
    lay = SlotLayout(1.0e3)
    blk = build(Seq([Writer(), rec]), lay, 0)
    st, p = State(lay), Params(lay)
    tables = []
    for _ in range(3):
        (run_compiled if compiled else run_reference)(blk, st, p, 10)
        tables.append(rec.harvest(st))
    out = rec.split(np.concatenate(tables), lay)
    n = np.arange(0, 30, 3)
    np.testing.assert_array_equal(out["v.w"], np.stack([n, -n], axis=1))
    np.testing.assert_array_equal(out["lock.w.state"], 2 * n)
    assert out["lock.w.state"].dtype == np.int64
    assert out["J.w"].shape == (10, 2, 2)
    np.testing.assert_array_equal(out["J.w"][:, 1, 1], n + 1j * n)
    if compiled:
        assert compile_loop(blk) is not None


def test_recorder_overflow_raises():
    rec = Recorder(["v.w"], decim=1, chunk_samples=4)
    lay = SlotLayout(1.0)
    blk = build(Seq([Writer(), rec]), lay, 0)
    st, p = State(lay), Params(lay)
    with pytest.raises(IndexError):
        run_reference(blk, st, p, 5)



@pytest.mark.parametrize("r", [2, 4])
def test_hold_delay_rule(r):
    """Phase 5 F1-3 (core.rates "Hold delay"): a block held at rate_div r delays by (r − 1)/2 base samples
    (not r/2); the remaining half sample of a physical D/A hold appears only when the downstream
    continuous-time block is discretised with zoh (a bilinear design adds nothing).

    ramp x[n] (integrator of 1) → hold (toy_gain, rate_div r) → 1st-order low-pass (zoh | bilinear);
    the steady-state lag of a unit-DC-gain LTI stage behind a ramp is its DC group delay."""
    import cranberry.blocks.ctrl.filters  # noqa: F401  (registers lpf)
    import cranberry.core.testing  # noqa: F401
    from cranberry.core.pipeline import spec_from_obj
    from cranberry.core.testing import toy

    fs, fc, n = 1.0e6, 1.0e3, 8000

    def run(method):
        obj = [toy("toy_integrator", "ramp", in_="v.one", out="v.x", unit="V", k=fs),
               toy("toy_gain", "hold", rate_div=r, in_="v.x", out="v.h", unit="V", gain=1.0),
               {"block": "lpf", "name": "lp", "params": {"in": "v.h", "out": "v.y", "unit": "V", "order": 1,
                                                         "f_c_hz": fc, "method": method,
                                                         "toggles": {"enabled": True}}}]
        lay = SlotLayout(fs)
        lay.add("v.one", (), np.float64, "V")
        blk = build(spec_from_obj(obj), lay, 0, external=["v.one"])
        st, pa = State(lay), Params(lay)
        st["v.one"] = 1.0
        x, h, y = (np.empty(n) for _ in range(3))
        for i in range(n):
            run_reference(blk, st, pa, 1)
            x[i], h[i], y[i] = float(st["v.x"]), float(st["v.h"]), float(st["v.y"])
        tail = slice(n - 4000, n)  # whole hold periods, transient gone (≫ 1/(2π fc))
        return np.mean((x - h)[tail]), np.mean((x - y)[tail])

    hold_zoh, lag_zoh = run("zoh")
    hold_bil, lag_bil = run("bilinear")
    assert hold_zoh == pytest.approx((r - 1) / 2, abs=1e-9) and hold_bil == pytest.approx((r - 1) / 2, abs=1e-9)
    tau_zoh = 1.0 / (1.0 - np.exp(-2 * np.pi * fc / fs)) - 0.5  # analogue 1/(a h) to O(a h)
    tau_bil = 1.0 / (2.0 * np.tan(np.pi * fc / fs))  # prewarped analogue prototype, no extra delay
    assert lag_zoh - (r - 1) / 2 - tau_zoh == pytest.approx(0.5, abs=2e-3)
    assert lag_bil - (r - 1) / 2 - tau_bil == pytest.approx(0.0, abs=2e-3)


@pytest.mark.parametrize("compiled", [False, True])
def test_recorder_antialias_boxcar_fir(compiled):
    """Phase 5 F1-7: boxcar records the mean of the last decim samples, fir a unit-DC-gain linear-phase
    low-pass (ramp → n − 8·decim after the transient); int slots stay point-sampled; compiled = reference."""
    from cranberry.runner.run import antialias_taps

    d = 3
    out = {}
    for mode in ("boxcar", "fir"):
        rec = Recorder(["J.w", "v.w", "lock.w.state"], decim=d, chunk_samples=50, antialias=mode)
        assert rec.rate_div == 1 and rec.decim == d
        lay = SlotLayout(1.0e3)
        blk = build(Seq([Writer(), rec]), lay, 0)
        st, p = State(lay), Params(lay)
        tables = []
        for _ in range(4):
            (run_compiled if compiled else run_reference)(blk, st, p, 50)
            tables.append(rec.harvest(st))
        out[mode] = rec.split(np.concatenate(tables), lay)
    n = np.arange(0, 200, d)
    bx = out["boxcar"]
    np.testing.assert_allclose(bx["v.w"][1:, 0], n[1:] - (d - 1) / 2, rtol=0, atol=1e-12)
    np.testing.assert_allclose(bx["v.w"][1:, 1], -(n[1:] - (d - 1) / 2), rtol=0, atol=1e-12)
    np.testing.assert_allclose(bx["J.w"][1:, 1, 1], (n[1:] - (d - 1) / 2) * (1 + 1j), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(bx["lock.w.state"], 2 * n)
    ntaps = len(antialias_taps("fir", d))
    assert ntaps == 16 * d + 1
    late = n >= ntaps
    np.testing.assert_allclose(out["fir"]["v.w"][late, 0], n[late] - 8 * d, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(out["fir"]["lock.w.state"], 2 * n)


def test_recorder_antialias_suppresses_folding():
    """A tone just below f_rec = f_sim/decim folds to near DC with full amplitude under point sampling;
    boxcar (sinc zero at f_rec) and fir attenuate it."""
    import cranberry.core.testing  # noqa: F401
    from cranberry.core.pipeline import from_spec, spec_from_obj
    from cranberry.core.testing import toy

    fs, d = 1.0e4, 10
    amp = {}
    for mode in ("none", "boxcar", "fir"):
        rec = Recorder(["v.s"], decim=d, chunk_samples=20000, antialias=mode)
        lay = SlotLayout(fs)
        sine = from_spec(spec_from_obj(toy("toy_sine", "s", out="v.s", unit="V", amp=1.0, freq_hz=0.98 * fs / d)))
        root = Seq([sine, rec])
        blk = build(root, lay, 0)
        st, p = State(lay), Params(lay)
        run_compiled(blk, st, p, 20000)
        y = rec.split(rec.harvest(st), lay)["v.s"][400:]
        amp[mode] = np.sqrt(2) * np.std(y)
    assert amp["none"] == pytest.approx(1.0, rel=0.02)
    assert amp["boxcar"] < 0.05 and amp["fir"] < 0.02


@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("mode", ["none", "boxcar", "fir"])
def test_recorder_cos_sin2_formed_before_decimation(mode, compiled):
    """Phase 5 M2C-01: ``cos_sin2`` records cos 2x / sin 2x formed every base sample and then filtered,
    so the mean of the recorded cos 2x is the full-rate ⟨cos 2x⟩ (J0(2a) for x = a sin) whatever the
    anti-alias mode, while cos 2(filtered x) would lose the fast jitter (boxcar: ≈ 1)."""
    import cranberry.core.testing  # noqa: F401
    from scipy.special import j0

    from cranberry.core.pipeline import from_spec, spec_from_obj
    from cranberry.core.testing import toy
    from cranberry.runner.run import antialias_taps

    fs, d, a, n = 1.0e4, 10, 0.3, 20000
    rec = Recorder(["v.s"], decim=d, chunk_samples=n, antialias=mode, cos_sin2=["v.s"])
    lay = SlotLayout(fs)
    sine = from_spec(spec_from_obj(toy("toy_sine", "s", out="v.s", unit="V", amp=a, freq_hz=0.98 * fs / d)))
    blk = build(Seq([sine, rec]), lay, 0)
    st, p = State(lay), Params(lay)
    (run_compiled if compiled else run_reference)(blk, st, p, n)
    out = rec.split(rec.harvest(st), lay)
    assert set(out) == {"v.s", "v.s.cos2", "v.s.sin2"}
    x = a * np.sin(2 * np.pi * 0.98 * fs / d * np.arange(n) / fs)
    taps = antialias_taps(mode, d)
    full_c, full_s = np.cos(2 * x), np.sin(2 * x)
    if taps is not None:  # causal FIR, zero history before n = 0
        full_c, full_s = (np.convolve(y, taps)[:n] for y in (full_c, full_s))
    np.testing.assert_allclose(out["v.s.cos2"], full_c[::d], atol=1e-9)
    np.testing.assert_allclose(out["v.s.sin2"], full_s[::d], atol=1e-9)
    assert np.mean(out["v.s.cos2"][100:]) == pytest.approx(j0(2 * a), abs=5e-3)
    with pytest.raises(ValueError, match="not recorded"):
        Recorder(["v.s"], decim=d, chunk_samples=n, cos_sin2=["v.t"])


def test_record_without_antialias_warns_once_when_decimating():
    import warnings

    import cranberry.core.testing  # noqa: F401
    import cranberry.runner.run as rr
    from cranberry.config.spec import RunSettings, RunSpec
    from cranberry.core.pipeline import spec_from_obj
    from cranberry.core.testing import toy

    def spec(rec):
        return RunSpec(run=RunSettings(duration_s=100 / 1e3, f_sim_hz=1e3, chunk_samples=64, seed=1),
                       pipeline=spec_from_obj([toy("toy_sine", "s", out="v.s", unit="V", amp=1.0, freq_hz=10.0)]),
                       measure={"record": rec})

    rr._ALIAS_WARNED = False
    with pytest.warns(rr.RecordAliasWarning):
        a = rr.run(spec({"slots": ["v.s"], "decim": 4}))
    with warnings.catch_warnings():
        warnings.simplefilter("error", rr.RecordAliasWarning)
        b = rr.run(spec({"slots": ["v.s"], "decim": 4}))  # once per process
        rr._ALIAS_WARNED = False
        c = rr.run(spec({"slots": ["v.s"], "decim": 4, "antialias": "none"}))  # explicit: no warning
        rr.run(spec({"slots": ["v.s"], "decim": 1}))  # no decimation: no warning
    np.testing.assert_array_equal(a.series["v.s"], c.series["v.s"])  # backward compatible
    np.testing.assert_array_equal(a.t_s, b.t_s)
    bx = rr.run(spec({"slots": ["v.s"], "decim": 4, "antialias": "boxcar"}))
    assert bx.series["v.s"].shape == a.series["v.s"].shape and not np.array_equal(bx.series["v.s"], a.series["v.s"])
