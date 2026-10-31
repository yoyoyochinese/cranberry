"""VT-33 ESA parts (WP-B5): T04 §12 V7, V8, V9 and the swept spectrum of §5.5.

Sources: notes/phase-noise-to-squeezing.md (T04) §5, §9.1, §12; 00_summary VT-33 and C24 (scanned-trace
minimum: only monotonicity in VBW is required, the depth of Ref A Fig. 3 is not). Ref A source:
L = 8 %, α = 9.06 /W, P = 0.64 W (prompt §7-1).

Finding (WP-B5, for the Phase 5 review): the log-average trace scatter of T04 V7 (0.093–0.099 dB) uses
the π²/6 = 1.645 variance ratio of independent samples (AN 1303). For a continuous RBW-filtered
signal averaged by a VBW ≪ RBW video filter the exact ratio is ∫Li₂(|ρ|²)dτ/∫|ρ|²dτ = 1.325 (rect RBW)
/ 1.341 (Gaussian), so log scatter is 0.083–0.089 dB (rect, ENBW factor 1.128–1.0); B3's radiometer MC
found the same 1.32 (``blocks/err/qnl.py``). The power-average range 0.072–0.077 dB holds for a
rectangular RBW (B_stat = ENBW); a Gaussian RBW has B_stat = √2·ENBW (0.0647 dB). The tests check the
Monte Carlo against the exact expressions and the power range against T04. Phase 5 (D65-2, M1-A-06):
targets per RBW shape — rect power 0.0725–0.0770 / log 0.0834–0.0886 dB, Gaussian (refA's shape) power
0.0609–0.0647 / log 0.0706–0.0750 dB (ENBW factor 1.128–1.0).
"""

import math

import numpy as np
import pytest

import cranberry.blocks.det.testing as T
from cranberry.blocks.meas import esa
from cranberry.blocks.meas.phase_to_quadrature import phase_to_quadrature
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.core.rng import stream
from cranberry.core.units import lin_to_db

pytestmark = pytest.mark.vt

L, R = 0.08, math.sqrt(9.06 * 0.64)
VM, VP = 0.5 * (L + (1 - L) * math.exp(-2 * R)), 0.5 * (L + (1 - L) * math.exp(2 * R))  # conv units


def mk(**over):
    return make_block(BlockSpec.from_obj(T.zero_span(**over)))


def _pooled_std(blk, n, seed):
    rng = stream(seed, blk.name)
    return math.sqrt(np.mean([blk.evaluate(0.5, rng=rng)["raw_db"].var() for _ in range(n)]))


@pytest.mark.parametrize("k", [1.0, 1.128])
def test_v7_radiometer_scatter_rect(k):
    """RBW 1 MHz / VBW 100 Hz: power 0.072–0.077 dB (T04 V7); MC = exact expression (±5 %, ≈3.5σ)."""
    for avg in ("power", "log"):
        want = esa.radiometer_std_db(1e6, 100.0, k, "rect", avg)
        if avg == "power":
            assert 0.072 <= round(want, 3) <= 0.077
        else:
            assert want == pytest.approx(0.0886 if k == 1.0 else 0.0834, abs=5e-4)  # not 0.093–0.099 (docstring)
        blk = mk(enbw_factor=k, average=avg, sweep_time_s=0.5, n_points=20000)
        assert _pooled_std(blk, 16, 11) == pytest.approx(want, rel=0.05)


def test_v7_radiometer_scatter_gaussian_rbw():
    blk = mk(rbw_shape="gaussian", oversample=4, sweep_time_s=0.5, n_points=20000)
    want = esa.radiometer_std_db(1e6, 100.0, 1.0, "gaussian", "power")
    assert want == pytest.approx(0.0770 / 2**0.25, rel=1e-3)
    assert _pooled_std(blk, 16, 12) == pytest.approx(want, rel=0.05)


#: D65-2 / M1-A-06 targets per RBW shape, RBW 1 MHz, VBW 100 Hz; the range is ENBW factor 1.128 → 1.0 (dB)
SCATTER_TARGETS = {("rect", "power"): (0.0725, 0.0770), ("rect", "log"): (0.0834, 0.0886),
                   ("gaussian", "power"): (0.0609, 0.0647), ("gaussian", "log"): (0.0706, 0.0750)}


@pytest.mark.parametrize("shape, avg", list(SCATTER_TARGETS))
def test_v7_scatter_targets_per_rbw_shape(shape, avg):
    """D58-4 as corrected by D65-2: the targets differ by RBW shape (refA uses the Gaussian one); the
    log-average range 0.083–0.089 dB is the rectangular value only. Gaussian with −3 dB width = RBW
    (ENBW = 1.0645·RBW, refA, M1-A-05): power 0.0627 dB, log 0.0727 dB."""
    lo, hi = SCATTER_TARGETS[shape, avg]
    assert esa.radiometer_std_db(1e6, 100.0, 1.128, shape, avg) == pytest.approx(lo, abs=5e-5)
    assert esa.radiometer_std_db(1e6, 100.0, 1.0, shape, avg) == pytest.approx(hi, abs=5e-5)
    if shape == "gaussian":
        k3db = math.sqrt(math.pi / (4 * math.log(2)))
        assert esa.radiometer_std_db(1e6, 100.0, k3db, shape, avg) == pytest.approx(
            {"power": 0.0627, "log": 0.0727}[avg], abs=5e-5)


def test_v7_radiometer_scatter_gaussian_rbw_log():
    """Gaussian-RBW log-average scatter: MC = exact expression (0.0750 dB at ENBW = RBW). Oversample 8:
    at 4 the coarse sampling of the Li₂ cusp biases the log variance by +1.3 % (D65-2)."""
    blk = mk(rbw_shape="gaussian", average="log", oversample=8, sweep_time_s=0.5, n_points=20000)
    want = esa.radiometer_std_db(1e6, 100.0, 1.0, "gaussian", "log")
    assert want == pytest.approx(0.0750, abs=5e-5)
    assert _pooled_std(blk, 16, 13) == pytest.approx(want, rel=0.05)


def test_v7_log_bias_and_ratio_independent_of_average():
    """Log average 2.51 ± 0.05 dB below power average; δ = 0: squeezed/shot ratio independent of the
    average type (< 0.01 dB; common random numbers for the two chains)."""
    raw = {a: mk(average=a, sweep_time_s=1.0).evaluate(0.5, rng=stream(3, "esa"))["raw_db"].mean()
           for a in ("power", "log", "voltage")}
    assert raw["log"] - raw["power"] == pytest.approx(-2.51, abs=0.05)
    assert raw["voltage"] - raw["power"] == pytest.approx(-1.05, abs=0.05)
    out = {a: mk(average=a, sweep_time_s=1.0).evaluate(VM, rng=stream(4, "esa"))["mean_db"] for a in ("power", "log")}
    assert out["log"] - out["power"] == pytest.approx(0.0, abs=0.01)
    assert out["power"] == pytest.approx(lin_to_db(VM / 0.5), abs=0.02)


def _slow_jitter(sigma, n, seed):
    """Gaussian δ held for 5 µs (slower than 1/RBW = 1 µs) on the 2 MHz ESA grid → V′(t) via
    phase_to_quadrature (F4). Held explicitly: the ESA interpolates V′ linearly between samples."""
    d = np.repeat(stream(seed, "vt33.delta").standard_normal(n) * sigma, 10)
    return phase_to_quadrature(d, VM, VP, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0)["v_prime_t"]


@pytest.mark.parametrize(("sigma", "diff", "noise"), [(0.009, 0.033, False), (0.030, 0.80, False),
                                                      (0.030, 0.80, True)])
def test_v8_log_bias_of_slow_jitter(sigma, diff, noise):
    """Ref A source, δ slower than RBW: log − power average = 0.033 / 0.80 dB at 9 / 30 mrad (±10 %)."""
    v = _slow_jitter(sigma, 200_000, 5)
    res = {a: mk(average=a, noise=noise, sweep_time_s=1.0).evaluate(v, dt_s=0.5e-6, rng=stream(6, "esa"))["mean_db"]
           for a in ("power", "log")}
    assert res["power"] - res["log"] == pytest.approx(diff, rel=0.10)
    if not noise and sigma == 0.009:
        assert res["power"] == pytest.approx(-10.149, abs=0.01) and res["log"] == pytest.approx(-10.181, abs=0.01)


def test_v9_phase_scanned_trace():
    """Δ = ω t, ω = π/44 ms (Ref A Fig. 3 read-off): maximum 20.5 dB; minimum deepens monotonically with
    VBW (T04 calc: +1.76 / −3.19 dB at 100 Hz, −9.9 at 1 kHz, −10.6 at 10 kHz); C24 open."""
    dt = 1e-6
    t = np.arange(int(0.15 / dt)) * dt
    d = math.pi / 44e-3 * t - 1.0
    v = VM * np.cos(d) ** 2 + VP * np.sin(d) ** 2
    mins = {}
    for avg in ("power", "log"):
        mins[avg] = []
        for vbw in (100.0, 1e3, 1e4):
            r = mk(average=avg, vbw_hz=vbw, noise=False, n_points=150000).evaluate(v, dt_s=dt)
            mins[avg].append(r["trace_db"].min())
            assert r["trace_db"].max() == pytest.approx(20.5, abs=0.1)
        assert mins[avg][0] > mins[avg][1] > mins[avg][2]
    np.testing.assert_allclose(mins["power"], [1.76, -9.9, -10.6], atol=0.1)
    assert mins["log"][0] == pytest.approx(-3.19, abs=0.1)
    assert mins["power"][2] == pytest.approx(lin_to_db(VM / 0.5), abs=0.05)


def test_swept_spectrum_with_frontend():
    """Ref A Fig. 4 analogue: BHD with dip and circuit noise C(f) = |H|²C_0; (b) dark-subtracted spectrum
    returns R, (a) unsubtracted = (C R + 1)/(C + 1) per bin; bin scatter = √2·radiometer σ (two traces)."""
    fe = make_block(BlockSpec.from_obj(T.frontend(clearance_lin=30.0, toggles={"imbalance_leak": False})))
    r = 0.25
    sig = lambda f: fe.psd(f, r)["signal"]  # noqa: E731
    shot = lambda f: fe.psd(f, 1.0)["shot"]  # noqa: E731
    dark = lambda f: fe.psd(f, 1.0)["dark"]  # noqa: E731
    sa = make_block(BlockSpec.from_obj(T.swept(n_sub=1, f_stop_hz=1.2e9)))  # n_sub = 1: no RBW smoothing
    out = sa.evaluate(sig, shot, dark)
    c = fe.clearance(out["f_hz"])
    np.testing.assert_allclose(10 ** (out["normalized_db"] / 10), (c * r + 1) / (c + 1), rtol=1e-12)
    np.testing.assert_allclose(10 ** (out["corrected_db"] / 10), r, rtol=1e-12)
    noisy = make_block(BlockSpec.from_obj(T.swept(noise=True, n_points=4000)))
    o = noisy.evaluate(sig, shot, dark, rng=stream(9, "sa"))
    sd = esa.radiometer_std_db(1e6, 100.0, 1.0, "rect", "power")
    resid = o["normalized_db"] - lin_to_db((c_n := fe.clearance(o["f_hz"])) * r / (c_n + 1) + 1 / (c_n + 1))
    assert resid.std() == pytest.approx(math.sqrt(2) * sd, rel=0.1)
    low = o["f_hz"] < 0.3e9  # (b) unbiased where the clearance is high
    assert np.mean(10 ** (o["corrected_db"][low] / 10)) == pytest.approx(r, rel=0.01)
