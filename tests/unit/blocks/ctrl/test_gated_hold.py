"""``gated_hold`` (WP6-G): zero / hold / noise_only outside the window, settle_s blanking, toggle, reference = kernel."""

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.ctrl.gated_hold import make_gated_hold
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba
FS = 1.0e6


def pipe(mode, settle_s=0.0, gating=True, **kw):
    p = {"in": "v.x", "gate": "gate.gw.lock", "out": "v.y", "unit": "V", "mode": mode, "settle_s": settle_s,
         "toggles": {"gating": gating}, **kw}
    g = {"block": "gate_schedule", "name": "gw",
         "params": {"detector_gate": "complement", "write_phase": False, "period_s": 20e-6, "lock_fraction": 0.4,
                    "phase": 0.0, "edge_s": 0.0, "timing_jitter_s": 0.0, "leak": 0.0,
                    "toggles": {"gating": True, "jitter": False, "edges": False, "leak": False}}}
    return [toy("toy_sine", "src", out="v.x", unit="V", amp=1.0, freq_hz=3.1e4), g,
            {"block": "gated_hold", "name": "gh", "params": p}, pt.record(["v.x", "v.y", "gate.gw.lock"], 200)]


def run(mode, compiled=False, **kw):
    import cranberry.core.testing  # noqa: F401  (registers toy_sine)

    rec, _, _ = pt.run_and_record(pipe(mode, **kw), FS, 200, compiled=compiled, seed=5)
    return rec["v.x"], rec["v.y"], rec["gate.gw.lock"] > 0.5


def test_zero_and_hold():
    x, y, g = run("zero")
    np.testing.assert_array_equal(y[g], x[g])
    np.testing.assert_array_equal(y[~g], 0.0)
    x, y, g = run("hold")
    np.testing.assert_array_equal(y[g], x[g])
    last = np.flatnonzero(g[:-1] & ~g[1:])  # last active sample of each window
    for i in last:
        j = i + 1
        while j < len(g) and not g[j]:
            assert y[j] == x[i]
            j += 1


def test_settle_blanking_after_each_edge():
    x, y, g = run("zero", settle_s=3e-6)
    act = np.zeros_like(g)
    for k in range(len(g)):
        start = k
        while start > 0 and g[start - 1] == g[k]:
            start -= 1
        act[k] = g[k] and (k - start >= 3 or (start == 0 and k >= 3))
    np.testing.assert_array_equal(y[act], x[act])
    np.testing.assert_array_equal(y[~act], 0.0)
    assert act.sum() == g.sum() - 3 * 10  # 8 active samples per 20-sample period instead of 8 + 3


def test_noise_only_and_gating_off():
    x, y, g = run("noise_only", noise_v_rms=0.5)
    np.testing.assert_array_equal(y[g], x[g])
    assert 0.3 < np.std(y[~g]) < 0.7 and np.all(y[~g] != 0.0)
    x, y, _ = run("noise_only", gating=False, noise_v_rms=0.5)
    np.testing.assert_array_equal(y, x)


@pytest.mark.parametrize("mode", ["zero", "hold", "noise_only"])
def test_reference_equals_kernel(mode):
    kw = {"noise_v_rms": 0.2} if mode == "noise_only" else {}
    a = run(mode, compiled=False, settle_s=2e-6, **kw)
    b = run(mode, compiled=True, settle_s=2e-6, **kw)
    np.testing.assert_allclose(a[1], b[1], rtol=0, atol=1e-12)


def test_validation():
    base = {"in": "v.x", "gate": "gate.g.lock", "out": "v.y", "unit": "V", "settle_s": 0.0, "toggles": {"gating": True}}
    with pytest.raises(ValueError, match="mode"):
        make_gated_hold({**base, "mode": "sample"})
    with pytest.raises(KeyError):
        make_gated_hold({**base, "mode": "noise_only"})  # noise_v_rms required (D34)
    with pytest.raises(ValueError, match="settle_s"):
        make_gated_hold({**base, "mode": "zero", "settle_s": -1.0})
    blk = make_gated_hold({**base, "mode": "hold"}, name="gh")
    assert blk.reads == ("v.x", "gate.g.lock") and blk.runtime_params == ("settle_s",)
