"""``liu_probe_interferometer`` optional ``in_gate`` (WP6-G): absent = unchanged, K ∝ √g, probe shot noise ∝ g."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
import cranberry.path.testing as pt
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba
FS = 1.0e6


def pipe(gate=None, noise=False, g_amp=0.0):
    p = {"mode": "homodyne_dc", "in_pr": "beam.pr.phi", "in_ref": "beam.ref.phi", "out": "v.liu", "p_pr_avg_w": 1e-6,
         "p_ref_w": 1e-3, "duty": 0.5, "visibility": 1.0, "coupler_sign": 1.0, "v_off_v": 0.0,
         "responsivity_a_per_w": 1.0, "gain_v_per_a": 1e3, "i_n_a_rthz": 1e-11,
         "toggles": {"shot_noise": noise, "electronic_noise": noise, "gating": True}}
    if gate is not None:
        p["in_gate"] = gate
    blocks = [toy("toy_sine", "x", out="beam.pr.phi", unit="rad", amp=3.0, freq_hz=2.1e4),
              toy("toy_sine", "r", out="beam.ref.phi", unit="rad", amp=0.0, freq_hz=1.0)]
    if gate == "v.g":
        blocks.append(toy("toy_sine", "g", out="v.g", unit="1", amp=g_amp, freq_hz=7e3))
    return blocks + [{"block": "liu_probe_interferometer", "name": "liu", "params": p},
                     pt.record(["v.liu", "beam.pr.phi", *(["v.g"] if gate == "v.g" else [])], 300)]


def test_absent_or_null_gate_is_unchanged():
    a, _, _ = pt.run_and_record(pipe(), FS, 300, seed=2)
    b, _, _ = pt.run_and_record(pipe(gate=None), FS, 300, seed=2)
    np.testing.assert_array_equal(a["v.liu"], b["v.liu"])
    k = 2 * 1e3 * math.sqrt(1e-6 / 0.5 * 1e-3)
    np.testing.assert_allclose(a["v.liu"], k * np.cos(a["beam.pr.phi"]), rtol=1e-12)


@pytest.mark.parametrize("compiled", [False, True])
def test_gate_scales_fringe_by_sqrt_g(compiled):
    rec, _, _ = pt.run_and_record(pipe(gate="v.g", g_amp=1.0), FS, 300, compiled=compiled, seed=2)
    g = np.clip(rec["v.g"], 0.0, None)  # a negative gate counts as 0 (probe off)
    k = 2 * 1e3 * np.sqrt(1e-6 / 0.5 * 1e-3 * g)
    np.testing.assert_allclose(rec["v.liu"], k * np.cos(rec["beam.pr.phi"]), rtol=1e-12, atol=1e-15)
    assert (g == 0.0).sum() > 50


def test_gate_zero_leaves_reference_and_electronic_noise():
    rec, _, _ = pt.run_and_record(pipe(gate="v.g", noise=True, g_amp=0.0), FS, 300, seed=2)
    # g = 0: no fringe; S_v = G²(i_n² + 2eℜP_ref)/D, per-sample std √(S_v f_b/2)
    s_v = 1e6 * (1e-22 + 2 * 1.602176634e-19 * 1e-3) / 0.5
    assert np.std(rec["v.liu"]) == pytest.approx(math.sqrt(s_v * FS / 2), rel=0.15)
    assert abs(np.mean(rec["v.liu"])) < 3 * math.sqrt(s_v * FS / 2) / math.sqrt(300)
