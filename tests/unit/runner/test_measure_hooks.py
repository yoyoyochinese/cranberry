"""Runner measurement-stage hooks ``measure.esa`` / ``measure.pulse_sample`` (WP-I1, D58-6).

Toy port-less loop (``blocks.meas.testing``): a sinusoidal LO jitter δ(t) = a sin(2πft) at the BHD, the
Ref A single-mode source, a 5 % tap on the squeezed segment. Checks:

* zero span: V′(t) from the recorded metric series gives the same mean variance as the metric's
  noise levels (full-rate moments) for the locked (0) and orthogonal (π/2) quadratures; with radiometer
  noise off the trace mean equals 10 log10 R (power average, flat V′ over the VBW);
* swept spectrum: R through ``bhd_frontend`` without circuit noise is flat at the metric's R′;
  ``exclude_losses`` re-computes η without the listed elements;
* pulse sample: ``r_slots`` items from metrics / numbers reach ``pulse_sample.evaluate`` unchanged;
* every hook result lands in ``results.npz`` (``measure.<hook>.<entry>.<field>``) when recording;
* D34 errors: missing keys, wrong block kind, record shorter than the sweep, unknown excluded loss.
"""

import math

import numpy as np
import pytest

import cranberry.blocks.det.testing as T
from cranberry.blocks.meas.pulse_sample import expected_r_meas
from cranberry.blocks.meas.testing import metric_cfg, path_blocks, refa_source, run_spec, single_port_topology, toy
from cranberry.core.units import lin_to_db
from cranberry.runner.run import MEASUREMENT_HOOKS, run

pytestmark = pytest.mark.numba
AMP = 0.02  # rad


def spec(esa=None, pulse=None, n=40_000, decim=4, record=True, fs=2e5):
    pipe = [toy("toy_sine", "jit", out="drv.lo_seg.stretcher", unit="1", amp=AMP, freq_hz=1e3),
            toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
            {"block": "tap", "name": "tap1", "params": {"tap_ratio": 0.05, "p_probe_w": 1e-3, "out": "v.tap",
                                                        "segment": "sq_seg", "applies_to": ["gbs", "bhd"],
                                                        "source": "[test]", "toggles": {"enabled": True}}},
            *path_blocks({"las": "laser.las.phi"}),
            {"block": "metric_i", "name": "mi", "params": {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                                                          "setpoint_rad": 0.0}}]
    offline = [refa_source(), T.zero_span("zs", noise=False, sweep_time_s=0.1, n_points=200),
               T.zero_span("zs_noise", noise=True, sweep_time_s=0.1, n_points=200),
               T.swept("sw", noise=False, subtract=True, f_stop_hz=1e8, n_points=50),
               T.frontend("fe", on=False), T.pulsed_kernel("pk"), T.pulse_sample("ps", on=False)]
    m = {"metrics": {"i": metric_cfg("refA.sqz", include_hf_phase_noise=False)}}
    if record:
        m["record"] = {"slots": ["beam.lo.phi"], "decim": decim}
    if esa is not None:
        m["esa"] = esa
    if pulse is not None:
        m["pulse_sample"] = pulse
    return run_spec(pipe, f_sim_hz=fs, n_samples=n, topology=single_port_topology(), offline={"blocks": offline},
                    measure=m)


def test_zero_span_matches_metric_levels():
    esa = {"sq": {"block": "zs", "metric": "i.port0", "quadrature_rad": 0.0},
           "anti": {"block": "zs", "metric": "i.port0", "quadrature_rad": math.pi / 2},
           "sq_noise": {"block": "zs_noise", "metric": "i.port0", "quadrature_rad": 0.0}}
    res = run(spec(esa=esa))
    nl = res.metrics.noise_levels["i.port0"]
    out = res.measurements["esa"]
    assert set(res.measurements) == {"esa"} and set(MEASUREMENT_HOOKS) == {"esa", "pulse_sample"}
    assert 10 * math.log10(out["sq"]["r_lin"]) == pytest.approx(nl["r_minus_db"], abs=2e-3)
    assert 10 * math.log10(out["anti"]["r_lin"]) == pytest.approx(nl["r_plus_db"], abs=2e-3)
    # radiometer noise off, power average, V′ flat over the 100 Hz VBW (1 kHz jitter averaged): mean = R
    assert out["sq"]["mean_db"] == pytest.approx(10 * math.log10(out["sq"]["r_lin"]), abs=0.02)
    assert out["anti"]["mean_db"] == pytest.approx(10 * math.log10(out["anti"]["r_lin"]), abs=0.02)
    assert out["sq_noise"]["std_db"] > 3 * out["sq"]["std_db"]  # radiometer scatter on
    # sinusoidal jitter: χ = J0(2a) (conv §6)
    from scipy.special import j0
    vm, vp = nl["v_minus"], nl["v_plus"]
    expect = 0.5 * (vm + vp) - 0.5 * (vp - vm) * j0(2 * AMP)
    assert out["sq"]["v_prime"] == pytest.approx(expect, rel=2e-3)


def test_swept_spectrum_and_excluded_losses():
    esa = {"flat": {"block": "sw", "metric": "i.port0", "quadrature_rad": 0.0, "frontend": "fe", "exclude_losses": []},
           "no_tap": {"block": "sw", "metric": "i.port0", "quadrature_rad": 0.0, "frontend": "fe",
                      "exclude_losses": ["tap1.tap"]}}
    res = run(spec(esa=esa))
    nl = res.metrics.noise_levels["i.port0"]
    flat, no_tap = res.measurements["esa"]["flat"], res.measurements["esa"]["no_tap"]
    assert flat["eta"] == pytest.approx(0.95) and no_tap["eta"] == 1.0
    np.testing.assert_allclose(flat["normalized_db"], nl["v_prime_db"], atol=1e-9)
    np.testing.assert_allclose(no_tap["normalized_db"], lin_to_db(no_tap["r_lin"]), atol=1e-9)
    # η = 1: the squeezed variance without the tap (pure loss algebra, conv §2)
    v_src = (nl["v_minus"] - 0.05 * 0.5) / 0.95
    assert no_tap["r_lin"] < flat["r_lin"] and no_tap["v_prime"] >= v_src


def test_pulse_sample_hook():
    pulse = {"frame": {"block": "ps", "kernel": "pk",
                       "r_slots": [{"metric": "i.port0", "quadrature_rad": 0.0}] + [1.0] * 31}}
    res = run(spec(pulse=pulse))
    r = res.measurements["pulse_sample"]["frame"]
    nl = res.metrics.noise_levels["i.port0"]
    assert r["r_slots"][0] == pytest.approx(10 ** (nl["v_prime_db"] / 10), rel=1e-12)
    np.testing.assert_allclose(r["expected"], expected_r_meas(r["r_slots"], "pulse_resolved", [0], [1.0], 0.0))


def test_hook_results_in_run_record(tmp_path):
    esa = {"sq": {"block": "zs", "metric": "i.port0", "quadrature_rad": 0.0}}
    res = run(spec(esa=esa), record_root=tmp_path)
    with np.load(res.run_dir / "results.npz") as z:
        keys = set(z.files)
        assert {"measure.esa.sq.trace_db", "measure.esa.sq.mean_db", "measure.esa.sq.r_lin"} <= keys
        assert float(z["measure.esa.sq.mean_db"]) == res.measurements["esa"]["sq"]["mean_db"]


@pytest.mark.parametrize("esa, kw, err", [
    ({"sq": {"block": "zs", "metric": "i.port0"}}, {}, KeyError),                       # quadrature_rad missing
    ({"sq": {"block": "fe", "metric": "i.port0", "quadrature_rad": 0.0}}, {}, KeyError),  # not an ESA block
    ({"sq": {"block": "zs", "metric": "i.port7", "quadrature_rad": 0.0}}, {}, KeyError),  # no such metric
    ({"sq": {"block": "zs", "metric": "i.port0", "quadrature_rad": 0.0}}, {"n": 10_000}, ValueError),  # record too short
    ({"sq": {"block": "zs", "metric": "i.port0", "quadrature_rad": 0.0}}, {"record": False}, ValueError),
    ({"s": {"block": "sw", "metric": "i.port0", "quadrature_rad": 0.0, "frontend": "fe",
            "exclude_losses": ["nope"]}}, {}, KeyError),
    ({"s": {"block": "sw", "metric": "i.port0", "quadrature_rad": 0.0, "frontend": "zs",
            "exclude_losses": []}}, {}, KeyError),
])
def test_hook_errors(esa, kw, err):
    with pytest.raises(err):
        run(spec(esa=esa, **kw))
