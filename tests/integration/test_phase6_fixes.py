"""Phase 6 fixes of reviews/phase6_wp6g_ab.md on the real presets: masked recorded series (G-1), points without a
measurement sample and per-point sweep errors (G-2, D75), decompose validity on an empty window, and D91-1 (DSC1 at
500 kHz on F_liu20 + a gate overlay: the gated stability check against the time loop)."""

from __future__ import annotations

import importlib.util
import types
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry.blocks.chk.gated_loop_stability import make_gated_loop_stability_check
from cranberry.runner.decompose import validity_violations
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand, sweep

_spec = importlib.util.spec_from_file_location("future_support", Path(__file__).with_name("future_support.py"))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

pytestmark = [pytest.mark.integration, pytest.mark.numba]


def gated(period, lf, record, duration=4e-4, overlays=("gate/borealis_window",), extra=None):
    s = fs.load(fs.nominal("F_liu20"), params={"t_settle_s": 0.0}, run={"duration_s": duration}, overlays=overlays,
                measure={"record": record})
    axes = {"offline.params.gate_period_s": [period], "offline.params.gate_lock_fraction": [lf], **(extra or {})}
    (pt,) = expand(s, axes)
    return pt.spec


def codes(r):
    return [e["code"] for e in r.events]


def go(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(spec, mode="reference")


def test_commensurate_decimation_warns_instead_of_empty_series():
    """G-1 case A: T = decim samples → the point-sampled mask is one window phase: CHK mask_series, moments exact."""
    r = go(gated(1e-5, 0.65, {"slots": ["lock.L1.err"], "decim": 10, "antialias": "none"}))
    assert "mask_series" in codes(r) and r.valid
    m = r.metrics.metric_ii[1]
    assert m.stats.n > 0 and m.theta_rad.size == 0


def test_fir_mask_never_mixes_windows():
    """G-1 case B: with the 161-tap FIR a recorded mask is 1 only if all contributing samples are in the window; the
    35-sample window is shorter than the filter, so nothing is kept — and that is reported, not hidden."""
    r = go(gated(1e-4, 0.65, {"slots": ["lock.L1.err"], "decim": 10, "antialias": "fir"}, duration=6e-4))
    m = np.asarray(r.series["gate.gw.meas"], float)
    assert set(np.unique(m)) <= {0.0, 1.0}
    assert "mask_series" in codes(r)
    box = go(gated(1e-4, 0.65, {"slots": ["lock.L1.err"], "decim": 5, "antialias": "boxcar"}, duration=6e-4))
    mb = np.asarray(box.series["gate.gw.meas"], float)
    assert set(np.unique(mb)) <= {0.0, 1.0} and mb.any()


def test_point_without_measurement_sample_is_kept_invalid():
    """G-2: T = 10 µs, D_lock 0.9 at 1 MHz: no sample in the measurement window → NaN moments, metric_no_sample
    (invalidating), the gated check's gate_window_short; no exception."""
    r = go(gated(1e-5, 0.9, {"slots": ["lock.L1.err"], "decim": 1, "antialias": "none"}))
    assert "metric_no_sample" in codes(r) and "gate_window_short" in codes(r) and not r.valid
    assert np.isnan(r.metrics.metric_ii[1].stats.std)


def test_sweep_keeps_failing_points():
    """G-2 / D75: runner.sweep records a point that raises as an invalid run_error result and continues."""
    s = fs.load(fs.nominal("F_liu20"), params={"t_settle_s": 0.0}, run={"duration_s": 3e-6, "seed_mode": "common"},
                overlays=("gate/borealis_window",))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = sweep(s, {"offline.params.gate_period_s": [1e-4, -1.0], "offline.params.gate_lock_fraction": [0.65]})
    assert len(res) == 2 and res[0].metrics is not None
    assert codes(res[1]) == ["run_error"] and not res[1].valid
    assert res[1].events[0]["context"]["sweep_point"]["offline.params.gate_period_s"] == -1.0


def test_decompose_validity_empty_window_is_invalid():
    r = types.SimpleNamespace(series={"m": np.zeros(10), "flag": np.zeros(10)}, final_state=None)
    bad = validity_violations(r, [{"slot": "flag", "abs_max": 0.5}], "m")
    assert bad and "cannot be tested" in bad[0]
    r.series["m"][3] = 1.0
    assert validity_violations(r, [{"slot": "flag", "abs_max": 0.5}], "m") == []


# --------------------------------------------------------------------------------------------- D91-1 verification

D91 = ("controller/dsc1", "gate/borealis_window")
QUIET = {f"{b}{k}.{t}": False for k in "01" for b, t in (("adc", "quantization"), ("liu", "shot_noise"),
                                                         ("liu", "electronic_noise"))}
QUIET.update({f"dac{k}f.quantization": False for k in "01"})


def _d91_spec(duration):
    return gated(1e-4, 0.65, {"slots": ["lock.L1.err", "gate.gw.meas"], "decim": 1, "antialias": "none"},
                 duration=duration, overlays=(*D91, "gate/borealis_si"), extra={"offline.params.dsc1_latency_s": [5.7e-6]})


def _chk_params(spec):
    (cp,) = [b.params for b in offline_block_specs(spec) if b.name == "chk.gated_stability"]
    return dict(cp)


def test_d91_dsc1_gate_check_models_the_time_loop():
    """F_liu20 + controller/dsc1 + gate/borealis_window(+si): pid / sup / dac at rate_div 2 (h = 2 µs), 120 / 100 kHz
    poles, τ_L 5.7 µs accounted as D61-2 (delay lines τ_L − 1/f_sim − 1 µs); the gated check models the same law, rate,
    poles, modes and delay (τ_L minus the hold it models itself) and the time loop locks."""
    spec = _d91_spec(3e-3)
    p = spec.offline["params"]
    assert p["pid_rate_div"] == 2 and p["pid_h_s"] == pytest.approx(2e-6)
    assert p["ctrl_explicit_delay_s"] == pytest.approx(1e-6)  # ADC average 0.5 µs + hold 0.5 µs
    cp = _chk_params(spec)
    assert cp["pid_discretisation"] == "dsc1" and cp["pid_rate_div"] == 2 and cp["poles_hz"][:2] == [1.2e5, 1.0e5]
    assert cp["delay_s"] == pytest.approx(5.7e-6 - 0.5e-6)
    chk = make_gated_loop_stability_check(cp, name="chk")
    assert chk.result["stable"]
    r = go_compiled(spec)
    assert r.valid and np.all(np.isfinite(np.asarray(r.series["lock.L1.err"], float)))


def go_compiled(spec, toggles=None):
    import dataclasses

    spec = dataclasses.replace(spec, toggles={**spec.toggles, **(toggles or {})})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(spec, mode="compiled")


@pytest.mark.slow
def test_d91_measurement_window_variance_matches_time_loop():
    """The check's var_meas with the free-running PSD of the same seed vs the time loop's masked variance (noise and
    quantisation off; reviews/scripts/p6_wp6g_timeloop.py §2 method). Measured at 0.1 s: ratio 1.023 (meas), 0.985
    (lock window)."""
    from scipy.signal import welch

    spec = _d91_spec(0.1)
    lk = go_compiled(spec, QUIET)
    fr = go_compiled(spec, {**QUIET, **{f"pid{k}{s}.feedback": False for k in "01" for s in "fs"}})
    keep = np.asarray(lk.t_s) >= 0.02
    f, s = welch(np.asarray(fr.series["lock.L1.err"], float)[keep], 1e6, nperseg=1 << 17, detrend="linear")
    res = make_gated_loop_stability_check({**_chk_params(spec), "disturbance": {"f_hz": f[1:].tolist(),
                                                                               "psd_rad2_per_hz": s[1:].tolist()}},
                                          name="c").result
    e = np.asarray(lk.series["lock.L1.err"], float)[keep]
    meas = np.asarray(lk.series["gate.gw.meas"], float)[keep] > 0.5
    e = e - e[meas].mean()
    assert np.var(e[meas]) == pytest.approx(res["var_meas_rad2"], rel=0.1)
    assert np.var(e[~meas]) == pytest.approx(res["var_lock_rad2"], rel=0.1)
