"""WP6-G gate overlays on the future presets: P2 Borealis window, P3 chopper (with the WDM probe overlays), P8 intermittent
slow loop. Every overlay builds on every base (the CHK blocks run); one F_* nominal per overlay runs the time loop
through several windows: the probe is off and the controller output frozen outside the lock window, the metrics count
only the measurement window, and the gated stability check (D80-2: loop_unstable → valid False) is present."""

from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_spec = importlib.util.spec_from_file_location("future_support", Path(__file__).with_name("future_support.py"))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

BASES = ("F_wg100", "F_wg6", "F_bor6", "F_liu20")
OVERLAYS = {
    "borealis": ("gate/borealis_window", "gate/borealis_si"),
    "chopper": ("probe/sep_wdm", "probe/freq_1THz", "gate/chopper_wdm", "gate/chopper_p3"),
    "intermittent": ("gate/intermittent_slow",),
}
REC = ["ctrl.sup1.u", "v.liu1", "v.adc1", "delay.gd1.out", "gate.gw.lock", "gate.gw.meas", "lock.L1.err"]


@pytest.mark.parametrize("base", BASES)
@pytest.mark.parametrize("which", sorted(OVERLAYS))
def test_overlay_builds_on_every_base(base, which):
    r = fs.build_check(fs.nominal(base), OVERLAYS[which])
    chk = fs.block(r, "chk.gated_stability")
    assert chk.result is not None and chk.result["n_active"] > 0
    assert r.valid == chk.result["stable"]
    assert ("loop_unstable" in fs.chk_codes(r)) == (not r.valid)


def _run(base, overlays, n):
    spec = fs.load(fs.nominal(base), params={"t_settle_s": 0.0}, run={"duration_s": n * 1e-6},
                   overlays=overlays, measure={"record": {"slots": REC, "decim": 1, "antialias": "none"}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(fs.first_point(spec), mode="reference")


@pytest.mark.parametrize("base, overlays, held", [
    ("F_liu20", OVERLAYS["borealis"], True),
    ("F_bor6", OVERLAYS["chopper"], False),                              # P3: low pass, no hold
    ("F_wg6", (*OVERLAYS["intermittent"], "gate/borealis_si"), True),   # a window short enough for a reference run
])
def test_nominal_time_loop_with_overlay(base, overlays, held):
    r = _run(base, overlays, 450)
    s = {k: np.asarray(v, dtype=float) for k, v in r.series.items()}
    p = fs.params(r)
    lock, meas, g_det = s["gate.gw.lock"], s["gate.gw.meas"] > 0.5, s["delay.gd1.out"]
    assert meas.any() and (~meas).any() and np.all(lock[meas] <= p["gate_leak"] + 1e-15)  # probe off in the window
    # the probe power reaches det_k with the flight delay; with the probe off liu_k carries no fringe
    k_v = p["liu_k_v"]
    off = g_det < 1e-3
    assert off.any() and np.max(np.abs(s["v.liu1"][off])) < 0.05 * k_v
    on = g_det > 0.999
    assert np.max(np.abs(s["v.liu1"][on])) > 0.5 * np.max(np.abs(s["v.liu1"][off]))
    if held:
        # zero-mode hold + I-only tustin pid: frozen output from the second sample after the controller gate closes
        frozen = off & np.r_[False, False, off[:-2]] & np.r_[False, off[:-1]]
        du = np.abs(np.diff(s["ctrl.sup1.u"]))[frozen[1:]]
        assert du.size and np.all(du == 0.0)
    mi = r.metrics.metric_ii[1]
    t = np.asarray(r.t_s)
    assert mi.stats.n == int(meas.sum()) == mi.theta_rad.size
    assert r.valid and fs.block(r, "chk.gated_stability").result["stable"]
    assert t.size == 450


def test_unstable_gate_point_is_flagged_not_pruned(tmp_path):
    """The bad configuration (D34) of holding the last ERROR (gated_hold hold mode: the integrator keeps ramping) over a
    long hold, on F_liu20: the gated check emits loop_unstable, RunResult.valid is False and the point still runs (D75,
    D80-2). λ ≈ e^{−ω_u T_on}(1 − ω_u T_off) = −16 for f_u = 3 kHz, T_on = 40 µs, T_off = 1.96 ms."""
    hold = tmp_path / "hold.yaml"
    hold.write_text("offline:\n  params: {gate_period_s: 2.0e-3, gate_lock_fraction: 0.02, f_u_hz: 3.0e+3}\n"
                    "  blocks:\n    - name: chk.gated_stability\n      params: {hold_mode: hold}\n"
                    "pipeline:\n  - name: adc1\n    seq:\n      - name: gh1\n        params: {mode: hold}\n")
    spec = fs.load(fs.nominal("F_liu20"), run={"duration_s": 3e-6}, overlays=(*OVERLAYS["borealis"][:1], str(hold)),
                   pipeline=[{"name": "mi", "params": {"t_settle_s": 0.0}}, {"name": "mii", "params": {"t_settle_s": 0.0}}])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run(fs.first_point(spec), mode="reference")
    chk = fs.block(r, "chk.gated_stability")
    assert not chk.result["stable"] and chk.result["multipliers"][0] > 10.0
    assert "loop_unstable" in fs.chk_codes(r) and r.valid is False
