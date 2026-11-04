"""Phase 6 review (WP6-G D90-2, D85-3): `mask_slot` with a decimated recorder.

The in-kernel moments of metric_i / metric_ii use every base sample with mask > 0.5 (correct). The recorded series
(MetricI/II series, the zero-span cos2/sin2 series, stray-light series, decompose validity) are masked with the
RECORDED mask, which goes through the same decimation / anti-alias filter as the data:

  A. F_common's own record (decim 10, antialias none) and a gate period of 10 base samples (gate_period_s = 1e-5, the
     lower end of gate/borealis_window's axis): the recorded mask is sampled at one fixed position of the window
     every period, so the masked series is either empty or the whole record, whatever the window.
  B. antialias boxcar (decim 10, T = 100 µs, D_lock 0.65): a recorded value averages 10 base samples; around the
     window edges it mixes lock-window samples in, and the boxcar-filtered mask (> 0.5) keeps such values.

Run: uv run python reviews/scripts/p6_wp6g_mask_decim.py
"""

from __future__ import annotations

import warnings

import numpy as np

from cranberry import config
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

NOMINAL = "presets/future/F_liu20_nominal.yaml"


def spec(period, lf, record, duration=4e-3):
    tree = {"extends": NOMINAL, "run": {"duration_s": duration}, "offline": {"params": {"t_settle_s": 0.0}}}
    if record:
        tree["measure"] = {"record": record}
    s = config.load(tree, overlays=["gate/borealis_window"])
    (pt,) = expand(s, {"offline.params.gate_period_s": [period], "offline.params.gate_lock_fraction": [lf]})
    return pt.spec


def report(label, r):
    mi = r.metrics.metric_ii[1]
    ser = np.asarray(r.series.get("gate.gw.meas", []), float)
    print(f"  {label}: recorded mask: {ser.size} values, fraction > 0.5 = {np.mean(ser > 0.5) if ser.size else float('nan'):.3f},"
          f" distinct {np.unique(np.round(ser, 3))[:6]} | in-kernel moments n = {mi.stats.n}, std {mi.stats.std:.4e} rad |"
          f" masked series n = {mi.theta_rad.size}, std {np.std(mi.theta_rad) if mi.theta_rad.size else float('nan'):.4e} rad")


def main():
    warnings.simplefilter("ignore")
    ensure_loaded()
    print("A. F_common record (decim 10, antialias none), T = 10 µs (= decim):")
    for lf in (0.3, 0.65, 0.9):
        try:
            report(f"D_lock {lf}", run(spec(1e-5, lf, None), mode="compiled"))
        except Exception as e:  # noqa: BLE001  (a grid point of gate/borealis_window: no sample in the window)
            print(f"  D_lock {lf}: RUN RAISES {type(e).__name__}: {e} (runner.sweep has no per-point handling)")
    print("   same, decim 1 (reference):")
    rec1 = {"slots": ["lock.L1.err"], "decim": 1, "antialias": "none"}
    report("D_lock 0.65", run(spec(1e-5, 0.65, rec1), mode="compiled"))
    print("B. T = 100 µs, D_lock 0.65: decim 10 boxcar vs decim 1")
    for rec in ({"slots": ["lock.L1.err"], "decim": 1, "antialias": "none"},
                {"slots": ["lock.L1.err"], "decim": 10, "antialias": "boxcar"},
                {"slots": ["lock.L1.err"], "decim": 10, "antialias": "fir"}):
        r = run(spec(1e-4, 0.65, rec), mode="compiled")
        report(f"decim {rec['decim']} {rec['antialias']}", r)
        m = np.asarray(r.series["gate.gw.meas"], float)
        part = (m > 0.5) & (m < 1.0 - 1e-12)
        print(f"     recorded mask values in (0.5, 1): {int(part.sum())} of {int((m > 0.5).sum())} kept samples "
              "(each averages lock-window samples into the measurement window)")


if __name__ == "__main__":
    main()
