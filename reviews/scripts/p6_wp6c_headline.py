"""Phase 5 review of WP6-C (D92-4): reproduce the F_liu20 headline numbers independently of the author's runs.

  uv run python reviews/scripts/p6_wp6c_headline.py --duration 30 --seeds 1 2 3

For each seed: runner.decompose on F_liu20_nominal (base) and on F_liu20_nominal + lock/dual_band_slow +
lock/points/dual_band_slow_nominal (dbs). Extra diagnostics (seed 1 only, --diag):
  * dbs comp-1 run with pcf1.counting_noise off: how much of the dbs component 1 is counting noise;
  * the seed fringe X_s − π/2 against Θ_1 (is the seed an exact copy of the squeezed path?);
  * the base comp-1 series: fraction of its variance below f_c = 2 Hz (what the slow loop can remove).
Read-only: no preset is modified.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import math
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.runner.decompose import decompose_metric_ii  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

OV = ("lock/dual_band_slow", "lock/points/dual_band_slow_nominal")
COMP = ("locked_path_residual", "common_points", "laser_decorrelation_oob", "squeezed_only")


def spec(duration, seed, overlays=(), toggles=None, record=None):
    meas = None
    s = fs.load(fs.nominal("F_liu20"), run={"duration_s": duration, "seed": seed}, overlays=overlays, toggles=toggles)
    if record:
        meas = dict(s.measure or {})
        meas["record"] = {"slots": list(record), "decim": 100}
        s = dataclasses.replace(s, measure=meas)
    return s


def comps(s):
    d = decompose_metric_ii(s)[1]
    return {c: 1e3 * math.sqrt(d.variances_rad2[c]) for c in COMP} | {"total": 1e3 * math.sqrt(d.total_var_rad2),
                                                                  "ratio": d.linearity_ratio, "status": d.status}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=30.0)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--save", default=None)
    ap.add_argument("--spectrum-only", action="store_true", help="only diagnostic (3): the base comp-1 spectrum")
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    out = {}
    for sd in ([] if a.spectrum_only else a.seeds):
        for tag, ov in (("base", ()), ("dbs", OV)):
            t0 = time.time()
            r = comps(spec(a.duration, sd, ov))
            out[f"{tag}_seed{sd}"] = r
            print(f"seed {sd} {tag:4s}: " + "  ".join(f"{k[:6]} {v:7.2f}" for k, v in r.items()
                                                     if isinstance(v, float)) + f"  {r['status']}  ({time.time()-t0:.0f} s)",
                  flush=True)
    if a.diag:
        sd = a.seeds[0]
        groups = spec(1.0, sd, OV).measure["decompose"]["groups"]
        only1 = {t: (c == "locked_path_residual") for c, ts in groups.items() for t in ts}
        # (1) dbs comp-1 with and without counting noise
        for cn in (True, False):
            tg = {**only1, "pcf1.counting_noise": cn}
            r = run(spec(a.duration, sd, OV, toggles=tg))
            print(f"dbs comp-1 only, counting_noise={cn}: std {1e3*float(r.metrics.metric_ii[1].stats.std):.2f} mrad",
                  flush=True)
        # (2) seed fringe vs Theta_1, all on
        rec = ["metric.mii.port1", "beam.dseed1.phi", "beam.dseed0.phi", "beam.sqz1.phi", "beam.sqz0.phi"]
        r = run(spec(min(a.duration, 2.0), sd, OV, record=rec))
        ser = r.series
        names = [k for k in ser if k.startswith("metric.") and "port1" in k and not k.endswith(".acc")]
        xs = np.asarray(ser["beam.dseed1.phi"]) - np.asarray(ser["beam.dseed0.phi"])
        th = np.asarray(ser["beam.sqz1.phi"]) - np.asarray(ser["beam.sqz0.phi"])
        print(f"max |(φ_dseed1 − φ_dseed0) − (θ_sqz1 − θ_sqz0)| = {np.max(np.abs(xs - th)):.3e} rad  (metric slots {names})")
    if a.diag or a.spectrum_only:
        sd = a.seeds[0]
        groups = spec(1.0, sd, OV).measure["decompose"]["groups"]
        only1 = {t: (c == "locked_path_residual") for c, ts in groups.items() for t in ts if not t.startswith("pcf1.")}
        # (3) base comp-1 spectrum: variance fraction below f_c
        r = run(spec(a.duration, sd, (), toggles=only1, record=["metric.mii.port1"]))
        x = np.asarray(r.series["metric.mii.port1"], dtype=float)
        dt = 100 / 1.0e6
        x = x[int(1.0 / dt):]                     # drop the first second (acquisition)
        x = x - x.mean()
        X = np.fft.rfft(x)
        f = np.fft.rfftfreq(x.size, dt)
        p = np.abs(X) ** 2
        for fc in (0.5, 2.0, 10.0):
            print(f"base comp-1: fraction of variance below {fc:4.1f} Hz = {p[(f > 0) & (f < fc)].sum()/p[f > 0].sum():.3f}")
        print(f"base comp-1 series std (after 1 s) {1e3*x.std():.2f} mrad")
    if a.save:
        Path(a.save).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
