"""Phase 5 (Ref D prediction follow-up): time-loop reproduction of the six stability mismatches.

Each point of predictions/refD/results/runs.json is rebuilt with the prediction's own ``build_spec`` (the
exact point values, preset expressions re-evaluated) and run for a short time with
  * ``full``:   as in the prediction (noise and quantisation on);
  * ``lin``:    every noise source, drift and both quantisers off (refD_support.NOISE_OFF): the loop is
                deterministic, so after pull-in the residual starts at the 1e-12 rad level and any
                oscillation that appears must grow out of the linear dynamics.
The 1 ms-window rms of lock.lockD.err (about the lock point) is printed against time with the dominant
frequency of the last 10 ms. A linear instability shows as exponential growth (a straight line in log rms)
from the 1e-12 level until the fringe nonlinearity bounds it near 1.6 rad.

Run: uv run python reviews/scripts/p5_refD_mismatch_timeloop.py lhs_002 lhs_042 lhs_024 lhs_034 lhs_044 lhs_057
"""

import dataclasses
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "predictions/refD"))
import run_prediction as rp  # noqa: E402  (read-only use of the frozen prediction's spec builder)

_sp = importlib.util.spec_from_file_location("rs", ROOT / "tests/integration/refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

from cranberry.runner.run import run  # noqa: E402

DUR = 0.06


def envelope(rid, mode):
    runs = {r["id"]: r for r in json.loads((ROOT / "predictions/refD/results/runs.json").read_text())["runs"]}
    spec = rp.build_spec(dict(runs[rid + "_q1"]), duration_s=DUR)
    if mode == "lin":
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **rs.NOISE_OFF})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(spec)
    e = np.asarray(res.series["lock.lockD.err"], float)
    t = np.asarray(res.t_s, float)
    fs = 1.0 / (t[1] - t[0])
    tail = e[t > DUR - 0.01]
    branch = round(float(np.mean(tail)) / np.pi) * np.pi  # nearest stable point (0 or ±π about +π/2)
    w = int(round(1e-3 * fs))
    n = len(e) // w
    blocks = (e[: n * w] - branch).reshape(n, w)
    rms = np.sqrt(np.mean((blocks - blocks.mean(axis=1, keepdims=True)) ** 2, axis=1))
    x = tail - tail.mean()
    sp = np.abs(np.fft.rfft(x * np.hanning(len(x))))
    f_dom = np.fft.rfftfreq(len(x), 1 / fs)[np.argmax(sp[1:]) + 1]
    return rms, f_dom, float(np.std(tail)), branch


def main(ids):
    for rid in ids:
        for mode in ("full", "lin"):
            rms, f_dom, std_tail, branch = envelope(rid, mode)
            pick = [1, 2, 3, 5, 8, 12, 20, 30, 45, 59]
            print(f"{rid} {mode:4s} branch {branch:+.2f} rad, std(last 10 ms) {std_tail:.3g} rad, "
                  f"dominant {f_dom / 1e3:.2f} kHz")
            print("    1-ms rms at t[ms] " + "  ".join(f"{k}:{rms[k]:.1e}" for k in pick if k < len(rms)))


if __name__ == "__main__":
    main(sys.argv[1:])
