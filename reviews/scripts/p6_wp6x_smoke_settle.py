"""Phase 5 review of WP6-X (D94 結果): the bench/refA_fiber smoke table (reviews/scripts/p6_bench_smoke.py) with the
lock-acquisition transient removed.

At L = 100 m (fiber_matched) the EOM of lock 2 saturates (|v.dac2| = 0.5 V) during the first ≈ 10 ms for some seeds,
until the 100 Hz EOM → stretcher offload has taken the initial fibre excursion; metric (i) is then elevated in that
window, and p6_bench_smoke.py averages it over the whole 50 ms run. This script reruns the smoke points and prints the
time-loop std of metric (i) over the whole run and over t ≥ T_SETTLE, the first-window max |v.dac2|, and σ_total with the
time-loop part replaced by the settled one (σ² − std²_whole + std²_settled; the ledger part unchanged).
Usage: uv run python reviews/scripts/p6_wp6x_smoke_settle.py [overlay ...]
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np

from cranberry import config
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
DUR, T_SETTLE, SEEDS = 0.05, 0.015, (1, 2, 3)
overlays = sys.argv[1:]


def spec_at(L, route, seed):
    base = config.load({"extends": str(ROOT / "presets/bench/refA_fiber_bench.yaml"),
                        "run": {"duration_s": DUR, "seed": seed}, "measure": {"esa": None}}, overlays=list(overlays))
    return expand(base, {"offline.params.L_sq_fiber_m": [L], "offline.params.lo_route": [route],
                         "offline.params.eta_fiber_out": [0.98]})[0].spec


print(f"refA_fiber_bench {' '.join(overlays)}: std of metric (i) [mrad], whole run / t ≥ {T_SETTLE * 1e3:.0f} ms")
print(f"{'L':>5s} {'route':>14s} {'std whole':>16s} {'std settled':>16s} {'σ_total':>16s} {'σ_total settled':>16s}"
      f" {'max|v.dac2| t<T':>16s}")
for route in ("fiber_matched", "free_space"):
    for L in (1.0, 10.0, 100.0):
        rows = []
        for seed in SEEDS:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = run(spec_at(L, route, seed))
            x = np.asarray(res.series["metric.mi.port0"])
            t = np.asarray(res.t_s)
            nl = res.metrics.noise_levels["i.port0"]
            s_all, s_set = float(x.std()), float(x[t >= T_SETTLE].std())
            sig = nl["sigma_total_rad"]
            sig_set = math.sqrt(max(sig**2 - nl["std_loop_rad"] ** 2 + s_set**2, 0.0))
            dac = float(np.abs(np.asarray(res.series["v.dac2"])[t < T_SETTLE]).max())
            rows.append((s_all, s_set, sig, sig_set, dac))
        a = np.array(rows)
        m, s = a.mean(axis=0), a.std(axis=0, ddof=1)
        print(f"{L:5.0f} {route:>14s} " + " ".join(f"{1e3 * m[j]:7.2f} ± {1e3 * s[j]:5.2f}" for j in range(4))
              + f" {m[4]:9.3f} (max {a[:, 4].max():.3f})", flush=True)
