"""Phase 6 WP6-X smoke result of bench configuration (a) (D85-12): presets/bench/refA_fiber_bench.yaml (refA +
bench/refA_fiber) at L = 1, 10, 100 m, both LO routes, 3 seeds each, 50 ms runs (ESA hooks off: they need 0.15 s).
Prints metric (i) of port 0 (refA is port-less: metric (ii) is not applicable, G13) — mean and std of Δ in the time loop,
total σ with the ledger (σ_total), the squeezing with the phase jitter V′ [dB], η at the BHD — and the wall time per
simulated second after a JIT warm-up. Optional overlays as arguments (e.g. bench/dsc1_refA).
Usage: uv run python reviews/scripts/p6_bench_smoke.py [overlay ...]
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from cranberry import config
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
DUR = 0.05
SEEDS = (1, 2, 3)


def spec_at(L, route, seed, overlays):
    base = config.load({"extends": str(ROOT / "presets/bench/refA_fiber_bench.yaml"),
                        "run": {"duration_s": DUR, "seed": seed}, "measure": {"esa": None}}, overlays=list(overlays))
    return expand(base, {"offline.params.L_sq_fiber_m": [L], "offline.params.lo_route": [route],
                         "offline.params.eta_fiber_out": [0.98]})[0].spec


def go(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = time.perf_counter()
        res = run(spec)
        return res, time.perf_counter() - t


overlays = sys.argv[1:]
go(spec_at(1.0, "fiber_matched", 0, overlays))  # JIT warm-up
print(f"refA_fiber_bench {' '.join(overlays)}  duration {DUR} s, f_sim {spec_at(1.0, 'fiber_matched', 0, overlays).run.f_sim_hz:.3g} Hz")
print(f"{'L [m]':>6s} {'LO route':>14s} {'mean Δ [mrad]':>16s} {'std Δ [mrad]':>14s} {'σ_total [mrad]':>16s} {'V′ [dB]':>14s}"
      f" {'η_bhd':>7s} {'valid':>6s} {'wall/s_sim':>11s}")
for route in ("fiber_matched", "free_space"):
    for L in (1.0, 10.0, 100.0):
        rows, walls, valid = [], [], True
        for seed in SEEDS:
            res, wall = go(spec_at(L, route, seed, overlays))
            nl = res.metrics.noise_levels["i.port0"]
            rows.append((nl["mean_loop_rad"], nl["std_loop_rad"], nl["sigma_total_rad"], nl["v_prime_db"], nl["eta"]))
            walls.append(wall / DUR)
            valid &= bool(res.valid) and not res.events
        a = np.array(rows)
        m, s = a.mean(axis=0), a.std(axis=0, ddof=1)
        print(f"{L:6.0f} {route:>14s} {1e3 * m[0]:+8.2f} ± {1e3 * s[0]:5.2f} {1e3 * m[1]:7.2f} ± {1e3 * s[1]:4.2f}"
              f" {1e3 * m[2]:8.2f} ± {1e3 * s[2]:4.2f} {m[3]:7.3f} ± {s[3]:4.3f} {m[4]:7.4f} {str(valid):>6s}"
              f" {np.mean(walls):9.1f} s")
