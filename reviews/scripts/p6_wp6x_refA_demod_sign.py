"""Phase 5 review of WP6-X (D94-2): sign of the lock-2 demodulation-phase absorption 2π f_sh Δτ in bench/refA_fiber.

At L = 100 m (the smoke table) 2π f_sh Δτ ≈ 3.06 rad ≈ π, so the absorption with the WRONG sign (error 2 × 3.06 rad ≈ 2π)
would lock just as well: the smoke table does not test the sign. At L = 50 m the probe term is 1.52 rad: the coded sign
gives full slope, no absorption leaves cos(1.52) ≈ 0.05 of it, the opposite sign cos(3.03) ≈ −0.99 (wrong loop sign).
Runs bench/refA_fiber_bench (fiber_matched, 20 ms, seed 1) with pr_bhd_static_extra_rad = +X (coded), 0, −X and prints
metric (i) and the lock-2 error.
Usage: uv run python reviews/scripts/p6_wp6x_refA_demod_sign.py [L_m]
"""
import math
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import yaml

from cranberry import config
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
L = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
X = 2 * math.pi * 1.0e6 * (1.468 * L - 1.0003) / 299792458.0
TMP = Path(tempfile.mkdtemp())

print(f"L = {L} m: 2π f_sh Δτ = {X:.3f} rad")
for label, extra in (("coded (+X)", X), ("none (0)", 0.0), ("flipped (−X)", -X)):
    pt = TMP / f"pt_{len(list(TMP.iterdir()))}.yaml"
    pt.write_text(yaml.safe_dump({"offline": {"params": {"pr_bhd_static_extra_rad": float(extra)}}}))
    base = config.load({"extends": str(ROOT / "presets/bench/refA_fiber_bench.yaml"),
                        "run": {"duration_s": 0.02, "seed": 1}, "measure": {"esa": None}}, overlays=[str(pt)])
    spec = expand(base, {"offline.params.L_sq_fiber_m": [L], "offline.params.lo_route": ["fiber_matched"],
                         "offline.params.eta_fiber_out": [0.98]})[0].spec
    assert spec.offline["params"]["pr_bhd_static_extra_rad"] == extra
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(spec)
    nl = res.metrics.noise_levels["i.port0"]
    half = len(res.t_s) // 2
    errs = {k: np.asarray(v)[half:] for k, v in res.series.items() if k.startswith("lock.") and k.endswith(".err")}
    e2 = {k: (float(np.mean(v)), float(np.std(v))) for k, v in errs.items()}
    print(f"  {label:>13s}: valid {res.valid}, events {len(res.events)}, mean Δ {1e3 * nl['mean_loop_rad']:+9.2f} mrad, "
          f"std Δ {1e3 * nl['std_loop_rad']:8.2f} mrad, lock errors (mean, std) [rad] "
          + ", ".join(f"{k}: {m:+.3f}, {s:.3f}" for k, (m, s) in e2.items()))
