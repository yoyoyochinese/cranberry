"""Phase 6 fix of reviews/phase6_wp6c_p.md (C-1..C-3, C-9): F_liu20 component table of the two-layer lock at the
physically sound nominal point and its calibration-slot options.

  uv run python reviews/scripts/p6_wp6c_fix_headline.py --duration 30 --seed 1 --option frames --save out.json

Options (all after lock/dual_band_slow):
  frames  lock/points/dual_band_slow_nominal: seed after the squeezer, pump on in the calibration slot (squeezed vacuum
          b' counted), calibration frames (20 % of the frames, N_b gates each)
  guard   ... + lock/points/dbs_cal_guard (n_guard 1: one seed per frame in a guard slot, L_k × 1.25)
  ideal   lock/points/dual_band_slow_ideal (the pre-review D92-4 point: ideal_copy, ideal_free, b' = 0)
  base    F_liu20_nominal without the slow layer
Prints the four components (mrad), the total, the decomposition's linearity ratio and status; a status other than ok
(linearity outside the 5 % tolerance) marks the decomposition invalid (review C-9).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.runner.decompose import decompose_metric_ii  # noqa: E402

DBS = "lock/dual_band_slow"
OPTIONS = {"frames": (DBS, "lock/points/dual_band_slow_nominal"),
           "guard": (DBS, "lock/points/dual_band_slow_nominal", "lock/points/dbs_cal_guard"),
           "ideal": (DBS, "lock/points/dual_band_slow_ideal"),
           "base": ()}
COMP = ("locked_path_residual", "common_points", "laser_decorrelation_oob", "squeezed_only")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--option", choices=sorted(OPTIONS), default="frames")
    ap.add_argument("--save", default=None)
    ap.add_argument("--extra", nargs="*", default=[], help="overlay files applied last (e.g. a point of another variant)")
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    t0 = time.time()
    s = fs.load(fs.nominal("F_liu20"), run={"duration_s": a.duration, "seed": a.seed},
                overlays=(*OPTIONS[a.option], *a.extra))
    d = decompose_metric_ii(s)[1]
    v = d.variances_rad2
    out = {c: 1e3 * math.sqrt(v[c]) for c in COMP} | {"total": 1e3 * math.sqrt(d.total_var_rad2),
                                                      "ratio": d.linearity_ratio, "status": d.status,
                                                      "L_k_m": float(s.offline["params"]["L_k_m"])}
    print(f"seed {a.seed} {a.option:6s}: " + "  ".join(f"{c[:8]} {out[c]:7.2f}" for c in COMP)
          + f"  total {out['total']:7.2f} mrad  linearity {out['ratio']:.3f}  {out['status']}"
          + ("" if out["status"] == "ok" else "  (INVALID decomposition)") + f"  L_k {out['L_k_m']:.0f} m"
          + f"  ({time.time() - t0:.0f} s)", flush=True)
    if a.save:
        Path(a.save).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
