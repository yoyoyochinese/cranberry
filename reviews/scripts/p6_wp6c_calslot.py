"""Phase 5 review of WP6-C (D92-4): the cost of the calibration slot. The nominal point counts one port-k calibration
gate per frame (dbs_f_gate_hz = 1/t_frame, cal_fraction 1) but F_liu20 has n_guard = 0: every slot of the frame is a
data slot, so the seed slot "next to the block" needs an extra slot per block (n_guard → 1, demux/guard_slot), which
lengthens L_k by (N_b + 1)/N_b and with it Δτ_k (component 3).

  uv run python reviews/scripts/p6_wp6c_calslot.py --duration 1

Component 3 does not depend on the run length (1 s = 30 s within 0.1 %, p6_wp6c_headline.py), so 1 s suffices for it.
Read-only.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.runner.decompose import decompose_metric_ii  # noqa: E402

DBS = ("lock/dual_band_slow", "lock/points/dual_band_slow_nominal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    for label, ov in (("dbs, n_guard 0 (as D92-4)", DBS), ("dbs + demux/guard_slot", ("demux/guard_slot", *DBS))):
        s = fs.load(fs.nominal("F_liu20"), run={"duration_s": a.duration, "seed": a.seed}, overlays=ov)
        d = decompose_metric_ii(s)[1]
        p = dict(s.offline["params"])
        v = d.variances_rad2
        print(f"{label:28s}: L_k {p['L_k_m']:.0f} m, t_frame {1e6*p['t_frame_s']:.2f} µs; comps " +
              "  ".join(f"{1e3*math.sqrt(v[c]):6.2f}" for c in ("locked_path_residual", "common_points",
                                                             "laser_decorrelation_oob", "squeezed_only")) +
              f"  total {1e3*math.sqrt(d.total_var_rad2):6.2f} mrad  {d.status}", flush=True)


if __name__ == "__main__":
    main()
