"""Phase 6 axis 2: artefact checks of the headline points (runner.convergence): step size (f_sim × 2, same root seed,
rate-tied expressions re-evaluated), chunk size (bit-identical), and the seed spread (taken from the 10-seed campaign
jobs, reported by p6_axis2_report). Also the run-length check 3 s vs 30 s (campaign S_/L_ jobs).

  uv run python reviews/scripts/p6_axis2_convergence.py [--labels base P3_wdm_chop6k ...] [--duration 3]

Quantities and tolerances: ii.port1.std_loop_rad and i.port1.sigma_total_rad, 5 mrad (pass: |Δ| ≤ tol/3 = 1.7 mrad,
PLAN §6.4), ii.port1.hf_variance_rad2 (ledger) 2.5e-5 rad². Writes runs/p6_axis2/convergence.json.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import warnings
from pathlib import Path

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

from cranberry.runner.convergence import convergence  # noqa: E402

Q = {"ii.port1.std_loop_rad": 0.005, "i.port1.sigma_total_rad": 0.005, "ii.port1.hf_variance_rad2": 2.5e-5}


def _plain(o):
    if isinstance(o, dict):
        return {str(k): _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    try:
        return float(o)
    except (TypeError, ValueError):
        return o if isinstance(o, (str, bool)) or o is None else str(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="F_liu20")
    ap.add_argument("--labels", nargs="*", default=["base", "P3_wdm_chop6k", "P2_borealis_si", "P5_dbs"])
    ap.add_argument("--duration", type=float, default=3.0)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    path = lib.OUT / "convergence.json"
    out = __import__("json").loads(path.read_text()) if path.exists() else {}
    for label in a.labels:
        spec = lib.build(a.base, label, duration_s=a.duration, seed=1)
        rep = convergence(spec, ("chunk", "f_sim"), quantities=Q, f_sim_factors=(2,))
        out[f"{a.base}/{label}"] = {"passed": rep.passed, "checks": _plain(rep.checks)}
        print(label, rep.passed, {k: v["passed"] for k, v in rep.checks.items()}, flush=True)
        for k, v in rep.checks["f_sim"]["values"]["x2"].items():
            print("   ", k, f"ref {v['reference']:.5g} x2 {v['value']:.5g} Δ {v['change']:.3g} (tol/3 {v['tolerance'] / 3:.3g})")
        path.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
