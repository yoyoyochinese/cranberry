"""Phase 6 axis 1 (b): refC_fiber_bench squeezing / anti-squeezing at the ESA sidebands (measure.esa hooks of refC,
0.6 s zero-span records, so full-length runs; seed 1 only — the phase-noise seed spread is in sweep_C).

Usage: uv run python reviews/scripts/p6_axis1_refC_esa.py
Writes docs/figures/phase6/axis1/b_refC_esa.json and prints a table.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

LENGTHS = [1.0, 10.0, 100.0]
CTRL = {"fpga": None, "moku": "bench/moku_delta_refC", "dsc1": "bench/dsc1_refC"}


def scalars(m):
    out = {}
    for k, v in (m or {}).items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                if isinstance(vv, (int, float)) and not isinstance(vv, bool):
                    out[f"{k}.{kk}"] = vv
        elif isinstance(v, (int, float)):
            out[k] = v
    return out


if __name__ == "__main__":
    jobs = [{"label": f"{c}|{Lm:g}", "measurements": True,
             "point": {**L.c_point(Lm, ctrl=ctrl, duration_s=0.625), "measure_esa_off": False}}
            for c, ctrl in CTRL.items() for Lm in LENGTHS]
    rows = L.pool(jobs, "refC_esa")
    out = {}
    for r in rows:
        esa = (r.get("measurements") or {}).get("esa", {})
        out[r["label"]] = {"error": r.get("error"), "esa": {k: scalars(v) if isinstance(v, dict) else v
                                                          for k, v in esa.items() if not k.startswith("spectrum")},
                           "std_rad": r.get("i.port0.std_rad"), "eta": r.get("i.port0.eta")}
        print(r["label"], r.get("error", ""), {k: {kk: (round(vv, 3) if isinstance(vv, float) else vv)
                                                  for kk, vv in v.items() if "db" in kk.lower()}
                                              for k, v in out[r["label"]]["esa"].items() if isinstance(v, dict)})
    (L.FIG / "b_refC_esa.json").write_text(json.dumps(out, indent=1, default=float))
    _ = np
