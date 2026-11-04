"""Quick text view of runs/p6_axis2 jobs (axis 2): uv run python reviews/scripts/p6_axis2_peek.py [substring ...]"""

import importlib.util
import sys
import warnings
from pathlib import Path

_sp = importlib.util.spec_from_file_location("r", Path(__file__).with_name("p6_axis2_report.py"))
r = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(r)
warnings.simplefilter("ignore")
COLS = ("ii_total", *r.lib.COMP, "mean_ii", "sigma_i", "eta_gbs", "stray_mode", "win_mode_mean", "fref_mode", "slips",
        "range_pzt1", "range_tec1", "valid")
for row in sorted(r.table(r.lib.load_results()), key=lambda x: x["tag"]):
    if sys.argv[1:] and not any(s in row["tag"] for s in sys.argv[1:]):
        continue
    parts = []
    for c in COLS:
        if c in row:
            m, s = row[c]["mean"], row[c]["std"]
            parts.append(f"{c[:6]}={m:.3g}±{s:.2g}")
    print(f"{row['tag']:48s} L={row['L_k_m']['mean']:.0f} " + " ".join(parts) + (f" [{row['events']}]" if row["events"] else ""))
