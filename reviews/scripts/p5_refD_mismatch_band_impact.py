"""Phase 5 (Ref D prediction follow-up): effect of the re-classification on the frozen prediction's §4 band
("LHS, analytically stable points, quantisation on"), from the frozen per-run scalars only (no re-run).

The stable set with the ``dsc1`` check (p5_refD_mismatch_reclassify.py): lhs_002, lhs_042, lhs_118 leave it,
lhs_024, lhs_034, lhs_044 join it (76 points either way). Prints old / new median [5 %, 95 %].

Run: uv run python reviews/scripts/p5_refD_mismatch_band_impact.py
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = {"lhs_002", "lhs_042", "lhs_118"}
IN = {"lhs_024", "lhs_034", "lhs_044"}
KEYS = ["err.std", "lock.lockD.err.rms_1_200000", "lock.lockD.err.rms_10000_100000", "err.offset", "v.pdb.std",
        "v.pdb.ptp_10ms_median", "v.dsc1.out.mean", "v.dsc1.out.ptp_10ms_median", "v.dsc1.out.rms_10000_100000",
        "t_loss_extrap_s"]


def q(x):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    return np.median(x), np.percentile(x, 5), np.percentile(x, 95)


def main():
    runs = json.loads((ROOT / "predictions/refD/results/runs.json").read_text())["runs"]
    for quant in (1, 0):
        lhs = {r["id"][:7]: r for r in runs if r["group"] == "lhs" and int(r["quant"]) == quant}
        old = [r for r in lhs.values() if r["stability"]["stable"]]
        new = [r for k, r in lhs.items() if (r["stability"]["stable"] and k not in OUT) or k in IN]
        print(f"quantisation {'on' if quant else 'off'}: n old {len(old)}, new {len(new)}")
        for k in KEYS:
            a, b = q([r["scalars"].get(k) for r in old]), q([r["scalars"].get(k) for r in new])
            print(f"  {k:34s} old {a[0]:10.4g} [{a[1]:10.4g}, {a[2]:10.4g}]   new {b[0]:10.4g} [{b[1]:10.4g}, {b[2]:10.4g}]")
        n_old = sum(r["locked"] for r in old)
        n_new = sum(r["locked"] for r in new)
        print(f"  locked in time loop: old {n_old}/{len(old)}, new {n_new}/{len(new)}")


if __name__ == "__main__":
    main()
