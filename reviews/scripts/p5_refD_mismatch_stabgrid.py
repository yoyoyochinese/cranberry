"""Phase 5 (Ref D prediction follow-up): how the D75 stable region (PREDICTION.md §0 / §6.2, stabgrid) moves when
``loop_stability_check`` uses the DSC1's discrete PID law instead of K_p + K_i/s + K_d s. Same 41 × 41
(f_res, Q) grid, nominal plant gain, the two T15 DSC1 delay models; prints the largest stable Q per f_res
(old continuous / new dsc1) at the columns of PREDICTION.md §6.2. Uses run_prediction.stability_grid's
recipe (chk params of a resolved refD point) read-only; nothing under predictions/ is written.

Run: uv run python reviews/scripts/p5_refD_mismatch_stabgrid.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "predictions/refD"))
import run_prediction as rp  # noqa: E402

from cranberry import config  # noqa: E402
from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop  # noqa: E402
from cranberry.runner.run import offline_block_specs  # noqa: E402
from cranberry.runner.sweep import expand  # noqa: E402

COLS_KHZ = [5, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60]


def main():
    warnings.simplefilter("ignore")
    axes = rp.axes_of_preset()
    nom = {a["name"]: a["nominal"] for a in axes}
    f_res = np.geomspace(5.0e3, 6.0e4, 41)
    q_res = np.geomspace(3.0, 100.0, 41)
    base = config.load(rp.REFD)
    for vname, vv in {"5.7us_poles": {"latency_s": 5.7e-6, "in_out_poles": True},
                      "8.3us_nopoles": {"latency_s": 8.3e-6, "in_out_poles": False}}.items():
        (pt,) = expand(base, {f"offline.params.{k}": [v] for k, v in {**nom, **vv}.items()})
        (cs,) = [b for b in offline_block_specs(pt.spec) if b.kind == "loop_stability_check"]
        p0 = dict(cs.params)
        lo, hi, n = p0["f_grid_hz"]
        fgrid = np.geomspace(float(lo), float(hi), int(n))
        qmax = {}
        for disc in ("continuous", "dsc1"):
            p = {**p0, "discretisation": disc}
            st = np.zeros((41, 41), bool)
            for i, fr in enumerate(f_res):
                for j, qq in enumerate(q_res):
                    p["modes"] = [{"f_hz": float(fr), "q": float(qq)}]
                    st[i, j] = nyquist_stable(fgrid, open_loop(fgrid, p), n_integrators(p))[0]
            qmax[disc] = [q_res[st[i]].max() if st[i].any() else None for i in range(41)]
            all_q = [f_res[i] for i in range(41) if st[i].all()]
            print(f"{vname} {disc:10s} stable fraction {st.mean():.3f}; f_res stable at every Q: "
                  + ", ".join(f"{f_res[i] / 1e3:.1f}" for i in range(41) if st[i].all()) + " kHz"
                  + f" (min {min(all_q) / 1e3:.1f} kHz)" if all_q else "")
        print(f"{vname}: largest stable Q per f_res (continuous / dsc1)")
        for c in COLS_KHZ:
            i = int(np.argmin(abs(f_res - c * 1e3)))
            a, b = qmax["continuous"][i], qmax["dsc1"][i]
            print(f"   {f_res[i] / 1e3:6.1f} kHz: {'-' if a is None else f'{a:5.1f}'} / {'-' if b is None else f'{b:5.1f}'}")


if __name__ == "__main__":
    main()
