"""Phase 5 (Ref D prediction follow-up): re-run ``loop_stability_check`` on every frozen prediction point with
the controller law the time loop runs (``discretisation: dsc1``, DSC1 p.19) and with the old continuous law
K_p + K_i/s + K_d s, and compare both verdicts with the frozen time-loop outcome (``locked`` in
predictions/refD/results/runs.json, quantisation-on runs; nothing is re-simulated).

Also: the closed-loop pole nearest the 915B resonance from the analytic model (root of 1 + L(s) by Newton,
s complex) for the points given on the command line, to compare with the growth rate / frequency measured
by p5_refD_mismatch_timeloop.py.

Run: uv run python reviews/scripts/p5_refD_mismatch_reclassify.py [lhs_002 lhs_042 ...]
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "predictions/refD"))
import run_prediction as rp  # noqa: E402  (read-only use of the frozen prediction's spec builder)

from cranberry.blocks.chk.loop_stability import _ctrl, _dynamics, n_integrators, nyquist_stable  # noqa: E402
from cranberry.core.registry import ensure_loaded, make_block  # noqa: E402
from cranberry.runner.run import offline_block_specs  # noqa: E402


def chk_params(job):
    spec = rp.build_spec(dict(job), duration_s=0.01)
    ensure_loaded()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (chk,) = [make_block(b) for b in offline_block_specs(spec) if b.kind == "loop_stability_check"]
    return dict(chk.params), chk.result


def verdict(p):
    lo, hi, n = p["f_grid_hz"]
    f = np.geomspace(float(lo), float(hi), int(n))
    s = 2j * np.pi * f
    lg = float(p["k0"]) * _ctrl(s, p) * np.exp(-s * float(p["delay_s"])) * _dynamics(s, p)
    return nyquist_stable(f, lg, n_integrators(p))


def pole_near_resonance(p):
    w = 2 * np.pi * float(p["modes"][0]["f_hz"])
    q = float(p["modes"][0]["q"])

    def g(s):  # (1 + L)·(s² + sω/Q + ω²): the open-loop resonance pole removed, analytic near it
        s = np.atleast_1d(s)
        lg = float(p["k0"]) * _ctrl(s, p) * np.exp(-s * float(p["delay_s"])) * _dynamics(s, p)
        return ((1.0 + lg) * (s * s + s * w / q + w * w))[0]
    s = complex(-w / (2 * q), w)
    for _ in range(100):
        d = 1e-7 * abs(s)
        ds = g(s) / ((g(s + d) - g(s - d)) / (2 * d))
        s -= ds
        if abs(ds) < 1e-9 * abs(s):
            break
    return s


def main(pole_ids):
    runs = json.loads((ROOT / "predictions/refD/results/runs.json").read_text())["runs"]
    if pole_ids and pole_ids[0] == "--poles-only":
        pole_ids, rows = pole_ids[1:], []
    else:
        rows = None
    if rows is None:
        rows = [r for r in runs if r["feedback"] and r["quant"] and r["group"] in ("lhs", "resgrid", "oat", "nominal")]
    seen = {}
    for r in rows:  # nominal: 8 seeds of one point
        key = json.dumps(r["values"], sort_keys=True)
        seen.setdefault(key, r)
    table = {}
    for r in seen.values():
        p, res = chk_params(r)
        new = verdict(p)
        old = verdict({**p, "discretisation": "continuous"})
        assert old[0] == r["stability"]["stable"], r["id"]  # the frozen verdict is the continuous one
        table[r["id"]] = (old, new, r["locked"], r["scalars"]["err.std"])
    for grp in ("lhs", "resgrid", "oat", "nominal"):
        ids = [k for k in table if k.startswith(grp)]
        agree_old = sum(table[k][0][0] == table[k][2] for k in ids)
        agree_new = sum(table[k][1][0] == table[k][2] for k in ids)
        n_st_old = sum(table[k][0][0] for k in ids)
        n_st_new = sum(table[k][1][0] for k in ids)
        print(f"{grp:8s} n={len(ids):3d}  stable old/new {n_st_old}/{n_st_new}  "
              f"agree with time loop old {agree_old}  new {agree_new}")
    print("points where old and new verdicts differ, or new disagrees with the time loop:")
    for k, (old, new, locked, std) in sorted(table.items()):
        if old[0] != new[0] or new[0] != locked:
            print(f"  {k:28s} old {'S' if old[0] else 'U'} (max|L|@-180 {old[1]:.3f} at {old[2] / 1e3:6.2f} kHz)  "
                  f"new {'S' if new[0] else 'U'} ({new[1]:.3f} at {new[2] / 1e3:6.2f} kHz)  "
                  f"time loop {'locked' if locked else 'not locked'} (err std {std:.3g} rad)")
    byid = {r["id"]: r for r in runs}
    for rid in pole_ids:
        p, _ = chk_params(byid[rid + "_q1"])
        for d in ("continuous", "dsc1"):
            s = pole_near_resonance({**p, "discretisation": d})
            print(f"{rid} {d:10s} pole near f_res: growth {s.real:+9.1f} /s ({np.exp(s.real * 1e-3):.3f} per ms), "
                  f"f {s.imag / 2 / np.pi / 1e3:.3f} kHz")


if __name__ == "__main__":
    main(sys.argv[1:])
