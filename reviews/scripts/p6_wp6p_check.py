"""Phase 5 review of WP6-P (D95): heterodyne programming, sign, step_response wrapping, metric (ii) with switching
transients. Read-only (overlay files are written to a temp dir).

  uv run python reviews/scripts/p6_wp6p_check.py            (≈ 5–10 min)

A  heterodyne RF program on F_wg6 with a table spanning ±1.5 rad: lock-coordinate and Θ bias at nominal K, K/2
   (visibility 0.5) and K/2 (power × 0.25); program sign −1 (preset) and +1.
B  step_response's Θ wrapping: synthetic segments (exact first-order step, no noise) with quadrature steps |Δ| < π/2
   and ≥ π/2; how often the shipped random tables (seed 0) produce |Δ(p_1 − p_0)| ≥ π/2.
C  metric (ii) std of a programmed run (π/4 walk, dwell 1 ms) against the same run without program: the set point
   follows the commanded target instantly, so every switching transient is counted in metric (ii).
"""

from __future__ import annotations

import importlib.util
import math
import tempfile
import warnings
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.blocks.ctrl.phase_program import random_table  # noqa: E402
from cranberry.blocks.meas.step_response import segment_response, _wrap  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

HET = ("probe/detect_heterodyne", "program/p6_steps")
QUIET = {"laser.fm_white": False, "laser.fm_flicker": False}
TMP = Path(tempfile.mkdtemp())


def overlay(name, table, dwell, sign=-1.0, liu=None, step_rad=None, source="table"):
    prog = {"sign": sign}
    if source == "table":
        prog |= {"source": "table", "times_s": [5e-3 + i * dwell for i in range(len(table))], "table_rad": table}
    else:
        prog |= {"step_rad": step_rad}
    pipe = [{"name": "squeezer", "seq": [{"name": "prog", "params": prog}]}]
    for k in (0, 1) if liu else ():
        pipe.append({"name": f"liu{k}", "seq": [{"name": f"liu{k}", "params": dict(liu)}]})
    p = TMP / f"{name}.yaml"
    p.write_text(yaml.safe_dump({"offline": {"params": {"prog_dwell_s": dwell}}, "pipeline": pipe}))
    return str(p)


def go(ovl, duration, toggles=QUIET, decim=10):
    spec = fs.load(fs.nominal("F_wg6"), run={"duration_s": duration}, toggles=toggles, overlays=ovl,
                   measure={"record": {"slots": ["lock.L0.err"], "decim": decim, "antialias": "none"}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(fs.first_point(spec))


def sec_a():
    table = [[1.5, -0.4], [-1.5, 0.7], [0.0, -1.5], [0.8, 1.5]]
    dwell = 4e-3
    print("A  heterodyne RF program, table ±1.5 rad (bias = mean over the last half of each 4 ms segment)")
    for label, liu, sign in (("nominal, sign −1", None, -1.0), ("visibility 0.5", {"visibility": 0.5}, -1.0),
                             ("power ×0.25", {"p_pr_avg_w": 2.5e-7}, -1.0), ("nominal, sign +1", None, 1.0)):
        r = go((*HET, overlay(f"a_{label[:5]}_{sign}", table, dwell, sign, liu)), 5e-3 + len(table) * dwell)
        sr = r.measurements["step_response"]["p6"]
        print(f"   {label:18s}: max|lock bias| {np.nanmax(np.abs(sr['bias_rad'])):.4f} rad, max|Θ bias| "
              f"{np.nanmax(np.abs(sr['theta_bias_rad'])):.4f} rad, slips {np.nansum(sr['slips']):.0f}, "
              f"Θ bias {np.round(sr['theta_bias_rad'][:, 0], 3)}; programmed Δ(p1−p0) "
              f"{np.round(np.diff(np.r_[[0.0], np.array(table)[:, 1] - np.array(table)[:, 0]]), 3)}")


def sec_b():
    print("B  step_response Θ wrapping (period π) on exact first-order steps (τ = 0.16 ms, tol 10 mrad)")
    t = np.arange(0, 3e-3, 1e-5)
    for d in (0.4, 1.2, 1.6, 2.4, math.pi - 0.05):
        g_prev, g = 0.0, d
        y = g_prev + d * (1 - np.exp(-t / 1.6e-4))            # the physical quadrature angle follows the shortest way?
        yq = g + _wrap(y - g, math.pi)                         # what step_response.table_of feeds segment_response
        r = segment_response(t, yq, 0.0, g_prev, g, 0.01)
        print(f"   Δ = {d:5.2f}: rise {1e3*r['rise_s']:.3f} ms (true {1e3*1.6e-4*math.log(9):.3f}), overshoot "
              f"{r['overshoot']:.3f}, settle {1e3*r['settle_s']:.3f} ms")
    for st in (None, math.pi / 8, math.pi / 4, math.pi / 2):
        tb = random_table(30, 2, st, 0)
        dth = np.diff(np.r_[[0.0], tb[:, 1] - tb[:, 0]])
        print(f"   random table step {st}: switches with |Δ(p1−p0)| ≥ π/2: {int(np.sum(np.abs(dth) >= math.pi/2 - 1e-12))}/30"
              f", max {np.max(np.abs(dth)):.2f}")


def sec_c():
    print("C  metric (ii) std with and without the program (F_wg6 + heterodyne, π/4 walk, dwell 1 ms, 30 switches)")
    dur = 5e-3 + 30 * 1e-3
    on = go((*HET, overlay("c_on", None, 1e-3, step_rad=math.pi / 4, source="random")), dur, toggles={})
    off_ovl = TMP / "c_off.yaml"
    off_ovl.write_text(yaml.safe_dump({"pipeline": [{"name": "squeezer", "seq": [{"name": "prog", "params": {
        "toggles": {"enabled": False, "quantization": True}}}]}], "offline": {"params": {"prog_dwell_s": 1e-3}}}))
    off = go((*HET, overlay("c_off0", None, 1e-3, step_rad=math.pi / 4, source="random"), str(off_ovl)), dur, toggles={})
    s_on, s_off = (float(x.metrics.metric_ii[1].stats.std) for x in (on, off))
    sr = on.measurements["step_response"]["p6"]
    print(f"   metric (ii) std: program on {1e3*s_on:.1f} mrad, off {1e3*s_off:.1f} mrad;"
          f" settle (port 0) median {1e3*np.nanmedian(sr['settle_s'][:, 0]):.2f} ms")


if __name__ == "__main__":
    sec_b()
    sec_a()
    sec_c()
