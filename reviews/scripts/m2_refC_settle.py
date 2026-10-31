"""Phase 5 M2 Ref C review: is the metric settle time t_settle_s = 0.02 s long enough over the f_12 axis?
Acquisition from the model's initial differential fibre phase: with the EOM clamped (±out_range_v), the
mirror stage slews at ki2·out_range·k_mirror = 2π f_12·(π/V_π)·out_range [rad/s] (independent of the
mirror range), i.e. 74 rad/s at f_12 = 30 Hz: an initial offset near π/2 needs ≈ 21 ms.
Records Δ for 0.1 s over several seeds and reports the last time |Δ| > 0.05 rad.
Run: uv run python reviews/scripts/m2_refC_settle.py
"""
import importlib.util
import math
from pathlib import Path

import numpy as np

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refC_support", ROOT / "tests" / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)
REC = {"record": {"slots": ["lock.lockC.err"], "decim": 4, "antialias": "none"}}
print("slew-limited acquisition time for |Δ0| = π/2:",
      {f12: round((math.pi / 2) / (2 * math.pi * f12 * math.pi / 4 * 0.5) * 1e3, 1) for f12 in (30, 300, 3000)}, "ms")
for f12 in (30.0, 300.0):
    out = []
    for seed in range(1, 9):
        spec = rs.load(params={"f_12_hz": f12}, run={"duration_s": 0.1, "seed": seed}, measure=REC)
        res = run(spec)
        e = np.asarray(res.series["lock.lockC.err"], dtype=float)
        d = np.pi / 2 - np.mod(np.pi / 2 - e, np.pi)
        t = np.arange(d.size) * 4 / spec.run.f_sim_hz
        big = np.nonzero(np.abs(d) > 0.05)[0]
        out.append((abs(d[0]), t[big[-1]] if big.size else 0.0))
    print(f"f_12 = {f12:6.1f} Hz: (|Δ(0)| rad, last |Δ| > 0.05 rad [ms]) per seed:",
          [(round(a, 2), round(b * 1e3, 1)) for a, b in out])
