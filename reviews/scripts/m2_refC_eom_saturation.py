"""Phase 5 M2 Ref C review: why does h_ref_conduit ×10 raise the non-dither residual 20× (4 → 82 mrad) instead
of √10×? Records the EOM drive (±out_range_v = ±0.39 rad) and the mirror drive, and repeats with a wider
EOM range to separate saturation from linear loop behaviour.

Run: uv run python reviews/scripts/m2_refC_eom_saturation.py
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
M = math.pi * 0.1 / 4.0
REC = {"record": {"slots": ["ctrl.pid1.u", "ctrl.pid2.u", "lock.lockC.err"], "decim": 1, "antialias": "none"}}

for lab, p in (("h_ref_conduit 4 (nominal)", {}), ("h_ref_conduit 40", {"h_ref_conduit_rad2_hz": 40.0}),
               ("h_ref_conduit 40, f_12 = 3 kHz", {"h_ref_conduit_rad2_hz": 40.0, "f_12_hz": 3e3}),
               ("f_12 = 30 Hz (axis bottom)", {"f_12_hz": 30.0})):
    spec = rs.load(params=p, run={"duration_s": 0.3}, measure=REC)
    res = run(spec)
    u1 = np.asarray(res.series["ctrl.pid1.u"])[int(0.02 * 4e6):]
    u2 = np.asarray(res.series["ctrl.pid2.u"])[int(0.02 * 4e6):]
    rng = rs.params_of(spec)["out_range_v"]
    nl = res.metrics.noise_levels["i.port0"]
    tot = rs.at(nl, "sigma_total_rad", 3e7)
    print(f"{lab:<32s} non-dither {math.sqrt(max(tot**2 - M**2/2, 0))*1e3:6.1f} mrad; EOM |u1| max {np.max(np.abs(u1)):.3f} V"
          f" (range {rng}), time at the clamp {np.mean(np.abs(u1) > 0.999 * rng)*100:.2f} %; mirror |u2| max {np.max(np.abs(u2)):.3f} V")
