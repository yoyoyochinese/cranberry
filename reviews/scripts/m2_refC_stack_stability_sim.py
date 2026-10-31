"""Time-domain confirmation of m2_refC_stack_stability.py: one grid point predicted unstable and its stable
neighbour (f_12 = 300 Hz), refC_nominal otherwise, all stochastic sources on, 0.1 s.
Run: uv run python reviews/scripts/m2_refC_stack_stability_sim.py
"""
import importlib.util
import math
from pathlib import Path

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refC_support", ROOT / "tests" / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)
M = math.pi * 0.1 / 4.0

for lab, p in (("predicted unstable: f_u1 30 kHz, f_12 3 kHz, f_res 10 kHz, Q 30",
                {"f_u1_hz": 3e4, "f_12_hz": 3e3, "mirror_f_res_hz": 1e4, "mirror_q": 30.0}),
               ("predicted unstable: f_u1 9.65 kHz, f_12 3 kHz, f_res 10 kHz, Q 30",
                {"f_u1_hz": 9655.0, "f_12_hz": 3e3, "mirror_f_res_hz": 1e4, "mirror_q": 30.0}),
               ("stable neighbour:   f_u1 9.65 kHz, f_12 300 Hz, f_res 10 kHz, Q 30",
                {"f_u1_hz": 9655.0, "f_12_hz": 300.0, "mirror_f_res_hz": 1e4, "mirror_q": 30.0})):
    res = run(rs.load(params=p, run={"duration_s": 0.1}))
    nl = res.metrics.noise_levels["i.port0"]
    tot = rs.at(nl, "sigma_total_rad", 3e7)
    print(f"{lab}: total {tot*1e3:8.1f} mrad, loop std {rs.at(nl, 'std_loop_rad', 3e7)*1e3:8.1f} mrad, "
          f"supervisor losses {int(res.final_state['lock.sup.n_loss'])}")
