"""Phase 5 (Ref D prediction follow-up): open loop around the 915B resonance at prediction points, three ways:
  chk      the frozen check (continuous PID K_p + K_i/s + K_d s),
  z        the exact multirate z-domain loop of the time-loop blocks, −g·C·(B·A + B'·A') (refD_support.loop_parts,
           image term at f + f_sim/2 included), sign normalised to negative feedback at 10 Hz,
  chk_dsc1 the check with the DSC1 discrete PID law (DSC1 p.19): K_p + K_i h/2 + K_i/s + K_d (1 − e^{−sh})/h.

Run: uv run python reviews/scripts/p5_refD_mismatch_loops.py lhs_002 lhs_042 lhs_024 lhs_034 lhs_044 lhs_057
"""
import importlib.util, json, sys, warnings  # noqa: E401
import numpy as np
sys.path.insert(0, "predictions/refD")
import run_prediction as rp
_sp = importlib.util.spec_from_file_location("rs", "tests/integration/refD_support.py")
rs = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(rs)
from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.runner.run import offline_block_specs
from cranberry.core.registry import ensure_loaded, make_block
warnings.simplefilter("ignore")
runs = json.load(open("predictions/refD/results/runs.json"))["runs"]
byid = {r["id"]: r for r in runs}
for rid in sys.argv[1:]:
    job = dict(byid[rid + "_q1"])
    spec = rp.build_spec(job, duration_s=0.01)
    ensure_loaded()
    (chk,) = [make_block(b) for b in offline_block_specs(spec) if b.kind == "loop_stability_check"]
    p = chk.params
    fr = p["modes"][0]["f_hz"]
    f = np.array([fr * x for x in (0.5, 0.8, 0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1, 1.2, 1.5)])
    lc = open_loop(f, {**p, "discretisation": "continuous"})  # the frozen prediction's check
    ld = open_loop(f, {**p, "discretisation": "dsc1", "h_s": float(p.get("h_s", 2e-6))})
    fz = rs.fringe(spec)
    kd = fz["v0"] * fz["vis"]
    a, b, c, g = rs.loop_parts(spec, f, k_d=kd)
    a2, b2, _, _ = rs.loop_parts(spec, f + spec.run.f_sim_hz / 2, k_d=kd)
    lz = -g * c * (b * a + b2 * a2)
    al, bl, cl, _ = rs.loop_parts(spec, np.array([10.0]), k_d=kd)
    sgn = np.sign((1j * (-g * cl * bl * al))[0].real)
    lz = sgn * lz
    print(rid, "k0_chk", p["k0"], "kd", kd, "check result (current preset):", chk.result)
    print("sign", sgn)
    for fi, x, y, w in zip(f, lc, lz, ld):
        print(f"{fi:10.1f}  chk |L|={abs(x):.4f} ph={np.degrees(np.angle(x)):8.2f}   z |L|={abs(y):.4f} ph={np.degrees(np.angle(y)):8.2f}  chk_dsc1 |L|={abs(w):.4f} ph={np.degrees(np.angle(w)):8.2f}")
