"""Phase 6 axis 1: convergence checks (runner.convergence, PLAN §6.4, prompt §8) on the headline points:
chunk (bit-identical, G11) and f_sim ×2 (and ×4 where affordable), tolerance per quantity (passed when the change
≤ tolerance/3). The seed spread comes from the 10-seed sweeps (p6_axis1_sweep.py).

Usage: uv run python reviews/scripts/p6_axis1_convergence.py [case ...]
Writes docs/figures/phase6/axis1/convergence.json.
"""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

F_KEEP = ["pid0s", "tec0", "pid1s", "tec1"]   # F_common's TEC stages run at a fixed rate_div (10 ms): keep it physical
CASES = {
    # (point, quantities {summary key: tolerance}, f_sim factors, keep_rate)
    "A_koh_10m": (L.a_point(10.0), {"i.port0.sigma_total_rad": 1.0e-3, "i.port0.std_loop_rad": 1.0e-3}, (2,), None),
    "A_dx1_lab_10m": (L.a_point(10.0, dx1=True, acoustic=True),
                      {"i.port0.sigma_total_rad": 5.0e-3, "i.port0.std_loop_rad": 5.0e-3}, (2,), None),
    "A_dx1_lab_fcopt_10m": (L.a_point(10.0, dx1=True, acoustic=True, offset="fc_opt"),
                            {"i.port0.sigma_total_rad": 5.0e-3}, (2,), None),
    "A_dx1_lab_dsc1_10m": (L.a_point(10.0, dx1=True, acoustic=True, ctrl="bench/dsc1_refA"),
                           {"i.port0.sigma_total_rad": 5.0e-3}, (2,), None),
    "C_10m": (L.c_point(10.0), {"i.port0.std_rad": 0.3e-3}, (2,), None),
    "C_dsc1_100m": (L.c_point(100.0, ctrl="bench/dsc1_refC"), {"i.port0.std_rad": 0.3e-3}, (2,), None),
    "P_liu_lab_10m": (L.p_point(10.0, lab=True), {"ii.port1.sigma_total_rad": 3.0e-3, "ii.port1.std_loop_rad": 3.0e-3},
                      (2, 4), F_KEEP),
    "P_liu_lab_dsc1_100m": (L.p_point(100.0, ctrl="controller/dsc1", lab=True),
                            {"ii.port1.sigma_total_rad": 3.0e-3, "ii.port1.std_loop_rad": 3.0e-3}, (2, 4), F_KEEP),
    "P_het_lab_100m": (L.p_point(100.0, ["probe/detect_heterodyne"], lab=True),
                       {"ii.port1.sigma_total_rad": 3.0e-3}, (2,), F_KEEP),
    "P_P3_lab_100m": (L.p_point(100.0, ["probe/sep_wdm", "probe/freq_1THz", "gate/chopper_wdm"], lab=True,
                                axes={"offline.params.gate_f_ch_hz": 6.0e3, "offline.params.gate_jitter_s": 0.0}),
                      {"ii.port1.sigma_total_rad": 3.0e-3}, (2,), F_KEEP),
    "D_lab": (L.d_point(), {"lock.lockD.err.std": 2.0e-3}, (2,), None),
}


def one(case):
    from cranberry.runner.convergence import convergence

    pt, qs, fac, keep = CASES[case]
    spec = L.point(**pt)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            rep = convergence(spec, ("chunk", "f_sim"), quantities=qs, f_sim_factors=fac, keep_rate=keep)
        except Exception as e:  # noqa: BLE001 - reported
            return case, {"error": f"{type(e).__name__}: {e}"[:400]}
    return case, rep.checks


def main(cases):
    with ProcessPoolExecutor(max_workers=min(L.N_PROC, len(cases))) as ex:
        out = dict(ex.map(one, cases))
    for c, chk in out.items():
        print(c, json.dumps(chk, default=float)[:900])
    p = L.FIG / "convergence.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(out)
    L.FIG.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(old, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1:] or list(CASES))
