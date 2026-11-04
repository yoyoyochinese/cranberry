"""Phase 6 axis 1: wall time and the summary keys of one short run per base (planning the seed budget).
Usage: uv run python reviews/scripts/p6_axis1_timing.py
"""
import sys

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from p6_axis1_lib import job  # noqa: E402

BASES = {
    "A": dict(base="presets/bench/refA_fiber_bench.yaml", duration_s=0.05,
              axes={"offline.params.L_sq_fiber_m": 10.0, "offline.params.lo_route": "fiber_matched",
                    "offline.params.eta_fiber_out": 0.98},
              extra={"pipeline": [{"name": "mi", "params": {"t_settle_s": 0.015}}]}),
    "C": dict(base="presets/bench/refC_fiber_bench.yaml", duration_s=0.1, axes={"offline.params.L_fiber_m": 10.0}),
    "P": dict(base="presets/bench/pulsed_wg_1port.yaml", duration_s=0.05, axes={"offline.params.L_port_m": 10.0}),
    "D": dict(base="presets/bench/refD_lab.yaml", duration_s=0.1,
              axes={"offline.params.acoustic_hf_slope": 2.0, "offline.params.acoustic_rho_arms": 0.0,
                    "offline.params.dac_output_noise_v": 0.0}),
}
if __name__ == "__main__":
    for k in sys.argv[1:] or BASES:
        r = job({"label": k, "point": BASES[k]})
        print(k, f"wall {r['wall_s']:.1f} s", r.get("error", ""))
        for kk, v in r.items():
            if kk not in ("psd", "measurements") and not kk.startswith("series."):
                print("   ", kk, v)
