"""M6 review: the F_liu20 metric (ii) decomposition of D74-3 (components, total, linearity, ledger part).
Usage: uv run python reviews/scripts/m6_decompose_liu20.py [base]
"""
import importlib.util
import math
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.decompose import decompose_metric_ii  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

base = sys.argv[1] if len(sys.argv) > 1 else "F_liu20"
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    spec = fs.load(fs.nominal(base), run={"duration_s": 0.03})
    c = decompose_metric_ii(spec)[1]
    full = run(spec)
print(f"{base}: status {c.status}, linearity {c.linearity_ratio:.4f}")
for k, v in c.variances_rad2.items():
    print(f"  {k:28s} {v:.4e} rad^2  {1e3*math.sqrt(max(v, 0)):7.2f} mrad")
s = sum(c.variances_rad2.values())
loop = full.metrics.metric_ii[1].stats.var
hf = full.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
print(f"  sum of components {s:.4e} ({1e3*math.sqrt(s):.1f} mrad); full run: loop var {loop:.4e}, ledger {hf:.4e}, "
      f"std(loop) {1e3*math.sqrt(loop):.1f} mrad, sqrt(loop+ledger) {1e3*math.sqrt(loop+hf):.1f} mrad")
print("  other attributes:", {k: v for k, v in vars(c).items() if k not in ("variances_rad2",)})

if "--totals" in sys.argv:   # metric (ii) totals of all four bases (D74-3 quotes sigma = 101-154 mrad)
    for b in ("F_wg100", "F_wg6", "F_bor6", "F_liu20"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = run(fs.load(fs.nominal(b), run={"duration_s": 0.03}))
        st = r.metrics.metric_ii[1].stats
        hf = r.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
        print(f"{b}: std {1e3*st.std:.1f} mrad, mean {1e3*st.mean:.1f} mrad, sqrt(var+ledger) "
              f"{1e3*math.sqrt(st.var+hf):.1f} mrad, rms incl. mean {1e3*math.sqrt(st.var+hf+st.mean**2):.1f} mrad, "
              f"V- {r.metrics.noise_levels['ii.port1']['v_minus_db']:.2f} dB")
