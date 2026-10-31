"""M6 review: per-port loss budget (eta_gbs / eta_bhd, every registered element) and stray photons per mode of the
four future base presets at their nominal points (port k = 31), plus the path_build overlays NF / MF / LL on F_liu20,
compared with T12 §4.2 / D38 / D51-3 (-1.5 / -4.7 / -7.3 dB one demux stage; LL with a 5-stage tree ~ -6.1 dB).
Usage: uv run python reviews/scripts/m6_loss_stray.py
"""
import importlib.util
import math
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)


def db(x):
    return 10 * math.log10(x)


def report(preset, overlays=(), params=None, detail=True):
    res = fs.build_check(preset, overlays, params=params)
    reg = res.offline.loss_registry
    eff = res.metrics.efficiency
    print(f"\n== {Path(preset).stem} {'+'.join(Path(o).stem for o in overlays)}")
    for k in (0, 1):
        print(f"  port {k}: eta_gbs = {eff[k].eta_gbs:.4f} ({db(eff[k].eta_gbs):+.2f} dB), eta_bhd = {eff[k].eta_bhd:.4f} "
              f"({db(eff[k].eta_bhd):+.2f} dB)")
    if detail:
        for path in sorted(reg.paths):
            for metric in ("gbs",):
                els = reg.elements_on(path, metric)
                if not els:
                    continue
                print(f"  {path} [{metric}]: " + ", ".join(f"{e.id}={e.eta:.4f}" for e in els))
    try:
        s = fs.stray_light(res)
        print("  stray (P̄ = %.1e W): " % fs.params(res)["p_pr_avg_w"] + ", ".join(f"{k}={v:.3e}" for k, v in s.items()
                                                                                  if isinstance(v, float)))
    except Exception as exc:  # noqa: BLE001
        print("  stray: n/a", exc)
    return res


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for b in ("F_wg100", "F_wg6", "F_bor6", "F_liu20"):
        report(fs.nominal(b))
    for pb in ("NF", "MF", "LL"):
        report(fs.nominal("F_liu20"), (f"presets/axes/path_build/{pb}.yaml",), detail=False)
