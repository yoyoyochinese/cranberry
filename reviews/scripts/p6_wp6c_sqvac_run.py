"""Phase 5 review of WP6-C (D92-4/-5): component 1 of metric (ii) on F_liu20 + lock/dual_band_slow (nominal point) when
the calibration gate also holds the squeezed vacuum of the pumped slot (seed before the squeezer).

  uv run python reviews/scripts/p6_wp6c_sqvac_run.py --duration 30 --seed 1

The squeezed vacuum is entered as its equivalent Poisson background b' (p6_wp6c_counting.py section E: r = 1.8,
η·qe = 0.59–0.68 → b' = 1.68–1.77) through dbs_dark_cps, i.e. as a KNOWN, subtracted background (best case: the
block subtracts dark but not bg). μ is set to the re-optimised 4.2. Only the component-1 sources are on (as
runner.decompose's component-1 run). Read-only.
"""

from __future__ import annotations

import argparse
import importlib.util
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.runner.run import run  # noqa: E402

OV = ("lock/dual_band_slow", "lock/points/dual_band_slow_nominal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    base = fs.load(fs.nominal("F_liu20"), run={"duration_s": 1.0}, overlays=OV)
    groups = base.measure["decompose"]["groups"]
    only1 = {t: (c == "locked_path_residual") for c, ts in groups.items() for t in ts}
    t_gate = 1.0e-8
    for label, params in (("nominal (no background)", {}),
                          ("b' = 1.73, μ = 3.43", {"dbs_dark_cps": 1.73 / t_gate}),
                          ("b' = 1.73, μ = 4.2", {"dbs_dark_cps": 1.73 / t_gate, "dbs_mu_per_gate": 4.2})):
        # the point file is applied after the base preset's params, so pass the override as an overlay-level param
        s = fs.load(fs.nominal("F_liu20"), run={"duration_s": a.duration, "seed": a.seed}, overlays=OV,
                    toggles=only1, params=None)
        if params:
            p = ROOT / "reviews/scripts/_tmp_sqvac_point.yaml"
            import yaml
            p.write_text(yaml.safe_dump({"offline": {"params": params}}))
            s = fs.load(fs.nominal("F_liu20"), run={"duration_s": a.duration, "seed": a.seed}, overlays=(*OV, str(p)),
                        toggles=only1)
            p.unlink()
        r = run(s)
        pp = fs.params(r)
        print(f"{label:26s}: comp-1 std {1e3*float(r.metrics.metric_ii[1].stats.std):6.2f} mrad "
              f"(dark/gate {pp['dbs_dark_per_gate']:.3g}, μ {pp['dbs_mu_per_gate']}, valid {r.valid})", flush=True)


if __name__ == "__main__":
    main()
