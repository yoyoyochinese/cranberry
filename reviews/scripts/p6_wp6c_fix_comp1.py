"""Phase 5 review of the WP6-C fix (D98-4): what raises component 1 of metric (ii) from 9.8 (ideal point) to 22.2 mrad
(nominal point) on F_liu20 + lock/dual_band_slow, 30 s. Only the component-1 sources are on (as runner.decompose's
component-1 run); one variant per line. Read-only.

  uv run python reviews/scripts/p6_wp6c_fix_comp1.py --duration 30 --seed 1 [--variants nominal sqvac_off ...]

Variants (overlays after lock/dual_band_slow):
  ideal        lock/points/dual_band_slow_ideal (D92-4 point)
  nominal      lock/points/dual_band_slow_nominal (cal_frames 0.2, b' counted, μ 4.2, seed after the squeezer)
  sqvac_off    nominal with pcf1.squeezed_vacuum off (b' = 0), μ 4.2
  count_off    nominal with pcf1.counting_noise off (expected counts: the non-counting floor)
  seed_off     nominal with env.mech.seed.mech off (no SEED strain noise)
  ideal_mu42   ideal with μ 4.2 (μ effect alone)
  gated        nominal with dbs_cal_pump gated (pump off in the calibration frames: b' = 0, no seed / squeezed-vacuum
               interference), μ 3.43 (the optimum without background)
"""

from __future__ import annotations

import argparse
import importlib.util
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

from cranberry.runner.run import run  # noqa: E402

DBS = "lock/dual_band_slow"
NOM = (DBS, "lock/points/dual_band_slow_nominal")
IDEAL = (DBS, "lock/points/dual_band_slow_ideal")
VARIANTS = {"ideal": (IDEAL, {}, {}), "nominal": (NOM, {}, {}),
            "sqvac_off": (NOM, {"pcf1.squeezed_vacuum": False}, {}),
            "count_off": (NOM, {"pcf1.counting_noise": False}, {}),
            "seed_off": (NOM, {"env.mech.seed.mech": False}, {}),
            "ideal_mu42": (IDEAL, {}, {"dbs_mu_per_gate": 4.2}),
            "gated": (NOM, {}, {"dbs_cal_pump": "gated", "dbs_mu_per_gate": 3.43})}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    base = fs.load(fs.nominal("F_liu20"), run={"duration_s": 1.0}, overlays=NOM)
    groups = base.measure["decompose"]["groups"]
    only1 = {t: (c == "locked_path_residual") for c, ts in groups.items() for t in ts}
    for v in a.variants:
        ov, tog, par = VARIANTS[v]
        extra = ()
        if par:
            import yaml
            p = ROOT / f"reviews/scripts/_tmp_comp1_{v}_{a.seed}.yaml"
            p.write_text(yaml.safe_dump({"offline": {"params": par}}))
            extra = (str(p),)
        t0 = time.time()
        s = fs.load(fs.nominal("F_liu20"), run={"duration_s": a.duration, "seed": a.seed}, overlays=(*ov, *extra),
                    toggles={**only1, **tog})
        r = run(s)
        if extra:
            Path(extra[0]).unlink()
        pp = fs.params(r)
        std = 1e3 * float(r.metrics.metric_ii[1].stats.std)
        print(f"seed {a.seed} {v:11s}: comp-1 std {std:6.2f} mrad  (μ {pp['dbs_mu_per_gate']}, b' "
              f"{pp['dbs_sqvac_per_gate']:.3f}, f_gate·cal {pp['dbs_f_gate_hz'] * pp['dbs_cal_fraction']:.4g} Hz, "
              f"gates/window {pp['dbs_f_gate_hz'] * pp['dbs_cal_fraction'] * pp['dbs_t_int_s']:.0f}, valid {r.valid}) "
              f"({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
