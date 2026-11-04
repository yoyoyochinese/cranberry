"""Phase 6 axis 2, Procedure 2 / prompt §8: ablation of every toggle of the four nominal future presets (Liu lock, port
k = 31) and the §8 combination of the effects that mattered (runner.explore.combine_effects), on metric (ii) σ_total
and on the GBS-input V₋ of metric (ii).

  uv run python reviews/scripts/p6_axis2_ablate.py [--bases F_liu20 ...] [--duration 3] [--seed 1]

Writes runs/p6_axis2/ablate_<base>.json: baseline, per-toggle Δ of the quantities, the mattering set (|Δσ| > 5 %),
the combined value and the cumulative ladder. One seed (common random numbers: the toggle is the only difference).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
import warnings
from pathlib import Path

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

import dataclasses  # noqa: E402

from cranberry.runner.ablate import list_toggles  # noqa: E402
from cranberry.runner.run import run  # noqa: E402
from cranberry.runner.explore import combine_effects  # noqa: E402

Q = ("ii.port1.sigma_total_rad", "ii.port1.v_minus_db", "ii.port1.r_minus_db", "i.port1.sigma_total_rad",
     "eff.port1.eta_gbs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bases", nargs="*", default=list(lib.BASES))
    ap.add_argument("--duration", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    for base in a.bases:
        t0 = time.time()
        spec = lib.build(base, "base", duration_s=a.duration, seed=a.seed)
        toggles = list_toggles(spec)
        # a toggle whose "off" state needs parameters the preset does not have (e.g. src.pulsed off on the OPO
        # presets: the CW source needs f_s_hz, D34) is reported and left out
        abl, skipped = {"baseline": run(spec)}, {}
        for t in toggles:
            try:
                abl[t] = run(dataclasses.replace(spec, toggles={**spec.toggles, t: False}))
            except Exception as e:  # noqa: BLE001
                skipped[t] = f"{type(e).__name__}: {str(e)[:200]}"
        comb = combine_effects(spec, Q[0], ablation=abl, threshold=0.05)
        b = abl["baseline"].summary()
        out = {"base": base, "duration_s": a.duration, "seed": a.seed, "baseline": {q: b.get(q) for q in Q},
               "impact": {t: {q: r.summary().get(q, float("nan")) - b.get(q, float("nan")) for q in Q} | {"valid": r.valid}
                          for t, r in abl.items() if t != "baseline"},
               "mattered": comb["mattered"], "combined": comb["combined"], "residual": comb["residual"],
               "ladder": comb["ladder"], "skipped": skipped, "wall_s": time.time() - t0}
        lib.OUT.mkdir(parents=True, exist_ok=True)
        (lib.OUT / f"ablate_{base}.json").write_text(json.dumps(out, indent=1, default=float))
        print(f"{base}: {len(toggles)} toggles, mattered {comb['mattered']}, combined {comb['combined']:.4f} vs "
              f"{comb['baseline']:.4f} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
