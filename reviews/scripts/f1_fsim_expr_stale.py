"""Phase 5 M1 framework review, finding F1-2: the f_sim convergence paths (runner.convergence('f_sim')
via sweep.expand, and core.rates.with_f_sim) do not re-evaluate expressions that depend on
run.f_sim_hz, so the D55-1 rule tau_s = τ_L − 1/f_sim is broken in every f_sim ×2 / ×4 run.

Run: uv run python reviews/scripts/f1_fsim_expr_stale.py
"""

import dataclasses

import numpy as np

from cranberry import config
from cranberry.config.spec import iter_block_specs
from cranberry.core.rates import with_f_sim
from cranberry.runner.sweep import expand


def delay_taus(spec):
    return {b.name: b.params["tau_s"] for _, b in iter_block_specs(spec.pipeline) if b.kind == "delay_buffer"}


def loop_total(tau_s, f_sim):
    return tau_s + 1.0 / f_sim  # D55-1 (explicitly modelled block delays aside)


if __name__ == "__main__":
    for preset, tau_L_expected in (("presets/reference/refD_nominal.yaml", None),
                                   ("presets/reference/refA.yaml", 1.0e-6)):
        rs = config.load(preset)
        f0 = rs.run.f_sim_hz
        base = delay_taus(rs)
        common = dataclasses.replace(rs, run=dataclasses.replace(rs.run, seed_mode="common"))
        print(f"{preset} (f_sim = {f0/1e6:g} MHz): tau_s = {base}")
        for k, p in zip((2, 4), expand(common, {"run.f_sim_hz": [2 * f0, 4 * f0]})):
            taus = delay_taus(p.spec)
            for name, t in taus.items():
                tot0, tot = loop_total(base[name], f0), loop_total(t, k * f0)
                print(f"  expand x{k}: {name}: tau_s = {t:.4g} s (unchanged: {t == base[name]}), "
                      f"tau_s + 1/f_sim = {tot*1e6:.4f} us vs {tot0*1e6:.4f} us at f_sim "
                      f"-> loop phase shift at 60 kHz {360*60e3*(tot0 - tot):.2f} deg")
        t2 = delay_taus(with_f_sim(rs, 2))
        print(f"  with_f_sim x2: tau_s = {t2} (unchanged: {t2 == base})")
    _ = np  # keep import for users extending the script
