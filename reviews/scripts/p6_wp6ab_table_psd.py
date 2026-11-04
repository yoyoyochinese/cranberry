"""Phase 6 review (WP6-AB D89-2/3): the table-shaped acoustic noise (`segment_mech_noise` mech_shape: table).

1. Per-arm split in the time loop: refD_nominal + lab/acoustic + lab/multimode915b, feedback off. One-sided Welch PSD
   (conv §7) of the two arm phases psi.a1.pr.obs, psi.a2.lo.obs and of the lock combination lock.lockD.err against
   the table T: S_a1 + S_a2 = T (γ = 1: S_i = (L_i/13.2 m) T; γ = 0: T/2 each) and S_combo = S_a1 + S_a2 (ρ = 0).
2. The extrapolated part: variance of T above 3 / 10 / 30 kHz and above f_N, for Ref D (13.2 m sum) and for the
   F_* spools of lab/acoustic_future (γ = 1: (L/13.2 m)·T per fibre), against an alternative steeper roll-off; and the
   part above the block Nyquist frequency that is both in the exactly-sampled OU series and in the ledger row.

Run: uv run python reviews/scripts/p6_wp6ab_table_psd.py [--duration 2.0]
"""

from __future__ import annotations

import argparse
import dataclasses
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.signal import welch

from cranberry import config
from cranberry.blocks.noise._lgs import loglog_interp, section_psd, table_fit
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

TABLE = [[10.0, 6.32e-4], [100.0, 2.0e-5], [300.0, 3.6e-6], [1000.0, 6.7e-7], [3000.0, 9.0e-8], [3.0e4, 9.0e-10]]
POINT = {"acoustic_scale": 1.0, "f_res_hz": 2.52e4, "q_res": 25.0, "mm2_f_hz": 5.4e3, "mm2_q": 15.0, "mm2_r": 0.11,
         "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07, "mm4_f_hz": 4.03e4, "mm4_q": 100.0, "mm4_r": 0.15}
SLOTS = ["psi.a1.pr.obs", "psi.a2.lo.obs", "lock.lockD.err"]


def split(duration, gamma):
    tree = {"extends": "presets/reference/refD_nominal.yaml", "overlays": ["lab/acoustic", "lab/multimode915b"],
            "run": {"duration_s": duration},
            "measure": {"record": {"slots": SLOTS, "decim": 2, "antialias": "fir"}}}
    (pt,) = expand(config.load(tree), {f"offline.params.{k}": [v] for k, v in {**POINT, "acoustic_gamma": gamma}.items()})
    spec = dataclasses.replace(pt.spec, toggles={**pt.spec.toggles, "dsc1.pid.feedback": False})
    r = run(spec)
    p = spec.offline["params"]
    m = np.asarray(r.t_s) >= 0.3
    T = loglog_interp(TABLE)
    out = {}
    for s in SLOTS:
        f, P = welch(np.asarray(r.series[s], float)[m], 5e5, nperseg=1 << 15, detrend="linear")
        out[s] = (f, P)
    f = out[SLOTS[0]][0]
    print(f"γ = {gamma}: L_a1 = {p.get('L_arm1_m', '?')}, ΔL = {p.get('dL_m', '?')}")
    for x in (100.0, 300.0, 1e3, 3e3, 1e4):
        k = (f > x / 1.1) & (f < x * 1.1)
        a1, a2, c = (float(np.median(out[s][1][k])) for s in SLOTS)
        print(f"  {x:7.0f} Hz: S_a1/T {a1 / T(x):.3f}  S_a2/T {a2 / T(x):.3f}  (S_a1+S_a2)/T {(a1 + a2) / T(x):.3f}  "
              f"S_combo/(S_a1+S_a2) {c / (a1 + a2):.3f}")


def extrapolation():
    T = loglog_interp(TABLE)
    print("table fit (analytic OU sum / T) at 4 per decade:",
          {x: round(float(sum(section_psd(s, np.array([x]))[0] for s in table_fit(T, 10.0, 3e4, 4)) / T(x)), 3)
           for x in (30.0, 100.0, 1e3, 1e4, 2e4, 3e4, 1e5)})

    def var_above(f1, steep_from=None):
        g = (lambda f: T(f)) if steep_from is None else (lambda f: T(f) * (1.0 if f < steep_from else (steep_from / f) ** 2))
        return quad(g, f1, 1e8, limit=400, points=[3e3, 3e4, 1e5])[0]

    print("variance of the Ref D sum table T above f1 [rad²] (f⁻² tail | tail steepened to f⁻⁴ above 10 kHz):")
    for f1 in (3e3, 1e4, 3e4, 5e5):
        print(f"  f1 = {f1:8.0f} Hz: {var_above(f1):.3e} | {var_above(f1, 1e4):.3e}")
    for name, L in (("F_liu20 port 1 (1268 m)", 1268.2), ("F_wg100 port 31 (3170 m)", 3170.0), ("F_wg6 / F_bor6 port 31 (1060 m)", 1060.0)):
        s = L / 13.2
        print(f"  lab/acoustic_future γ = 1, {name}: ×{s:.0f} → above 3 kHz {s * var_above(3e3):.3e} rad² "
              f"({1e3 * np.sqrt(s * var_above(3e3)):.0f} mrad), above 10 kHz {1e3 * np.sqrt(s * var_above(1e4)):.0f} mrad, "
              f"above f_N = 500 kHz (counted twice: OU samples + ledger) {1e3 * np.sqrt(s * var_above(5e5)):.1f} mrad; "
              f"f⁻⁴ above 10 kHz: above 3 kHz {1e3 * np.sqrt(s * var_above(3e3, 1e4)):.0f} mrad")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=2.0)
    ap.add_argument("--no-time-loop", action="store_true")
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    ensure_loaded()
    extrapolation()
    if not a.no_time_loop:
        split(a.duration, 1.0)
        split(a.duration, 0.0)


if __name__ == "__main__":
    main()
