"""Phase 6 review (WP6-AB D89-4, item 5): is the 5.3 kHz side mode identifiable from |S|²?

At refD_nominal + lab/multimode915b (D89-4 fit), replace the 5.4 kHz side mode (r = 0.11, Q 15) by
  (a) an extra loop delay Δτ,  (b) an extra real pole f_p  — scanned for the best |S|² peak near 5.3 kHz —
  (c) a series biquad (zero pair / pole pair) placed anywhere in the loop (e.g. the 914 drive or the detection),
and compare |S|² with the measured values (COMPARISON.md §2) at S = 0 / 0.5 / 1.5 / 2 V (K_d ∝ √(V0² − S²)).
(c) is built to equal the parallel sum's factor near 5.4 kHz: with the main mode ≈ flat there,
Σ r_i H_i ≈ (1 − r_2) + r_2 H_2 = biquad with poles at (f_2, Q_2) and zeros at f_2/√(1 − r_2).

Run: uv run python reviews/scripts/p6_wp6ab_multimode_alt.py
"""

from __future__ import annotations

import math
import warnings

import numpy as np

from cranberry import config
from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import offline_block_specs
from cranberry.runner.sweep import expand

POINT = {"acoustic_scale": 1.0, "acoustic_gamma": 1.0, "f_res_hz": 2.52e4, "q_res": 25.0, "mm2_f_hz": 5.4e3,
         "mm2_q": 15.0, "mm2_r": 0.11, "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07, "mm4_f_hz": 4.03e4,
         "mm4_q": 100.0, "mm4_r": 0.15}
MEAS = {"300": [0.020, 0.013, 0.028, 0.032], "1k": [0.16, 0.20, 0.28, 0.43], "3k": [0.62, 0.62, 0.80, 1.0],
        "5.3k": [13.6, 9.2, 7.0, 3.8], "8k": [1.58, 1.77, 1.59, 1.39], "12k": [1.60, 1.57, 1.48, 1.55]}
SP = [0.0, 0.5, 1.5, 2.0]
V0 = 2.86
F = np.geomspace(200.0, 6e4, 40000)


def params():
    tree = {"extends": "presets/reference/refD_nominal.yaml", "overlays": ["lab/acoustic", "lab/multimode915b"]}
    (pt,) = expand(config.load(tree), {f"offline.params.{k}": [v] for k, v in POINT.items()})
    (b,) = [x for x in offline_block_specs(pt.spec) if x.name == "chk.lockD_stability"]
    return dict(b.params)


def summary(L):
    out = {}
    for i, s in enumerate(SP):
        s2 = np.abs(1 / (1 + L * math.sqrt(V0**2 - s**2) / V0)) ** 2
        m = (F > 4.5e3) & (F < 6.5e3)
        k = np.argmax(s2[m])
        half = s2[m] > s2[m][k] / 2
        fw = F[m][half].max() - F[m][half].min()
        out[s] = {"pk": s2[m][k], "fpk": F[m][k], "fwhm": fw,
                  **{n: s2[np.argmin(abs(F - x))] for n, x in (("300", 300), ("1k", 1e3), ("3k", 3e3), ("8k", 8e3), ("12k", 1.2e4))}}
    return out


def show(lab, L):
    o = summary(L)
    print(f"  {lab}")
    for n in ("5.3k", "300", "1k", "3k", "8k", "12k"):
        key = "pk" if n == "5.3k" else n
        print(f"     {n:>5}: " + "  ".join(f"{o[s][key]:6.3g} ({MEAS[n][i]:.3g})" for i, s in enumerate(SP))
              + (f"   peak f {o[0.0]['fpk']:.0f} Hz → {o[2.0]['fpk']:.0f} Hz (S 0 → 2 V), FWHM {o[0.0]['fwhm']:.0f} Hz"
                 if n == "5.3k" else ""))


def main():
    warnings.simplefilter("ignore")
    ensure_loaded()
    p = params()
    show("D89-4 fit (parallel-sum side mode 5.4 kHz, r 0.11, Q 15)", open_loop(F, p))
    m = [dict(x) for x in p["modes"]]
    no2 = [dict(m[0], r=m[0]["r"] + m[1]["r"]), m[2], m[3]]  # side mode removed, Σr = 1 kept
    q0 = {**p, "modes": no2}
    show("side mode removed", open_loop(F, q0))
    best = None
    for dt in np.linspace(0, 30e-6, 301):
        o = summary(open_loop(F, {**q0, "delay_s": p["delay_s"] + dt}))
        err = abs(math.log(o[0.0]["pk"] / 13.6)) + 3 * abs(math.log(o[0.0]["fpk"] / 5.35e3))
        best = min(best or (1e9,), (err, dt))
    show(f"(a) extra delay Δτ = {best[1] * 1e6:.1f} µs (best peak near 5.3 kHz)",
         open_loop(F, {**q0, "delay_s": p["delay_s"] + best[1]}))
    best = None
    for fp in np.geomspace(1e3, 2e5, 400):
        o = summary(open_loop(F, {**q0, "poles_hz": [*p["poles_hz"], fp]}))
        err = abs(math.log(o[0.0]["pk"] / 13.6)) + 3 * abs(math.log(o[0.0]["fpk"] / 5.35e3))
        best = min(best or (1e9,), (err, fp))
    show(f"(b) extra pole f_p = {best[1]:.0f} Hz (best peak near 5.3 kHz)",
         open_loop(F, {**q0, "poles_hz": [*p["poles_hz"], best[1]]}))
    w2, q2, r2 = 2 * math.pi * 5.4e3, 15.0, 0.11
    wz = w2 / math.sqrt(1 - r2)
    s = 2j * math.pi * F
    biq = (1 - r2) * (s * s + s * w2 / q2 * 1.0 + wz * wz) / (s * s + s * w2 / q2 + w2 * w2) * (w2 * w2) / (wz * wz) / (1 - r2)
    biq = ((1 - r2) * (s * s + s * w2 / q2) + w2 * w2) / (s * s + s * w2 / q2 + w2 * w2)
    show("(c) series biquad ((1−r)(s² + sω/Q) + ω²)/(s² + sω/Q + ω²) anywhere in the loop, main mode r 0.78",
         open_loop(F, q0) * biq)
    from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable
    ff = np.geomspace(1e-2, 1e7, 400001)
    print("  sensitivity to mm2_r (main r = 1 − Σ side r), S = 0: r → (stable, max|L|@−180°, |S|² peak, FWHM)")
    for r in (0.06, 0.08, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15):
        mm = [dict(m[0], r=1 - r - m[2]["r"] - m[3]["r"]), dict(m[1], r=r), m[2], m[3]]
        q = {**p, "modes": mm}
        st, mx, _ = nyquist_stable(ff, open_loop(ff, q), n_integrators(q))
        o = summary(open_loop(F, q))[0.0]
        print(f"     r {r:.2f}: stable {st}, max|L| {mx:.3f}, |S|² {o['pk']:.1f}, FWHM {o['fwhm']:.0f} Hz")
    print("  height 13.6 AND FWHM ≈ 300 Hz (COMPARISON §2) over mm2_q × mm2_r (S = 0): best few")
    rows = []
    for q2 in (3.0, 5.0, 8.0, 10.0, 15.0, 25.0, 40.0):
        for r in np.linspace(0.02, 0.4, 77):
            mm = [dict(m[0], r=1 - r - m[2]["r"] - m[3]["r"]), dict(m[1], r=r, q=q2), m[2], m[3]]
            o = summary(open_loop(F, {**p, "modes": mm}))
            rows.append((abs(math.log(o[0.0]["pk"] / 13.6)) + abs(math.log(o[0.0]["fwhm"] / 300.0)), q2, r, o[0.0]["pk"],
                         o[0.0]["fwhm"], o[2.0]["pk"]))
    for row in sorted(rows)[:4]:
        print(f"     Q2 {row[1]:4.0f}, r {row[2]:.3f}: |S|² {row[3]:.1f}, FWHM {row[4]:.0f} Hz, at S = 2 V {row[5]:.1f} (meas 3.8)")


if __name__ == "__main__":
    main()
