"""Phase 6 review (WP6-AB D89-5, item 6): what `controller/dsc1`'s h = 1/f_sim (PID every base sample, 1 µs) does
to the stability verdict and |S|², against the device (DSC1 PID update 500 kHz, h = 2 µs, DSC1 p.2) — input to the
pending user decision "rate_div 2 / h = 2 µs".

A. Future presets (F_* + controller/dsc1; the normalised L-Liu loop is the same for every base): chk.lock_stability
   over f_u × f_res × Q (27 points; latency fixed at 5.7 µs by the overlay; TEC axes at the nominal point) for
     V0   as shipped: dsc1 law h = 1 µs, delay latency + 0.5/f_sim, no DSC1 I/O poles;
     V1   rate_div 2, τ compensated: h = 2 µs, delay 6.2 µs (tau_s − hold, total τ_L kept);
     V2   rate_div 2, τ not compensated: h = 2 µs, delay 6.7 µs (+ (r − 1)/(2 f_sim) hold, D68-2);
     Vdev the DSC1 model validated on Ref D: h = 2 µs, delay = τ_L = 5.7 µs, I/O poles 120 kHz / 100 kHz (refD.yaml:
          "5.7 µs with the 120 kHz / 100 kHz poles, 8.3 µs pure delay"; devices/dsc1.yaml adc/dac bandwidth);
     V0+p as shipped plus the I/O poles (isolates the missing-pole effect).
B. Ref D (refD_nominal + lab/acoustic + lab/multimode915b at the D89-4 fit point, h_s = 2 µs as validated): the same
   chk.lockD_stability with h = 1 µs, (i) ki_per_s, kd_s fixed (only the law's h changes) and (ii) the DSC1 entries
   dsc1_i, dsc1_d fixed (what the device would do if it really updated at 1 MHz).

Run: uv run python reviews/scripts/p6_wp6ab_dsc1_h.py
"""

from __future__ import annotations

import itertools
import math
import warnings

import numpy as np

from cranberry import config
from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import offline_block_specs
from cranberry.runner.sweep import expand

FF = np.geomspace(1e-6, 1e8, 400001)
FS = np.geomspace(300.0, 1e5, 20000)


def chk(spec, name):
    (b,) = [x for x in offline_block_specs(spec) if x.name == name]
    return dict(b.params)


def verdict(p):
    st, _, _ = nyquist_stable(FF, open_loop(FF, p), n_integrators(p))
    # high-frequency gain margin: max |L| at the −180° (mod 360°) crossings above 100 Hz (the TEC double integrator
    # makes a harmless low-frequency crossing with |L| ≫ 1 that nyquist_stable's max|L| reports)
    fh = FF[FF > 100.0]
    lh = open_loop(fh, p)
    ph = np.unwrap(np.angle(lh))
    k = np.flatnonzero(np.diff(np.floor((ph + np.pi) / (2 * np.pi))) != 0)
    mx, f180 = (float(np.abs(lh[k]).max()), float(fh[k][np.argmax(np.abs(lh[k]))])) if k.size else (0.0, math.nan)
    s2 = np.abs(1 / (1 + open_loop(FS, p))) ** 2
    j = int(np.argmax(s2))
    return st, mx, f180, s2[j], FS[j]


def future():
    spec = config.load({"extends": "presets/future/F_liu20_nominal.yaml"}, overlays=["controller/dsc1"])
    base = chk(spec, "chk.lock_stability")
    pp = spec.offline["params"]
    k0 = base["k0"]
    lat = float(pp["latency_s"])
    print(f"A. F_liu20_nominal + controller/dsc1: latency {lat * 1e6:.1f} µs, chk h_s {base['h_s'] * 1e6:.1f} µs, delay "
          f"{base['delay_s'] * 1e6:.2f} µs, poles {base['poles_hz']}")
    variants = {
        "V0 shipped (h 1 µs, 6.2 µs)": dict(h_s=1e-6, delay_s=lat + 0.5e-6, poles_hz=[]),
        "V0+p (+120/100 kHz poles)": dict(h_s=1e-6, delay_s=lat + 0.5e-6, poles_hz=[1.2e5, 1.0e5]),
        "V1 rate_div 2, τ comp. (6.2 µs)": dict(h_s=2e-6, delay_s=lat + 0.5e-6, poles_hz=[]),
        "V2 rate_div 2, no comp. (6.7 µs)": dict(h_s=2e-6, delay_s=lat + 1.0e-6, poles_hz=[]),
        "Vdev DSC1 (h 2 µs, 5.7 µs+poles)": dict(h_s=2e-6, delay_s=lat, poles_hz=[1.2e5, 1.0e5]),
    }
    grid = list(itertools.product([1e3, 5.477e3, 3e4], [5e3, 3.25e4, 6e4], [3.0, 17.32, 100.0]))
    res = {v: [] for v in variants}
    for fu, fr, q in grid:
        for v, upd in variants.items():
            p = {**base, **upd, "ki_per_s": 2 * math.pi * fu / k0, "modes": [{"f_hz": fr, "q": q}]}
            res[v].append(verdict(p))
    ref = [r[0] for r in res["Vdev DSC1 (h 2 µs, 5.7 µs+poles)"]]
    for v in variants:
        st = [r[0] for r in res[v]]
        flips = [(g, a) for g, a, b in zip(grid, st, ref) if a != b]
        print(f"  {v:34s} unstable {27 - sum(st):2d}/27; differs from Vdev at {len(flips)} points "
              f"({sum(1 for _, a in flips if a)} optimistic): {[(int(g[0]), int(g[1]), g[2]) for g, _ in flips][:8]}")
    print("  nominal-type points (f_u, f_res, Q): gain margin [dB] / max|S|² (f) per variant")
    for fu, fr, q in [(1e3, 2e4, 10.0), (5.477e3, 3.25e4, 17.32), (1e4, 2e4, 10.0), (3e4, 6e4, 3.0)]:
        cells = []
        for v, upd in variants.items():
            p = {**base, **upd, "ki_per_s": 2 * math.pi * fu / k0, "modes": [{"f_hz": fr, "q": q}]}
            st, mx, f180, s2, fs2 = verdict(p)
            cells.append(f"{'S' if st else 'U'} {-20 * math.log10(mx) if mx > 0 else float('inf'):5.1f} dB/{s2:6.2f}@{fs2 / 1e3:4.1f}k")
        print(f"    ({fu:6.0f}, {fr:6.0f}, {q:5.1f}): " + " | ".join(cells))


def refd():
    tree = {"extends": "presets/reference/refD_nominal.yaml", "overlays": ["lab/acoustic", "lab/multimode915b"]}
    point = {"acoustic_scale": 1.0, "acoustic_gamma": 1.0, "f_res_hz": 2.52e4, "q_res": 25.0, "mm2_f_hz": 5.4e3,
             "mm2_q": 15.0, "mm2_r": 0.11, "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07, "mm4_f_hz": 4.03e4,
             "mm4_q": 100.0, "mm4_r": 0.15}
    (pt,) = expand(config.load(tree), {f"offline.params.{k}": [v] for k, v in point.items()})
    p = chk(pt.spec, "chk.lockD_stability")
    h = p["h_s"]
    print(f"B. refD_nominal + lab/* (fit point): h_s {h * 1e6:.1f} µs, kp {p['kp']:.3e}, ki {p['ki_per_s']:.3e} /s, "
          f"kd {p['kd_s']:.3e} s, K_i h/2 = {p['ki_per_s'] * h / 2:.3e}")
    bands = {"5.3k": (5e3, 5.8e3), "10k": (9.3e3, 1.07e4), "24.9k": (2.2e4, 2.75e4), "40k": (3.7e4, 4.3e4)}
    f = np.geomspace(300.0, 6e4, 30000)
    for lab, q in [("validated h = 2 µs", p),
                   ("h = 1 µs, ki/kd fixed", {**p, "h_s": 1e-6}),
                   ("h = 1 µs, dsc1_i/dsc1_d fixed", {**p, "h_s": 1e-6, "ki_per_s": p["ki_per_s"] * h / 1e-6,
                                                     "kd_s": p["kd_s"] * 1e-6 / h})]:
        st, mx, f180 = nyquist_stable(FF, open_loop(FF, q), n_integrators(q))
        s2 = np.abs(1 / (1 + open_loop(f, q))) ** 2
        pk = {k: s2[(f > a) & (f < b)].max() for k, (a, b) in bands.items()}
        print(f"  {lab:32s} stable {st}, gain margin {-20 * math.log10(mx):5.1f} dB at {f180:.0f} Hz | |S|² "
              + ", ".join(f"{k} {v:.2f}" for k, v in pk.items()))


def main():
    warnings.simplefilter("ignore")
    ensure_loaded()
    future()
    refd()


if __name__ == "__main__":
    main()
