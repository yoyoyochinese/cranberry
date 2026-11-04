"""Phase 6 review (WP6-AB D89-4, item 5): are the 25 / 40 kHz bumps of the time-loop "|S|²" the DAC quantisation?

|S|²_time = PSD(locked err)/PSD(free err) also counts any in-loop noise the free run does not have. Same seed, S = 0,
refD_nominal + lab/acoustic + lab/multimode915b at the D89-4 point: quantisation of dsc1.dac on / off (and the ADC
quantisation / noise off as well), against the analytic |S|² and a white-quantisation estimate
|S·P|² (Δ²/12)/(f_u/2) / S_free with Δ = 8 V/2^12 (12-bit ±4 V, 500 kHz), P = k_act·M(f) (actuator path).

Run: uv run python reviews/scripts/p6_wp6ab_quant.py [--duration 1.5]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import warnings
from pathlib import Path

import numpy as np

from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import run

_s = importlib.util.spec_from_file_location("ab", Path(__file__).with_name("p6_ab_check.py"))
ab = importlib.util.module_from_spec(_s)
_s.loader.exec_module(ab)


def one(tg, duration, seed=20261031):
    out = []
    for fb in (True, False):
        spec = ab.spec_at(0.0, seed, duration, fb)
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **tg})
        r = run(spec)
        m = r.t_s >= 0.5
        f, p = ab.psd(np.asarray(r.series["lock.lockD.err"])[m])
        out.append(p)
    return f, out[0] / out[1], out[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=1.5)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    ensure_loaded()
    p = ab.check_params(ab.spec_at())
    bands = {"24.9k": (2.2e4, 2.75e4), "40k": (3.7e4, 4.3e4), "5.3k": (5.0e3, 5.8e3)}
    for lab, tg in [("all on", {}), ("DAC quantisation off", {"dsc1.dac.quantization": False}),
                    ("DAC + ADC quantisation and ADC noise off", {"dsc1.dac.quantization": False,
                                                                  "dsc1.adc.quantization": False, "dsc1.adc.noise": False})]:
        f, s2, free = one(tg, a.duration)
        cells = []
        for k, (lo, hi) in bands.items():
            mm = (f > lo) & (f < hi)
            cells.append(f"{k} {s2[mm].max():.2f}@{f[mm][np.argmax(s2[mm])] / 1e3:.2f}k")
        print(f"  {lab:42s} time-loop |S|²: " + ", ".join(cells))
    # analytic |S|² and white DAC quantisation estimate at the bump frequencies
    k_act = p["k0"] / 2.86  # rad per V at the DSC1 output (k0 = K_d · k_act · g, K_d = 2.86 V/rad)
    lsb = 8.0 / 2**12
    sq = lsb**2 / 12 / 2.5e5
    for x in (2.47e4, 4.01e4):
        L = open_loop(np.array([x]), p)[0]
        S = 1 / (1 + L)
        q = dict(p, kp=1.0, ki_per_s=0.0, kd_s=0.0, k0=k_act)  # actuator path k_act·M(f)·delay/poles
        P = open_loop(np.array([x]), q)[0]
        mm = (f > x / 1.03) & (f < x * 1.03)
        print(f"  {x / 1e3:.1f} kHz: analytic |S|² {abs(S)**2:.2f}; white DAC estimate |S P|² S_q / S_free = "
              f"{abs(S * P)**2 * sq / np.median(free[mm]):.2f} (S_free {np.median(free[mm]):.2e} rad²/Hz)")


if __name__ == "__main__":
    main()
