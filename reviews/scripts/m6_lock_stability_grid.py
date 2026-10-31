"""M6 review: linear stability of the F_common L-Liu fine loop over its sweep axes (no CHK in the future presets).

F_common: ki = 2 pi f_u / (K k_pzt) with plant K k_pzt  =>  L(s) = (2 pi f_u / s) * w_r^2/(s^2 + s w_r/Q + w_r^2) * e^{-s tau},
tau = latency_s (tau_s + 1/f_sim, D55-1) + 0.5/f_sim (DAC ZOH, D68-2) at f_sim = 1 MHz. The TEC offload stage
(f_12 <= 0.3 Hz) does not matter at the resonance. Axes (F_common.yaml): f_u {1e3, 5.48e3, 3e4}, pzt_f_res
{5e3, 32.5e3, 6e4}, pzt_q {3, 17.3, 100}, latency {1, 2} us.  Uses the project's own loop_stability_check helpers.
Also: the closed-loop runs of m6_laser_decorrelation.py (f_u = 10 / 30 kHz at the nominal 20 kHz, Q 10 resonance)
return metric (ii) 80-100x the open-loop value with valid = True (validity watches only saturation / lock loss).
Usage: uv run python reviews/scripts/m6_lock_stability_grid.py
"""
import itertools

import numpy as np

from cranberry.blocks.chk.loop_stability import nyquist_stable, open_loop

f = np.geomspace(1e-2, 1e8, 400001)
rows, n_unst = [], 0
grid = list(itertools.product([1e3, np.sqrt(1e3 * 3e4), 3e4], [5e3, 32.5e3, 6e4], [3.0, np.sqrt(300.0), 100.0], [1e-6, 2e-6]))
for fu, fr, q, lat in grid:
    p = {"k0": 1.0, "kp": 0.0, "ki_per_s": 2 * np.pi * fu, "kd_s": 0.0, "delay_s": lat + 0.5e-6, "poles_hz": [],
         "modes": [{"f_hz": fr, "q": q}]}
    st, mx, f180 = nyquist_stable(f, open_loop(f, p), 1)
    n_unst += not st
    rows.append((fu, fr, q, lat, st, mx, f180))
print(f"unstable: {n_unst} / {len(grid)} grid points of (f_u, f_res, Q, latency)")
for fu in sorted({r[0] for r in rows}):
    sub = [r for r in rows if r[0] == fu]
    print(f"  f_u = {fu:7.0f} Hz: unstable {sum(not r[4] for r in sub):2d} / {len(sub)}")
nom = [r for r in rows if r[0] == 1e3 and r[1] == 32.5e3]
p_nom = {"k0": 1.0, "kp": 0.0, "ki_per_s": 2 * np.pi * 1e3, "kd_s": 0.0, "delay_s": 1.5e-6, "poles_hz": [],
         "modes": [{"f_hz": 2e4, "q": 10.0}]}
st, mx, f180 = nyquist_stable(f, open_loop(f, p_nom), 1)
print(f"nominal (f_u 1 kHz, 20 kHz, Q 10, 1 us): stable={st}, max|L| at -180 = {mx:.3f} at {f180:.0f} Hz (GM {1/mx:.2f})")
for fu in (1e4, 3e4):
    p_nom["ki_per_s"] = 2 * np.pi * fu
    st, mx, f180 = nyquist_stable(f, open_loop(f, p_nom), 1)
    print(f"nominal resonance, f_u = {fu:.0f}: stable={st}, max|L| at -180 = {mx:.2f}")
    assert not st
assert n_unst > len(grid) // 3
print("self-checks passed")
