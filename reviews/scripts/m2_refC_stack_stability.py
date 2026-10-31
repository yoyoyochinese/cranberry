"""Phase 5 M2 Ref C review: linear stability of the three-stage stack over the refC sweep grid.

Continuous model of VT-20(c) (tests/vt/test_vt_20_refC_total_residual.py::_loop_noise_bandwidth):
L(s) = ω_u1/s·(1 + ω_12/s·M(s) + ω_12 ω_23/s²)·LPF(s)·e^{−sτ_L}, M = mirror resonance (r = 1), LPF = 2-pole
100 kHz demodulator. Delay by a (6,6) Padé approximant; closed-loop poles = roots of num + den.
Axes: f_u1, f_12, f_23, mirror_f_res, mirror_q, latency (refC.yaml sweep values).
Run: uv run python reviews/scripts/m2_refC_stack_stability.py
"""
import itertools
import math

import numpy as np
from numpy.polynomial import polynomial as P
from scipy.interpolate import pade

F_U1 = np.geomspace(1e3, 3e4, 4)
F_12 = np.geomspace(30, 3e3, 3)
F_23 = np.geomspace(0.1, 10, 3)
F_RES = np.linspace(1e4, 2e4, 3)
Q = np.geomspace(3, 30, 3)
TAU = (2.5e-7, 1e-6, 2e-6)
F_LPF = 1e5


def poly(c):  # coefficients highest power first -> numpy.polynomial (lowest first)
    return np.asarray(c, dtype=float)[::-1]


def closed_loop_unstable(fu1, f12, f23, fr, q, tau):
    wu, w12, w23, wr, wl = (2 * math.pi * v for v in (fu1, f12, f23, fr, F_LPF))
    # L = wu * [s²·Mden + w12·wr²·s + w12·w23·Mden] / (s³ · Mden) · wl²/(s + wl)² · e^{−sτ}
    mden = poly([1, wr / q, wr**2])
    s2 = poly([1, 0, 0])
    num = P.polyadd(P.polyadd(P.polymul(s2, mden), poly([w12 * wr**2, 0])), w12 * w23 * mden)
    num = wu * wl**2 * num
    den = P.polymul(P.polymul(poly([1, 0, 0, 0]), mden), P.polypow(poly([1, wl]), 2))
    # Padé of e^{-x}, x = sτ (scaled variable x to keep conditioning)
    pn, pd = pade([(-1) ** k / math.factorial(k) for k in range(13)], 6)
    pn = np.array([c * tau**k for k, c in enumerate(pn.coeffs[::-1])])
    pd = np.array([c * tau**k for k, c in enumerate(pd.coeffs[::-1])])
    char = P.polyadd(P.polymul(num, pn), P.polymul(den, pd))
    roots = P.polyroots(char)
    return np.max(roots.real) > 0, np.max(roots.real)


bad = []
n = 0
for pt in itertools.product(F_U1, F_12, F_23, F_RES, Q, TAU):
    n += 1
    uns, _ = closed_loop_unstable(*pt)
    if uns:
        bad.append(pt)
print(f"unstable grid points: {len(bad)} / {n}")
by = {}
for fu1, f12, f23, fr, q, tau in bad:
    by.setdefault((round(fu1), round(f12), round(q, 2)), set()).add((round(fr), tau))
for k in sorted(by):
    print(f"  f_u1={k[0]:>6} Hz  f_12={k[1]:>5} Hz  Q={k[2]:>5}: {len(by[k])} (f_res, τ_L) combos x 3 f_23")
print("nominal point (10 kHz, 300 Hz, 1 Hz, 15 kHz, Q 10, 1 µs):", closed_loop_unstable(1e4, 300, 1, 1.5e4, 10, 1e-6))
