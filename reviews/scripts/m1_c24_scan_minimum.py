"""Phase 5 M1 review: C24 -- depth of the Ref A Fig. 3 phase-scan minima vs VBW (independent of esa.py).

V(D)/V_SN = R- cos^2 D + R+ sin^2 D (L = 8 %, r = sqrt(9.06*0.64), theta = 0), D = w t, w = 71 rad/s
(minima 44 ms apart = pi rad, T04 §5.4). Video filter: one-pole (-3 dB at VBW) in the power or log
domain (no radiometer noise: its mean effect is the -2.51 dB log bias, which cancels against the shot
trace). Also: the scan speed needed for a -9.5 dB minimum at VBW 100 Hz, and the VBW needed at 71 rad/s.
Run: uv run python reviews/scripts/m1_c24_scan_minimum.py
"""

import math

import numpy as np
from scipy.optimize import brentq
from scipy.signal import lfilter

L, R = 0.08, math.sqrt(9.06 * 0.64)
RM, RP = L + (1 - L) * math.exp(-2 * R), L + (1 - L) * math.exp(2 * R)


def min_depth_db(vbw, w=71.0, domain="power", fs=2e6, periods=3):
    t = np.arange(int(periods * math.pi / w * fs)) / fs
    v = RM * np.cos(w * t) ** 2 + RP * np.sin(w * t) ** 2
    a = math.exp(-2 * math.pi * vbw / fs)
    x = v if domain == "power" else 10 * np.log10(v)
    y = lfilter([1 - a], [1, -a], x, zi=[x[0] * a])[0]
    y = y[len(y) // 3:]  # settled
    return float(10 * np.log10(y.min()) if domain == "power" else y.min())


print(f"locked R- = {10 * math.log10(RM):.2f} dB, R+ = {10 * math.log10(RP):.2f} dB")
for vbw in (100.0, 300.0, 1e3, 3e3, 1e4):
    print(f"  VBW {vbw:7.0f} Hz: minimum {min_depth_db(vbw):+6.2f} dB (power) {min_depth_db(vbw, domain='log'):+6.2f} dB (log)")
w_need = brentq(lambda w: min_depth_db(100.0, w=w) + 9.5, 0.1, 71.0)
print(f"  at VBW 100 Hz a -9.5 dB minimum (power) needs w <= {w_need:.2f} rad/s, i.e. a pi period >= {math.pi / w_need:.2f} s"
      f" (Fig. 3: 0.044 s)")
d100 = min_depth_db(100.0)
assert 1.0 < d100 < 2.5, d100  # T04 §5.4: +1.76 dB
