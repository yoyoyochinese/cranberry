"""Phase 5 M1 review: C23 / D63-7 -- does lock 2's gain (prop. to A+(P)) raise theta through the Bode integral?

Linear frequency-domain model of refA's lock 2 (no simulator): open-loop Delta noise
S_open(f) = |H_FC e^{-i w dtau} - 1|^2 h0/f^2 (laser through the LO filter cavity vs the squeezed light,
dtau = LO - sqz delay), loop L(f) = (2 pi f_u / (i w)) * e^{-i w tau_L} * LPF_demod(f) * ZOH(f)
with f_u = 20 kHz * A+(P)/A+(0.64 W) (preset: integral gain fixed, plant slope prop. to A+),
tau_L = 1 us (preset dl2), demod LPF 250 kHz first order, ZOH at f_sim = 20 MHz (negligible).
sigma^2(P) = int S_open |1/(1+L)|^2 df. Compared with the open-loop value (theta of the constant-theta curve).
Run: uv run python reviews/scripts/m1_lock2_bode.py
"""

import math

import numpy as np

H0, DNU, C = 200.0 / math.pi, 3.3e6, 299792458.0
TAU_LO = (5.0 * 1.468 + 0.3 * 1.0003) / C
TAU_SQ = (1.0003 * (1.0 + 0.3 + 0.2 + 1.0) + 0.045 * 2.2) / C
DTAU = TAU_LO - TAU_SQ
f = np.geomspace(1.0, 1e12, 2_000_001)
w = 2 * np.pi * f


def s_open(dnu):
    hfc = 1.0 / (1.0 + 2j * f / dnu)
    return np.abs(hfc * np.exp(-1j * w * DTAU) - 1.0) ** 2 * H0 / f**2


def sigma2(fu, dnu, tau_l=1e-6, f_lpf=2.5e5):
    loop = (2 * np.pi * fu / (1j * w)) * np.exp(-1j * w * tau_l) / (1 + 1j * f / f_lpf)
    s = s_open(dnu) / np.abs(1 + loop) ** 2
    return float(np.trapezoid(s * f, np.log(f)))


def a_plus(p):
    return math.exp(math.sqrt(9.06 * p))


for label, dnu in (("9 mrad preset (dnu_FC 3.3 MHz)", 3.3e6), ("14 mrad overlay (dnu_FC 1.175 MHz)", 1.175e6)):
    s_ol = float(np.trapezoid(s_open(dnu) * f, np.log(f)))
    print(f"{label}: open loop theta = {1e3 * math.sqrt(s_ol):.3f} mrad")
    base = None
    for p in (0.04, 0.2, 0.4, 0.64, 1.0, 1.2):
        fu = 2e4 * a_plus(p) / a_plus(0.64)
        s2 = sigma2(fu, dnu)
        base = base or s2
        print(f"   P = {1e3 * p:5.0f} mW  f_u = {fu / 1e3:6.1f} kHz  theta = {1e3 * math.sqrt(s2):.4f} mrad"
              f"  (theta/theta_open - 1 = {100 * (math.sqrt(s2 / s_ol) - 1):+.2f} %)")
