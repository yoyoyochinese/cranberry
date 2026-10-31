"""Phase 5 M2 Ref C review: (1) the preset's closed-form QNL band integrals (refC.yaml qnl_*) against direct
quadrature of S_x = A + C cos 2Δ (T14 §6.1) and against qnl_band_integrals; (2) the effect of the BHD_B
dip (r = 0.5, f_dip = 0.83 GHz, |H|² on the optical noise, dark noise flat as in T13 C(f) = |H|²C0) on the
sensor-limited phase PSD ∝ Q/ΔP² (square-law form; the log law scales the same at Δ = 0) — the preset
assumes a flat |H| in 8–300 MHz.
Run: uv run python reviews/scripts/m2_refC_qnl_integrals.py
"""
import math

import numpy as np
from scipy.integrate import quad

from cranberry.blocks.err.qnl import qnl_band_integrals

eta, x, fwhm, f1, f2, cl = 0.79, 0.56, 1.75e9, 8e6, 3e8, 11.0
fh = fwhm / 2
sel = 10 ** (-cl / 10)
rm = lambda f: 1 - eta * 4 * x / ((1 + x) ** 2 + (f / fh) ** 2)  # noqa: E731
rp = lambda f: 1 + eta * 4 * x / ((1 - x) ** 2 + (f / fh) ** 2)  # noqa: E731

# closed forms copied from refC.yaml
at = math.atan
I1 = (at(f2 / fh / (1 - x)) - at(f1 / fh / (1 - x))) / (1 - x)
I2 = (at(f2 / fh / (1 + x)) - at(f1 / fh / (1 + x))) / (1 + x)
J = lambda a: (f2 / fh) / (2 * a**2 * (a**2 + (f2 / fh) ** 2)) - (f1 / fh) / (2 * a**2 * (a**2 + (f1 / fh) ** 2)) + (at(f2 / fh / a) - at(f1 / fh / a)) / (2 * a**3)  # noqa: E731
J11, J22 = J(1 - x), J(1 + x)
J12 = (I1 - I2) / ((1 + x) ** 2 - (1 - x) ** 2)
k, c = 2 * eta * x, 1 + sel
closed = {"p_mean": c * (f2 - f1) + k * fh * (I1 - I2), "p_delta": k * fh * (I1 + I2), "p_dark": sel * (f2 - f1),
          "q0": c**2 * (f2 - f1) + 2 * c * k * fh * (I1 - I2) + k**2 * fh * (J11 - 2 * J12 + J22),
          "q1": -2 * c * k * fh * (I1 + I2) - 2 * k**2 * fh * (J11 - J22), "q2": k**2 * fh * (J11 + 2 * J12 + J22)}


def direct(h2=lambda f: 1.0):
    A = lambda f: h2(f) * (rp(f) + rm(f)) / 2 + sel  # noqa: E731
    C = lambda f: -h2(f) * (rp(f) - rm(f)) / 2  # noqa: E731
    q = lambda g: quad(g, f1, f2, limit=200, epsrel=1e-12)[0]  # noqa: E731
    return {"p_mean": q(A), "p_delta": -q(C), "p_dark": q(lambda f: sel), "q0": q(lambda f: A(f) ** 2),
            "q1": q(lambda f: 2 * A(f) * C(f)), "q2": q(lambda f: C(f) ** 2)}


d = direct()
f = np.linspace(f1, f2, 200001)
tb = qnl_band_integrals(f, np.ones_like(f), rm(f), rp(f), np.full_like(f, sel), np.ones_like(f))
print("closed form vs quad vs qnl_band_integrals (relative):")
for key in closed:
    print(f"  {key:7s} {closed[key]: .6e}  quad rel {closed[key]/d[key]-1: .1e}  trapz rel {closed[key]/tb[key]-1: .1e}")

tau = 0.5 / 0.83e9
for r in (0.3, 0.5, 0.9):
    hd = lambda f, r=r: abs((1 + r * np.exp(-2j * np.pi * f * tau)) / (1 + r)) ** 2  # noqa: E731
    dd = direct(hd)
    q_at0 = dd["q0"] + dd["q1"] + dd["q2"]
    q_flat = d["q0"] + d["q1"] + d["q2"]
    s_ratio = (q_at0 / dd["p_delta"] ** 2) / (q_flat / d["p_delta"] ** 2)
    print(f"dip r={r}: |H|²(300 MHz) = {hd(3e8):.3f}; sensor-limited phase PSD (Δ=0) × {s_ratio:.3f} "
          f"→ rms × {math.sqrt(s_ratio):.3f}; ΔP × {dd['p_delta']/d['p_delta']:.3f}")
