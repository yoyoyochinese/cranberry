"""Phase 5 M2 Ref C review: what total phase jitter σ do Ref C's own zero-span numbers (Fig. 4) imply?

Closed-form OPO spectra (T02, conv §4): R± = 1 ± η·4x/((1 ∓ x)² + (f/f_h)²), f_h = FWHM/2 = 875 MHz,
jitter: V′ = (R− + R+)/2 ± (R− − R+)/2·e^{−2σ²} (conv §6, Gaussian).
Data [Ref C Fig. 4]: anti/sqz at 30 MHz, 0.25, 0.69, 1.21 GHz = 10.0/−4.7, 8.5/−4.7, 5.2/−3.0, 2.8/−2.2 dB.
Run: uv run python reviews/scripts/m2_refC_fig4_fit.py
"""
import numpy as np
from scipy.optimize import brentq, least_squares

F = np.array([3e7, 2.5e8, 6.9e8, 1.21e9])
ANTI = np.array([10.0, 8.5, 5.2, 2.8])
SQZ = np.array([-4.7, -4.7, -3.0, -2.2])
FH = 0.875e9
M = np.pi * 0.1 / 4  # dither amplitude [rad]; rms-equivalent m/√2
X_T04 = 0.5585017864288694


def rpm(f, eta, x, fh=FH):
    u2 = (f / fh) ** 2
    return 1 - eta * 4 * x / ((1 + x) ** 2 + u2), 1 + eta * 4 * x / ((1 - x) ** 2 + u2)


def jitter(rm, rp, s):
    chi = np.exp(-2 * s * s)
    return 0.5 * (rm + rp) + 0.5 * (rm - rp) * chi, 0.5 * (rm + rp) - 0.5 * (rm - rp) * chi


def db(v):
    return 10 * np.log10(v)


print("1) preset point eta = 0.79, x = X_T04: model vs Fig. 4 for sigma = 0, m/sqrt2, 0.078, 0.1 rad")
for s in (0.0, M / np.sqrt(2), 0.078, 0.1):
    vm, vp = jitter(*rpm(F, 0.79, X_T04), s)
    print(f"  sigma={s:.4f}: sqz " + " ".join(f"{a:+.2f}" for a in db(vm)) + "   anti " + " ".join(f"{a:+.2f}" for a in db(vp)))
print("  paper       : sqz " + " ".join(f"{a:+.2f}" for a in SQZ) + "   anti " + " ".join(f"{a:+.2f}" for a in ANTI))


def res30(p):
    vm, vp = jitter(*rpm(F[:1], 0.79, p[0]), p[1])
    return [db(vm)[0] - SQZ[0], db(vp)[0] - ANTI[0]]


r = least_squares(res30, [0.56, 0.08], bounds=([0.1, 0], [0.95, 0.5]))
print("2) 30 MHz pair alone, eta = 0.79: solve (x, sigma)")
print(f"  x = {r.x[0]:.4f}, sigma = {r.x[1]:.4f} rad, non-dither part sqrt(s^2 - m^2/2) = "
      f"{np.sqrt(max(r.x[1]**2 - M**2 / 2, 0)):.4f} rad")


def res_all(p, fh=FH, skip=()):
    vm, vp = jitter(*rpm(F, p[0], p[1], fh), p[2])
    out = np.r_[db(vm) - SQZ, db(vp) - ANTI]
    return out[[i for i in range(8) if i % 4 not in skip]]


print("3) least squares over the Fig. 4 numbers (fh fixed 875 MHz): free (eta, x, sigma)")
for skip, lab in (((), "all four"), ((2,), "without 0.69 GHz")):
    r = least_squares(res_all, [0.79, 0.56, 0.08], bounds=([0.3, 0.1, 0], [1, 0.95, 0.5]), kwargs={"skip": skip})
    print(f"  {lab}: eta={r.x[0]:.3f} x={r.x[1]:.3f} sigma={r.x[2]:.4f} rad, rms residual "
          f"{np.sqrt(np.mean(r.fun**2)):.3f} dB, residuals sqz(4)/anti(4) = {np.round(res_all(r.x), 2)}")

print("4) 0.69 GHz: sigma needed to reach -3.0 dB at eta=0.79, x=X_T04 (phase noise only)")
rm, rp = rpm(np.array([6.9e8]), 0.79, X_T04)
s = brentq(lambda s: db(jitter(rm, rp, s)[0])[0] + 3.0, 0, 0.5)
print(f"  sigma = {s:.3f} rad -> anti {db(jitter(rm, rp, s)[1])[0]:.2f} dB (paper 5.2); 30 MHz would then be "
      f"{db(jitter(*rpm(F[:1], 0.79, X_T04), s)[0])[0]:.2f} dB (paper -4.7)")
print("  dither-only (m/sqrt2) degradation of the squeezed level [dB] vs f:")
for f in (3e7, 2.5e8, 3e8, 5e8, 6.9e8, 1.21e9):
    rm, rp = rpm(np.array([f]), 0.79, X_T04)
    d = db(jitter(rm, rp, M / np.sqrt(2))[0])[0] - db(rm)[0]
    print(f"    {f/1e9:5.2f} GHz: +{d:.3f} dB")
print("  phase-insensitive excess N (both quadratures + N) that reproduces 0.69 GHz [shot-noise units]:")
vm, vp = jitter(rm := rpm(np.array([6.9e8]), 0.79, X_T04)[0], rpm(np.array([6.9e8]), 0.79, X_T04)[1], M / np.sqrt(2))
print(f"    sqz needs +{10**(-0.3) - vm[0]:.3f}, anti needs +{10**0.52 - vp[0]:.3f}")

print("5) linewidth free too (eta, x, sigma, fh)")
r = least_squares(lambda p: res_all(p[:3], p[3] * 1e9), [0.79, 0.56, 0.08, 0.875],
                  bounds=([0.3, 0.1, 0, 0.3], [1, 0.95, 0.5, 2]))
print(f"  eta={r.x[0]:.3f} x={r.x[1]:.3f} sigma={r.x[2]:.4f} fh={r.x[3]:.3f} GHz (FWHM {2*r.x[3]:.2f}), "
      f"residuals {np.round(r.fun, 2)}")
