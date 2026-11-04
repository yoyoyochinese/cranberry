"""Phase 5 review of WP6-X (D94-2, D94 結果): LO delay mismatch and the filter-cavity (FC) delay compensation of
bench/refA_fiber, independent of the simulator (closed-form white-FM integral, conv §7).

theta(t) = phi_LO − theta_sq carries the laser phase as H_FC·phi_L(t − tau_lo) − phi_L(t − tau_sq), so with
Δ = tau_lo − tau_sq and the white FM S_phi = h0/f²:

    sigma²(Δ) = ∫_0^∞ |H_FC(f) e^{−i2πfΔ} − 1|² h0/f² df,   H_FC = 1/(1 + i f/f_c), f_c = Δν_FC/2

(the additive form 2π²h0|Δ| + π h0/Δν_FC used in the refA header is exact only at Δ = 0 or without the FC).
The FC delays the LO's laser phase by tau_FC = 1/(2π f_c) at low f, so a squeezed path LONGER than the LO by about
tau_FC compensates it. Geometry of bench/refA_fiber (header): free_space Δ = (7.64 − 1.60 − n_g L)/c,
fiber_matched Δ ≈ 20.1 ns for every L.
Prints: sigma_laser per L for both routes (exact and additive), the optimum L of the free_space route, the prediction
of the smoke table's free_space column from its L = 1 m value (other contributions held fixed), and the DX1 numbers.
Usage: uv run python reviews/scripts/p6_wp6x_fc_delay.py
"""
import math

import numpy as np
from scipy import integrate, optimize

C = 299792458.0
H0_KOH = 200.0 / math.pi          # refA koheras stand-in [Hz²/Hz]
H0_DX1 = 2.55e4                   # D35 DX1 white FM
DNU_FC = 3.3e6                    # refA fc.dnu_fc_hz
N_G = 1.468


def sigma2(delta_s, h0=H0_KOH, dnu_fc=DNU_FC):
    fc = dnu_fc / 2.0

    def g(f):
        h = 1.0 / (1.0 + 1j * f / fc)
        return abs(h * np.exp(-2j * math.pi * f * delta_s) - 1.0) ** 2 * h0 / f**2

    tot = 0.0
    edges = [1e-3, 1e3, 1e5, 1e6, 1e7, 1e8, 1e10]
    for a, b in zip(edges[:-1], edges[1:]):
        tot += integrate.quad(g, a, b, limit=2000)[0]
    return tot


def additive(delta_s, h0=H0_KOH, dnu_fc=DNU_FC):
    return 2 * math.pi**2 * h0 * abs(delta_s) + math.pi * h0 / dnu_fc


def d_free(L):
    return (7.64 - 1.60 - N_G * L) / C


assert math.isclose(sigma2(0.0), math.pi * H0_KOH / DNU_FC, rel_tol=1e-3)            # FC term alone
assert math.isclose(sigma2(50e-9, dnu_fc=1e15), 2 * math.pi**2 * H0_KOH * 50e-9, rel_tol=2e-3)   # delay term alone

tau_fc = 1.0 / (math.pi * DNU_FC)
print(f"tau_FC = 1/(π Δν_FC) = {tau_fc * 1e9:.1f} ns")
print(f"{'L [m]':>6s} {'route':>14s} {'Δ [ns]':>9s} {'exact [mrad]':>13s} {'additive [mrad]':>16s}")
for route in ("fiber_matched", "free_space"):
    for L in (1.0, 3.0, 10.0, 20.0, 30.0, 100.0):
        d = 20.1e-9 if route == "fiber_matched" else d_free(L)
        print(f"{L:6.0f} {route:>14s} {d * 1e9:9.1f} {1e3 * math.sqrt(sigma2(d)):13.2f} {1e3 * math.sqrt(additive(d)):16.2f}")

res = optimize.minimize_scalar(lambda L: sigma2(d_free(L)), bounds=(1.0, 100.0), method="bounded")
d_opt = d_free(res.x)
print(f"\nfree_space optimum: L = {res.x:.1f} m (Δ = {d_opt * 1e9:.1f} ns, −Δ/tau_FC = {-d_opt / tau_fc:.2f}), "
      f"sigma_laser = {1e3 * math.sqrt(res.fun):.2f} mrad")

smoke = {1.0: 9.03, 10.0: 5.97, 100.0: 20.57}       # p6_bench_smoke.py free_space sigma_total [mrad] (this review's run)
c0 = (smoke[1.0] * 1e-3) ** 2 - sigma2(d_free(1.0))
print(f"other contributions from L = 1 m: {1e3 * math.sqrt(max(c0, 0.0)):.2f} mrad (c0 = {c0:.3e} rad²)")
for L, s in smoke.items():
    pred = math.sqrt(c0 + sigma2(d_free(L)))
    print(f"  L = {L:5.0f} m: predicted {1e3 * pred:6.2f} mrad, smoke {s:6.2f} mrad")

print("\nDX1 (h0 = 2.55e4):")
print(f"  ratio h0 = {H0_DX1 / H0_KOH:.1f}")
print(f"  FC term π h0/Δν_FC = {math.pi * H0_DX1 / DNU_FC:.4f} rad² ({1e3 * math.sqrt(math.pi * H0_DX1 / DNU_FC):.0f} mrad)")
for d, what in ((20.1e-9, "fiber_matched 20.1 ns"), (d_free(100.0), "free_space 100 m")):
    print(f"  {what}: delay term 2π²h0|Δ| = {2 * math.pi**2 * H0_DX1 * abs(d):.4f} rad²; exact total "
          f"{math.sqrt(sigma2(d, H0_DX1)):.3f} rad, additive {math.sqrt(additive(d, H0_DX1)):.3f} rad")
dmax = optimize.brentq(lambda d: sigma2(d, H0_DX1) - 2 * sigma2(0.0, H0_DX1), 0.0, 1e-6)
print(f"  |Δ| (LO later) that doubles the FC term: {dmax * 1e9:.1f} ns")
# sigma² ∝ h0, so the optimum Δ is the Koheras one (d_opt above)
print(f"  optimum Δ = {d_opt * 1e9:.1f} ns (squeezed longer): {math.sqrt(sigma2(d_opt, H0_DX1)):.3f} rad; "
      f"Δ = 0: {math.sqrt(sigma2(0.0, H0_DX1)):.3f} rad")
for off in (-1e-9, 1e-9, -1e-8, 1e-8):
    print(f"    Δ_opt {off * 1e9:+5.0f} ns: {math.sqrt(sigma2(d_opt + off, H0_DX1)):.4f} rad")
print("  without a filter cavity (Δν_FC → ∞), delay term only: Δ = 1 ns → "
      f"{1e3 * math.sqrt(sigma2(1e-9, H0_DX1, 1e15)):.1f} mrad, 0.2 ns → {1e3 * math.sqrt(sigma2(2e-10, H0_DX1, 1e15)):.1f} mrad")
