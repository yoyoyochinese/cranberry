"""M6 review: laser linewidth requirement from metric (ii) component 3 (D31-3, D74-4, D77-1).

White FM S_nu = h0 (Lorentzian FWHM dnu = pi h0), worst port k = K - 1, dtau = (K-1)(N_b + n_guard)/f_rep:
    sigma_Theta^2 = R(f_u) * 2 pi^2 h0 dtau = R * 2 pi dnu dtau   =>   dnu_max = sigma^2 / (2 pi R dtau)
R(f_u) >= 1 is the closed-loop factor of the Liu lock (m6_liu_lock_long_fiber.py; 1.02-1.05 at f_u = 1 kHz, the
nominal point). The table uses R = 1 (the lock never helps white FM, so these are upper bounds on dnu).
Also the slow-drift budget of the orthogonal-PM + TEC-offload configuration (m6_laser_drift_tec_offload.py):
    |Theta| = (|K_T|/(k_T+|K_T|)) 2 pi dnu_drift tau_k   =>   dnu_drift,max = sigma / (0.0445 * 2 pi tau_k).
Usage: uv run python reviews/scripts/m6_linewidth_table.py
"""
import math

NG, C, L0 = 1.468, 299792458.0, 2.0
SIG = (5e-3, 10e-3, 20e-3, 50e-3)
FRAC = 2.05 / (44.0 + 2.05)


def dnu_max(sig, dtau, R=1.0):
    return sig**2 / (2 * math.pi * R * dtau)


def fmt(x):
    return f"{x:.3g}"


print("Table A: K = 32 ports, worst port 31; dnu_max [Hz] (h0_max = dnu/pi) for sigma_Theta from the laser alone")
print("| configuration | f_rep | N_b | L_31 [m] | dtau_31 [us] | " + " | ".join(f"{1e3*s:g} mrad" for s in SIG) +
      " | drift_max @10 mrad (PM orth + TEC) |")
print("|---|---|---|---|---|" + "---|" * len(SIG) + "---|")
for name, fr, nb in (("D77-1 例 (N_b = 1)", 1e8, 1), ("100 MHz, N_b = 10", 1e8, 10), ("F_wg100", 1e8, 50),
                     ("F_liu20", 2e7, 4), ("20 MHz, N_b = 1", 2e7, 1), ("F_wg6 / F_bor6", 6e6, 1)):
    dt = 31 * nb / fr
    lk = L0 + C * dt / NG
    tau_k = NG * lk / C
    print(f"| {name} | {fr/1e6:g} MHz | {nb} | {lk:.0f} | {dt*1e6:.3f} | " +
          " | ".join(fmt(dnu_max(s, dt)) for s in SIG) + f" | {fmt(1e-2 / (FRAC * 2 * math.pi * tau_k))} Hz |")

print("\nTable B: F_liu20 timing (20 MHz, N_b = 4), dnu_max [Hz] at 10 mrad vs number of ports K (worst port K - 1)")
print("| K | dtau [us] | dnu_max (R = 1) | dnu_max (R = 1.02, f_u = 1 kHz) | dnu_max (R = 1.20, f_u = 10 kHz) |")
print("|---|---|---|---|---|")
for K in (2, 4, 8, 16, 32, 64):
    dt = (K - 1) * 4 / 2e7
    print(f"| {K} | {dt*1e6:.2f} | {fmt(dnu_max(1e-2, dt))} | {fmt(dnu_max(1e-2, dt, 1.02))} | {fmt(dnu_max(1e-2, dt, 1.20))} |")

# self-checks against D74-4 / D77-1
assert abs(dnu_max(1e-2, 31 / 1e8) - 51.3) < 0.5                  # D77-1: 51 Hz at 63 m
assert abs(dnu_max(1e-2, 31 * 4 / 2e7) - 2.57) < 0.01             # D74-4: Hz class at 6.2 us
assert abs(math.sqrt(2 * math.pi**2 * 64 * 31 * 4 / 2e7) - 0.0885) < 1e-3   # F_liu20 component 3: 88 mrad at h0 = 64
print("self-checks passed")
