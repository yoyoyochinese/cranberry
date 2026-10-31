"""M6 review: does the Liu-type lock remove the in-band laser decorrelation for km-scale port fibres?

Closed-loop path model (conv §7 emission time, conv §9; the general form of reviews/phase5_M4 §2.1 with the actuator at
the upstream end x_a = 0 and T_px = 0):  per port j,
    theta_j(w) = e^{-iw(T_a+tau_j)} - T(w) * [e^{-iw(T_b+tau_j)} - e^{-iw T_ref}] * e^{-iw tau_j}
(the stretcher correction a_j = -T*P_j reaches the GBS tau_j after it is applied), T = L/(1+L),
L(s) = (2 pi f_u / s) e^{-s tau_L}  (integrator + loop delay; the 915B-like resonance is left out here).
Theta_k = theta_k - theta_0.  White FM h0:  S_phi = h0/f^2.
Ratio R(f_u) = int |Theta_k|^2 S_phi df / int |e^{-iw tau_k} - e^{-iw tau_0}|^2 S_phi df  (locked / unlocked).
D77-1 predicts R ~ 1 - 2 dtau f_u (in-band residual = |1/(1+L)|^2 x out-of-band). That holds only for 2 pi f_u tau_k << 1:
for tau_k of the F_* presets (5-16 us) the ideal-lock term |H_ideal| ~ w |Y| |e^{-iw tau_k}-e^{-iw tau_0}|, |Y| ~ tau_k,
is O(1) inside the loop band and the lock ADDS decorrelation.  Self-check at the end (asserts).
Usage: uv run python reviews/scripts/m6_liu_lock_long_fiber.py
"""
import numpy as np

C = 299792458.0
NG = 1.468


def ratio(tau_k, tau_0, f_u, T_a=8.4e-9, T_b=34e-9, T_ref=9.8e-9, tau_L=1.5e-6, f_max=5e5):
    f = np.logspace(0, np.log10(f_max), 400001)
    w = 2 * np.pi * f
    L = (2 * np.pi * f_u / (1j * w)) * np.exp(-1j * w * tau_L)
    T = L / (1 + L) if f_u > 0 else 0.0
    th = lambda tau: np.exp(-1j * w * (T_a + tau)) - T * (np.exp(-1j * w * (T_b + tau)) - np.exp(-1j * w * T_ref)) * np.exp(-1j * w * tau)  # noqa: E731
    s = 1.0 / f**2
    locked = np.trapezoid(np.abs(th(tau_k) - th(tau_0)) ** 2 * s, f)
    free = np.trapezoid(np.abs(np.exp(-1j * w * tau_k) - np.exp(-1j * w * tau_0)) ** 2 * s, f)
    return locked / free


cases = {   # name: (L_k, L_0) in m; F_* nominal port 31
    "D77-1 63 m (100 MHz, N_b=1)": (2.0 + 31 * C / (1e8 * NG), 2.0),
    "F_wg6 / F_bor6 1.06 km": (2.0 + 31 * C / (6e6 * NG), 2.0),
    "F_liu20 1.27 km": (2.0 + 31 * 4 * C / (2e7 * NG), 2.0),
    "F_wg100 3.17 km": (2.0 + 31 * 50 * C / (1e8 * NG), 2.0),
}
print(f"{'case':30s} {'tau_k[us]':>9s} " + " ".join(f"f_u={fu:>6.0f}" for fu in (300, 1e3, 3e3, 1e4, 3e4)) + "   (1-2 dtau f_u @1e4)")
res = {}
for name, (lk, l0) in cases.items():
    tk, t0 = NG * lk / C, NG * l0 / C
    r = [ratio(tk, t0, fu) for fu in (300, 1e3, 3e3, 1e4, 3e4)]
    res[name] = r
    print(f"{name:30s} {tk*1e6:9.3f} " + " ".join(f"{x:10.4f}" for x in r) + f"   {1-2*(tk-t0)*1e4:.4f}")

# self-checks: D77-1 predicts R = 1 - 2 dtau f_u < 1; the closed loop gives R > 1 in every case (Bode waterbed of the
# delayed loop on a spectrum that is flat in f up to ~1/(2 dtau), plus the ideal-lock term H_ideal ~ w Y for km fibres)
for name, r in res.items():
    assert all(x > 1.0 for x in r), name
tk, t0 = NG * cases["D77-1 63 m (100 MHz, N_b=1)"][0] / C, NG * 2.0 / C
assert res["D77-1 63 m (100 MHz, N_b=1)"][3] > 1 - 2 * (tk - t0) * 1e4 + 0.03
assert res["F_liu20 1.27 km"][3] > 1.15 and res["F_liu20 1.27 km"][4] > 1.5
assert res["F_wg100 3.17 km"][3] > 1.35
# small-(w Y) regime: R ~ 1 - pi f_u dtau (in-band |S|^2) + 4 pi f_u Y (cross term -2 Re[S open conj(T H_ideal)], ratio
# 2 w_u Y at every f > f_u up to ~1/Y), Y ~ dtau for a long port: R ~ 1 + 3 pi f_u dtau.  Check on the 63 m case.
for i, fu in enumerate((300, 1e3, 3e3, 1e4)):
    approx = 1 + 3 * np.pi * fu * (tk - t0)
    print(f"63 m, f_u = {fu:6.0f}: R = {res['D77-1 63 m (100 MHz, N_b=1)'][i]:.4f}, 1 + 3 pi f_u dtau = {approx:.4f}")
    assert abs(res["D77-1 63 m (100 MHz, N_b=1)"][i] - approx) < 0.15 * (approx - 1) + 1e-4
print("self-checks passed")
