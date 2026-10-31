"""VT-07 (prompt §7-5, T06 V1): fundamental soliton invariance over 5 z0.

N = 1 sech, α = β3 = f_R = 0, shock off, 200 fixed RK4IP steps, 2^12 points, T_win = 80 T0.
Tolerance: relative L2 error vs √P0 sech(T/T0) e^{iz/(2L_D)} (phase included) ≤ 1e−5 (prototype 2.9e−6).
"""

import numpy as np

from cranberry.blocks.fib.gnlse import GnlseModel, pulse_field, solve_pulsed, time_axis

B2, GAMMA, T0 = -21.7e-27, 1.3e-3, 100e-15
LD = T0**2 / abs(B2)
P1 = abs(B2) / (GAMMA * T0**2)
Z0 = np.pi / 2 * LD


def test_fundamental_soliton_invariance():
    t, dt = time_axis(4096, 80 * T0)
    m = GnlseModel(4096, dt, lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1.0,
                   c_raman=0.0, raman_model=None, self_steepening=False)
    a0 = pulse_field("sech", P1, T0, t)
    r = solve_pulsed(m, a0, t, 5 * Z0, n_steps=200, n_save=1, control="fixed", rtol=None)
    exact = a0 * np.exp(1j * 5 * Z0 / (2 * LD))
    err = np.linalg.norm(r.a_out - exact) / np.linalg.norm(exact)
    assert err <= 1e-5, err
