"""VT-08 (prompt §7-5, T06 V2): second-order soliton recurrence at the soliton period z0.

N = 2 sech, α = β3 = f_R = 0, shock off, 400 fixed steps, 2^12 points, T_win = 80 T0.
Tolerance: max_T |I(z0) − I(0)| / max I(0) ≤ 1e−5 (prototype 2.6e−6).
"""

import numpy as np

from cranberry.blocks.fib.gnlse import GnlseModel, pulse_field, solve_pulsed, time_axis

B2, GAMMA, T0 = -21.7e-27, 1.3e-3, 100e-15
LD = T0**2 / abs(B2)
P1 = abs(B2) / (GAMMA * T0**2)
Z0 = np.pi / 2 * LD


def test_second_order_soliton_recurrence():
    t, dt = time_axis(4096, 80 * T0)
    m = GnlseModel(4096, dt, lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1.0,
                   c_raman=0.0, raman_model=None, self_steepening=False)
    a0 = pulse_field("sech", 4 * P1, T0, t)
    r = solve_pulsed(m, a0, t, Z0, n_steps=400, n_save=4, control="fixed", rtol=None)
    i0 = np.abs(a0) ** 2
    err = np.max(np.abs(np.abs(r.a_out) ** 2 - i0)) / i0.max()
    assert err <= 1e-5, err
    # it does change in between (peak compression at z0/2), so the recurrence is not trivial
    assert np.max(np.abs(r.a_tz[2]) ** 2) > 2 * i0.max()
