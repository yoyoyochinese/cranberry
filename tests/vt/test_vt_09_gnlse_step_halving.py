"""VT-09 (prompt §7-5, T06 V3, PLAN §6.4): global error ~1/16 per step halving (RK4IP is 4th order).

N = 2 soliton to z0, fixed steps without extrapolation, 200/400/800 steps vs a 6400-step reference.
Tolerance: successive error ratios in [13, 19] in the asymptotic range (prototype 14.2, 15.2); the
coarser 100 → 200 ratio is reported but not judged (T06 §3.2: < 16 before convergence is expected).
"""

import numpy as np

from cranberry.blocks.fib.gnlse import GnlseModel, pulse_field, solve_pulsed, time_axis

B2, GAMMA, T0 = -21.7e-27, 1.3e-3, 100e-15
LD = T0**2 / abs(B2)
P1 = abs(B2) / (GAMMA * T0**2)
Z0 = np.pi / 2 * LD


def test_step_halving_error_ratio():
    t, dt = time_axis(4096, 80 * T0)
    m = GnlseModel(4096, dt, lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1.0,
                   c_raman=0.0, raman_model=None, self_steepening=False)
    a0 = pulse_field("sech", 4 * P1, T0, t)

    def run(n):
        return solve_pulsed(m, a0, t, Z0, n_steps=n, n_save=1, control="fixed", rtol=None).a_out

    ref = run(6400)
    errs = [np.linalg.norm(run(n) - ref) / np.linalg.norm(ref) for n in (200, 400, 800)]
    ratios = [errs[0] / errs[1], errs[1] / errs[2]]
    assert all(13.0 <= r <= 19.0 for r in ratios), (errs, ratios)
