"""M6 review: OPO detuning -> theta_sq coefficients of F_liu20 / F_bor6 vs the x_ss sweep axis.

The presets hard-code c_len_sq_per_kappa (-1.154 Liu @ x_ss=4.129, -1.063 Borealis @ 1.3) and
c_pt_per_kappa while x_ss is a sweep axis. Recompute the WP-B11 tables (full b, pump-cavity phase included)
at the axis points: c_len = cs + 2 cp (cavity-length change, Delta_p = 2 Delta_s), per kappa_c.
"""
from cranberry.blocks.source.pulsed_testing import BOREALIS, LIU, pulsed_opo

for name, base, grid in (("Liu (F_liu20)", LIU, [1.0, 2.5645, 4.129]), ("Borealis (F_bor6)", BOREALIS, [1.3, 1.85, 2.4])):
    for x in grid:
        b = pulsed_opo(base, x_ss=x, discretization="van_loan")
        kc = b.kappa_c
        cs, cp, cb = b.angle_coefficients(x)
        c_len = (cs + 2 * cp) * kc
        c_pt = 0.5 * 2 * cb * kc   # transmitted pump, theta-units, length change (Delta_p = 2 Delta_s)
        print(f"{name:18s} x_ss={x:6.3f}  dth/dDs={cs*kc:+.4f}  dth/dDp={cp*kc:+.4f}  "
              f"length c_len={c_len:+.4f}  pump(theta,len)={c_pt:+.4f}  per kappa_c")
