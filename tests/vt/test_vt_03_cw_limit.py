"""VT-03 (prompt §7-1 CW limit / §5.1 LASER_MODE = CW, T06 V5, PLAN §2.2): GNLSE CW shortcut — GNLSE part.

Owner of this file: WP-B9 (the pulsed-source CW-limit parts belong to WP-B11 and may be added here).

* CW shortcut (T06 eq. (5)) == RK4IP of the same CW field on a time grid with every term on
  (loss, β2, β3, β4, Kerr, Raman, self-steepening, external XPM): relative difference ≤ 1e−12.
* The CW nonlinear phase is (γP + bγP̄) L_eff and does not change when Raman is switched off (D49).
* Sideband rotation θ_disp(Ω) = L Σ_even β_k Ω^k/k! equals the even part of the GNLSE linear operator.
* All toggles off ⇒ identity.
"""

import numpy as np
import pytest

from cranberry.blocks.fib.gnlse import TOGGLES, db_per_km_to_alpha_per_m, l_eff, model_from_params
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block

FIBER = dict(L_m=1000.0, alpha_per_m=db_per_km_to_alpha_per_m(0.2), beta2_s2_per_m=-21.7e-27,
             beta3_s3_per_m=0.126e-39, beta_higher=[-1e-55], gamma_per_w_m=1.3e-3, lambda0_m=1550e-9,
             raman_model="blow_wood", f_R=0.18, xpm_b=2.0, xpm_power_w=0.1)
ALL_ON = dict.fromkeys(TOGGLES, True)


def _block(toggles, **extra):
    params = dict(FIBER, toggles=toggles, laser_mode="cw", power_w=1.0, grid={"n_t": 64, "t_window_s": 1e-9},
                  steps={"n_steps": 2000, "n_save": 1, "control": "fixed"}, **extra)
    return make_block(BlockSpec("fiber_gnlse", "fib", params))


def test_cw_shortcut_equals_rk4ip_on_grid():
    cw = _block(ALL_ON).result
    grid = _block(dict(ALL_ON, cw_shortcut=False)).result
    out = grid.a_out
    assert np.max(np.abs(out - out[0])) <= 1e-12 * abs(out[0])  # stays CW
    assert abs(out[0] - cw.a_out) <= 1e-12 * abs(cw.a_out)
    le = l_eff(FIBER["alpha_per_m"], FIBER["L_m"])
    assert cw.phi_nl_rad == pytest.approx((1.3e-3 * 1.0 + 2.0 * 1.3e-3 * 0.1) * le, rel=1e-14)
    assert cw.eta == pytest.approx(10 ** (-0.02), rel=1e-12)


def test_raman_off_keeps_total_gamma():
    on, off = _block(ALL_ON).result, _block(dict(ALL_ON, raman=False)).result
    assert off.phi_nl_rad == pytest.approx(on.phi_nl_rad, rel=1e-15)


def test_sideband_rotation_is_even_part_of_linear_operator():
    cw = _block(ALL_ON).result
    w = np.linspace(-2 * np.pi * 5e9, 2 * np.pi * 5e9, 11)
    tg = dict(ALL_ON)
    m = model_from_params(FIBER, tg, 64, 1e-11)
    even = lambda x: FIBER["L_m"] * sum(b * x ** (k + 2) / [2, 6, 24][k] for k, b in enumerate(  # noqa: E731
        (FIBER["beta2_s2_per_m"], FIBER["beta3_s3_per_m"], FIBER["beta_higher"][0])) if k % 2 == 0)
    np.testing.assert_allclose(cw.sideband_rotation_rad(w), even(w), rtol=1e-14, atol=0)
    # θ_disp(Ω) = L·Im[D(Ω) + D(−Ω)]/2 on the model's own grid (Nyquist bin 32 has no −Ω partner)
    d, k = m.D, np.r_[0:32, 33:64]
    np.testing.assert_allclose(cw.sideband_rotation_rad(m.omega[k]), FIBER["L_m"] * (d + d[(-np.arange(64)) % 64]).imag[k] / 2,
                               rtol=1e-12, atol=1e-18)
    # Ref C order of magnitude (T06 §8.7): 1 km, 1 GHz → |θ| = 4.3e-4 rad (negative for β2 < 0)
    assert cw.sideband_rotation_rad(2 * np.pi * 1e9) == pytest.approx(-4.28e-4, rel=0.01)


def test_all_toggles_off_identity():
    r = _block(dict.fromkeys(TOGGLES, False) | {"cw_shortcut": True}).result
    assert r.a_out == r.a_in and r.phi_nl_rad == 0.0 and r.eta == 1.0
    assert np.all(r.sideband_rotation_rad(np.array([1e9, 1e10])) == 0.0)
