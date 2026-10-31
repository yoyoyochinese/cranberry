"""VT-18 (WP-B6 unit part; WP-I4 adds the preset run): Ref B definition discrimination and back-solve
(prompt §7-9 per D48; 00_summary §4 VT-18; T02 V-B4, V-B5, V-B6).

* V-B4: x = 0.999, θ̃ = 1.5°: R′+ = +21.4 dB with Ω_B = 2πf/γ_B; the paper's literal Ω = f/γ (2π missing,
  emulated by evaluating at f/2π) gives ≈ +37 dB — the discriminator of the Ω definition.
* V-B5: L_ic = 0.0038 (measured at 100 mW), θ̃ = 1.5°: R′∓ = −9.44 / +14.95 dB (±0.02); with the
  circuit noise −21.7 dB (pure loss η_el = 1/(1 + C)) −9.21 / +14.92 dB, within 0.25 dB of the observed
  −9.01 / +15.12 dB. Back-solve from the reported generated R∓ = −10.12 / +15.15 dB: x = 0.753 ± 0.002,
  total loss 0.0708 (paper 0.0709).
* V-B6: lossless, Ω = 0 (toggle frequency_dependence off): min_x R′− = 10 log10(sin 2θ̃) = −8.67 / −12.81 dB
  for 3.9° / 1.5° (paper −8.66 / −12.81, ±0.02).
"""

import numpy as np
import pytest
from scipy.optimize import brentq, minimize_scalar

from cranberry.blocks.source.opo_cw import closed_form
from cranberry.blocks.source.opo_testing import REF_B, REF_B_ETA_D, db, measured_db, opo

TH15 = np.radians(1.5)


def test_v_b4_two_pi_in_omega():
    good = measured_db(opo(REF_B, pump_mode="x", x=0.999), REF_B_ETA_D, TH15)[1][0]
    bad = measured_db(opo(REF_B, pump_mode="x", x=0.999, f_s_hz=[1e6 / (2 * np.pi)]), REF_B_ETA_D, TH15)[1][0]
    assert good == pytest.approx(21.4, abs=0.05)
    assert bad == pytest.approx(37.4, abs=0.2)


def test_v_b5_100mw_point_and_circuit_noise():
    rm, rp = measured_db(opo(REF_B, loss_ic=0.0038), REF_B_ETA_D, TH15)
    assert (float(rm[0]), float(rp[0])) == pytest.approx((-9.44, 14.95), abs=0.02)
    c = 10 ** (-21.7 / 10)  # circuit noise / shot noise
    with_c = [db((10 ** (v / 10) + c) / (1 + c)) for v in (float(rm[0]), float(rp[0]))]
    assert with_c == pytest.approx([-9.21, 14.92], abs=0.02)
    assert with_c == pytest.approx([-9.01, 15.12], abs=0.25)  # observed [Ref B §4]


def test_v_b5_back_solve():
    b = opo(REF_B, loss_ic=0.0038)
    w = 2 * np.pi * 1e6 / b.kappa(b.x)  # κ independent of x for a constant L_ic
    rm, rp = 10 ** (-10.12 / 10), 10 ** (15.15 / 10)

    def ratio(x):  # (R+ − 1)/(1 − R−) does not depend on η_tot
        m, p = closed_form(x, 1.0, w)
        return (p - 1) / (1 - m)

    x = brentq(lambda x: ratio(x) - (rp - 1) / (1 - rm), 0.3, 0.95, xtol=1e-12)
    eta = (rp - 1) / (closed_form(x, 1.0, w)[1] - 1)
    assert x == pytest.approx(0.753, abs=0.002)
    assert 1 - eta == pytest.approx(0.0708, abs=0.0005)


@pytest.mark.parametrize("th_deg, target", [(3.9, -8.67), (1.5, -12.81)])
def test_v_b6_ideal_limit_sin_2theta(th_deg, target):
    th = np.radians(th_deg)

    def rmin(x):
        b = opo(REF_B, loss_ic=0.0, pump_mode="x", x=x, toggles={"frequency_dependence": False})
        return float(measured_db(b, 1.0, th)[0][0])

    res = minimize_scalar(rmin, bounds=(0.0, 0.9999), method="bounded", options={"xatol": 1e-9})
    assert res.fun == pytest.approx(db(np.sin(2 * th)), abs=1e-4)
    assert res.fun == pytest.approx(target, abs=0.02)
