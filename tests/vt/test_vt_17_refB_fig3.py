"""VT-17 (WP-B6 unit part; WP-I4 adds the preset run): Ref B Fig. 3 (prompt §7-9 read per D48;
00_summary §4 VT-17; T02 V-B3).

x = √(P/P_th), P_th = 180 mW, Ω_B = 2πf/γ_B (f = 1 MHz), L_ic(x) = 0.00249 + 0.00222x, T = 0.123,
l = 0.5 m, η = 0.998, ξ = 0.988, ζ = 0.99 (downstream efficiency ζξ²η), Ref B Eq. (1) jitter.

* T02 V-B3 table: R′∓ at each P for θ̃ = 1.5° and 3.9°, ±0.02 dB (same equations, numerical error only).
* Minima: θ̃ = 1.5° x_opt = 0.758 (±0.01), R′−,min = −9.35 dB; θ̃ = 3.9° x_opt = 0.607 (±0.01),
  R′−,min = −7.17 dB (±0.02 dB). The 0.1 dB band is x = 0.708–0.806 / 0.558–0.654 (±0.005).
* Fig. 3 read values (±0.3 dB): x → 1 solid ≈ −7.5 / +21.5 dB, dotted ≈ −1.5 dB. "Degrades above
  x ≈ 0.82" is a statement about measured points and is not a target (D48).
"""

import numpy as np
import pytest
from scipy.optimize import brentq, minimize_scalar

from cranberry.blocks.source.opo_testing import REF_B, REF_B_ETA_D, measured_db, opo

TH = {1.5: np.radians(1.5), 3.9: np.radians(3.9)}

# T02 V-B3: P [mW] → {θ̃: (R′−, R′+) [dB]}
TABLE = {
    10: {1.5: (-3.66, 3.89), 3.9: (-3.58, 3.88)},
    20: {1.5: (-5.11, 5.63), 3.9: (-4.93, 5.61)},
    50: {1.5: (-7.72, 9.47), 3.9: (-6.92, 9.46)},
    65: {1.5: (-8.51, 11.15), 3.9: (-7.17, 11.14)},
    90: {1.5: (-9.25, 13.87), 3.9: (-6.69, 13.85)},
    100: {1.5: (-9.35, 14.94), 3.9: (-6.22, 14.93)},
    120: {1.5: (-9.20, 17.08), 3.9: (-4.93, 17.06)},
    130: {1.5: (-8.98, 18.11), 3.9: (-4.19, 18.10)},
    150: {1.5: (-8.35, 19.96), 3.9: (-2.70, 19.94)},
}


def at_x(x, th):
    rm, rp = measured_db(opo(REF_B, pump_mode="x", x=x), REF_B_ETA_D, TH[th])
    return float(rm[0]), float(rp[0])


@pytest.mark.parametrize("p_mw", sorted(TABLE))
@pytest.mark.parametrize("th", [1.5, 3.9])
def test_v_b3_table(p_mw, th):
    rm, rp = measured_db(opo(REF_B, pump_w=p_mw * 1e-3), REF_B_ETA_D, TH[th])
    assert (float(rm[0]), float(rp[0])) == pytest.approx(TABLE[p_mw][th], abs=0.02)


def test_v_b3_near_threshold():
    assert at_x(0.999, 1.5) == pytest.approx((-7.63, 21.43), abs=0.02)
    assert at_x(0.999, 3.9) == pytest.approx((-1.43, 21.41), abs=0.02)
    # Fig. 3 read values (±0.3 dB)
    assert at_x(0.999, 1.5) == pytest.approx((-7.5, 21.5), abs=0.3)
    assert at_x(0.999, 3.9)[0] == pytest.approx(-1.5, abs=0.3)


@pytest.mark.parametrize("th, x_opt, r_min, band", [(1.5, 0.758, -9.35, (0.708, 0.806)),
                                                     (3.9, 0.607, -7.17, (0.558, 0.654))])
def test_v_b3_minimum_and_flat_band(th, x_opt, r_min, band):
    res = minimize_scalar(lambda x: at_x(x, th)[0], bounds=(0.3, 0.95), method="bounded", options={"xatol": 1e-7})
    assert res.x == pytest.approx(x_opt, abs=0.01)
    assert res.fun == pytest.approx(r_min, abs=0.02)
    edges = (brentq(lambda x: at_x(x, th)[0] - res.fun - 0.1, 0.3, res.x),
             brentq(lambda x: at_x(x, th)[0] - res.fun - 0.1, res.x, 0.99))
    assert edges == pytest.approx(band, abs=0.005)
