"""VT-04 (prompt §7-2, WP-B1 part): Ref A Fig. 5, four curves × 8 pump powers and the optima.

Sources: S §4 VT-04; T01 V-T01-2, §7.2 table; T04 V2; D48 (§7-2 reading). Single-mode reference
model, α = 9.06 /W, lumped loss, small-angle average (the §7-1 form the figure uses).
Tolerances: ±0.01 dB to the T01 §7.2 table (values given to 0.01 dB); optima P_opt = 612 mW at
−10.152 dB (θ 9 mrad, L 8 %) and 503 mW at −8.397 dB (14 mrad, 12 %) to ±1 mW / ±1e−3 dB.
Excluded (D48): the right end (1200 mW) of the "previous work" dashed curve — the figure shows
≈ −4.8 dB against the computed −6.07 dB, so it is not a shape target. The ±0.35 dB check against the
figure's data points needs digitised points that are not in the notes (left to WP-I1).
"""

import numpy as np
import pytest

from cranberry.blocks.source.single_mode import make_single_mode_reference, optimal_pump_w

P_MW = [40, 100, 200, 400, 640, 800, 1000, 1200]
CURVES = {  # (θ [rad], L): [(R′− dB, R′+ dB)] per P, T01 §7.2
    "this_work": (0.009, 0.08, [(-4.48, 4.98), (-6.62, 7.96), (-8.43, 11.36), (-9.84, 16.18), (-10.15, 20.56),
                                (-9.98, 23.02), (-9.47, 25.78), (-8.67, 28.28)]),
    "previous_work": (0.014, 0.12, [(-4.15, 4.85), (-5.98, 7.80), (-7.40, 11.18), (-8.32, 15.99), (-8.29, 20.36),
                                    (-7.91, 22.83), (-7.14, 25.59), None]),  # right end excluded (D48)
    "blue": (0.0, 0.08, [(-4.49, 4.98), (-6.63, 7.96), (-8.47, 11.36), (-9.98, 16.18), (-10.58, 20.56),
                         (-10.75, 23.02), (-10.85, 25.78), (-10.90, 28.28)]),
    "green": (0.014, 0.08, [(-4.48, 4.98), (-6.61, 7.96), (-8.39, 11.35), (-9.64, 16.18), (-9.60, 20.56),
                            (-9.08, 23.02), (-8.06, 25.78), (-6.71, 28.28)]),
}


def curve_db(theta, loss, p_w):
    b = make_single_mode_reference(dict(loss=loss, alpha_per_w=9.06, pump_w=p_w, theta_rms_rad=theta,
                                        phi_opa0_rad=0.0, average="small_angle"))
    return b.r_prime_db()


@pytest.mark.parametrize("name", list(CURVES))
def test_fig5_curves(name):
    theta, loss, pts = CURVES[name]
    for p_mw, target in zip(P_MW, pts, strict=True):
        if target is None:
            continue
        rm, rp = curve_db(theta, loss, p_mw * 1e-3)
        assert rm == pytest.approx(target[0], abs=0.01), (name, p_mw)
        assert rp == pytest.approx(target[1], abs=0.01), (name, p_mw)


@pytest.mark.parametrize("theta, loss, p_opt_mw, rmin_db", [(0.009, 0.08, 612, -10.152), (0.014, 0.12, 503, -8.397)])
def test_optima(theta, loss, p_opt_mw, rmin_db):
    p = optimal_pump_w(9.06, theta)
    assert p * 1e3 == pytest.approx(p_opt_mw, abs=1.0)
    assert curve_db(theta, loss, p)[0] == pytest.approx(rmin_db, abs=1e-3)
    assert 10 * np.log10(loss + (1 - loss) * np.sin(2 * theta)) == pytest.approx(rmin_db, abs=1e-3)
