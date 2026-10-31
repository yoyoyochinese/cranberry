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


# ------------------------------------------------------------------------------------------------------------------
# WP-I4: the same targets through the refB preset and the runner (prompt §7-9 per D48; refB.yaml header).
#
# * calculated curves (±0.02 dB): the runner's offline stage (opo_cw_spectrum with the L(x) callable, LossRegistry
#   ζξ²η) at each Fig. 3 pump power, mapped by Ref B Eq. (1) at θ̃ = 1.5° / 3.9°, gives the V-B3 table; x_opt and
#   R′−,min from the preset's own opo block;
# * closed loop (lock tuned to θ̃ ≈ 1.5°: refB_nominal; ≈ 3.9°: overlay presets/axes/refB_theta/single_layer_pzt.yaml,
#   the 2 kHz single-layer PZTs): θ̃ of the run at 100 mW within 1.5 ± 0.1° / 3.9 ± 0.3° (seed spread 0.02° at
#   0.5 s / 0.09° at 1 s), the noise levels are Eq. (1) at the run's own θ̃ (0.02 dB) and lie within the Fig. 3
#   read accuracy ±0.3 dB of the curves. θ̃ of the closed loop grows towards low x (≈ 1.8° / 4.5° at 10 mW):
#   lock-1's residual reaches Δ through A−/A+ = (1 − x)/(1 + x) (T14 §3.4), a model prediction the constant-θ̃
#   curves of Fig. 3 do not contain; it stays inside ±0.3 dB because R′ is insensitive to θ̃ at low x.
# ------------------------------------------------------------------------------------------------------------------

import importlib.util  # noqa: E402
from pathlib import Path  # noqa: E402

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parents[1] / "integration" / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)


@pytest.fixture(scope="module")
def solid_sweep():
    """refB_nominal closed loop at the nine Fig. 3 pump powers (0.3 s each, the preset seed)."""
    from cranberry.runner.run import run

    return {p: run(rs.load(params={"pump_w": p * 1e-3}, run={"duration_s": 0.3})) for p in sorted(TABLE)}


@pytest.mark.numba
@pytest.mark.integration
@pytest.mark.parametrize("th", [1.5, 3.9])
def test_refB_preset_calculated_curves(solid_sweep, th):
    for p, res in solid_sweep.items():
        nl = res.metrics.noise_levels["i.port0"]
        got = rs.eq1_db(nl["v_minus_db"][rs.I_1MHZ], nl["v_plus_db"][rs.I_1MHZ], TH[th])
        assert got == pytest.approx(TABLE[p][th], abs=0.02), p


@pytest.mark.numba
@pytest.mark.integration
@pytest.mark.parametrize("th, x_opt, r_min", [(1.5, 0.758, -9.35), (3.9, 0.607, -7.17)])
def test_refB_preset_minimum(solid_sweep, th, x_opt, r_min):
    """x_opt from the preset's opo block (pump axis P = x² P_th) and the run's η."""
    res = solid_sweep[100]
    eta = res.offline.loss_registry.eta("sqz", "bhd")
    params = rs.offline_block_params(res.spec, "opo")
    p_th = rs.params_of(res.spec)["p_th_w"]

    def r_minus(x):
        return float(measured_db(opo(params, pump_w=x * x * p_th), eta, TH[th])[0][rs.I_1MHZ])

    m = minimize_scalar(r_minus, bounds=(0.3, 0.95), method="bounded", options={"xatol": 1e-7})
    assert m.x == pytest.approx(x_opt, abs=0.01)
    assert m.fun == pytest.approx(r_min, abs=0.02)


@pytest.mark.numba
@pytest.mark.integration
def test_refB_closed_loop_solid_curve(solid_sweep, record_property):
    th = {}
    for p, res in solid_sweep.items():
        nl = res.metrics.noise_levels["i.port0"]
        t = float(nl["sigma_total_rad"][rs.I_1MHZ])
        th[p] = float(np.degrees(t))
        rm, rp = float(nl["r_minus_db"][rs.I_1MHZ]), float(nl["r_plus_db"][rs.I_1MHZ])
        assert (rm, rp) == pytest.approx(rs.eq1_db(nl["v_minus_db"][rs.I_1MHZ], nl["v_plus_db"][rs.I_1MHZ], t), abs=0.02)
        assert (rm, rp) == pytest.approx(TABLE[p][1.5], abs=0.3), p     # Fig. 3 solid, read accuracy
    record_property("theta_deg_vs_pump_mw", {k: round(v, 3) for k, v in th.items()})
    assert th[100] == pytest.approx(1.5, abs=0.1)
    assert all(1.3 < v < 1.95 for v in th.values())
    assert th[10] > th[100] > th[150]                                     # A−/A+ coupling of lock 1 (header)


@pytest.mark.slow
@pytest.mark.numba
@pytest.mark.integration
@pytest.mark.parametrize("p_mw", [65, 100, 130])
def test_refB_closed_loop_dotted_curve(p_mw):
    """2 s records (the 2 kHz loops leave slow residual fluctuations: seed spread 0.09° at 1 s)."""
    from cranberry.runner.run import run

    res = run(rs.load(params={"pump_w": p_mw * 1e-3}, run={"duration_s": 2.0}, overlays=(rs.SINGLE_LAYER,)))
    nl = res.metrics.noise_levels["i.port0"]
    t = float(nl["sigma_total_rad"][rs.I_1MHZ])
    rm, rp = float(nl["r_minus_db"][rs.I_1MHZ]), float(nl["r_plus_db"][rs.I_1MHZ])
    assert (rm, rp) == pytest.approx(rs.eq1_db(nl["v_minus_db"][rs.I_1MHZ], nl["v_plus_db"][rs.I_1MHZ], t), abs=0.02)
    assert (rm, rp) == pytest.approx(TABLE[p_mw][3.9], abs=0.3)          # Fig. 3 dotted, read accuracy
    if p_mw == 100:
        assert np.degrees(t) == pytest.approx(3.9, abs=0.3)
