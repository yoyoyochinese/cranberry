"""VT-19 (Ref C preset part, WP-I3): the dither contribution in the closed Ref C loop (prompt §7-10 first
half; 00_summary §4 VT-19 ← T04 V5 / §9.1, T14 V8, T16 V15). Block part: test_vt_19_dither_contribution.py.

``refC_nominal`` at the T04 source point (x with 10.0 dB jitter-free anti-squeezing at 30 MHz, η_tot = 0.79;
e^{2r} = 12.4 in the single-mode form), all noise sources off, loop closed:

(a) analytic dither (preset default): the time-loop residual is 0, the metric's dither part is
    m/√2 = 0.0555 rad (m = π·0.1 V/4 V = 0.0785 rad, ±4.5°) with χ = J0(2m) = 0.99384, and R′−(30 MHz)
    = −5.18 dB; without the dither entry in the metric (σ = 0) −5.63 dB (T04 §9.1, ±0.01 dB);
(b) time-series dither (presets/axes/refC_dither/timeseries.yaml: dither on the EOM drive, RF-rate QNL and
    demodulator): the loop still locks and the residual itself is the dither: rms 0.0555 rad ±2 % (the
    demodulator's finite LPF lets a little 1 MHz ripple through the loop, +1.4 %), R′−(30 MHz) equal to
    (a) within ±0.05 dB (T04 V5's time-domain tolerance);
(c) 0.1 rad total with Ref C noise sources: the preset's own laser-phase-noise path (white FM h0 with a
    2 m LO/squeezed path mismatch, T05 §3.5: σ² ≈ 2π² h0 τ, mostly above the loop band → BandLedger)
    sized to 0.1 rad total gives R′−(30 MHz) = −4.32 dB (T04 §9.1) ±0.05 dB (σ within ±5 %).
"""

import importlib.util
import math
from pathlib import Path

import pytest
from scipy.special import j0

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parents[1] / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

M = math.pi * 0.1 / 4.0
RUN = {"duration_s": 0.06}


@pytest.fixture(scope="module")
def analytic():
    return run(rs.load(params={"x": rs.X_T04}, run=RUN, toggles=rs.NOISE_OFF))


def test_vt19a_analytic_dither_levels(analytic):
    nl = analytic.metrics.noise_levels["i.port0"]
    assert M == pytest.approx(0.0785, abs=5e-5) and M / math.sqrt(2) == pytest.approx(0.0555, abs=5e-5)
    assert rs.at(nl, "std_loop_rad", 3e7) < 1e-5                       # nothing else moves Δ
    assert nl["sigma_dither_rad"] == pytest.approx(M / math.sqrt(2), rel=1e-12)
    assert nl["chi_dither"] == pytest.approx(j0(2 * M), rel=1e-12)
    assert nl["chi_dither"] == pytest.approx(0.99384, abs=1e-5)
    assert rs.at(nl, "v_plus_db", 3e7) == pytest.approx(10.0, abs=1e-4)          # the T04 source point
    assert rs.at(nl, "r_minus_db", 3e7) == pytest.approx(-5.18, abs=0.01)
    assert rs.at(nl, "v_minus_db", 3e7) == pytest.approx(-5.63, abs=0.01)   # no jitter
    # the same run without the dither entry in the metric: σ = 0 → −5.63 dB
    no_d = run(rs.load(params={"x": rs.X_T04}, run=RUN, toggles=rs.NOISE_OFF,
                       measure={"metrics": {"i": {"dither_m_rad": []}}}))
    assert rs.at(no_d.metrics.noise_levels["i.port0"], "r_minus_db", 3e7) == pytest.approx(-5.63, abs=0.01)


def test_vt19b_time_series_dither_equals_analytic(analytic):
    ts = run(rs.load(params={"x": rs.X_T04}, run=RUN, toggles=rs.NOISE_OFF, overlays=(rs.TIMESERIES,)))
    nl, na = ts.metrics.noise_levels["i.port0"], analytic.metrics.noise_levels["i.port0"]
    assert nl["chi_dither"] == 1.0                                     # not counted twice
    assert rs.at(nl, "sigma_loop_rad", 3e7) == pytest.approx(M / math.sqrt(2), rel=0.02)
    assert abs(rs.at(nl, "mean_loop_rad", 3e7)) < 1e-3                 # locked at Δ = 0
    for f in (3e7, 2.5e8, 6.9e8, 1.21e9):
        assert rs.at(nl, "r_minus_db", f) == pytest.approx(rs.at(na, "r_minus_db", f), abs=0.05)


def test_vt19c_total_0p1_rad_gives_minus_4p32_db():
    tau = 1.468 * 2.0 / 299792458.0                                    # 2 m LO-path excess (n_g of the SMF)
    h0 = (0.1**2 - M**2 / 2) / (2 * math.pi**2 * tau)                  # 2π² h0 τ = 0.1² − m²/2
    tog = {k: v for k, v in rs.NOISE_OFF.items() if k != "laser.fm_white"}
    res = run(rs.load(params={"x": rs.X_T04, "dL_m": 2.0, "h0_hz2_per_hz": h0}, run=RUN, toggles=tog))
    nl = res.metrics.noise_levels["i.port0"]
    sigma = rs.at(nl, "sigma_total_rad", 3e7)
    assert sigma == pytest.approx(0.1, rel=0.05)
    assert rs.at(nl, "r_minus_db", 3e7) == pytest.approx(-4.32, abs=0.05)
