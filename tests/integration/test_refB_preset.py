"""WP-I4 (M3): the refB preset (Takeno et al., Opt. Express 15, 4321 (2007); D03, D34, D48 §7-9, D51-10, D56-9;
PLAN §1.2 M3, 00_summary §2.1 `refB` / `L-B`).

* every quantity Ref B leaves open is a named sweep axis under ``offline.params`` (D34) and ``refB_nominal``
  resolves all of them;
* offline products: x = √(P/P_th) = 0.7454 at 100 mW, L(x) callable [Ref B Eq. (3)], κ = c(T + L)/(2l) (FWHM
  12.1 MHz, T02 §3), η_det = ζξ²η [Ref B §4] in the LossRegistry, the OPO-segment detuning coefficients equal
  ``opo_cw_spectrum``'s (T02 §2.4) and A± = (1/(1 − x), 1/(1 + x)) (D51-10);
* the error slopes written in the preset (dither J1 series, A±) equal the blocks' open-loop slopes (1e−6), so
  the unity-gain frequencies hold along a pump sweep;
* D56-9: a lock-1 set-point offset δu = ±0.1 rad moves the lock-2 zero (and Δ) by ±14.64 mrad at x = 0.745
  (tan ε = (A−/A+) tan δu, T14 §3.4) through the closed loops;
* the time-series dither overlay reproduces the analytic dither (θ̃ within 2 %); the MCC low-pass overlay adds
  the LO–pump laser decorrelation (θ̃ grows);
* Fig. 2 through the ESA hook: zero-span means = the metric's R′∓ (radiometer scatter), within 0.3 dB of the
  circuit-noise-subtracted −9.22 / +15.15 dB [Ref B §4];
* YAML round trip (VT-F1 form) of refB / refB_nominal and a recorded CLI run.
"""

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.source.opo_cw import probe_phase_coef, squeezing_angle_coef
from cranberry.blocks.source.opo_testing import REF_B_ETA_D
from cranberry.config.yaml_io import Expr, loads_yaml
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parent / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

AXES = {"pump_w", "h0_hz2_per_hz", "peak_area_hz2", "peak_quantity", "h_ref_rad2_hz", "mech_f_c_hz", "opo_h_m2_hz",
        "opo_f_ugf_hz", "dnu_mcc_hz", "m_dither_rad", "p_pd_w", "p_pr_w", "latency_s", "dac_bits", "loop_gain_frac",
        "pzt1_f_res_hz", "pzt3_f_res_hz", "pzt_q", "pzt_k_rad_per_v"}
X_100 = math.sqrt(0.1 / 0.18)


def test_refB_unknowns_are_named_sweep_axes():
    markers = config.sweep_markers(config.load(rs.PRESET))
    axes = {k.removeprefix("offline.params.") for k in markers if k.startswith("offline.params.")
            and not isinstance(markers[k], Expr)}
    assert axes == AXES
    assert all(isinstance(v, Expr) for k, v in markers.items() if k.removeprefix("offline.params.") not in AXES)
    assert config.sweep_markers(config.load(rs.NOMINAL)) == {}
    # the Fig. 3 pump powers are the pump axis [Ref B §4]
    assert list(markers["offline.params.pump_w"].options) == [0.010, 0.020, 0.050, 0.065, 0.090, 0.100, 0.120, 0.130, 0.150]


def test_refB_offline_products():
    res = run(rs.load(run={"duration_s": 0.012}))
    p = rs.params_of(res.spec)
    assert p["x"] == pytest.approx(X_100, rel=1e-12) and p["x"] == pytest.approx(0.7454, abs=1e-4)   # D56-9
    l_ic = 0.00249 + 0.00222 * X_100                                     # [Ref B Eq. (3)] (T02 §3.4: 0.00414)
    src = res.offline.source_outputs["opo"]
    prov = src.provenance
    assert prov["loss_ic"] == pytest.approx(l_ic, rel=1e-12)
    kappa = 299792458.0 * (0.123 + l_ic) / (2 * 0.5)
    assert prov["kappa_rad_s"] == pytest.approx(kappa, rel=1e-12) == p["kappa_rad_s"]
    assert prov["fwhm_hz"] == pytest.approx(12.1e6, abs=0.05e6)          # T02 §3: γ_B/2π ≈ 12.1 MHz
    assert prov["escape_efficiency"] == pytest.approx(0.123 / (0.123 + l_ic), rel=1e-12)
    reg = res.offline.loss_registry
    assert reg.eta("sqz", "bhd") == pytest.approx(REF_B_ETA_D, rel=1e-12)   # ζ ξ² η [Ref B §4], dark-subtracted
    ids = {e.id: e.eta for e in reg.elements_on("sqz", "bhd")}
    assert ids["optics.propagation"] == 0.99 and ids["mm.overlap"] == pytest.approx(0.988**2) and ids["hd_det.qe"] == 0.998
    assert p["eta_tot"] == pytest.approx(prov["escape_efficiency"] * REF_B_ETA_D, rel=1e-12)
    # OPO-segment detuning coefficients = opo_cw_spectrum's (T02 §2.4, D58-3)
    coef = {(c["beam"], c["effect"]): c["value"] for c in res.spec.topology["coef"]
            if c["segment"] == "OPO" and c["driver"] == "detuning"}
    assert coef[("sqz", "phase")] == pytest.approx(squeezing_angle_coef(X_100, kappa), rel=1e-12)
    assert coef[("pr_bhd", "phase")] == pytest.approx(probe_phase_coef(kappa), rel=1e-12)
    assert coef[("prf", "phase")] == coef[("pr_bhd", "phase")]
    tabs = src.coef_tables
    assert float(np.ravel(tabs["OPO.detuning.squeezed"].values)[0]) == pytest.approx(coef[("sqz", "phase")], rel=1e-12)
    assert float(np.ravel(tabs["OPO.detuning.probe"].values)[0]) == pytest.approx(coef[("pr_bhd", "phase")], rel=1e-12)
    assert (p["a_plus"], p["a_minus"]) == pytest.approx((1 / (1 - X_100), 1 / (1 + X_100)), rel=1e-12)
    assert p["a_minus"] / p["a_plus"] == pytest.approx(0.1459, abs=1e-4)   # D56-9
    assert res.metrics.metric_ii_status == "not_applicable"


@pytest.mark.parametrize("pump_w", [0.020, 0.100, 0.150])
def test_refB_error_slopes_match_blocks(pump_w):
    """Open loops (feedback off, u = offset_v), every noise off, one PZT moved at a time: I1 = k1 sin 2δψ (PZT1 on
    the probe, Φ1 = 2δψ) and I2 = k2 sin δψ (PZT3 on the LO, Δ = δψ) with the preset's slope expressions."""
    dv = 2.0e-4
    toggles = {**rs.NOISE_OFF, "pid1.feedback": False, "pid2.feedback": False}
    rec = {"record": {"slots": ["dem.dem1.i", "dem.dem2.i", "metric.mi.port0"], "decim": 1}}
    out = []
    for v1, v2 in ((dv, 0.0), (0.0, dv)):
        pipe = [{"name": "pid1", "params": {"offset_v": v1}}, {"name": "pid2", "params": {"offset_v": v2}}]
        out.append(run(rs.load(params={"pump_w": pump_w}, run={"duration_s": 0.012}, toggles=toggles, pipeline=pipe,
                               measure=rec)))
    p = rs.params_of(out[0].spec)
    dpsi = p["pzt_k_rad_per_v"] * dv
    assert out[0].series["dem.dem1.i"][-1] == pytest.approx(p["k1_v_per_rad"] * math.sin(2 * dpsi), rel=1e-6)
    assert out[1].series["dem.dem2.i"][-1] == pytest.approx(p["k2_v_per_rad"] * math.sin(dpsi), rel=1e-6)
    assert out[1].series["metric.mi.port0"][-1] == pytest.approx(dpsi, rel=1e-9)   # Δ = PZT3's LO phase


@pytest.mark.parametrize("du", [0.1, -0.1])
def test_refB_lock2_shift_d56_9(du):
    """D56-9: lock 1 held at Φ1 = 2δu (pd1 frame rotated by 2δu) → lock 2's zero, hence Δ, moves by
    ε = atan((A−/A+) tan δu) = ±14.64 mrad at x = 0.7454 (T14 §3.4); the probe sits δu off θ_sq (`cop`)."""
    res = run(rs.load(run={"duration_s": 0.03}, toggles=rs.NOISE_OFF,
                      pipeline=[{"name": "pd1", "params": {"phi_opa_rad": math.pi + 2 * du}}],
                      measure={"record": {"slots": ["lock.lock1.err", "lock.cop.err", "metric.mi.port0"], "decim": 1}}))
    eps = math.atan((1 - X_100) / (1 + X_100) * math.tan(du))
    assert abs(eps) == pytest.approx(14.64e-3, abs=0.005e-3)            # D56-9
    s = res.series
    assert s["lock.lock1.err"][-1] == pytest.approx(2 * du, rel=1e-9)
    assert s["lock.cop.err"][-1] == pytest.approx(du, rel=1e-9)
    assert s["metric.mi.port0"][-1] == pytest.approx(eps, rel=1e-6)
    assert res.metrics.noise_levels["i.port0"]["mean_loop_rad"][0] == pytest.approx(eps, rel=1e-4)


def test_refB_dither_timeseries_and_mcc_overlays():
    """The 295 kHz dither in the time series (RF-rate PD / HD and demodulators) gives the analytic dither's
    residual; the MCC low-pass (Δν 1 MHz) adds |1 − H_MCC|² S_φ of the Ti:Sa noise to Δ."""
    base = run(rs.load(run={"duration_s": 0.15}))
    ts = run(rs.load(run={"duration_s": 0.15}, overlays=(rs.TIMESERIES,)))
    assert ts.spec.toggles["pd1.rf_rate_mode"] and ts.spec.toggles["dither.enabled"]
    th = [r.metrics.noise_levels["i.port0"]["sigma_total_rad"][rs.I_1MHZ] for r in (base, ts)]
    assert th[1] == pytest.approx(th[0], rel=0.02)
    with pytest.warns(UserWarning, match="FC pair on 'MCC'"):   # D66-1 (open): discrete > continuous at 2 MHz
        mcc = run(rs.load(run={"duration_s": 0.15}, overlays=(rs.MCC,)))
    assert mcc.metrics.noise_levels["i.port0"]["sigma_total_rad"][rs.I_1MHZ] > 1.1 * th[0]


def test_refB_fig2_esa_zero_span():
    """Fig. 2 traces (ii)/(iii) through measure.esa: 1 MHz, RBW 30 kHz, VBW 300 Hz, 100 ms, 20 averages."""
    res = run(rs.load(esa=True))
    nl = res.metrics.noise_levels["i.port0"]
    sq, anti = res.measurements["esa"]["squeezed"], res.measurements["esa"]["anti_squeezed"]
    assert sq["mean_db"] == pytest.approx(nl["r_minus_db"][rs.I_1MHZ], abs=0.1)
    assert anti["mean_db"] == pytest.approx(nl["r_plus_db"][rs.I_1MHZ], abs=0.05)
    # circuit-noise-subtracted levels of Ref B §4 (the model uses L(0.745) = 0.00414, not the measured 0.0038)
    assert sq["mean_db"] == pytest.approx(-9.22, abs=0.3)
    assert anti["mean_db"] == pytest.approx(15.15, abs=0.3)


def test_refB_round_trip_and_cli_record(tmp_path):
    for p in (rs.PRESET, rs.NOMINAL):
        s = config.load(p)
        assert config.load(loads_yaml(config.dump_yaml(s))) == s
    from cranberry.__main__ import main

    nominal_short = tmp_path / "refB_short.yaml"
    nominal_short.write_text(f"extends: {rs.ROOT / rs.NOMINAL}\nrun: {{duration_s: 0.02}}\nmeasure: {{esa: null}}\n")
    out = tmp_path / "runs"
    assert main(["run", str(nominal_short), "--out", str(out)]) in (0, None)
    (run_dir,) = list(out.iterdir())
    m = json.loads((run_dir / "metrics.json").read_text())
    assert m["noise_levels"]["i.port0"]["f_s_hz"] == [2.95e5, 1.0e6]
    reg = json.loads((run_dir / "loss_registry.json").read_text())
    assert any(e["id"] == "optics.propagation" and "Ref B" in e["source"] for e in reg["elements"])
