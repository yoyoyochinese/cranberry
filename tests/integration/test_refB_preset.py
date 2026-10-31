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
* the time-series dither overlay reproduces the analytic dither (θ̃ within 2 %); the cavity low-pass overlay
  (MCC on the LO, SHG on the pump, M3B-06) adds |H_SHG − H_MCC|² S_φ: nothing at equal linewidths, θ̃ grows when
  the SHG is much wider than the MCC;
* Phase 5 M3 review fixes (reviews/phase5_M3_refB.md, D76): M3B-01 `gain_handling` (normalized / fixed) and the
  lock stability checks (fixed gains: lock 1 unstable at 150 mW only); M3B-02 laser frequency → OPO detuning
  (θ̃ 1.5° → ≈ 2.2° at the top of the h0 axis, the FM-only metric = 2π√(h0 f_max)·(τ_geo + τ_det), negligible at
  the nominal h0); M3B-03 lock-1 PD after the OPO (Δ above the lock-1 band moves the lock-2 zero by
  (A−/A+)·Δ/(2κ)); M3B-04 dither-filtered slopes (×0.965 / ×0.763 at 100 / 150 mW); M3B-05 `lock2f` = ε;
* Fig. 2 through the ESA hook: zero-span means = the metric's R′∓ (radiometer scatter), within 0.3 dB of the
  circuit-noise-subtracted −9.22 / +15.15 dB [Ref B §4];
* YAML round trip (VT-F1 form) of refB / refB_nominal and a recorded CLI run.
"""

import importlib.util
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.source.opo_cw import probe_phase_coef, probe_transit_coef, squeezing_angle_coef
from cranberry.blocks.source.opo_testing import REF_B_ETA_D
from cranberry.config.yaml_io import Expr, loads_yaml
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import ensure_loaded, make_block
from cranberry.runner.run import offline_block_specs, run

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parent / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

AXES = {"pump_w", "h0_hz2_per_hz", "peak_area_hz2", "peak_quantity", "h_ref_rad2_hz", "mech_f_c_hz", "opo_h_m2_hz",
        "opo_f_ugf_hz", "dnu_mcc_hz", "m_dither_rad", "p_pd_w", "p_pr_w", "latency_s", "dac_bits", "loop_gain_frac",
        "pzt1_f_res_hz", "pzt3_f_res_hz", "pzt_q", "pzt_k_rad_per_v", "gain_handling", "dnu_shg_hz"}
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
    # M3B-03: transmitted probe −1/κ (transit), lock-1 PD probe and lock-2 frame −1/(2κ) (gain frame)
    assert coef[("pr_bhd", "phase")] == pytest.approx(probe_transit_coef(kappa), rel=1e-12) == pytest.approx(-1 / kappa)
    assert coef[("pr_pd", "phase")] == pytest.approx(probe_phase_coef(kappa), rel=1e-12) == pytest.approx(-0.5 / kappa)
    assert coef[("prf", "phase")] == coef[("pr_pd", "phase")]
    tabs = src.coef_tables
    assert float(np.ravel(tabs["OPO.detuning.squeezed"].values)[0]) == pytest.approx(coef[("sqz", "phase")], rel=1e-12)
    assert float(np.ravel(tabs["OPO.detuning.probe"].values)[0]) == pytest.approx(coef[("prf", "phase")], rel=1e-12)
    assert float(np.ravel(tabs["OPO.detuning.probe_transit"].values)[0]) == pytest.approx(coef[("pr_bhd", "phase")],
                                                                                        rel=1e-12)
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
    ε = atan((A−/A+) tan δu) (T14 §3.4): 14.64 mrad at x = 0.7454 with the quasi-static A± (D56-9), 15.2 mrad with
    the dither-filtered Re A±(Ω_m) the blocks use (M3B-04; review: 15.18). The probe sits δu off θ_sq (`cop`).
    M3B-05: lock 2 holds `lock2f` = φ_LO − φ_prf = ε; `lock2` = φ_LO − φ_pr_bhd = ε − δu is not the lock-2 error."""
    res = run(rs.load(run={"duration_s": 0.03}, toggles=rs.NOISE_OFF,
                      pipeline=[{"name": "pd1", "params": {"phi_opa_rad": math.pi + 2 * du}}],
                      measure={"record": {"slots": ["lock.lock1.err", "lock.cop.err", "metric.mi.port0",
                                                    "lock.lock2f.err", "lock.lock2.err"], "decim": 1}}))
    eps_qs = math.atan((1 - X_100) / (1 + X_100) * math.tan(du))
    assert abs(eps_qs) == pytest.approx(14.64e-3, abs=0.005e-3)         # D56-9 (quasi-static)
    p = rs.params_of(res.spec)
    eps = math.atan(p["re_a_minus_m"] / p["re_a_plus_m"] * math.tan(du))
    assert abs(eps) == pytest.approx(15.18e-3, abs=0.05e-3)             # M3B-04 (review §2 B)
    s = res.series
    assert s["lock.lock1.err"][-1] == pytest.approx(2 * du, rel=1e-9)
    assert s["lock.cop.err"][-1] == pytest.approx(du, rel=1e-9)
    assert s["metric.mi.port0"][-1] == pytest.approx(eps, rel=1e-6)
    assert s["lock.lock2f.err"][-1] == pytest.approx(eps, rel=1e-6)      # M3B-05: the quantity lock 2 holds
    assert s["lock.lock2.err"][-1] == pytest.approx(eps - du, rel=1e-6)  # diagnostic: −85 mrad, not a lock error
    assert res.metrics.noise_levels["i.port0"]["mean_loop_rad"][0] == pytest.approx(eps, rel=1e-4)


def test_refB_dither_timeseries_and_cavity_overlays():
    """The 295 kHz dither in the time series (RF-rate PD / HD and demodulators) gives the analytic dither's
    residual. M3B-06: the cavity overlay low-passes the LO (MCC) AND the pump (SHG bow tie), so the Ti:Sa noise
    reaches Δ as |H_SHG − H_MCC|² S_φ: at equal linewidths (refB_nominal, 1 MHz both) θ̃ is the base value; with the
    SHG 10× wider than the MCC (top of its axis) the LO low pass dominates and θ̃ grows."""
    th = lambda r: r.metrics.noise_levels["i.port0"]["sigma_total_rad"][rs.I_1MHZ]  # noqa: E731
    base = run(rs.load(run={"duration_s": 0.15}))
    ts = run(rs.load(run={"duration_s": 0.15}, overlays=(rs.TIMESERIES,)))
    assert ts.spec.toggles["pd1.rf_rate_mode"] and ts.spec.toggles["dither.enabled"]
    assert th(ts) == pytest.approx(th(base), rel=0.02)
    with warnings.catch_warnings():   # D66-1 (open): a clipped FC ledger row warns when discrete > continuous
        warnings.simplefilter("ignore", UserWarning)
        equal = run(rs.load(run={"duration_s": 0.15}, overlays=(rs.MCC,)))
        wide = run(rs.load(params={"dnu_shg_hz": 1.0e7}, run={"duration_s": 0.15}, overlays=(rs.MCC,)))
    assert all(equal.spec.toggles[f"{b}.laser_lowpass"] for b in ("mcc", "shg_sqz", "shg_prf", "shg_ppd"))
    assert th(equal) == pytest.approx(th(base), rel=0.02)              # the two low passes cancel
    assert th(wide) > 1.1 * th(base)


@pytest.mark.parametrize("pump_w, factor", [(0.100, 0.965), (0.150, 0.763)])
def test_refB_m3b04_dither_filtered_slopes(pump_w, factor):
    """M3B-04: the 295 kHz dither sidebands pass the OPO with A±(Ω_m) = 1/(1 ∓ x − iΩ_m/κ); the demodulated slopes
    use Re A±(Ω_m): lock 2 ∝ Re A+(Ω_m), lock 1 ∝ A+(0)Re A+(Ω_m) − A−(0)Re A−(Ω_m), ×0.965 / ×0.763 of the
    quasi-static values at 100 / 150 mW (review §2 B). The ERR blocks get the same effective gains, so
    test_refB_error_slopes_match_blocks holds them to these slopes."""
    spec = rs.load(params={"pump_w": pump_w})
    p = rs.params_of(spec)
    x, w = p["x"], 2 * math.pi * 2.95e5 / p["kappa_rad_s"]
    ap, am = 1 / (1 - x - 1j * w), 1 / (1 + x - 1j * w)
    assert (p["re_a_plus_m"], p["re_a_minus_m"]) == pytest.approx((ap.real, am.real), rel=1e-12)
    assert p["re_a_plus_m"] / p["a_plus"] == pytest.approx(factor, abs=1e-3)
    assert p["dg2_eff"] / (p["a_plus"] ** 2 - p["a_minus"] ** 2) == pytest.approx(factor, abs=2e-3)
    k2_qs = 2 * p["j1_m"] * p["hd_gain_v_per_a"] * 2 * p["responsivity_a_per_w"] * math.sqrt(p["p_lo_w"] * p["p_pr_w"])
    assert p["k2_v_per_rad"] == pytest.approx(k2_qs * p["re_a_plus_m"], rel=1e-12)
    blocks = {b.name: b.params for b in spec.pipeline.children}          # flat pipeline
    assert (blocks["bhdbeat"]["a_plus"], blocks["bhdbeat"]["a_minus"]) == (p["re_a_plus_m"], p["re_a_minus_m"])
    assert blocks["pd1"]["a_plus"] ** 2 == pytest.approx(p["a_plus"] * p["re_a_plus_m"], rel=1e-12)


PUMPS_W = [0.010, 0.020, 0.050, 0.065, 0.090, 0.100, 0.120, 0.130, 0.150]


def test_refB_m3b01_gain_handling_and_stability(record_property):
    """M3B-01: `gain_handling` normalized keeps f_u (K_i ∝ 1/slope(x)); fixed keeps K_i at its 100 mW value, so the
    loop gain follows the error slopes. The offline stability checks (actual slopes, PZT resonance, demodulator
    LPF, τ_L) flag exactly one point of the pump axis × gain_handling: lock 1 with fixed gains at 150 mW (|L| > 1 at
    the −180° crossing near the 14 kHz resonance); warned (ChkWarning), never pruned."""
    ensure_loaded()
    ki, flagged = {}, {}
    for gh in ("normalized", "fixed"):
        for pw in PUMPS_W:
            spec = rs.load(params={"pump_w": pw, "gain_handling": gh})
            p = rs.params_of(spec)
            ki[gh, pw] = (p["ki1_per_s"], p["ki2_per_s"], p["k1_v_per_rad"], p["k2_v_per_rad"])
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                chks = [make_block(b) for b in offline_block_specs(spec) if b.kind == "loop_stability_check"]
            assert sorted(c.name for c in chks) == ["chk.lock1_stability", "chk.lock2_stability"]
            bad = sorted(c.name for c in chks if not c.result["stable"])
            assert sum(issubclass(x.category, ChkWarning) for x in w) == len(bad)
            flagged[gh, pw] = (bad, max(c.result["max_l_at_180"] for c in chks))
    record_property("max_l_at_180", {f"{k[0]}@{k[1]}": round(v[1], 3) for k, v in flagged.items()})
    assert {k: v[0] for k, v in flagged.items() if v[0]} == {("fixed", 0.150): ["chk.lock1_stability"]}
    assert ki["fixed", 0.100] == ki["normalized", 0.100]
    for pw in PUMPS_W:
        k1n, k2n, s1, s2 = ki["normalized", pw]
        assert k1n * s1 == pytest.approx(ki["normalized", 0.100][0] * ki["normalized", 0.100][2], rel=1e-12)   # f_u const
        assert k2n * s2 == pytest.approx(ki["normalized", 0.100][1] * ki["normalized", 0.100][3], rel=1e-12)
        assert ki["fixed", pw][:2] == pytest.approx(ki["fixed", 0.100][:2], rel=1e-12)                       # K_i const


def test_refB_m3b02_laser_frequency_detunes_the_opo():
    """M3B-02: Δ = ω_cav − ω_L, so −2π·δν_L enters the length lock's input (adder `opo_lf`); the cavity follows
    the laser below f_ugf only. Through the OPO coefficients θ_sq lags φ_L by τ_det = 1/(κ(1 + x²)) = 16.9 ns.
    White FM alone: Δ_metric is white up to the metric's f_max = κ/2π, σ = 2π√(h0 f_max)·τ with τ = τ_geo (LO −
    sqz path, the WP-I4 model) or τ_geo + τ_det (loop below f_N, ledger pair with the combined delay above).
    With every noise on, θ̃ rises from 1.55° to ≈ 2.2° at the top of the h0 axis (review §2.3: 1.5° → 2.1°, which
    ignores the coherent geometric part) and by < 0.02° at the nominal h0."""
    tog = {k: v for k, v in rs.NOISE_OFF.items() if k != "tisa.fm_white"}
    sig = {}
    for ld in (False, True):
        res = run(rs.load(params={"h0_hz2_per_hz": 1.0e4, "laser_detuning": ld}, run={"duration_s": 0.1}, toggles=tog))
        p = rs.params_of(res.spec)
        sig[ld] = float(res.metrics.noise_levels["i.port0"]["sigma_total_rad"][rs.I_1MHZ])
        assert p["tau_det_s"] == (pytest.approx(1 / (p["kappa_rad_s"] * (1 + p["x"] ** 2)), rel=1e-12) if ld else 0.0)
    assert p["tau_det_s"] == pytest.approx(16.9e-9, abs=0.05e-9)
    tau_geo = 0.5 / 299792458.0             # (PUMP + SQ)·1.0003 + OPO·1.0 − (LO + MCC)·1.0003 = 0.5 m
    scale = 2 * math.pi * math.sqrt(1.0e4 * p["kappa_rad_s"] / (2 * math.pi))
    assert sig[False] == pytest.approx(scale * tau_geo, rel=0.05)
    assert sig[True] == pytest.approx(scale * (tau_geo + p["tau_det_s"]), rel=0.05)
    th = {}
    for ld in (False, True):
        for h0 in (100.0, 1.0e4):
            res = run(rs.load(params={"h0_hz2_per_hz": h0, "laser_detuning": ld}, run={"duration_s": 0.3}))
            th[ld, h0] = math.degrees(float(res.metrics.noise_levels["i.port0"]["sigma_total_rad"][rs.I_1MHZ]))
    assert th[True, 100.0] - th[False, 100.0] < 0.02                  # nominal h0: negligible
    assert th[False, 1.0e4] == pytest.approx(th[False, 100.0], abs=0.02)
    assert 2.0 < th[True, 1.0e4] < 2.4                                 # top of the h0 axis


def test_refB_m3b03_lock1_pd_after_the_opo():
    """M3B-03: the lock-1 PD sits after the OPO (Fig. 1), its probe pr_pd carries the gain-frame coefficient
    −1/(2κ), so Φ1 contains −Δ/κ: lock 1 holds the probe at the TRUE (detuned) extremum; pr_bhd carries the transit
    coefficient −1/κ. Static Δ with lock 1 closed: unchanged result (Δ/κ)(1 − x²)/(2(1 + x²)) (VT-23). Lock 1 open
    (feedback off: Δ above the lock-1 band, inside the lock-2 band): Φ1 = −Δ/κ, the probe is δ′ = −Δ/(2κ) off the
    extremum and lock 2's zero moves by atan((A−/A+) tan δ′) (review §4: −(A−/A+)·S1·T2·Δ/(2κ), S1 = T2 = 1)."""
    d = 1e-3
    p = rs.params_of(rs.load())
    x, kappa = p["x"], p["kappa_rad_s"]
    want = d * (1 - x * x) / (2 * (1 + x * x))
    out = {}
    for fb in (True, False):
        res = run(rs.load(run={"duration_s": 0.03}, toggles={**rs.NOISE_OFF, "pid1.feedback": fb},
                          pipeline=[{"name": "opo_ll", "params": {"offset_rad_s": d * kappa}}],
                          measure={"record": {"slots": ["lock.lock1.err", "metric.mi.port0", "lock.lock2f.err"],
                                              "decim": 1}}))
        out[fb] = {k: v[-1] for k, v in res.series.items()}
    assert out[True]["lock.lock1.err"] == pytest.approx(0.0, abs=1e-12)
    assert out[True]["metric.mi.port0"] == pytest.approx(want, rel=1e-6)
    assert out[False]["lock.lock1.err"] == pytest.approx(-d, rel=1e-9)
    eps = math.atan(p["re_a_minus_m"] / p["re_a_plus_m"] * math.tan(-d / 2))
    assert out[False]["lock.lock2f.err"] == pytest.approx(eps, rel=1e-6)
    assert out[False]["metric.mi.port0"] == pytest.approx(want + eps, rel=1e-6)
    assert out[False]["metric.mi.port0"] / want == pytest.approx(0.47, abs=0.01)   # 0.143 − 0.076 (Δ/κ units)


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
