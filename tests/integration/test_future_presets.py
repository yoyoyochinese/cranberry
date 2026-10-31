"""WP-I5 (M4–M6): the future 32-port presets F_wg100 / F_wg6 / F_bor6 / F_liu20, their overlay axes and lock topologies
(PLAN §5.4, §7.1; 00_summary §2.2–§2.4; D08–D13, D30–D32, D34, D38–D47, D50, D53; G13 / D54).

* every open quantity of a base preset is a named sweep axis under ``offline.params`` and ``<name>_nominal`` resolves
  all of them (D34);
* every preset / overlay combination of ``presets/future/manifest.yaml`` builds and validates through the runner
  (full matrix ``slow``; a smoke subset runs by default);
* a short compiled run of each base preset gives metric (i) per port and metric (ii) of port k vs reference port 0,
  noise levels, η at GBS and BHD; the four-component decomposition of metric (ii) is linear (T04 V12, 5 %);
* CHK warnings appear exactly for the infeasible settings: rep-rate tail (bad_bor_100MHz), demux toggle rate (N_b = 1
  at 100 MHz) and rise time, gate shorter than the mode, SBS (CW LO in a 3 km delivery fibre; 10 mW probe peak);
* D31-3 / D47 / D78: with the locks off the laser decorrelation equals 2π²h0Δτ (time loop below f_N + ledger above
  it, within the D65-4 interpolation deficit); the Liu lock multiplies it by R(f_u) ≥ 1, equal to the closed-loop path
  model within 3 % (Phase 5 M6-05); the squeezed-only component contains the demux switch phase (2S(κπδV/V_π)²,
  T11 §3.5).
* Phase 5 M6 fixes: the L-Liu loop stability CHK (M6-01: unstable points recorded, valid = False), component 2 over
  the run vs a 30 s window (M6-02), the *_nominal presets keep their own values and the stray-light CHK (M6-03), laser
  drift through the TEC offload (M6-04), components with their ledger rows (M6-06), the OPO detuning coefficient vs
  x_ss (M6-07).

The deliberately bad presets are checked in test_bad_presets.py; the 32-port scale estimate in
test_future_scale_estimate.py.
"""

import importlib.util
import itertools
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.dmx.demux import stage_counts
from cranberry.config.yaml_io import Expr, load_yaml
from cranberry.runner.decompose import decompose_metric_ii
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

BASES = ("F_wg100", "F_wg6", "F_bor6", "F_liu20")
OFF = {"laser.fm_white": False, "laser.fm_flicker": False, "laser.drift": False, "env.mech.ports.mech": False,
       "env.mech.lo.mech": False, "env.mech.dist.mech": False, "env.mech.fs_lock.mech": False,
       "env.mech.fs_sq.mech": False, "dmxsw.drive_phase_noise": False, "dmxsw.bias_drift": False,
       "env.thermal.F0.drift": False, "env.thermal.F1.drift": False, "env.thermal.LO.drift": False,
       "env.thermal.xtal.drift": False, "env.thermal.dist.drift": False, "env.thermal.dist0.drift": False,
       "env.thermal.dist1.drift": False, "env.dist.window.enabled": False, "liu0.shot_noise": False,
       "liu1.shot_noise": False, "liu0.electronic_noise": False, "liu1.electronic_noise": False}
LOCK_OFF = {"pid0f.feedback": False, "pid1f.feedback": False, "pid0s.feedback": False, "pid1s.feedback": False}


# ------------------------------------------------------------------------------------------------ manifest / axes


def test_manifest_covers_every_new_file():
    m = fs.manifest()
    listed = {f"presets/future/{b}.yaml" for b in m["base"]} | {f"presets/future/{b}_nominal.yaml" for b in m["base"]}
    listed |= {f"presets/bad/{b}.yaml" for b in m["bad"]}
    for o in [*m["overlays_all"], *m["overlays_opo"], *m["overlays_waveguide"],
              *(x for c in m["combinations"] for x in c["overlays"])]:
        listed.add(o if o.startswith("presets/") else f"presets/axes/{o}.yaml")
    listed |= {"presets/future/F_common.yaml", "presets/future/points/F_common_nominal.yaml"}
    axes = {"source", "demux", "path_build", "fiber", "probe", "largerange", "lo", "bhd", "opo_lock", "controller",
            "ports"}
    files = {str(p.relative_to(fs.ROOT)) for d in ("future", "bad", "locks") for p in (fs.ROOT / "presets" / d).rglob("*.yaml")}
    files |= {str(p.relative_to(fs.ROOT)) for a in axes for p in (fs.ROOT / "presets" / "axes" / a).glob("*.yaml")}
    files.discard("presets/future/manifest.yaml")
    assert files - listed == set()
    assert listed - files == set()
    assert len(m["bad"]) == 16                       # 00_summary §2.5


@pytest.mark.parametrize("base", BASES)
def test_open_quantities_are_named_sweep_axes(base):
    markers = config.sweep_markers(config.load(fs.ROOT / f"presets/future/{base}.yaml"))
    axes = {k for k, v in markers.items() if not isinstance(v, Expr)}
    assert len(axes) >= 30
    assert all(k.startswith("offline.params.") for k in axes)          # D34: axes live in one place
    for a in ("port_k", "h0_hz2_per_hz", "p_pr_avg_w", "probe_duty", "rho_ports", "K_T_pm_rad_per_m_K",
              "dmx_t_stage", "f_u_hz", "tec_tau_th_s", "D_pol_rad2_per_s"):
        assert f"offline.params.{a}" in axes
    assert config.sweep_markers(config.load(fs.ROOT / fs.nominal(base))) == {}


def _smoke():
    m = fs.manifest()
    out = [(fs.nominal(b), ()) for b in m["base"]]
    out += [(fs.nominal("F_liu20"), (o,)) for o in m["overlays_all"] if o.startswith("presets/locks")]
    out += [(fs.nominal(c["base"]), tuple(c["overlays"])) for c in m["combinations"][:5]]
    out += [(f"presets/bad/{b}.yaml", ()) for b in m["bad"]]
    return out


@pytest.mark.parametrize("preset,overlays", _smoke(), ids=lambda x: Path(x).stem if isinstance(x, str) else "+".join(
    Path(o).stem for o in x))
def test_manifest_smoke_builds(preset, overlays):
    res = fs.build_check(preset, overlays)
    assert set(res.metrics.metric_ii) == {1}


@pytest.mark.slow
@pytest.mark.parametrize("preset,overlays", fs.combinations(), ids=lambda x: Path(x).stem if isinstance(x, str) else
                         "+".join(Path(o).stem for o in x))
def test_manifest_every_combination_builds(preset, overlays):
    res = fs.build_check(preset, overlays)
    assert 0 in res.metrics.metric_i and 1 in res.metrics.metric_ii


# ------------------------------------------------------------------------------------------------ short runs


@pytest.fixture(scope="module")
def nominal_runs():
    return {b: run(fs.load(fs.nominal(b), run={"duration_s": 0.03})) for b in BASES}


@pytest.mark.parametrize("base", BASES)
def test_base_preset_short_run_metrics(nominal_runs, base):
    res = nominal_runs[base]
    m = res.metrics
    assert set(m.metric_i) == {0, 1} and set(m.metric_ii) == {1} and m.metric_ii_status == "not_run"
    # the Liu locks acquired X = π/2 and hold it (port 0: short fibre, mrad; port k: within the laser self-beat)
    assert abs(res.series["lock.L0.err"][-1]) < 0.02
    assert np.std(res.series["lock.L0.err"][len(res.series["lock.L0.err"]) // 2:]) < 0.01
    assert np.std(res.series["lock.L1.err"][len(res.series["lock.L1.err"]) // 2:]) < 0.3
    # metric (ii): finite, below the 0.3 rad class at the nominal point; metric (i) of port 0 is mrad-level
    s_ii = m.metric_ii[1].stats
    assert 0.0 < s_ii.std < 0.3 and abs(s_ii.mean) < 0.2
    assert m.metric_i[0].stats.std < 0.05
    # η at the GBS input and at the BHD reported separately (D46): BHD < GBS, port k < port 0 (longer fibre)
    e0, e1 = m.efficiency[0], m.efficiency[1]
    assert 0.0 < e1.eta_bhd < e1.eta_gbs < e0.eta_gbs < 1.0
    nl = m.noise_levels["ii.port1"]
    assert nl["loss_target"] == "gbs" and nl["v_minus_db"] < 0.0 and nl["hf_variance_rad2"] > 0.0
    assert m.noise_levels["i.port1"]["loss_target"] == "bhd"
    assert fs.chk_codes(res) == [] and res.valid        # nominal points are feasible, the loop stable (M6-01)
    # stray light at the nominal P̄ within T17's 5.8e−4 per mode (M6-03: 1 nW for the OPO presets, 1 µW for the others)
    assert fs.stray_light(res)["photons_per_mode"] < fs.params(res)["stray_n_mode_max"]


def test_decomposition_F_liu20_run_window_vs_30s_window():
    """Components include their ledger rows (M6-06) and component 2 depends on the metric window (M6-02): inside the
    30 ms run the laser dominates (D74-4); over a 30 s window the common points' f⁻² band adds ≈ 55 mrad (M6 review)."""
    out = {}
    for t_cal in (0.02, 30.0):                                     # 0.02 s ≤ T_rec = 25 ms: no window row
        c = decompose_metric_ii(fs.load(fs.nominal("F_liu20"), params={"t_cal_s": t_cal}, run={"duration_s": 0.03}))[1]
        assert c.status == "ok" and abs(c.linearity_ratio - 1.0) <= 0.05
        out[t_cal] = c.variances_rad2
    run_w, win = out[0.02], out[30.0]
    print("\n[metric (ii) components, F_liu20 nominal, rad] 30 ms: " + ", ".join(f"{k} {math.sqrt(v):.4f}" for k, v in
          run_w.items()) + " | 30 s: " + ", ".join(f"{k} {math.sqrt(v):.4f}" for k, v in win.items()))
    assert all(x > 0 for x in run_w.values())
    assert run_w["laser_decorrelation_oob"] > 10 * max(run_w["locked_path_residual"], run_w["common_points"],
                                                       run_w["squeezed_only"])
    # laser component = time loop + ledger above f_N = R·2π²h0Δτ (1.02–1.05 at f_u = 1 kHz, M6 review §2.1)
    p = fs.params(fs.build_check(fs.nominal("F_liu20"), ()))
    assert run_w["laser_decorrelation_oob"] / (2 * math.pi**2 * p["h0_hz2_per_hz"] * (p["tau_k_s"] - p["tau_0_s"])) == \
        pytest.approx(1.03, abs=0.03)
    assert math.sqrt(win["common_points"]) == pytest.approx(0.055, rel=0.1)
    assert win["common_points"] > 0.3 * win["laser_decorrelation_oob"]       # "laser dominant" only inside the run
    for k in ("locked_path_residual", "laser_decorrelation_oob", "squeezed_only"):
        assert win[k] == pytest.approx(run_w[k], rel=1e-9)                    # same seeds; only the window row differs


def test_common_point_window_row_and_temperature_coefficient():
    """M6-02: PRS_k carries k_T·L (2 m ≈ 88 rad/K); the window row is the f⁻² band [1/t_cal, 1/T_rec] of PRS0/1, REF0/1
    plus the branch drift (rT)²/12, and is 0 when t_cal ≤ T_rec."""
    from scipy.integrate import quad

    res = fs.build_check(fs.nominal("F_liu20"), ())
    p = fs.params(res)
    coef = {c["id"]: c["value"] for c in res.spec.topology["coef"]}
    assert coef["PRS1.dT.all.phase"] == pytest.approx(88.0, rel=1e-12)
    t_rec = p["t_rec_s"]

    def s(f, L):
        return (L / 1e3) * p["h_ref_dist_rad2_hz"] / f**2 / ((1 + (p["mech_f_lo_hz_dist"] / f) ** 2)
                                                           * (1 + (f / p["mech_f_c_hz"]) ** 2))

    band = sum(2 * quad(lambda u, L: s(math.exp(u), L) * math.exp(u), -math.log(p["t_cal_s"]), -math.log(t_rec),
                        args=(L,), limit=400)[0] for L in (p["L_prs_m"], p["L_refk_m"]))      # ∫ over ln f
    drift = (p["k_T_rad_per_m_K"] * (p["L_prs_m"] - p["L_refk_m"]) * p["dr_T_dist_k_per_s"]) ** 2 \
        * (p["t_cal_s"] ** 2 - t_rec**2) / 12
    assert p["cp_window_var_rad2"] == pytest.approx(band + drift, rel=1e-6)
    assert fs.params(fs.build_check(fs.nominal("F_liu20"), (), params={"t_cal_s": 1e-6}))["cp_window_var_rad2"] == 0.0


@pytest.mark.slow
@pytest.mark.parametrize("base", ("F_wg100", "F_wg6", "F_bor6"))
def test_decomposition_other_bases(base):
    c = decompose_metric_ii(fs.load(fs.nominal(base), run={"duration_s": 0.05}))[1]
    assert c.status == "ok" and abs(c.linearity_ratio - 1.0) <= 0.05
    assert c.variances_rad2["laser_decorrelation_oob"] == max(c.variances_rad2.values())


# ------------------------------------------------------------------------------------------------ CHK (G16-j)


@pytest.mark.parametrize("preset,overlays,params,expected", [
    ("presets/bad/bad_bor_100MHz.yaml", (), None, {"rep_rate_tail", "gate_shorter_than_mode", "demux_rise_time"}),
    (fs.nominal("F_wg100"), ("demux/pockels_tree",), None, {"demux_toggle_rate"}),
    (fs.nominal("F_wg100"), (), {"dmx_t_rise_s": 1.0e-8}, {"demux_rise_time"}),
    (fs.nominal("F_wg100"), ("lo/cw_digital",), None, {"sbs_gain"}),       # 10 mW CW LO, 3.2 km: G_B ≈ 17 > 5
    (fs.nominal("F_liu20"), (), {"dmx_f_tog_max_hz": 1.0e+6}, {"demux_toggle_rate"}),
    ("presets/bad/bad_mw_probe_samepol.yaml", (), None, {"sbs_gain", "stray_light_limit"}),
    (fs.nominal("F_liu20"), (), {"p_pr_avg_w": 1.0e-6}, {"stray_light_limit"}),   # 1 µW orthogonal, 1.27 km: 8e−3 (M6-03)
    (fs.nominal("F_liu20"), (), {"f_u_hz": 1.0e+4}, {"loop_unstable"}),          # M6-01: f_u·Q/f_res = 5 at 20 kHz, Q 10
], ids=["rep_rate_tail", "toggle_rate_Nb1_100MHz", "rise_time_10ns", "sbs_cw_lo", "toggle_rate_liu", "sbs_probe",
        "stray_1uW_liu", "loop_unstable"])
def test_chk_warnings_for_infeasible_settings(preset, overlays, params, expected):
    res = fs.build_check(preset, overlays, params=params)
    assert set(fs.chk_codes(res)) == expected


def test_bor100_tail_matches_t17():
    res = fs.build_check("presets/bad/bad_bor_100MHz.yaml", ())
    tail = res.offline.source_outputs["src"].tail_into_next_slot
    assert 0.10 < tail < 0.16                           # T17 V-P11: ε_tail(10 ns) ≈ 0.138 (x_ss = 2.4); 1.3 here
    assert fs.block(res, "chk_tail").result["f_rep_max_hz"] == pytest.approx(2 * math.pi * 66e6 / math.log(8.8e3), rel=1e-3)


# ------------------------------------------------------------------------------------------------ D31 / D47 effects


def _laser_only(params, toggles=None):
    res = run(fs.load(fs.nominal("F_liu20"), params=params, run={"duration_s": 0.03},
                      toggles={**OFF, "laser.fm_white": True, **(toggles or {})}))
    p = fs.params(res)
    v = res.metrics.metric_ii[1].stats.var + res.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
    return res, p, v


@pytest.mark.parametrize("port_k,h0", [(31, 64.0), (31, 1.0), (1, 64.0)])
def test_laser_decorrelation_open_loop_equals_white_fm_requirement(port_k, h0):
    """Locks off, only the laser: time-loop variance (≤ f_N) + ledger variance above f_N = 2π²h0·Δτ_k (conv §7 delayed
    self-interference), Δτ_k = n_g(L_k − L_0)/c, less the linear-interpolation deficit of the loop (D65-4: 1–4 %)."""
    res, p, v = _laser_only({"h0_hz2_per_hz": h0, "port_k": port_k}, LOCK_OFF)
    ratio = v / (2 * math.pi**2 * h0 * (p["tau_k_s"] - p["tau_0_s"]))
    assert 0.95 < ratio <= 1.005
    assert res.metrics.noise_levels["i.port1"]["hf_variance_rad2"] < 1e-9 * v   # EM1 marker: cancels in metric (i)


@pytest.mark.parametrize("f_u", [1.0e3, 1.0e4])
def test_liu_lock_increases_white_fm_decorrelation(f_u):
    """D78-1 / M6-05: the same-laser Liu lock does not cancel the white-FM decorrelation, it multiplies it by
    R(f_u) = locked/unlocked ≥ 1, which equals the closed-loop path model (fs.liu_lock_ratio) within 3 % (stretcher
    resonance off, as the model; review: 1.02 at 1 kHz, 1.20 at 10 kHz for port 31 of F_liu20)."""
    prm = {"f_u_hz": f_u}
    no_res = {"pzt0.resonance": False, "pzt1.resonance": False}
    _, p, on = _laser_only(prm, no_res)
    _, _, off = _laser_only(prm, {**no_res, **LOCK_OFF})
    r_sim, r_model = on / off, fs.liu_lock_ratio(p, f_u, 1.0e6)
    assert r_model > 1.0 and r_sim > 1.0
    assert r_sim == pytest.approx(r_model, rel=0.03)


def test_squeezed_only_component_contains_demux_switch_phase():
    tg = {**OFF, "dmxsw.drive_phase_noise": True}
    res = run(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.03}, toggles=tg))
    p = fs.params(res)
    n = stage_counts(int(p["K_ports"]), p["dmx_layout"])
    per_stage = (p["dmx_kappa_pm"] * math.pi * p["dmx_drive_noise_v_rms"] / p["dmx_v_pi_v"]) ** 2
    # T11 §3.5: drive jitter independent per port (own pulse times) even on shared switches (M4-01)
    assert res.metrics.metric_ii[1].stats.var == pytest.approx((n[0] + n[int(p["port_k"])]) * per_stage, rel=0.08)
    # not seen by the Liu probe: what is left in its error is the 16-bit DAC limit cycle (D65-3: mrad class)
    assert np.std(res.series["lock.L1.err"][-1000:]) < 0.01


# ------------------------------------------------------------------------------------------------ Phase 5 M6 fixes


def _presets_with_overlays():
    out = [fs.nominal(b) for b in fs.manifest()["base"]]
    return out + [f"presets/bad/{b}.yaml" for b in fs.manifest()["bad"]
                  if load_yaml(fs.ROOT / f"presets/bad/{b}.yaml").get("overlays")]


@pytest.mark.parametrize("preset", _presets_with_overlays(), ids=lambda x: Path(x).stem)
def test_presets_keep_their_own_values(preset):
    """M6-03 / D74-7: a preset's own offline.params survive its own ``overlays:`` (F_bor6_nominal's 1 nW was overwritten
    by points/F_common_nominal.yaml's 1 µW)."""
    own = (load_yaml(fs.ROOT / preset).get("offline") or {}).get("params") or {}
    res = config.load(fs.ROOT / preset).offline["params"]
    for k, v in own.items():
        if isinstance(v, (Expr, str, bool)):
            continue
        if k.endswith("_db"):                                    # D56-1: <x>_db → <x>_lin at load
            k, v = k[:-3] + "_lin", 10 ** (float(v) / 10)
        assert float(res[k]) == pytest.approx(float(v), rel=1e-12), k


def test_nominal_probe_powers():
    p = {b: config.load(fs.ROOT / fs.nominal(b)).offline["params"]["p_pr_avg_w"] for b in BASES}
    assert p == {"F_wg100": 1e-6, "F_wg6": 1e-6, "F_bor6": 1e-9, "F_liu20": 1e-9}   # T14 §8.2 / D50 (OPO: ≲ 1–3 nW)


def _chk_stability(params):
    from cranberry.core.registry import ensure_loaded, make_block
    from cranberry.runner.run import offline_block_specs

    ensure_loaded()
    spec = fs.load(fs.nominal("F_liu20"), params=params)
    (b,) = [make_block(x) for x in offline_block_specs(spec) if x.name == "chk.lock_stability"]
    return b


@pytest.mark.slow
def test_lock_stability_grid_counts_unstable_points():
    """M6-01: over f_u × f_res × Q × latency of F_common, 40 / 54 points are linearly unstable (M6 review); each is
    recorded by the CHK, none pruned (D75)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid = itertools.product([1e3, math.sqrt(3e7), 3e4], [5e3, 3.25e4, 6e4], [3.0, math.sqrt(300.0), 100.0],
                                 [1e-6, 2e-6])
        res = [_chk_stability({"f_u_hz": a, "pzt_f_res_hz": b, "pzt_q": c, "latency_s": d}).result["stable"]
               for a, b, c, d in grid]
    assert sum(not x for x in res) == 40


def test_unstable_lock_point_is_flagged_invalid():
    """M6-01: f_u = 10 kHz at the nominal 20 kHz, Q 10 resonance: the CHK warns, RunResult.valid is False and the metric
    (ii) of the oscillating loop is far above the locked value (review: 87× the open-loop laser term)."""
    assert _chk_stability({}).result["stable"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bad = run(fs.load(fs.nominal("F_liu20"), params={"f_u_hz": 1.0e4}, run={"duration_s": 0.03}))
        good = run(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.03}))
    assert [e["code"] for e in bad.events] == ["loop_unstable"] and not bad.valid
    assert good.valid
    assert bad.metrics.metric_ii[1].stats.std > 3 * good.metrics.metric_ii[1].stats.std


@pytest.mark.slow
def test_laser_drift_leaks_through_tec_offload_orthogonal_pm():
    """M6-04: a laser drift ν̇ puts 2πδν·τ_k into the Liu error; the TEC takes it over and moves the fast axis (probe) by
    (k_T + |K_T|)LΔT but the slow axis (squeezed) by k_T·LΔT, so dΘ_k/dt = (|K_T|/(k_T + |K_T|))·2πν̇τ_k ≈ 0.045·2πν̇τ_k
    (review: 1 %); a same-polarization probe has no leak."""
    quiet = {**OFF, "laser.fm_white": False, "opo.noise.length_noise": False, "laser.drift": True}
    prm = {"latency_s": 1.0e-5, "tec_rate_div": 1000, "tau_th_spool_s": 1.0, "tec_tau_th_s": 1.0, "f_12_hz": 0.1,
           "laser_nu_dot_hz_per_s": 3.0e4}
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for lab, ov in (("orth", ()), ("same", ("probe/pol_same",))):
            out[lab] = run(fs.load(fs.nominal("F_liu20"), params=prm, run={"f_sim_hz": 1.0e5, "duration_s": 10.0},
                                   overlays=ov, toggles=quiet))
    p = fs.params(out["orth"])
    kt, kk = p["k_T_rad_per_m_K"], abs(p["K_T_pm_rad_per_m_K"])
    expect = kk / (kt + kk) * 2 * math.pi * p["laser_nu_dot_hz_per_s"] * p["tau_k_s"]
    s_o = fs.slope(out["orth"].t_s, out["orth"].series["metric.mii.port1"], 6.0)
    s_s = fs.slope(out["same"].t_s, out["same"].series["metric.mii.port1"], 6.0)
    assert kk / (kt + kk) == pytest.approx(0.045, abs=0.002)
    assert s_o == pytest.approx(expect, rel=0.05)
    assert abs(s_s) < 0.01 * abs(s_o)


@pytest.mark.parametrize("base,x_ss", [("F_liu20", 1.0), ("F_liu20", 2.5645), ("F_liu20", 4.129), ("F_bor6", 1.3),
                                       ("F_bor6", 2.4)])
def test_opo_detuning_coefficient_follows_x_ss(base, x_ss):
    """M6-07: c_len_sq_per_kappa (θ_sq per Δ_s/κ_c, length change Δ_p = 2Δ_s) follows the x_ss axis: equal to the
    source's own WP-B11 angle_coefficients at that x_ss within 5e−4 (review: −1.040 … −1.154 Liu, −1.063 … −1.091
    Borealis)."""
    res = fs.build_check(fs.nominal(base), (), params={"x_ss": x_ss})
    src = fs.block(res, "src")
    cs, cp, _ = src.angle_coefficients(x_ss)
    assert fs.params(res)["c_len_sq_per_kappa"] == pytest.approx((cs + 2 * cp) * src.kappa_c, abs=5e-4)
