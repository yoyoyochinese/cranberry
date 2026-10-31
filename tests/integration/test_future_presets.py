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
* D31-3 / D47: the laser decorrelation component equals 2π²h0Δτ (time loop below f_N + ledger above it), the
  squeezed-only component contains the demux switch phase (2S(κπδV/V_π)², T11 §3.5).

The deliberately bad presets are checked in test_bad_presets.py; the 32-port scale estimate in
test_future_scale_estimate.py.
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import Expr
from cranberry.runner.decompose import decompose_metric_ii
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

BASES = ("F_wg100", "F_wg6", "F_bor6", "F_liu20")
OFF = {"laser.fm_white": False, "env.mech.ports.mech": False, "env.mech.lo.mech": False, "env.mech.dist.mech": False,
       "env.mech.fs_lock.mech": False, "env.mech.fs_sq.mech": False, "dmx0.drive_phase_noise": False,
       "dmx1.drive_phase_noise": False, "dmx0.bias_drift": False, "dmx1.bias_drift": False,
       "env.thermal.F0.drift": False, "env.thermal.F1.drift": False, "env.thermal.LO.drift": False,
       "env.thermal.xtal.drift": False, "liu0.shot_noise": False, "liu1.shot_noise": False,
       "liu0.electronic_noise": False, "liu1.electronic_noise": False}


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
    assert fs.chk_codes(res) == []                      # nominal points are feasible


def test_decomposition_F_liu20_linear_and_laser_dominated():
    comps = decompose_metric_ii(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.03}))
    c = comps[1]
    assert c.status == "ok" and abs(c.linearity_ratio - 1.0) <= 0.05
    v = c.variances_rad2
    assert all(x > 0 for x in v.values())
    # D31-3 / D47: at h0 = 64 Hz²/Hz and Δτ = 6.2 µs (port 31 at 20 MHz, N_b = 4) the laser term dominates
    assert v["laser_decorrelation_oob"] > 10 * max(v["locked_path_residual"], v["common_points"], v["squeezed_only"])


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
    ("presets/bad/bad_mw_probe_samepol.yaml", (), None, {"sbs_gain"}),
], ids=["rep_rate_tail", "toggle_rate_Nb1_100MHz", "rise_time_10ns", "sbs_cw_lo", "toggle_rate_liu", "sbs_probe"])
def test_chk_warnings_for_infeasible_settings(preset, overlays, params, expected):
    res = fs.build_check(preset, overlays, params=params)
    assert set(fs.chk_codes(res)) == expected


def test_bor100_tail_matches_t17():
    res = fs.build_check("presets/bad/bad_bor_100MHz.yaml", ())
    tail = res.offline.source_outputs["src"].tail_into_next_slot
    assert 0.10 < tail < 0.16                           # T17 V-P11: ε_tail(10 ns) ≈ 0.138 (x_ss = 2.4); 1.3 here
    assert fs.block(res, "chk_tail").result["f_rep_max_hz"] == pytest.approx(2 * math.pi * 66e6 / math.log(8.8e3), rel=1e-3)


# ------------------------------------------------------------------------------------------------ D31 / D47 effects


@pytest.mark.parametrize("port_k,h0", [(31, 64.0), (31, 1.0), (1, 64.0)])
def test_laser_decorrelation_equals_white_fm_requirement(port_k, h0):
    """Metric (ii) with only the laser on: time-loop variance (≤ f_N, minus the lock-band cancellation of the same-laser
    probe) + ledger variance above f_N = 2π²h0·Δτ_k (conv §7 delayed self-interference), Δτ_k = n_g(L_k − L_0)/c."""
    tg = {**OFF, "laser.fm_white": True}
    res = run(fs.load(fs.nominal("F_liu20"), params={"h0_hz2_per_hz": h0, "port_k": port_k}, run={"duration_s": 0.03},
                      toggles=tg))
    p = fs.params(res)
    dtau = p["tau_k_s"] - p["tau_0_s"]
    loop = res.metrics.metric_ii[1].stats.var
    hf = res.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
    assert (loop + hf) / (2 * math.pi**2 * h0 * dtau) == pytest.approx(1.0, abs=0.15)
    assert res.metrics.noise_levels["i.port1"]["hf_variance_rad2"] < 1e-9 * hf   # EM1 marker: cancels in metric (i)


def test_squeezed_only_component_contains_demux_switch_phase():
    tg = {**OFF, "dmx0.drive_phase_noise": True, "dmx1.drive_phase_noise": True}
    res = run(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.03}, toggles=tg))
    p = fs.params(res)
    per_port = p["dmx_n_stages"] * (p["dmx_kappa_pm"] * math.pi * p["dmx_drive_noise_v_rms"] / p["dmx_v_pi_v"]) ** 2
    assert res.metrics.metric_ii[1].stats.var == pytest.approx(2 * per_port, rel=0.08)   # T11 §3.5, two ports
    # not seen by the Liu probe: what is left in its error is the 16-bit DAC limit cycle (D65-3: mrad class)
    assert np.std(res.series["lock.L1.err"][-1000:]) < 0.01
