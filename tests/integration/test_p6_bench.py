"""Phase 6 WP6-X (D85-12, D85-9, D86, D85-15): the axis-1 bench configurations (presets/bench/*, presets/axes/bench/*) and
the probe overlays actuator_pathlength / two_tone / detect_heterodyne_refshift build and run 10–50 ms; the resolved
values that carry the physics (fibre length, delay mismatch, loss, controller rates and checks, the two-tone bias and
the heterodyne offset) are checked against their formulas."""

import dataclasses
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import load_yaml
from cranberry.core.registry import ensure_loaded, make_block
from cranberry.path.topology import as_topology
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

pytestmark = [pytest.mark.numba, pytest.mark.integration]
ROOT = Path(__file__).resolve().parents[2]
C = 299792458.0


def load(preset, overlays=(), duration=0.01, drop_esa=True, **point):
    tree = {"extends": str(ROOT / preset), "run": {"duration_s": duration}}
    if drop_esa:
        tree["measure"] = {"esa": None}   # zero-span hooks need ≥ 0.15 s records
    spec = config.load(tree, overlays=list(overlays))
    if config.sweep_markers(spec):
        axes = {f"offline.params.{k}": [v] for k, v in point.items()}
        spec = expand(spec, axes)[0].spec
    return spec


def go(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(spec)
    assert res.valid and not res.events
    return res


def chk(spec, name):
    ensure_loaded()
    (b,) = [make_block(x) for x in offline_block_specs(spec) if x.name == name]
    return b.result


def find(node, name):
    if getattr(node, "name", None) == name and hasattr(node, "params"):
        return node
    for f in dataclasses.fields(node):
        for c in (v if isinstance(v := getattr(node, f.name), tuple) else (v,)):
            if dataclasses.is_dataclass(c) and not isinstance(c, type) and (r := find(c, name)) is not None:
                return r
    return None


# ------------------------------------------------------------------------------------------------ (a) Ref A + fibre

REFA = "presets/bench/refA_fiber_bench.yaml"


@pytest.mark.parametrize("route", ["fiber_matched", "free_space"])
def test_refA_fiber_geometry_delay_and_loss(route):
    L = 31.6
    spec = load(REFA, L_sq_fiber_m=L, lo_route=route, eta_fiber_out=0.95)
    topo = as_topology(spec.topology)
    assert topo.segments["SQ"].kind == "pm" and topo.segments["SQ"].length_m == pytest.approx(L)
    assert topo.segments["LO"].length_m == pytest.approx(5.0 + (L if route == "fiber_matched" else 0.0))
    # delay mismatch LO − squeezed from the topology, and the laser's FC ledger pair reads the same lengths
    d = topo.timing("lo").tau_s - topo.timing("sqz").tau_s
    want = (7.6400 - 1.5994 - (0.0 if route == "fiber_matched" else 1.468 * L)) / C
    assert d == pytest.approx(want, abs=2e-11)
    pair = find(spec.pipeline, "koheras").params["oob_pairs"][0]
    assert pair["tau_ref_s"] == pytest.approx(topo.timing("sqz").tau_s, rel=1e-12)
    assert pair["tau_fc_s"] == pytest.approx(topo.timing("lo").tau_s, rel=1e-12)   # LO + LO_fs (0 m) + FC
    # fibre loss on the squeezed path: 0.98 × 0.95 coupling, 0.5 dB/km
    el = {e["id"]: e for e in next(b for b in spec.offline["blocks"] if b["name"] == "optics")["params"]["elements"]}
    assert el["sq_fiber_coupling"]["eta"] == pytest.approx(0.98 * 0.95)
    assert el["sq_fiber_attenuation"]["eta"] == pytest.approx(10 ** (-0.05 * L / 1000))
    assert el["components"]["eta"] == 0.99   # refA's element kept
    # lock-2 demodulation phase absorbs the probe's 2π f_sh Δτ
    dpr = 2 * math.pi * 1e6 * (1.468 * L - 1.0003) / C
    assert find(spec.pipeline, "dem2").params["phi_d_rad"] == pytest.approx(-math.pi / 2 - dpr)


def test_refA_fiber_100m_runs_and_offload_holds_eom():
    res = go(load(REFA, duration=0.02, L_sq_fiber_m=100.0, lo_route="fiber_matched", eta_fiber_out=0.98))
    nl = res.metrics.noise_levels["i.port0"]
    assert abs(nl["mean_loop_rad"]) < 3e-3 and nl["std_loop_rad"] < 0.03
    assert np.max(np.abs(res.series["v.dac2"])) < 0.45   # EOM inside ±0.5 V: the stretcher offload takes the fibre drift
    assert 0.85 < nl["eta"] < 0.9


def test_refA_fiber_acoustic_table_equals_lab_acoustic():
    lab = load_yaml(ROOT / "presets/axes/lab/acoustic.yaml")
    ours = load_yaml(ROOT / "presets/axes/bench/refA_fiber_acoustic.yaml")
    t_lab = lab["pipeline"][0]["params"]["psd_table"]
    t_ours = ours["pipeline"][0]["seq"][0]["params"]["psd_table"]
    assert [r[0] for r in t_lab] == [r[0] for r in t_ours]
    assert [r[1].source for r in t_lab] == [r[1].source for r in t_ours]
    spec = load(REFA, ("bench/refA_fiber_acoustic",), L_sq_fiber_m=10.0, lo_route="fiber_matched", eta_fiber_out=0.98,
                acoustic_scale=1.0, acoustic_gamma=1.0)
    mech = find(spec.pipeline, "mech").params
    assert mech["mech_shape"] == "table" and mech["segments"] == ["PR", "LO", "SQ"]
    assert [r[1] for r in mech["psd_table"][1:5]] == pytest.approx([2e-5, 3.6e-6, 6.7e-7, 9e-8])
    go(spec)


def test_refA_laser_dx1_overlay():
    spec = load("presets/reference/refA.yaml", ("bench/refA_laser_dx1",))
    assert find(spec.pipeline, "koheras").params["h0_hz2_per_hz"] == 2.546e4
    res = go(spec)
    # π·h0/Δν_FC alone is 0.024 rad² at the DX1 linewidth: the LO filter-cavity term grows ×400 (refA header)
    assert res.metrics.noise_levels["i.port0"]["sigma_total_rad"] > 0.1


def test_refA_dsc1_overlay_rates_checks_and_run():
    spec = load(REFA, ("bench/dsc1_refA",), L_sq_fiber_m=10.0, lo_route="fiber_matched", eta_fiber_out=0.98)
    f_sim = spec.run.f_sim_hz
    for p, d, dl in (("pid1", "dac1", "dl1"), ("pid2", "dac2", "dl2")):
        assert find(spec.pipeline, p).rate_div == 40 and find(spec.pipeline, d).rate_div == 40
        assert find(spec.pipeline, p).params["discretisation"] == "dsc1"
        assert find(spec.pipeline, d).params["bits"] == 12
        assert find(spec.pipeline, dl).params["tau_s"] == pytest.approx(8.3e-6 - 1 / f_sim - 39 / (2 * f_sim))
    for name in ("chk.lock1_stability", "chk.lock2_stability"):
        r = chk(spec, name)
        assert r["stable"] and 0.5 < r["max_l_at_180"] < 0.8 and 2.5e4 < r["f_180_hz"] < 4e4   # ≈ 4 dB margin at 31–33 kHz
    go(spec)


def test_refA_moku_overlay_checks_and_run():
    spec = load(REFA, ("bench/moku_delta_refA",), L_sq_fiber_m=10.0, lo_route="free_space", eta_fiber_out=0.98)
    assert find(spec.pipeline, "dac2").params["bits"] == 14
    assert find(spec.pipeline, "dl1").params["tau_s"] == pytest.approx(2.1e-7 - 1 / spec.run.f_sim_hz)
    for name in ("chk.lock1_stability", "chk.lock2_stability"):
        r = chk(spec, name)
        assert r["stable"] and r["max_l_at_180"] < 0.05
    go(spec)


# ------------------------------------------------------------------------------------------------ (b) Ref C + bench fibre

REFC = "presets/bench/refC_fiber_bench.yaml"


@pytest.mark.parametrize("overlays", [(), ("bench/dsc1_refC",), ("bench/moku_delta_refC",)], ids=["fpga", "dsc1", "moku"])
def test_refC_fiber_bench(overlays):
    L = 10.0
    spec = load(REFC, overlays, duration=0.03, L_fiber_m=L)
    p = spec.offline["params"]
    topo = as_topology(spec.topology)
    assert topo.segments["SQF_spool"].length_m == pytest.approx(L) and topo.segments["SQF_dep"].length_m == 0.0
    assert p["eta_fiber_link"] == pytest.approx(0.9294118 / 10 ** (-0.02) * 10 ** (-0.2e-3 * L / 10))
    r = chk(spec, "chk.lockC_stability")
    assert r["stable"]
    if overlays == ("bench/dsc1_refC",):
        assert find(spec.pipeline, "pid1").rate_div == 8 and find(spec.pipeline, "dac3").params["range_v"] == [-4.0, 4.0]
        c = next(b for b in spec.offline["blocks"] if b["name"] == "chk.lockC_stability")["params"]
        assert c["discretisation"] == "dsc1" and all(s["discretisation"] == "dsc1" for s in c["stages"])
    if overlays == ("bench/moku_delta_refC",):
        assert p["latency_s"] == pytest.approx(2.5e-7) and p["dac_bits"] == 14
    res = go(spec)
    assert res.metrics.metric_i[0].stats.std < 0.02


# ------------------------------------------------------------------------------------------------ (c) pulsed WG, 1 port


@pytest.mark.parametrize("overlays", [(), ("controller/dsc1",), ("controller/moku_delta",)], ids=["fpga", "dsc1", "moku"])
def test_pulsed_wg_1port(overlays):
    L = 10.0
    spec = load("presets/bench/pulsed_wg_1port.yaml", overlays, drop_esa=False, L_port_m=L)
    p = spec.offline["params"]
    assert p["L_k_m"] == pytest.approx(L)
    assert as_topology(spec.topology).segments["F1"].length_m == pytest.approx(L)
    d1 = find(spec.pipeline, "delay1").params
    assert d1["port"] == 0 and d1["L0_m"] == pytest.approx(L)
    assert find(spec.pipeline, "dmxsw").params["toggles"]["enabled"] is False
    assert chk(spec, "chk.lock_stability")["stable"]
    res = go(spec)
    assert set(res.metrics.metric_i) == {0, 1} and 1 in res.metrics.metric_ii


# ------------------------------------------------------------------------------------------------ (d) Ref D lab


def test_refD_lab_point():
    spec = load("presets/bench/refD_lab.yaml", duration=0.02, drop_esa=False)
    assert config.sweep_markers(spec) == {}
    r = chk(spec, "chk.lockD_stability")
    assert r["stable"] and r["max_l_at_180"] == pytest.approx(10 ** (-3.3 / 20), rel=0.03)   # D89-4: 3.3 dB at 5.41 kHz
    assert r["f_180_hz"] == pytest.approx(5.41e3, rel=0.01)
    st = find(spec.pipeline, "stretcher915b").params["modes"]
    assert [m["r"] for m in st] == pytest.approx([0.67, 0.11, 0.07, 0.15])
    assert find(spec.pipeline, "env.mech").params["mech_shape"] == "table"
    assert config.load(ROOT / "presets/reference/refD_nominal.yaml").offline["params"]["f_res_hz"] == 2.0e4   # unchanged
    go(spec)


# ------------------------------------------------------------------------------------------------ probe overlays


def test_actuator_pathlength_runs_with_fringe_hop():
    spec = load("presets/future/F_liu20_nominal.yaml",
                ("probe/pol_same", "probe/freq_1THz", "probe/actuator_pathlength", "largerange/fringe_hop"),
                drop_esa=False, bias_model="R_prime", R_s_per_m_K=3.0e-16, R_prime_s_per_m_K=5.6e-16)
    p = spec.offline["params"]
    assert p["R_bias_s_per_m_K"] == 5.6e-16 and p["pl_g"] == pytest.approx(1 + 1.007 * 1e12 * p["lambda0_m"] / C)
    assert p["hop_bias_rad"] == pytest.approx(-2 * math.pi * 1.007 * 1e12 * p["lambda0_m"] / C, rel=0.01)
    go(spec)


def test_bad_kT_kappa_independent_axis():
    spec = config.load(ROOT / "presets/future/F_liu20_nominal.yaml",
                       overlays=["probe/freq_1THz", "probe/actuator_pathlength", "probe/bad_kT_kappa_independent"])
    assert "offline.params.kappa_tau_s_per_m_K" in config.sweep_markers(spec)
    pts = expand(spec, {"offline.params.R_s_per_m_K": [3e-16], "offline.params.R_prime_s_per_m_K": [5.6e-16]})
    Rs = sorted({pt.spec.offline["params"]["R_bias_s_per_m_K"] for pt in pts})
    assert len(Rs) == 2 and Rs[0] < Rs[1] and Rs[1] > 1.5e-15   # outside the glass-data range (D86)


@pytest.mark.parametrize("r,detune", [(0.1, 0.3), (0.2, 0.785)])
def test_two_tone_bias_matches_formula(r, detune):
    spec = load("presets/future/F_liu20_nominal.yaml", ("probe/two_tone",), drop_esa=False,
                tt_r=r, tt_dL0_detune_rad=detune, tt_on_grid=True, tt_dnu_target_hz=7.0e9)
    p = spec.offline["params"]
    assert p["tt_dnu_hz"] == 7.0e9 and p["tt_m"] % 2 == 0
    b = -math.atan(r * math.tan(detune))
    assert p["tt_bias_rad"] == pytest.approx(b)
    res = go(spec)
    e0 = res.series["lock.L0.err"]
    assert np.mean(e0[len(e0) // 2:]) == pytest.approx(b, abs=2e-3)       # short port: within 2 mrad of −atan(r tan E)
    e1 = res.series["lock.L1.err"]
    assert np.mean(e1[len(e1) // 2:]) == pytest.approx(b, abs=0.03)       # port 31: same E on the grid (laser self-beat ≈ 0.09 rad)


def test_two_tone_grid_condition_off_grid():
    """Off the grid, port k sees E_k = E_0 + 2πδν·port_k(N_b + n_guard)/f_rep (F_wg6: 7 GHz is not a multiple of 3 MHz)."""
    on = load("presets/future/F_wg6_nominal.yaml", ("probe/two_tone",), drop_esa=False, tt_r=0.0, tt_dL0_detune_rad=0.0,
              tt_on_grid=True, tt_dnu_target_hz=7.0e9).offline["params"]
    off = load("presets/future/F_wg6_nominal.yaml", ("probe/two_tone",), drop_esa=False, tt_r=0.0, tt_dL0_detune_rad=0.0,
               tt_on_grid=False, tt_dnu_target_hz=7.0e9).offline["params"]
    for p, want_cos in ((on, 1.0), (off, math.cos(2 * math.pi * 7e9 * 31 / 6e6))):
        dk = 2 * math.pi * p["tt_dnu_hz"] * p["port_k"] * (p["N_b"] + p["n_guard"]) / p["f_rep_hz"]
        assert math.cos(dk) == pytest.approx(want_cos, abs=1e-6)
    assert on["tt_dnu_hz"] == pytest.approx(7.002e9) and abs(math.cos(2 * math.pi * 7e9 * 31 / 6e6)) < 0.6


def test_detect_heterodyne_refshift_keeps_probe_offset():
    for ovs in (("probe/freq_1THz", "probe/detect_heterodyne_refshift"), ("probe/detect_heterodyne_refshift", "probe/freq_1THz")):
        spec = config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=list(ovs))
        beams = as_topology(spec.topology).beams
        assert spec.offline["params"]["probe_freq_offset_hz"] == 1e12
        assert beams["pr1"].freq_offset_hz == 1e12 and beams["ref1"].freq_offset_hz == 1e12 - 8e7
        assert find(spec.pipeline, "liu1").params["mode"] == "heterodyne"
    # the issue it fixes: probe/detect_heterodyne after probe/freq_1THz turns the probe into an 80 MHz probe
    old = config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=["probe/freq_1THz", "probe/detect_heterodyne"])
    assert old.offline["params"]["probe_freq_offset_hz"] == 8e7
    go(load("presets/future/F_liu20_nominal.yaml", ovs, drop_esa=False))
