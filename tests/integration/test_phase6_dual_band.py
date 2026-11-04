"""Phase 6 WP6-C (D85-5): the two-layer lock overlay ``lock/dual_band_slow`` (P5) builds on the four future presets at
its first sweep point and its choice axes, and on F_liu20_nominal its slow layer (photon_counting_fringe on in-band
seeds → set point of the L-Liu loop, or a REF_k stretcher) closes with the right sign: it removes the static metric (ii)
offset the fast proxy leaves, windows close every T_int, the slow-loop stability check passes and the common-point
window row is suppressed by the slow loop's sensitivity."""

import importlib.util
import math
import warnings
from pathlib import Path

import pytest
import yaml

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.integration]

OV = "lock/dual_band_slow"
NOM = "lock/points/dual_band_slow_nominal"


def point(tmp, name, values):
    path = tmp / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"offline": {"params": values}}))
    return str(path)


@pytest.mark.parametrize("base", ["F_wg100", "F_wg6", "F_bor6", "F_liu20"])
def test_builds_on_every_future_preset(base):
    res = fs.build_check(fs.nominal(base), (OV,))
    assert set(res.metrics.metric_ii) == {1}
    assert "cnt.pcf1.window_index" in res.final_state.layout
    p = fs.params(res)
    assert p["dbs_cal_slot"] == "cal_frames" and p["dbs_seed_inject"] == "after_squeezer"   # first point = primary
    assert p["dbs_f_gate_hz"] == pytest.approx(p["N_b"] / p["t_frame_s"])


@pytest.mark.parametrize("values", [{"dbs_target": "ref_stretcher"}, {"dbs_seed_inject": "before_squeezer"},
                                    {"dbs_seed_inject": "before_squeezer", "dbs_seed_model": "ideal_copy"},
                                    {"dbs_seed_inject": "before_squeezer", "dbs_seed_prep": "pm_pi4"}])
def test_choice_axes_build(tmp_path, values):
    res = fs.build_check(fs.nominal("F_liu20"), (OV, NOM, point(tmp_path, "pt", values)))
    p = fs.params(res)
    if values.get("dbs_seed_model") == "ideal_copy":
        assert p["dbs_seed_k"] == "beam.dseed1.phi" and p["dbs_psa_r"] == 0.0
    elif values.get("dbs_seed_prep") == "pm_pi4":
        assert p["dbs_phi_offset_rad"] == pytest.approx(math.pi / 2 - 2 * math.atan(math.exp(-3.6)))
        assert p["dbs_p_in_k"] == pytest.approx(p["dbs_p_in_0"])
    elif "dbs_seed_inject" in values:
        assert p["dbs_seed_k"] == "beam.dsin1.phi" and p["dbs_psa_r"] == 1.8 and p["dbs_phi_offset_rad"] == 0.0
        assert p["dbs_eta_coupler"] == 1.0 and p["dbs_sqvac_per_gate"] > 1.5
    else:
        assert p["dbs_g_setpoint"] == 0.0 and p["dbs_ref_k_rad_per_v"] > 0.0


def test_nominal_point_and_cal_slot_options():
    """Nominal (review C-1..C-4): after_squeezer with the tap loss registered, b' counted, calibration frames (L_k
    unchanged); the guard point adds the slot (L_k × 1.25, CHK dbs_cal_guard_slot, SNSPD recovery CHK); the ideal point
    keeps the pre-review assumptions and is flagged (dbs_cal_slot_ideal)."""
    nom = fs.build_check(fs.nominal("F_liu20"), (OV, NOM))
    p = fs.params(nom)
    assert p["dbs_seed_k"] == "beam.dseedu1.phi" and p["dbs_eta_coupler"] == pytest.approx(0.99)
    assert p["dbs_sqvac_per_gate"] == pytest.approx(1.72, abs=0.01)
    assert p["dbs_phi_offset_rad"] == pytest.approx(math.pi / 2)
    base_lk = p["L_k_m"]
    g = fs.build_check(fs.nominal("F_liu20"), (OV, NOM, "lock/points/dbs_cal_guard"))
    assert fs.params(g)["L_k_m"] - p["L0_m"] == pytest.approx(1.25 * (base_lk - p["L0_m"]), rel=1e-9)
    ideal = fs.params(fs.build_check(fs.nominal("F_liu20"), (OV, "lock/points/dual_band_slow_ideal")))
    assert ideal["dbs_seed_k"] == "beam.dseed1.phi" and ideal["dbs_sqvac_per_gate"] == 0.0


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("dbs")
    d = {"duration_s": 0.6}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {"base": run(fs.load(fs.nominal("F_liu20"), run=d)),
                "sp": run(fs.load(fs.nominal("F_liu20"), run=d, overlays=(OV, NOM))),
                "ref": run(fs.load(fs.nominal("F_liu20"), run=d,
                                   overlays=(OV, NOM, point(tmp, "ref", {"dbs_target": "ref_stretcher"}))))}


def test_slow_layer_removes_the_proxy_offset(runs):
    base = runs["base"].metrics.metric_ii[1].stats.mean
    assert abs(base) > 0.02                                    # the fast proxy alone leaves ≈ 35 mrad (nominal seed)
    for k in ("sp", "ref"):
        r = runs[k]
        assert r.valid and not [e for e in r.events if e["source"] == "chk.dbs_stability"]
        assert abs(r.metrics.metric_ii[1].stats.mean) < 0.2 * abs(base)
        st = r.final_state
        assert int(st["cnt.pcf1.window_index"]) == 12          # 0.6 s / 50 ms
        assert float(st["cnt.pcf1.x_hat"]) == pytest.approx(math.pi / 2, abs=0.05)
    p = fs.params(runs["sp"])
    assert p["dbs_cp_window_var_rad2"] < 0.05 * p["cp_window_var_rad2"]   # window row under the slow loop
    # the gates per window follow f_gate·cal_fraction·T_int (N_b gates per calibration frame, 20 % of the frames)
    assert float(runs["sp"].final_state["cnt.pcf1.n_gates"]) == pytest.approx(
        0.2 * 0.05 * p["N_b"] / p["t_frame_s"], abs=1.0)


@pytest.mark.parametrize("fct,unstable", [(0.30, False), (0.33, True)])
def test_slow_loop_check_at_the_exact_discrete_limit(tmp_path, fct, unstable):
    """Review C-6: window average → hold → integrator is stable for f_c·T_int < 1/π = 0.318 (Jury); chk.dbs_stability
    (τ_eq = (π/4)·T_int + τ_L) flags 0.33 and passes 0.30 (the old e^{−s(T + τ_L)} model flagged everything above 0.25)."""
    res = fs.build_check(fs.nominal("F_liu20"), (OV, NOM, point(tmp_path, "fc", {"dbs_fc_x_tint": fct})))
    hit = [e for e in res.events if e["source"] == "chk.dbs_stability" and e["code"] == "loop_unstable"]
    assert bool(hit) == unstable
