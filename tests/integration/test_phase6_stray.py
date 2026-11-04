"""Phase 6 (D85-4, D85-11) stray-light overlays on future nominal presets (F_liu20, F_wg6): ``ref/aurora_calframe``,
``ref/aurora_guardslot`` (forward_ref_stray_table) and ``stray/window_tail`` (window_stray_light) build and run through
the full runner; the runner calls their ``check()`` (events in ``RunResult.events``); the guard slot doubles / stretches
L_k and reports the D78 linewidth consequence; the window tail equals the analytic 2τ_k box."""

import importlib.util
from pathlib import Path

import pytest
import yaml

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.integration]

REF_POINT = dict(fref_p_peak_w=1.0e-6, fref_t_pulse_s=1.0e-9, fref_eps_mod=1.0e-5, fref_edge_tau_s=1.0e-9,
                 det_p_afterpulse=1.0e-3, det_tau_afterpulse_s=1.0e-7, det_recovery_s=4.3e-8)
WIN_POINT = dict(win_period_s=1.0e-4, win_duty_lock=0.65, win_edge_s=1.0e-7, win_eps_off=1.0e-5,
                 win_limit_statistic="mean")


def point(tmp, name, values):
    """A point overlay applied AFTER the axis overlay (a file's own keys lose to its overlays, D74-7)."""
    path = tmp / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"offline": {"params": values}}))
    return str(path)


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p6stray")
    pre, ref = fs.nominal("F_liu20"), point(tmp, "ref", REF_POINT)
    return {"cal": fs.build_check(pre, ("ref/aurora_calframe", ref)),
            "guard": fs.build_check(pre, ("ref/aurora_guardslot", ref)),
            # F_wg6: stray1 weights by the τ_mode Gaussian there, as the overlay does (the OPO presets weight stray1 by
            # the OPO mode; the overlay keeps the Gaussian stand-in, see window_stray module docstring)
            "win": fs.build_check(fs.nominal("F_wg6"), ("stray/window_tail", point(tmp, "win", WIN_POINT))),
            "base": fs.build_check(pre, ())}


def test_forward_reference_placements(runs):
    cal, guard = fs.block(runs["cal"], "frefstray1").result, fs.block(runs["guard"], "frefstray1").result
    assert cal["lk_factor"] == 1.0 and guard["lk_factor"] == pytest.approx(5 / 4)   # F_liu20: N_b = 4
    assert "forward_ref_guard_slot" in fs.chk_codes(runs["guard"])
    assert "forward_ref_guard_slot" not in fs.chk_codes(runs["cal"])
    # the guard-slot pulse sits one slot from the data every frame; calibration frames only next to t_cal_s frames
    assert guard["edges_mode_mean"] > 1e6 * cal["edges_mode_mean"]
    assert guard["photons_per_mode_mean"] > cal["photons_per_mode_mean"]
    # the guard overlay really lengthens the port fibre of the run (D62-3)
    p_guard, p_base = fs.params(runs["guard"]), fs.params(runs["base"])
    L0 = p_base["L0_m"]
    assert (p_guard["L_k_m"] - L0) == pytest.approx(1.25 * (p_base["L_k_m"] - L0))


def test_window_tail_on_nominal(runs):
    p = fs.params(runs["win"])
    r = fs.block(runs["win"], "winstray1").result
    w = 2 * p["n_g"] * p["L_k_m"] / fs.C_LIGHT
    assert r["tail_s"] == pytest.approx(w)
    t_meas, dt, e = 0.35e-4 - 1e-7, p["t_frame_s"], 1e-7
    # GBS slot centres once per frame from the end of the ramp; h = 1 − τ/W (τ from the 50 % switch-off point)
    tau = [e / 2 + (j + 0.5) * dt for j in range(int(t_meas // dt))]
    h = sum(max(0.0, 1 - x / w) for x in tau) / len(tau)
    assert r["n_slots"] == len(tau)
    assert r["mean_factor"] == pytest.approx(1e-5 + (1 - 1e-5) * h, rel=1e-6)
    # the window counts are the stray1 (stray_light_table) numbers at the run's P̄ times the window factor
    assert r["photons_per_mode_mean"] == pytest.approx(r["mean_factor"] * r["photons_per_mode_on"], rel=1e-9)
    assert r["photons_per_mode_on"] == pytest.approx(fs.stray_light(runs["win"])["photons_per_mode"], rel=1e-9)
