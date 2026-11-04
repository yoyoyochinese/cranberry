"""``gated_loop_stability_check`` (WP6-G) against the time loop: Floquet multiplier at stable / marginal / unstable points,
the cyclostationary residual under a coloured disturbance, the ungated limit, CHK events."""

import math
import warnings

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401  (toy_noise_ar1)
import cranberry.path.testing as pt
from cranberry.blocks.chk.gated_loop_stability import make_gated_loop_stability_check
from cranberry.blocks.chk.loop_stability import make_loop_stability_check
from cranberry.core.compile import run_compiled
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy
from cranberry.runner.run import INVALIDATING_CHK

pytestmark = pytest.mark.numba
FS = 1.0e6
WU = 2 * math.pi * 1e3          # integrator loop gain [rad/s] (k0 = 1)
DELAY = 2e-6                    # total loop delay τ_s + 1/f_sim
T_ON = 65.5e-6                  # lock window (66 samples, the Borealis 65 µs class)


def loop_pipe(period, lf, mode, ki=WU, source=None):
    """d → e = d − y → gated_hold → pid (tustin, I only) → delay_buffer → y (the plant is a unit gain)."""
    src = source or {"block": "ramp", "name": "d", "params": {"out": "v.d", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                               "period_s": 1.0, "offset_v": 1.0,
                                                               "toggles": {"enabled": True}}}
    return [
        src,
        {"block": "gate_schedule", "name": "gw",
         "params": {"detector_gate": "complement", "write_phase": False, "period_s": period, "lock_fraction": lf,
                    "phase": 0.0, "edge_s": 0.0, "timing_jitter_s": 0.0, "leak": 0.0,
                    "toggles": {"gating": True, "jitter": False, "edges": False, "leak": False}}},
        {"block": "adder", "name": "plant", "params": {"ins": ["v.d", "delay.dl.out"], "gains": [1.0, -1.0],
                                                       "out": "v.e", "unit": "V"}},
        {"block": "gated_hold", "name": "gh", "params": {"in": "v.e", "gate": "gate.gw.lock", "out": "v.c", "unit": "V",
                                                         "mode": "hold" if mode == "hold" else "zero", "settle_s": 0.0,
                                                         "toggles": {"gating": True}}},
        {"block": "pid", "name": "pid",
         "params": {"in": "v.c", "out": "ctrl.pid.u", "device": None, "discretisation": "tustin", "anti_windup": "none",
                    "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
                    "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1e12, "u_max_v": 1e12,
                    "toggles": {"feedback": True, "p_term": False, "i_term": True, "d_term": False}}},
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out", "method": "linear",
                                                           "tau_s": DELAY - 1 / FS, "unit": "V"}},
    ]


def check_params(period, lf, mode, ki=WU, disturbance=None, **kw):
    p = {"k0": 1.0, "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "delay_s": DELAY, "poles_hz": [], "modes": [],
         "f_grid_hz": [1.0, 1.0e7, 20001], "pid_discretisation": "tustin", "f_sim_hz": FS, "period_s": period,
         "lock_fraction": lf, "gate_delay_s": 0.0, "settle_s": 0.0, "hold_mode": mode, "off_gain": 0.0, "meas_delay_s": 0.0,
         "meas_duty": None, "disturbance": disturbance, "round_trip": None, "toggles": {"enabled": True}}
    p.update(kw)
    return p


def check(*a, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_gated_loop_stability_check(check_params(*a, **kw), name="chk.gated")


def per_period_ratio(period, lf, mode, n_per=14):
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(loop_pipe(period, lf, mode)), lay, 1)
    st, pr = State(lay), Params(lay)
    n = int(round(period * FS))
    e = []
    for _ in range(n_per):
        run_compiled(blk, st, pr, n)
        e.append(float(st["v.e"]))
    return abs(e[-1] / e[-2])


@pytest.mark.parametrize("t_off, verdict", [(200e-6, "stable"), (395.5e-6, "marginal"), (600e-6, "unstable")])
def test_floquet_multiplier_matches_time_loop(t_off, verdict):
    """hold mode (the held error keeps the integrator ramping): λ ≈ e^{−ω_u T_on}(1 − ω_u T_off) → unstable for long
    holds. ρ of the check is the per-period growth of the time loop's error (dominant multiplier)."""
    period = T_ON + t_off
    blk = check(period, T_ON / period, "hold")
    rho = blk.result["rho"]
    ratio = per_period_ratio(period, T_ON / period, "hold")
    lam = math.exp(-WU * 66e-6) * (1 - WU * (round(period * FS) - 66) * 1e-6)   # continuous, no delay
    assert rho == pytest.approx(abs(lam), abs=0.035)
    if verdict == "stable":
        assert rho < 0.25 and ratio < 0.25 and blk.result["stable"]
        assert [e.code for e in blk.chk_log.events] == ["gate_period_not_integer"]  # 265.5 samples (G-7)
    else:
        assert ratio == pytest.approx(rho, rel=1e-6)
        assert blk.result["stable"] == (verdict == "marginal" and rho < 1.0)
    if verdict == "marginal":
        assert abs(rho - 1.0) < 0.01
    if verdict == "unstable":
        assert [e.code for e in blk.chk_log.events if e.code != "gate_period_not_integer"] == ["loop_unstable"] and "loop_unstable" in INVALIDATING_CHK


def test_zero_mode_holds_and_ungated_limit():
    """zero mode with an I-only controller is an exact hold: stable for any hold length here; lock_fraction 1 gives the
    ungated discrete loop (ρ = ρ_sample^N) and the continuous Nyquist verdict of ``loop_stability_check``."""
    r = check(T_ON + 5e-3, T_ON / (T_ON + 5e-3), "zero").result
    assert r["stable"] and r["rho"] < math.exp(-WU * 60e-6)
    one = check(100e-6, 1.0, "zero").result
    assert one["rho"] == pytest.approx(one["rho_ungated_per_sample"] ** 100, rel=1e-6)
    lsc = make_loop_stability_check({k: v for k, v in check_params(100e-6, 1.0, "zero").items()
                                     if k in ("k0", "kp", "ki_per_s", "kd_s", "delay_s", "poles_hz", "modes",
                                              "f_grid_hz", "toggles")}, name="lsc").result
    assert one["ungated_stable"] == lsc["stable"] is True
    bad = check(100e-6, 1.0, "zero", ki=2 * math.pi * 150e3).result   # ω_u τ > π/2
    assert not bad["stable"] and not bad["ungated_stable"]


@pytest.mark.parametrize("mode", ["zero", "hold"])
def test_residual_matches_time_loop(mode):
    """AR(1) disturbance (corner ≈ 80 Hz, in-band): the exact cyclostationary variance over the measurement / lock
    window equals the time loop's to the statistical error (~5 %)."""
    a = 0.9995
    period, lf = 300e-6, 100.5e-6 / 300e-6
    n, n_per = 300, 3000
    src = toy("toy_noise_ar1", "d", out="v.d", unit="V", a=a, b=1.0)
    pipe = loop_pipe(period, lf, mode, ki=2 * math.pi * 3e3, source=src) + [pt.record(["v.e", "gate.gw.meas"], n * n_per)]
    rec, _, _ = pt.run_and_record(pipe, FS, n * n_per, compiled=True, seed=11)
    e = rec["v.e"].reshape(n_per, n)[100:]  # drop the start-up
    meas = rec["gate.gw.meas"][:n] > 0.5
    f = np.geomspace(1.0, 0.4999 * FS, 3000)
    s = 2 / FS / np.abs(1 - a * np.exp(-2j * np.pi * f / FS)) ** 2      # one-sided PSD of the AR(1) sequence
    r = check(period, lf, mode, ki=2 * math.pi * 3e3, disturbance={"f_hz": f.tolist(), "psd_rad2_per_hz": s.tolist()}).result
    assert r["var_meas_rad2"] == pytest.approx(e[:, meas].var(), rel=0.08)
    assert r["var_lock_rad2"] == pytest.approx(e[:, :101].var(), rel=0.08)
    assert r["var_period_rad2"] == pytest.approx(e.var(), rel=0.08)
    assert r["var_meas_rad2"] > r["var_lock_rad2"] > r["var_ungated_rad2"]   # the hold costs; gating costs
    assert 0.0 < r["var_conversion_rad2"] < r["var_period_rad2"]
    assert r["var_hold_drift_rad2"] > 0.0


def test_conversion_vanishes_without_gating_and_round_trip_warning():
    f = np.geomspace(1.0, 4e5, 400)
    dist = {"f_hz": f.tolist(), "psd_rad2_per_hz": (1e-6 / f**2).tolist()}
    r = check(100e-6, 1.0, "zero", disturbance=dist).result
    assert abs(r["var_conversion_rad2"]) < 1e-9 * r["var_period_rad2"]
    assert r["var_period_rad2"] == pytest.approx(r["var_ungated_rad2"], rel=0.05)  # discrete ≈ continuous loop
    # P3 numbers: 6 kHz chopper, D_pr 0.5, D_det 0.25, 2.1 km → f_max ≈ 6.1 kHz (lit-nardelli2026 §3.1-3)
    rt = {"n_g": 1.468, "L_m": 2100.0, "t_edge_s": 0.0, "duty_det": 0.25}
    ok = check(1 / 6e3, 0.5, "zero", round_trip=rt)
    assert ok.result["round_trip_f_max_hz"] == pytest.approx(0.25 / (4 * 1.468 * 2100.0 / 299792458.0), rel=1e-12)
    assert not [e for e in ok.chk_log.events if e.code == "gate_round_trip"]
    fast = check(1 / 10e3, 0.5, "zero", round_trip=rt)
    assert [e.code for e in fast.chk_log.events if e.code != "gate_period_not_integer"] == ["gate_round_trip"]


def test_validation_and_toggle():
    with pytest.raises(ValueError, match="hold_mode"):
        check(100e-6, 0.5, "sample")
    with pytest.raises(ValueError, match="delay_s"):
        check(100e-6, 0.5, "zero", delay_s=0.5e-6)
    with pytest.raises(KeyError):
        make_gated_loop_stability_check({k: v for k, v in check_params(1e-4, 0.5, "zero").items() if k != "off_gain"})
    off = make_gated_loop_stability_check(check_params(1e-4, 0.5, "zero", toggles={"enabled": False}))
    assert off.result is None and not off.chk_log.events


# ---- Phase 6 fixes (reviews/phase6_wp6g_ab.md G-2/G-4/G-7, D91-1) ------------------------------------------------


def per_period_ratio_rd(period, lf, mode, rate_div, n_per=14):
    pipe = loop_pipe(period, lf, mode)
    pipe[4]["rate_div"] = rate_div  # the pid (D91: DSC1 at 500 kHz = rate_div 2 at 1 MHz)
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(pipe), lay, 1)
    st, pr = State(lay), Params(lay)
    n = int(round(period * FS))
    e = []
    for _ in range(n_per):
        run_compiled(blk, st, pr, n)
        e.append(float(st["v.e"]))
    return abs(e[-1] / e[-2])


@pytest.mark.parametrize("t_off", [600e-6, 801e-6])
def test_pid_rate_div_matches_time_loop(t_off):
    """D91-1: a PID updating every 2nd sample (h = 2 µs, output held) — the check's per-period multiplier equals the
    time loop's per-period growth, also for an odd period (monodromy over lcm(N, 2))."""
    period = T_ON + t_off
    blk = check(period, T_ON / period, "hold", pid_rate_div=2)
    ratio = per_period_ratio_rd(period, T_ON / period, "hold", 2, n_per=16)
    assert not blk.result["stable"]
    assert ratio == pytest.approx(blk.result["rho"], rel=1e-5)
    r1 = check(period, T_ON / period, "hold").result["rho"]
    assert blk.result["rho"] != pytest.approx(r1, rel=1e-6)  # the hold is modelled, not ignored


def test_modes_with_r_are_a_parallel_sum():
    """G-4: every mode with r → Σ r ω²/(...) (loop_stability rule, D89-1); a single r = 1 mode equals the series form."""
    m1 = [{"f_hz": 2.0e4, "q": 10.0}]
    one = check(100e-6, 1.0, "zero", modes=[{**m1[0], "r": 1.0}]).result
    ser = check(100e-6, 1.0, "zero", modes=m1).result
    assert one["rho"] == pytest.approx(ser["rho"], rel=1e-9)
    two = [{"f_hz": 2.0e4, "q": 10.0, "r": 0.8}, {"f_hz": 5.0e3, "q": 20.0, "r": 0.2}]
    par = check(100e-6, 1.0, "zero", modes=two).result
    prod = check(100e-6, 1.0, "zero", modes=[{k: v for k, v in m.items() if k != "r"} for m in two]).result
    lsc = make_loop_stability_check({k: v for k, v in check_params(100e-6, 1.0, "zero", modes=two).items()
                                     if k in ("k0", "kp", "ki_per_s", "kd_s", "delay_s", "poles_hz", "modes",
                                              "f_grid_hz", "toggles")}, name="lsc").result
    assert par["rho"] != pytest.approx(prod["rho"], rel=1e-3)
    assert par["ungated_stable"] == lsc["stable"]


def test_window_shorter_than_a_sample_warns():
    """G-2: T = 10 µs, D_lock 0.9 at 1 MHz → the measurement window (1 µs − 0.1 µs edge) holds no sample."""
    blk = check(10e-6, 0.9, "zero", meas_delay_s=1e-7)
    codes = [e.code for e in blk.chk_log.events]
    assert "gate_window_short" in codes and blk.result["n_meas"] == 0
    assert "gate_window_short" not in INVALIDATING_CHK  # the metrics' metric_no_sample invalidates the point
