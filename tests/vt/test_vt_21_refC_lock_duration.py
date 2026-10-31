"""VT-21 (Ref C part, WP-I3): the lock duration is limited by the range of the slow actuator when the drift
exceeds it (prompt §7-10 second half; [Ref C §4] "stability for several minutes, limited by the range of the
fiber-based phase shifter"; 00_summary §4 VT-21; T16 §3.6, §4.4, V12; T08 §2.4).

``refC_nominal`` (analytic dither, f_sim = 1 MHz for speed: the quantities are slow), every noise source off
except the differential temperature drift of the two fibre spools (T08 eq. 2.6 lag τ_th): fast values
dr_T = 4.4e−4 K/s, τ_th = 0.1 s so the test runs in seconds (the physics is the same at the nominal
5e−6 K/s, only slower). The differential fibre phase d(t) = k_T·L_spool·dr_T·(t − τ(1 − e^{−t/τ})) (rate
ρ = k_T L dr_T ≈ 15.9 rad/s) is followed by the stack: in the steady ramp the type-3 loop holds Δ ≈ 0, the
EOM stage returns to 0, the mirror holds the constant ρ/(2π f_23) (its integrator feeds the stretcher stage)
and the stretcher the rest. The supervisor (range criterion on the stretcher stage) therefore reports the
loss when

    d(t_loss) = R3·(1 − 2^(1−bits)) + ρ/(2π f_23)          (R3: stretcher half range)

within 1 % (T16 V12's "± 1 sample" is for an ideal integrator; here the stack's static shares are part of the
prediction). Until then |Δ| stays < 10 mrad; afterwards the mirror (range R2) and the EOM (±π/4 V/V_π) take
the rest and Δ leaves the lock (|Δ| > π/4) at d = R3 + R2 + R1 within 5 %. Halving R3 shortens the lock
accordingly. At the nominal point (32.5 rad, 5e−6 K/s, τ_th = 300 s) the same relation gives ≈ 6.7 min, the
"several minutes" of Ref C (analytic, not simulated).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parents[1] / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

DR, TAU = 4.4e-4, 0.1
TOG = {k: v for k, v in rs.NOISE_OFF.items() if k not in ("env.thermal.spool_sq.drift", "env.thermal.spool_lo.drift")}
REC = {"record": {"slots": ["lock.lockC.err", "ctrl.pid3.u", "act.stretcher.dpsi", "act.mirror.dpsi"], "decim": 10}}


def t_reach(level, rho, tau):
    """Time at which ρ(t − τ(1 − e^{−t/τ})) = level."""
    return brentq(lambda t: rho * (t - tau * (1 - math.exp(-t / tau))) - level, 0.0, 1e6)


def simulate(r3):
    spec = rs.load(params={"dr_T_spool_k_per_s": DR, "tau_th_spool_s": TAU, "stretcher_range_rad": r3},
                   run={"duration_s": 3.2, "f_sim_hz": 1.0e6}, toggles=TOG, measure=REC)
    return spec, run(spec)


@pytest.mark.parametrize("r3", [32.5, 16.25])
def test_vt21_refC_lock_lost_at_stretcher_range(r3):
    spec, res = simulate(r3)
    p = rs.params_of(spec)
    rho = p["k_T_rad_per_m_K"] * (p["L_fiber_m"] - p["L_deployed_m"]) * DR
    share2 = rho / (2 * math.pi * p["f_23_hz"])                          # mirror's static share
    r3_eff = r3 * (1 - 2.0 ** (1 - p["dac_bits"]))
    t_pred = t_reach(r3_eff + share2, rho, TAU)
    fs = res.final_state
    assert int(fs["lock.sup.n_loss"]) == 1
    t_loss = float(fs["lock.sup.t_loss"])
    assert t_loss == pytest.approx(t_pred, rel=0.01)
    # locked until then; the stretcher is at its range at the loss
    t = np.arange(res.series["lock.lockC.err"].size) * 10 / spec.run.f_sim_hz
    err = res.series["lock.lockC.err"]
    delta = np.pi / 2 - np.mod(np.pi / 2 - err, np.pi)                   # quadrature angle, (−π/2, π/2]
    settled = (t > 0.05) & (t < t_loss)
    assert np.max(np.abs(delta[settled])) < 0.01
    assert float(np.interp(t_loss, t, res.series["act.stretcher.dpsi"])) == pytest.approx(r3_eff, rel=1e-3)
    # the faster stages hold it a little longer, then Δ leaves the lock point
    r1 = math.pi / p["v_pi_v"] * p["out_range_v"]
    t_break = t[np.argmax(np.abs(delta) > np.pi / 4)]
    assert t_break == pytest.approx(t_reach(r3 + p["mirror_range_rad"] + r1, rho, TAU), rel=0.05)


def test_vt21_refC_no_loss_without_drift():
    tog = dict(TOG, **{"env.thermal.spool_sq.drift": False, "env.thermal.spool_lo.drift": False})
    res = run(rs.load(run={"duration_s": 1.0, "f_sim_hz": 1.0e6}, toggles=tog))
    assert int(res.final_state["lock.sup.n_loss"]) == 0


def test_vt21_refC_nominal_lock_time_is_minutes():
    p = rs.params_of(rs.load())
    rho = p["k_T_rad_per_m_K"] * (p["L_fiber_m"] - p["L_deployed_m"]) * p["dr_T_spool_k_per_s"]
    t = t_reach(p["stretcher_range_rad"] + rho / (2 * math.pi * p["f_23_hz"]), rho, p["tau_th_spool_s"])
    assert 2 * 60 < t < 10 * 60                                          # "several minutes" [Ref C §4]
