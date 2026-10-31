"""VT-21 (Ref D thermal-ramp part, WP-I2): the lock is lost when the thermal drift of ΔL exceeds the
915B range driven by the DSC1 (prompt §7-10 second half; T15 §2.11, §7 test 8; T08 §2.4;
notes/00_summary.md §4 VT-21). The block-level ramp / random-walk parts are test_vt_21_lock_duration.py.

``refD_nominal`` with every noise source and quantiser off except the ambient ramp of the thermal
zone (both arms; T08 eq. 2.6 lag τ_th): r_T = 2e−2 K/s and τ_th = 0.5 s (fast values so the test runs
in seconds; the physics is the same at 1e−4 K/s, only slower). With the loop holding φ_D at the lock
point, the DSC1 output follows

    k_act·g·u(t) = φ_0 − π/2 − k_T·ΔL·ΔT(t),   ΔT(t) = r_T·(t − τ_th(1 − e^{−t/τ_th}))

(φ_0: φ_D with the 915B at rest; g: DAC load divider), so the CV limit −4 V is reached when
k_T·ΔL·ΔT = 4·k_act·g + φ_0 − π/2 (Ref D's own R/r: R = 4 V × 6.49 rad/V ≈ 26 rad). The supervisor's
loss time must agree within 1 % (T15 test 8), and after the loss the phase slides over the fringe
(fringe slips counted on the true phase).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parents[1] / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

R_T, TAU_TH = 2.0e-2, 0.5


def test_vt21_refD_lock_lost_when_thermal_drift_exceeds_range():
    tog = {k: v for k, v in rs.NOISE_OFF.items() if k != "env.thermal.drift"}
    spec = rs.load(params={"r_T_k_per_s": R_T, "tau_th_s": TAU_TH}, run={"duration_s": 4.8},
                   measure={"record": {"slots": ["lock.lockD.err", "v.dsc1.cv"], "decim": 100}}, toggles=tog)
    res = run(spec)
    st = res.final_state
    pr = spec.offline["params"]
    a1, a2 = rs.segment(spec, "a1"), rs.segment(spec, "a2")
    dl = a1["L_m"] - a2["L_m"]
    leaves, params = rs.built_leaves(spec)
    k_act = leaves["stretcher915b"].gain_rad_per_v(params)
    g = leaves["dsc1.dac"].gain()
    assert k_act == pytest.approx(6.491, abs=5e-4)  # T16 V1
    phi0 = float(res.series["lock.lockD.err"][0]) + math.pi / 2  # first record: 915B still at rest
    reach = 4.0 * k_act * g + phi0 - math.pi / 2
    dT = lambda t: R_T * (t - TAU_TH * (1 - math.exp(-t / TAU_TH)))  # noqa: E731
    t_exp = brentq(lambda t: pr["k_T_rad_per_m_K"] * dl * dT(t) - reach, 0.1, 100.0)
    t_loss = float(st["lock.dsc1.sup.t_loss"])
    assert int(st["lock.dsc1.sup.n_loss"]) == 1
    assert t_loss == pytest.approx(t_exp, rel=0.01), (t_loss, t_exp)
    # before the loss the lock held (true phase error ≪ 1 rad), after it the phase slid over the fringe
    t, err = res.t_s, res.series["lock.lockD.err"]
    held = (t > 0.05) & (t < 0.99 * t_loss)
    assert np.max(np.abs(err[held])) < 0.05
    assert int(st["lock.dsc1.sup.n_slip"]) >= 1
    assert float(res.series["v.dsc1.cv"][-1]) == pytest.approx(-4.0, abs=1e-12)  # parked at the CV limit
