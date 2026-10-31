"""VT-15 (block part, WP-B4): integrator-loop residual σ² = π S_0 / (2 f_c) within 5 %
(prompt §7-7; notes/00_summary.md §4 VT-15; T15 §7 test 7, T15 §2.10).

Loop: random-walk disturbance φ_n (one-sided S_φ = S_0/f², S_0 = D/π², Var = 2Dt) → e = φ_n − φ_a →
``pid`` (integral only, DSC1 discretisation) → ``delay_buffer`` → static actuator (K = 1 rad/V).
L ≈ f_c/(i f) with f_c = K·K_i/(2π); residual variance π S_0/(2 f_c) = D/(2π f_c). The discrete loop
gives σ²/(2g − g²) with g = 2π f_c/f_sim (0.6 % above the continuous value here). Sample size: the
3σ error of the variance estimate of the AR(1) residual (ρ = 1 − g), 3·√(2(1 + ρ²)/((1 − ρ²) N)), is
computed below and kept under 2.5 %.

The lock-2 reduction part of VT-15 (residual of Δ scaled by A−/A+, Ref A) needs the OPA / error-signal
blocks and is WP-I1's.
"""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.pipeline import spec_from_obj
from cranberry.core.testing import toy
from cranberry.runner.run import run

pytestmark = pytest.mark.numba
FS, F_C, SIGMA = 1.0e6, 2.0e3, 1e-3


def loop_obj(ki):
    return [
        toy("toy_noise_ar1", "phin", out="v.phin", unit="rad", a=1.0, b=SIGMA),
        {"block": "adder", "name": "err", "params": {"ins": ["v.phin", "act.a.dpsi"], "gains": [1.0, -1.0],
                                                     "out": "v.e", "unit": "rad"}},
        {"block": "pid", "name": "pid", "params": {
            "in": "v.e", "out": "ctrl.pid.u", "device": None, "discretisation": "dsc1", "anti_windup": "none",
            "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
            "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1e9, "u_max_v": 1e9,
            "toggles": {"feedback": True, "p_term": False, "i_term": True, "d_term": False}}},
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out", "tau_s": 0.0,
                                                           "method": "linear"}},
        {"block": "pzt_fiber_stretcher", "name": "a", "params": {
            "in": "delay.dl.out", "gain_mode": "direct", "k_rad_per_v": 1.0, "polarity": 1.0, "modes": [],
            "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -1e9, "drive_max_v": 1e9,
            "hysteresis_frac": 0.0, "segment": None, "insertion_loss_db": None,
            "toggles": {"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": False}}},
    ]


def test_vt15_integrator_loop_residual():
    ki = 2 * math.pi * F_C  # K = 1 rad/V
    g = 2 * math.pi * F_C / FS
    n_keep = 3_000_000
    rho = 1.0 - g
    rel_3sigma = 3 * math.sqrt(2.0 * (1 + rho**2) / ((1 - rho**2) * n_keep))
    assert rel_3sigma < 0.025
    spec = RunSpec(run=RunSettings(duration_s=(n_keep + 20_000) / FS, f_sim_hz=FS, chunk_samples=1 << 18, seed=15),
                   pipeline=spec_from_obj(loop_obj(ki)), measure={"record": {"slots": ["v.e"], "decim": 1}})
    e = run(spec, keep_state=False).series["v.e"][20_000:]
    d = SIGMA**2 * FS / 2  # Var(φ_n) = σ² n = 2 D t
    s0 = d / math.pi**2
    target = math.pi * s0 / (2 * F_C)
    assert np.var(e) == pytest.approx(target, rel=0.05)
    assert np.var(e) == pytest.approx(SIGMA**2 / (2 * g - g * g), rel=rel_3sigma)
