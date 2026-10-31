"""VT-15 (block part, WP-B4): integrator-loop residual σ² = π S_0 / (2 f_c) within 5 %
(prompt §7-7; notes/00_summary.md §4 VT-15; T15 §7 test 7, T15 §2.10).

Loop: random-walk disturbance φ_n (one-sided S_φ = S_0/f², S_0 = D/π², Var = 2Dt) → e = φ_n − φ_a →
``pid`` (integral only, DSC1 discretisation) → ``delay_buffer`` → static actuator (K = 1 rad/V).
L ≈ f_c/(i f) with f_c = K·K_i/(2π); residual variance π S_0/(2 f_c) = D/(2π f_c). The discrete loop
gives σ²/(2g − g²) with g = 2π f_c/f_sim (0.6 % above the continuous value here). Sample size: the
3σ error of the variance estimate of the AR(1) residual (ρ = 1 − g), 3·√(2(1 + ρ²)/((1 − ρ²) N)), is
computed below and kept under 2.5 %.

Lock-2 reduction part (WP-I1; T14 §3.4, V3; S §4 VT-15 "lock 1 のセンサ雑音だけ入れたとき Δ の残差が
A−/A+ 倍（10 % 以内）"), through the runner on presets/reference/refA.yaml with every other noise source
and the converter quantisation switched off:

* static: lock 1 held at u′ = δu (demodulation phase φ_d1 = π/2 − 2δu, so Φ1^det = 2δu) moves the
  lock-2 zero by ε with tan ε = (A−/A+) tan(δu + u0) (u0 = 21 mrad non-common probe phase, from the
  baseline offset); δu = 0.1 / 0.3 / 0.6 rad, ±10 % (achieved ~0.1 %). Beyond ≈ 0.6 rad lock 1's
  actuator range (DAC ±0.5 V, V_π 3.5 V) ends;
* noise: lock-1 sensor noise only (i_n = 1 µA/√Hz to make it visible): std(Δ)/std(u) = A−/A+ = e^{−2r}
  ±10 %, with lock 1 slowed to ≈ 2 kHz so that lock 2 (20 kHz) follows u (with equal bandwidths the
  ratio is lower still, ≈ 0.74 A−/A+, since lock 2 tracks only part of u).
"""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry import config
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.config.yaml_io import load_yaml
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


# ------------------------------------------------------------------------------- lock 2 (Ref A, WP-I1)

QUIET = {"koheras.fm_white": False, "mech.mech": False, "pd1.shot_noise": False, "pd1.electronic_noise": False,
         "bhdbeat.shot_noise": False, "adc1.quantization": False, "adc2.quantization": False,
         "dac1.quantization": False, "dac2.quantization": False}
A_RATIO = math.exp(-2 * math.sqrt(9.06 * 0.64))  # A−/A+ = e^{−2r} at 640 mW (T14 §3.1)
KI1 = 2 * math.pi * 2.0e4 / (0.08 * 4.95 * 2 * math.pi / 3.5)  # refA pid1 (20 kHz unity gain)


def refa(duration_s, toggles, params=None):
    t = load_yaml("presets/reference/refA.yaml")
    t["run"]["duration_s"] = duration_s
    t["measure"].pop("esa")  # short runs: no zero-span traces
    t["toggles"] = dict(toggles)
    for b in t["pipeline"]:
        for (name, key), v in (params or {}).items():
            if b["name"] == name:
                b["params"][key] = v
    return run(config.load(t), keep_state=False)


def test_vt15_lock2_static_offset_scaled_by_gain_ratio():
    eps0 = refa(0.02, QUIET).metrics.noise_levels["i.port0"]["mean_loop_rad"]
    u0 = math.atan(math.tan(eps0) / A_RATIO)
    assert u0 == pytest.approx(2 * math.pi * 1e6 * 1.0003 / 299792458.0, rel=0.02)  # SQ: f_sh·τ difference
    for du in (0.1, 0.3, 0.6):
        r = refa(0.02, QUIET, {("dem1", "phi_d_rad"): math.pi / 2 - 2 * du})
        assert r.series["lock.lock1.err"][-200:].mean() / 2 == pytest.approx(du, rel=1e-3)  # u at the PD-OPA
        eps = r.metrics.noise_levels["i.port0"]["mean_loop_rad"]
        assert eps == pytest.approx(math.atan(A_RATIO * math.tan(du + u0)), rel=0.10)


def test_vt15_lock1_sensor_noise_reaches_delta_scaled():
    tg = dict(QUIET, **{"pd1.electronic_noise": True})
    r = refa(0.1, tg, {("pd1", "i_n_a_rthz"): 1e-6, ("pid1", "ki_per_s"): 0.1 * KI1})
    u = r.series["lock.lock1.err"] / 2  # u = Φ1^det/2 (true lock-1 error, lock_error)
    d = r.series["metric.mi.port0"]     # Δ (metric (i))
    assert u.std() > 5e-3               # lock 1 visibly noisy
    assert d.std() / u.std() == pytest.approx(A_RATIO, rel=0.10)
