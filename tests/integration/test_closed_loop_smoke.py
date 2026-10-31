"""Closed-loop smoke test (WP-F1; groundwork for prompt §7-7, VT-13/14): integrator plant + PI
controller + delay buffer, driven by a sinusoidal output disturbance through ``runner.run``. The
measured sensitivity S = y/d matches

* the exact z-domain transfer function of the pipeline (relative 1e−8), and
* the continuous-time model S(s) = 1/(1 + L(s)), L(s) = K (k_p + k_i/s) e^{−sτ}/s with
  τ = tau_s + 1/f_sim (loop-closure sample, ``core.delay``), within 2 % below the 10 kHz crossover
  and 6 % up to f_sim/30.

Pipeline (execution order; the plant reads the actuator output of the previous sample):
    dist: d = A sin(2πft)  →  plant: x ← x + hK·a  →  sens: y = x + d  →  pid: u = −(k_p y + k_i h Σ y)
    →  delay_buffer (k samples)  →  act: a = u
z-domain: L(z) = hK z⁻¹/(1 − z⁻¹) · z⁻ᵏ · (k_p + k_i h/(1 − z⁻¹)),  S = 1/(1 + L).
"""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.compile import loop_cache_size
from cranberry.core.pipeline import spec_from_obj
from cranberry.core.testing import toy
from cranberry.runner.run import run

pytestmark = [pytest.mark.numba, pytest.mark.integration]

FS = 1.0e6
H = 1.0 / FS
K = 1.0
KP = 2 * np.pi * 1.0e4 / K  # unity-gain near 10 kHz
KI = KP * 2 * np.pi * 1.0e3  # PI corner at 1 kHz
K_DELAY = 4  # samples in the delay buffer
N_TOTAL, N_WIN = 40_000, 10_000  # 40 ms run, last 10 ms analysed (integer periods for f = m·100 Hz)


def loop_spec(f_hz):
    obj = [
        toy("toy_sine", "dist", out="v.d", unit="rad", amp=1e-3, freq_hz=f_hz),
        toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="rad", k=K),
        toy("toy_sum", "sens", ins=["psi.plant.x", "v.d"], gains=[1.0, 1.0], out="v.y", unit="rad"),
        toy("toy_pid", "pid", in_="v.y", out="ctrl.pid.u", unit="V", kp=KP, ki=KI, kd=0.0, setpoint=0.0),
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out",
                                                           "tau_s": K_DELAY * H, "method": "linear"}},
        toy("toy_gain", "act", in_="delay.dl.out", out="act.u", unit="V", gain=1.0),
    ]
    return RunSpec(run=RunSettings(duration_s=N_TOTAL * H, f_sim_hz=FS, chunk_samples=8192, seed=0),
                   pipeline=spec_from_obj(obj), measure={"record": {"slots": ["v.y", "v.d"], "decim": 1}})


def s_z(f_hz):
    z1 = np.exp(-2j * np.pi * f_hz * H)
    L = H * K * z1 / (1 - z1) * z1**K_DELAY * (KP + KI * H / (1 - z1))
    return 1.0 / (1.0 + L)


def s_cont(f_hz):
    s = 2j * np.pi * f_hz
    L = K * (KP + KI / s) * np.exp(-s * (K_DELAY + 1) * H) / s
    return 1.0 / (1.0 + L)


def measured_s(f_hz):
    res = run(loop_spec(f_hz))
    n = np.arange(N_TOTAL - N_WIN, N_TOTAL)
    ph = np.exp(-2j * np.pi * f_hz * n * H)
    y, d = res.series["v.y"][n], res.series["v.d"][n]
    return np.sum(y * ph) / np.sum(d * ph)


@pytest.mark.parametrize("f_hz", [300.0, 1_000.0, 3_000.0, 10_000.0, 30_000.0])
def test_sensitivity_matches_z_domain_exactly(f_hz):
    sm = measured_s(f_hz)
    assert abs(sm - s_z(f_hz)) < 1e-8 * abs(s_z(f_hz))


@pytest.mark.parametrize("f_hz", [300.0, 1_000.0, 3_000.0, 10_000.0, 30_000.0])
def test_sensitivity_matches_continuous_model(f_hz):
    sm, sc = measured_s(f_hz), s_cont(f_hz)
    tol = 0.02 if f_hz < 10_000 else 0.06  # discretisation O(f/f_sim) near the 10 kHz crossover
    assert abs(sm - sc) < tol * abs(sc)


def test_loop_reused_across_frequencies_and_rejection():
    measured_s(500.0)
    size = loop_cache_size()
    s_lo = measured_s(200.0)  # runtime parameter (freq_hz) only: same compiled loop
    assert loop_cache_size() == size
    assert abs(s_lo) < 0.01  # strong rejection well inside the loop bandwidth (|S| ≈ f²/(f_u f_i))
