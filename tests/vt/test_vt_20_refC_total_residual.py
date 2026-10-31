"""VT-20 (WP-I3): Ref C total residual ≈ 0.1 rad and insensitivity to the controller model (prompt §7-10;
00_summary §4 VT-20 ← T08 V8, T16 V19, T14 V7; D25 / D51-7: no Moku:Pro preset).

(a) ``refC_nominal`` (all Ref C noise sources at their nominal axis points, analytic dither): the total
    residual σ_total = √(σ_loop² + σ_HF² + m²/2) is "of order 0.1 rad" — between the dither alone (0.0555 rad)
    and the paper's ≈ 0.1 rad fit — and the non-dither part is ≤ 0.083 rad = √(0.1² − 0.0555²) (T08 V8: a
    consistency check of the parameter ranges with the report, not a fit). The paper's 0.1 rad itself lies
    inside the axis ranges: a slower EOM loop (f_u1 = 1 kHz) with a 2 m path mismatch and a 31 kHz laser gives
    0.07–0.1 rad.
(b) Controller model (T16 V19): τ_L = 0.1, 0.5, 1, 2 µs (fpga_generic sweep; f_sim = 10 MHz so 0.1 µs is one
    sample, D55-1) and a 14-bit DAC change the non-dither residual and σ_total by < 10 % of the 1 µs / 16-bit
    value (same seed: common random numbers, so the differences are systematic).
(c) Sensor noise alone (T14 V7 with the D56-8 cyclostationary weighting, log detector κ_det = 1.32, the OPO's
    band shape): the closed-loop residual equals √(S_Δ ∫|T|² df) within 10 %, with S_Δ = S_b/K_d² from the
    quantum_noise_power parameters of the preset (b-weight 2⟨2 sin²θ S_v(θ)⟩, K_d = ∂b/∂Δ0 of the log law) and
    T = L/(1 + L) of the three-stage stack (integrators, mirror resonance, 2-pole demodulator LPF, τ_L).
    T14's 3.2 mrad (square law, flat band, 10 kHz noise bandwidth) is the same model with other inputs.
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parents[1] / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

M = math.pi * 0.1 / 4.0
F30 = 3e7


def residuals(res):
    nl = res.metrics.noise_levels["i.port0"]
    tot = rs.at(nl, "sigma_total_rad", F30)
    return tot, math.sqrt(max(tot**2 - M**2 / 2, 0.0))


def test_vt20a_nominal_total_of_order_0p1_rad():
    tot, rest = residuals(run(rs.load(run={"duration_s": 0.2})))
    assert M / math.sqrt(2) < tot <= 0.1
    assert rest <= 0.083
    # the paper's ≈ 0.1 rad is inside the axis ranges (slow EOM loop, path mismatch, broad laser)
    tot2, _ = residuals(run(rs.load(params={"f_u1_hz": 1.0e3, "dL_m": 2.0, "h0_hz2_per_hz": 1.0e4},
                                    run={"duration_s": 0.2})))
    assert 0.07 < tot2 <= 0.1


def test_vt20b_insensitive_to_controller_model():
    out = {}
    for tau, bits in ((1e-7, 16), (5e-7, 16), (1e-6, 16), (2e-6, 16), (1e-6, 14)):
        res = run(rs.load(params={"latency_s": tau, "dac_bits": bits}, run={"duration_s": 0.1, "f_sim_hz": 1.0e7}))
        out[(tau, bits)] = residuals(res)
    ref_tot, ref_rest = out[(1e-6, 16)]
    for k, (tot, rest) in out.items():
        assert tot == pytest.approx(ref_tot, rel=0.10), k
        assert rest == pytest.approx(ref_rest, rel=0.10), k


def _loop_noise_bandwidth(p, kappa, f_sim):
    """∫_0^{f_N} |T|² df of the stack's continuous model (module docstring (c))."""
    f = np.geomspace(1e-2, f_sim / 2, 400_000)
    s = 2j * np.pi * f
    w_u1, w12, w23 = (2 * np.pi * p[k] for k in ("f_u1_hz", "f_12_hz", "f_23_hz"))
    wr, q = 2 * np.pi * p["mirror_f_res_hz"], p["mirror_q"]
    mirror = wr**2 / (s**2 + wr / q * s + wr**2)
    lpf = 1.0 / (1.0 + s / (2 * np.pi * p["demod_lpf_hz"])) ** 2
    L = kappa * w_u1 / s * (1 + w12 / s * mirror + w12 * w23 / s**2) * lpf * np.exp(-s * p["latency_s"])
    t2 = np.abs(L / (1 + L)) ** 2
    return float(np.trapezoid(t2, f))


def test_vt20c_sensor_noise_alone_matches_radiometer_model():
    tog = {k: v for k, v in rs.NOISE_OFF.items() if k != "qnl.radiometer_noise"}
    spec = rs.load(params={"x": 0.56}, run={"duration_s": 0.1}, toggles=tog)
    res = run(spec)
    sim = rs.at(res.metrics.noise_levels["i.port0"], "std_loop_rad", F30)
    p = rs.params_of(spec)
    k, kap = p["log_amp_slope_v"], p["kappa_det"]
    pm, pd, q0, q1, q2 = (p[f"qnl_{n}"] for n in ("p_mean", "p_delta", "q0", "q1", "q2"))
    m = p["m_dither_rad"]
    th = 2 * np.pi * np.arange(4096) / 4096

    def b_of(d0):
        pw = pm - pd * np.cos(2 * d0 + 2 * m * np.sin(th))
        return np.mean(2 * np.sin(th) * k * 10 * np.log10(pw))

    kd = (b_of(1e-5) - b_of(-1e-5)) / 2e-5
    c2 = np.cos(2 * m * np.sin(th))
    pw = pm - pd * c2
    s_v = kap * (10 / math.log(10) * k / pw) ** 2 * 2 * (q0 + q1 * c2 + q2 * c2**2)
    s_b = 2 * np.mean(2 * np.sin(th) ** 2 * s_v)
    assert kd == pytest.approx(p["qnl_slope_v_per_rad"], rel=0.01)      # the preset's loop-gain normalisation
    pred = math.sqrt(s_b / kd**2 * _loop_noise_bandwidth(p, kd / p["qnl_slope_v_per_rad"], spec.run.f_sim_hz))
    assert sim == pytest.approx(pred, rel=0.10)
    assert 2e-3 < sim < 1e-2                                             # mrad level: not sensor-limited [Ref C §4]
