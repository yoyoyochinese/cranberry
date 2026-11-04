"""``loop_identification`` (D85-7) on a synthetic discrete loop with a known L(z): lock-in and spectral estimators,
the minus_l and plant relations, coherence, margins; the Moku:Delta summing-point reduction."""

import math

import numpy as np
import pytest
from scipy import signal

from cranberry.blocks.ctrl.excitation import excitation_schedule, schroeder_phases
from cranberry.blocks.meas.loop_identification import identify, loop_margins
from cranberry.hardware.moku import MokuDeltaController, reduce_summing_point

FS = 1.0e5
K = 0.2  # L(z) = K z⁻¹ / (1 − z⁻¹): integrator + one sample delay, f_u ≈ K·FS/2π ≈ 3.2 kHz


def l_true(f):
    zi = np.exp(-2j * np.pi * f / FS)
    return K * zi / (1 - zi)


def stepped(sc, n):
    t0 = round(sc["t_start_s"] * FS)
    dn = round(sc["dwell_s"] * FS)
    x = np.zeros(n)
    for k, (f, g) in enumerate(zip(sc["freqs_hz"], sc["gains"])):
        idx = np.arange(t0 + k * dn, min(n, t0 + (k + 1) * dn))
        x[idx] = sc["amplitude"] * g * np.sin(2 * np.pi * f * (idx - t0 - k * dn) / FS)
    return x


def multisine(sc, n):
    m = np.round(sc["freqs_hz"] * sc["period_s"])
    t = (np.arange(n) - round(sc["t_start_s"] * FS)) / FS
    y = (np.cos(2 * np.pi * np.outer(t, m / sc["period_s"]) + schroeder_phases(m.size))).sum(axis=1) / math.sqrt(m.size)
    return np.where(t >= 0, sc["amplitude"] * y, 0.0)


def closed_loop(x, noise=0.0, seed=0):
    """e = −L(e + x) + d with L = K z⁻¹/(1 − z⁻¹): e = −T x + S d, T = K z⁻¹/(1 − (1 − K) z⁻¹), S = 1 − T."""
    b, a = [0.0, K], [1.0, -(1.0 - K)]
    e = -signal.lfilter(b, a, x)
    if noise:
        d = noise * np.random.default_rng(seed).standard_normal(x.size)
        e += d - signal.lfilter(b, a, d)
    return e


SC_STEP = dict(waveform="stepped_sine", freqs_hz=[300.0, 1.0e3, 3.0e3, 1.0e4, 2.5e4], dwell_s=0.03, settle_s=0.01,
               amp_f_corner_hz=math.inf, t_start_s=0.005, amplitude=0.1)
SC_MS = dict(waveform="multisine", f_start_hz=200.0, f_stop_hz=3.0e4, n_freqs=20, spacing="log", period_s=0.01,
             settle_s=0.01, amp_f_corner_hz=math.inf, t_start_s=0.005, amplitude=0.1)


@pytest.mark.parametrize("cfg", [SC_STEP, SC_MS], ids=["stepped", "multisine"])
def test_minus_l_exact_without_noise(cfg):
    sc = excitation_schedule(cfg)
    n = int(round((0.005 + 5 * 0.03) * FS)) if cfg["waveform"] == "stepped_sine" else int(round(0.3 * FS))
    t = np.arange(n) / FS
    x = stepped(sc, n) if cfg["waveform"] == "stepped_sine" else multisine(sc, n)
    e = closed_loop(x)
    out = identify({"x": x, "reference": e + x, "response": e}, t, sc, FS, relation="minus_l", n_blocks=3)
    lt = l_true(out["f_hz"])
    np.testing.assert_allclose(out["l_hat"], lt, rtol=1e-9)
    np.testing.assert_allclose(out["s_hat"], 1 / (1 + lt), rtol=1e-9)
    np.testing.assert_allclose(out["t_hat"], lt / (1 + lt), rtol=1e-9)
    np.testing.assert_allclose(out["coherence"], 1.0, atol=1e-5)   # sub-blocks cut at the nearest sample
    if cfg["waveform"] == "multisine":
        np.testing.assert_allclose(out["f_hz"], sc["freqs_hz"])


def test_noise_lowers_coherence_but_not_bias():
    sc = excitation_schedule(SC_STEP)
    n = int(round((0.005 + 5 * 0.03) * FS))
    t = np.arange(n) / FS
    x = stepped(sc, n)
    rel = []
    for seed in range(6):
        e = closed_loop(x, noise=0.02, seed=seed)
        out = identify({"x": x, "reference": e + x, "response": e}, t, sc, FS, relation="minus_l", n_blocks=4)
        assert np.all(out["coherence"] < 1.0)
        rel.append(out["l_hat"] / l_true(out["f_hz"]) - 1)
    rel = np.array(rel)
    assert np.all(np.abs(rel.mean(axis=0)) < 4 * rel.std(axis=0) / math.sqrt(6) + 1e-3)   # unbiased (IV)


def test_plant_relation_with_controller():
    """x at the actuator input with the controller replaced: resp = P x, P = −K z⁻¹/(kp) (C = kp)."""
    kp = 2.0
    sc = excitation_schedule(SC_STEP)
    n = int(round((0.005 + 5 * 0.03) * FS))
    t = np.arange(n) / FS
    x = stepped(sc, n)
    resp = -(K / kp) * np.concatenate([[0.0], x[:-1]])
    out = identify({"x": x, "response": resp}, t, sc, FS, relation="plant",
                   controller={"kp": kp, "ki_per_s": 0.0, "kd_s": 0.0}, n_blocks=3)
    zi = np.exp(-2j * np.pi * out["f_hz"] / FS)
    np.testing.assert_allclose(out["p_hat"], -(K / kp) * zi, rtol=1e-9)
    np.testing.assert_allclose(out["l_hat"], K * zi, rtol=1e-9)


def test_model_comparison_and_margins():
    """L = k/s · e^{−sτ}: f_u = k/2π, PM = 90° − 360° f_u τ, f_180 = 1/(4τ), GM = −20 log10(k/(2π f_180))."""
    k, tau = 2 * np.pi * 2.0e3, 20e-6
    model = {"k0": 1.0, "kp": 0.0, "ki_per_s": k, "kd_s": 0.0, "delay_s": tau, "poles_hz": [], "modes": []}
    f = np.geomspace(100.0, 1.0e5, 3000)
    lf = k / (2j * np.pi * f) * np.exp(-2j * np.pi * f * tau)
    m = loop_margins(f, lf)
    assert m["f_u_hz"] == pytest.approx(2.0e3, rel=1e-4)
    assert m["phase_margin_deg"] == pytest.approx(90 - 360 * 2.0e3 * tau, abs=0.01)
    assert m["f_180_hz"] == pytest.approx(1 / (4 * tau), rel=1e-4)
    assert m["gain_margin_db"] == pytest.approx(-20 * math.log10(k / (2 * math.pi / (4 * tau))), abs=0.01)
    assert math.isnan(loop_margins(f[:10], lf[:10])["f_180_hz"])
    # identify's model comparison: a measured L equal to the model gives zero error
    sc = excitation_schedule(SC_STEP)
    n = int(round((0.005 + 5 * 0.03) * FS))
    t = np.arange(n) / FS
    x = stepped(sc, n)
    e = closed_loop(x)
    out = identify({"x": x, "reference": e + x, "response": e}, t, sc, FS, relation="minus_l", n_blocks=3,
                   model_params={**model, "ki_per_s": K * FS, "delay_s": 0.5 / FS})
    low = out["f_hz"] <= 1.0e4   # continuous model of a discrete integrator: (ωh)² errors above
    assert np.all(np.abs(out["mag_err"][low]) < 0.05) and np.all(np.abs(out["phase_err_deg"][low]) < 5.0)
    assert "model_f_u_hz" in out and np.all(np.isfinite(out["l_phase_deg"]))


def test_moku_summing_point_reduction_and_plan():
    lt = np.array([5.0 - 2j, 1.0 + 0.5j, 0.1 - 0.3j])
    for g in (1.0, 0.5):
        h = -g * lt / (1 + lt)                                            # FRA ratio e/x with e = −L(e + g x)
        r = reduce_summing_point(h, g)
        np.testing.assert_allclose(r["l_hat"], lt)
        np.testing.assert_allclose(r["s_hat"], 1 / (1 + lt))
        np.testing.assert_allclose(r["t_hat"], lt / (1 + lt))
    ctl = MokuDeltaController({"v.adc0": "In1", "ctrl.pid0f.u": "Out1"})
    plan = ctl.loop_gain_plan("lock0", "stepped_sine", [1e3, 1e4, 1e5], amplitude=0.05, settle_s=5e-3, dwell_s=25e-3,
                              error_in="v.adc0", output="ctrl.pid0f.u")
    assert plan["fra"]["points"] == 3 and "UNVERIFIED" in plan["pid_control_matrix"]["gain_unit"]
    with pytest.raises(NotImplementedError):
        ctl.acquire(1.0, ["v.adc0"])
    with pytest.raises(NotImplementedError):
        ctl.loop_gain_plan("lock0", "multisine", [1e3, 1e4], amplitude=0.05, settle_s=0.0, dwell_s=1.0,
                           error_in="v.adc0", output="ctrl.pid0f.u")
