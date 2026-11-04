"""WP6-P (D85-6, P6 Zhong et al. 2021): program the GBS port phases through the lock set point.

* every program overlay builds on every future base (heterodyne RF phase, its dither bad config, the DC set-point bad
  config);
* validation on F_wg6 nominal: under heterodyne detection the lock point is X = π/2 − φ_d, so the RF-phase program sets
  X_k − π/2 = p_k (and Θ_1 = p_1 − p_0) one to one whatever the probe power and visibility, while the DC set point
  calibrated at the nominal K moves to arcsin((K_cal/K) sin p) after a power / visibility drift; programming the dither
  demodulation phase does not move the lock point at all;
* settling time against the loop bandwidth f_u (integrator loop: t_s ≈ ln(|Δ|/tol)/(2π f_u));
* reference loop = compiled loop with the program (metric_ii setpoint_slot, slip_center_slot, in_program).
"""

from __future__ import annotations

import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_spec = importlib.util.spec_from_file_location("future_support", Path(__file__).with_name("future_support.py"))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

pytestmark = pytest.mark.numba
BASES = ("F_wg100", "F_wg6", "F_bor6", "F_liu20")
HET = ("probe/detect_heterodyne", "program/p6_steps")
DITHER = ("probe/detect_dither", "program/p6_steps")
DC = ("program/p6_steps_dc",)
TABLE = [[0.3, -0.2], [-0.2, 0.3]]
QUIET = {"laser.fm_white": False, "laser.fm_flicker": False}  # the lock point, not the laser decorrelation noise


@pytest.mark.parametrize("base", BASES)
@pytest.mark.parametrize("ovl", [HET, DITHER, DC], ids=["het", "dither", "dc"])
def test_overlay_builds_on_every_base(base, ovl):
    r = fs.build_check(fs.nominal(base), ovl)
    assert "p6" in r.measurements["step_response"]
    assert fs.block(r, "sr").kind == "step_response" and "ctrl.prog.target1" in r.series


def _overlay(tmp: Path, scheme: str, dwell_s: float, f_u_hz: float, liu: dict) -> str:
    """Last overlay: explicit 2-setting table, f_u, and the probe power / visibility of both liu blocks changed AFTER
    calibration (liu_k_v, i.e. the loop gain normalisation and K_cal, stay nominal). It must come after the detection
    overlay (detect_heterodyne / detect_dither replace liu0/liu1 wholesale)."""
    import yaml

    pipe = [{"name": "squeezer", "seq": [{"name": "prog", "params": {
        "source": "table", "times_s": [5.0e-3, 5.0e-3 + dwell_s], "table_rad": TABLE}}]}]
    for k in (0, 1) if liu else ():
        pipe.append({"name": f"liu{k}", "params": dict(liu)} if scheme == "dc"
                    else {"name": f"liu{k}", "seq": [{"name": f"liu{k}", "params": dict(liu)}]})
    p = tmp / f"o_{scheme}_{dwell_s}_{f_u_hz}_{'_'.join(f'{a}{b}' for a, b in liu.items())}.yaml"
    p.write_text(yaml.safe_dump({"offline": {"params": {"prog_dwell_s": dwell_s, "f_u_hz": f_u_hz}}, "pipeline": pipe}))
    return str(p)


def _run(tmp, scheme, dwell_s=3e-3, f_u_hz=1e3, mode="compiled", decim=10, **liu):
    ovl = {"het": HET, "dither": DITHER, "dc": DC}[scheme]
    spec = fs.load(fs.nominal("F_wg6"), run={"duration_s": 5e-3 + 2 * dwell_s}, toggles=QUIET,
                   overlays=(*ovl, _overlay(tmp, scheme, dwell_s, f_u_hz, liu)),
                   measure={"record": {"slots": ["lock.L0.err"], "decim": decim, "antialias": "none"}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(fs.first_point(spec), mode=mode)


CASES = {"nominal": {}, "visibility_0.5": {"visibility": 0.5}, "power_x0.25": {"p_pr_avg_w": 2.5e-7}}


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p6")
    return {(s, c): _run(tmp, s, **kw) for s in ("het", "dc") for c, kw in CASES.items()} | {
        ("dither", "nominal"): _run(tmp, "dither")}


def test_rf_program_is_one_to_one_whatever_the_probe_power_and_visibility(runs):
    for case in CASES:
        r = runs[("het", case)]
        sr = r.measurements["step_response"]["p6"]
        assert r.valid and sr["t_switch_s"].size == 2 and np.all(sr["slips"] == 0)
        np.testing.assert_allclose(sr["bias_rad"], 0.0, atol=3e-3)          # X_k − π/2 = p_k (both ports, both settings)
        np.testing.assert_allclose(sr["theta_bias_rad"], 0.0, atol=1e-2)    # Θ_1 = p_1 − p_0 on the squeezed light
        np.testing.assert_allclose(sr["pair_bias_rad"], 0.0, atol=3e-3)     # cross-port: (X_1 − X_0) vs p_1 − p_0
        assert sr["cross_d_settle_s"].size == 0  # TABLE has no equal steps on the two ports
    # K/2 really reached the loop: half the loop gain, about twice the settling time — but the same lock points
    nom = runs[("het", "nominal")].measurements["step_response"]["p6"]["settle_s"][:, 0]
    for case in ("visibility_0.5", "power_x0.25"):
        assert np.all(runs[("het", case)].measurements["step_response"]["p6"]["settle_s"][:, 0] > 1.5 * nom)


def test_dc_setpoint_program_moves_with_K(runs):
    sr = runs[("dc", "nominal")].measurements["step_response"]["p6"]
    np.testing.assert_allclose(sr["bias_rad"], 0.0, atol=5e-3)              # calibrated: correct
    for case, ratio in (("visibility_0.5", 2.0), ("power_x0.25", 2.0)):     # K_cal/K = 2 (K ∝ V √P̄)
        sr = runs[("dc", case)].measurements["step_response"]["p6"]
        g = np.array(TABLE)
        want = np.arcsin(np.clip(ratio * np.sin(g), -1, 1)) - g             # X − π/2 = arcsin((K_cal/K) sin p)
        np.testing.assert_allclose(sr["bias_rad"], want, atol=0.02)
        assert np.min(np.abs(sr["bias_rad"])) > 0.15


def test_dither_demodulation_phase_does_not_program(runs):
    sr = runs[("dither", "nominal")].measurements["step_response"]["p6"]
    g = np.array(TABLE)
    # dither lock: X stays at the extremum (lock.L<k>.err = X − 0 with the dither set points), so the error is −p
    np.testing.assert_allclose(sr["bias_rad"], -g, atol=5e-3)
    np.testing.assert_allclose(sr["theta_bias_rad"][:, 0], -(g[:, 1] - g[:, 0]), atol=1e-2)


def test_settling_time_scales_with_the_loop_bandwidth(tmp_path):
    """Port 0, |Δ| = 0.3 and 0.5 rad, tol 10 mrad: an integrator loop rises (10–90 %) in ln 9/(2π f_u) and settles in
    ln(|Δ|/tol)/(2π f_u). At f_u = 0.5 and 1 kHz the stretcher resonance (20 kHz, Q 10: |L(f_res)| = Q f_u/f_res ≤ 0.5)
    and the ~2 µs delay add little (measured 1.10/1.25 ms and 0.53/0.67 ms against 1.08/1.25 and 0.54/0.62 ms; at
    1.5 kHz the ringing makes it 0.52/0.77 ms against 0.36/0.42 ms, and 3 kHz is loop_unstable)."""
    for f_u in (5e2, 1e3):
        sr = _run(tmp_path, "het", dwell_s=3e-3, f_u_hz=f_u, decim=1).measurements["step_response"]["p6"]
        steps, ts = np.abs(sr["step_rad"][:, 0]), sr["settle_s"][:, 0]
        ratio = ts / (np.log(steps / 0.01) / (2 * math.pi * f_u))
        assert np.all((ratio > 0.85) & (ratio < 1.25)), (f_u, ratio)
        np.testing.assert_allclose(sr["rise_s"][:, 0], math.log(9) / (2 * math.pi * f_u), rtol=0.1)
        assert np.all(sr["overshoot"][:, 0] < 0.03)


def test_reference_equals_compiled(tmp_path):
    a = _run(tmp_path, "het", dwell_s=2e-4, mode="reference", decim=1)
    spec_dur = a.t_s.size
    b = _run(tmp_path, "het", dwell_s=2e-4, mode="compiled", decim=1)
    assert spec_dur == b.t_s.size
    for k in ("lock.L0.err", "metric.mii.port1", "ctrl.prog.target1"):
        np.testing.assert_allclose(a.series[k], b.series[k], rtol=0, atol=1e-12)
