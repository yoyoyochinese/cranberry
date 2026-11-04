"""Phase 6 WP6-AB (D85-1): the lab/acoustic and lab/multimode915b overlays (Ref D and future variants) build, feed the
same mode list to the stretcher and the stability check, and the controller overlays keep the stability check on the
controller law the time loop runs (D82 lesson: controller/dsc1 → discretisation dsc1; Moku:Delta stays continuous)."""

import dataclasses
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.blocks.err.testing import rig
from cranberry.core.registry import ensure_loaded, make_block
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]
REFD_NOMINAL = "presets/reference/refD_nominal.yaml"
MM = {"mm2_f_hz": 5.4e3, "mm2_q": 15.0, "mm2_r": 0.11, "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07,
      "mm4_f_hz": 4.03e4, "mm4_q": 100.0, "mm4_r": 0.15}
# Phase 6 fix axes (B-1 slope, B-6 arm correlation, B-2 DAC dither) at their D89 values
NEW_AXES = {"acoustic_hf_slope": 2.0, "acoustic_rho_arms": 0.0, "dac_output_noise_v": 0.0}
REFD_POINT = {"acoustic_scale": 1.0, "acoustic_gamma": 1.0, "f_res_hz": 2.52e4, "q_res": 25.0, **MM, **NEW_AXES}


def find(spec, name):
    """The BlockSpec ``name`` anywhere in a pipeline spec tree."""
    if getattr(spec, "name", None) == name and hasattr(spec, "params"):
        return spec
    for f in dataclasses.fields(spec):
        for c in (v if isinstance(v := getattr(spec, f.name), tuple) else (v,)):
            if dataclasses.is_dataclass(c) and not isinstance(c, type) and (r := find(c, name)) is not None:
                return r
    return None


def chk(spec, name):
    ensure_loaded()
    (b,) = [make_block(x) for x in offline_block_specs(spec) if x.name == name]
    return b


def refd_point(**extra):
    base = config.load(REFD_NOMINAL, overlays=["lab/acoustic", "lab/multimode915b"])
    (pt,) = expand(base, {f"offline.params.{k}": [v] for k, v in {**REFD_POINT, **extra}.items()})
    return pt.spec


def test_refd_overlays_same_modes_and_table():
    spec = refd_point()
    c = chk(spec, "chk.lockD_stability")
    st = find(spec.pipeline, "stretcher915b").params
    assert c.params["modes"] == st["modes"] and len(st["modes"]) == 4
    assert sum(m["r"] for m in st["modes"]) == pytest.approx(1.0, abs=1e-12)
    assert c.result["stable"]  # the fitted point keeps the D53 loop stable (gain margin ≈ 3.3 dB at 5.4 kHz)
    mech = find(spec.pipeline, "env.mech").params
    assert mech["mech_shape"] == "table" and mech["gamma"] == 1.0
    assert [r[1] for r in mech["psd_table"][1:5]] == pytest.approx([2e-5, 3.6e-6, 6.7e-7, 9e-8])  # COMPARISON §3 (A)
    # γ = 0: every arm carries half of the measured (two-arm) PSD
    g0 = find(refd_point(acoustic_gamma=0.0).pipeline, "env.mech").params
    assert g0["psd_table"][1][1] == pytest.approx(1e-5)


def test_refd_overlays_short_run_locks():
    spec = dataclasses.replace(refd_point(), run=dataclasses.replace(refd_point().run, duration_s=0.02))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(spec)
    assert res.valid and not res.events


@pytest.mark.parametrize("overlays", [("lab/acoustic_future",), ("lab/multimode915b_future",),
                                      ("lab/acoustic_future", "lab/multimode915b_future", "controller/dsc1")],
                         ids=["acoustic", "multimode", "both+dsc1"])
def test_future_overlays_build(overlays):
    """F_liu20_nominal + the future overlays at one point of their axes (pinned with expand: overlays are applied after
    a config's own offline.params, so fs.load(params=...) would not pin them), three reference-mode samples."""
    settle0 = [{"name": "mi", "params": {"t_settle_s": 0.0}}, {"name": "mii", "params": {"t_settle_s": 0.0}}]
    base = fs.load(fs.nominal("F_liu20"), run={"duration_s": 3.0e-6}, overlays=overlays, pipeline=settle0)
    point = {"acoustic_scale": 1.0, "acoustic_gamma": 0.0, "pzt_f_res_hz": 2.52e4, "pzt_q": 25.0, **MM, **NEW_AXES,
             "dsc1_latency_s": 5.7e-6}
    axes = {k: [v] for k, v in ((f"offline.params.{k}", v) for k, v in point.items()) if k in config.sweep_markers(base)}
    (pt,) = expand(base, axes)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(pt.spec, mode="reference")
    assert 1 in res.metrics.metric_ii
    if "lab/multimode915b_future" in overlays:
        c = fs.block(res, "chk.lock_stability")
        assert len(c.params["modes"]) == 4 and all("r" in m for m in c.params["modes"])
        assert find(pt.spec.pipeline, "pzt1").params["modes"] == c.params["modes"]
    if "lab/acoustic_future" in overlays:
        assert find(pt.spec.pipeline, "env.mech.lo").params["psd_table"][1][1] == pytest.approx(1e-5)  # γ = 0: T/2


def _pid_response(disc, ki, f, f_sim):
    obj = {"block": "pid", "name": "p", "params": {
        "in": "v.x", "out": "ctrl.p.u", "device": None, "anti_windup": "none", "discretisation": disc,
        "kp": 0.0, "ki_per_s": ki, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
        "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1.0, "u_max_v": 1.0,
        "toggles": {"feedback": True, "p_term": True, "i_term": True, "d_term": True}}}
    root, _, pr = rig([obj], f_sim, external={"v.x": "V"})
    return root.freq_response(f, pr) if hasattr(root, "freq_response") else root.children[0].freq_response(f, pr)


def _check_ctrl(c, f):
    """The check's controller C(i2πf) alone (k0 = 1, no delay / poles / modes / stages)."""
    p = {**c.params, "k0": 1.0, "delay_s": 0.0, "poles_hz": [], "modes": [], "stages": []}
    return open_loop(f, p)


def test_dsc1_overlay_discretises_future_stability_check():
    """D82 lesson + D91-1 on the future presets: with controller/dsc1 the fine-stage PIDs run the DSC1 law at the
    device's 500 kHz (rate_div 2, h = 2 µs), so chk.lock_stability uses discretisation dsc1 with that h; its controller
    then follows the time-loop block's discrete law (I-term approximation (ωh)²/12), where the continuous law is off by
    the K_i·h/2 lead (ωh/2)."""
    base = fs.load(fs.nominal("F_liu20"), overlays=("controller/dsc1",))
    (pt,) = expand(base, {"offline.params.dsc1_latency_s": [5.7e-6]})
    spec = pt.spec
    c = chk(spec, "chk.lock_stability")
    f_sim = spec.run.f_sim_hz
    assert c.params["discretisation"] == "dsc1" and c.params["h_s"] == pytest.approx(2.0 / f_sim)
    pid = find(spec.pipeline, "pid0f")
    assert pid.params["device"] == "dsc1" and pid.params["discretisation"] == "dsc1" and pid.rate_div == 2
    assert float(pid.params["ki_per_s"]) == pytest.approx(float(c.params["ki_per_s"]))
    f = np.geomspace(1e2, 3e4, 20)
    want = _pid_response("dsc1", float(c.params["ki_per_s"]), f, f_sim / 2)
    np.testing.assert_allclose(_check_ctrl(c, f), want, rtol=1.5e-2)  # (ωh)²/12 = 1.2 % at 30 kHz, h = 2 µs
    cont = chk(fs.load(fs.nominal("F_liu20")), "chk.lock_stability")
    assert cont.params.get("discretisation", "continuous") == "continuous"
    lead = np.angle(want[-1] / _check_ctrl(cont, f[-1:])[0])
    assert lead == pytest.approx(np.pi * f[-1] * 2 / f_sim, rel=0.05)  # ωh/2 = 10.8° at 30 kHz (h = 2 µs)


def test_moku_overlay_check_stays_continuous():
    """controller/moku_delta: Tustin at h = 1/f_sim (1 µs); for F_common's pure-I fine stage it has the exact −90° phase
    and |C| within (ωh)²/12 = 0.3 % at 30 kHz of the continuous law, so the check needs no discrete form."""
    spec = fs.load(fs.nominal("F_liu20"), overlays=("controller/moku_delta",))
    c = chk(spec, "chk.lock_stability")
    assert c.params.get("discretisation", "continuous") == "continuous"
    assert find(spec.pipeline, "pid0f").params["discretisation"] == "tustin"
    f = np.geomspace(1e2, 3e4, 20)
    np.testing.assert_allclose(_check_ctrl(c, f), _pid_response("tustin", float(c.params["ki_per_s"]), f,
                                                                 spec.run.f_sim_hz), rtol=4e-3)
