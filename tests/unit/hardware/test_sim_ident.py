"""D85-7 additions to the controller interface on SimController: per-loop engage / disengage, ``acquire(...,
injections_on=)``, ``measure_loop_gain``. The runner is a stand-in that closes a known discrete loop around the
excitation the controller configured; the old call signatures keep their behaviour."""

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import pytest
from scipy import signal

from cranberry.blocks.ctrl.excitation import excitation_schedule
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core import registry
from cranberry.core.block import BaseBlock, BlockSpec
from cranberry.core.pipeline import SeqSpec
from cranberry.hardware.interface import ControllerInterface
from cranberry.hardware.sim import SimController
from cranberry.runner.run import RunResult

FS = 1.0e5
K = 0.2


class _Stub(BaseBlock):
    def __call__(self, state, params):
        return state


@registry.register_block("test_ident_pid")
def make_test_ident_pid(params, *, name=None):
    return _Stub("test_ident_pid", name, params, runtime_params=("kp", "ki_per_s", "kd_s"), toggle_names=("feedback",))


def make_spec():
    exc = BlockSpec("excitation", "exc0", {"in": None, "out": "v.ref", "unit": "V", "waveform": "stepped_sine",
                                           "freqs_hz": [1.0e3], "dwell_s": 0.01, "settle_s": 0.002,
                                           "amp_f_corner_hz": math.inf, "amplitude": 0.0, "t_start_s": 0.0,
                                           "toggles": {"enabled": False}})
    fast = BlockSpec("test_ident_pid", "pid0f", {"kp": 0.0, "ki_per_s": 1.0, "kd_s": 0.0, "toggles": {"feedback": False}})
    slow = BlockSpec("test_ident_pid", "pid0s", {"kp": 0.0, "ki_per_s": 1.0, "kd_s": 0.0, "toggles": {"feedback": False}})
    return RunSpec(run=RunSettings(1e-3, FS, 4096, 7), pipeline=SeqSpec((exc, fast, slow)), preset="unit_test",
                   measure={"record": {"slots": [], "decim": 1}})


def leaf(spec, name):
    return next(b for b in spec.pipeline.children if b.name == name)


class LoopRunner:
    """x from the excitation block of the spec (its schedule), e = −T x with L = K z⁻¹/(1 − z⁻¹)."""

    def __init__(self):
        self.specs = []

    def __call__(self, spec):
        self.specs.append(spec)
        n = spec.run.n_samples
        t = np.arange(n) / FS
        p = leaf(spec, "exc0").params
        x = np.zeros(n)
        if spec.toggles.get("exc0.enabled"):
            sc = excitation_schedule(p)
            t0, dn = round(sc["t_start_s"] * FS), round(sc["dwell_s"] * FS)
            for k, (f, g) in enumerate(zip(sc["freqs_hz"], sc["gains"])):
                idx = np.arange(t0 + k * dn, min(n, t0 + (k + 1) * dn))
                x[idx] = p["amplitude"] * g * np.sin(2 * np.pi * f * (idx - t0 - k * dn) / FS)
        e = -signal.lfilter([0.0, K], [1.0, -(1.0 - K)], x)
        series = {"ctrl.exc0.x": x, "v.ref": e + x, "v.err": e, "v.other": np.zeros(n)}
        return RunResult(spec=spec, run_dir=None, t_s=t, series=series)


CONFIG = {
    "loops": {"lock0.fast": {"pid_block": "pid0f", "device": None}, "lock0.slow": {"pid_block": "pid0s", "device": None}},
    "excitations": {"err0": {"block": "exc0", "loop": "lock0.fast", "relation": "minus_l", "reference": "v.ref",
                             "response": "v.err", "controller": None, "model": None}},
}


@pytest.fixture
def ctl():
    runner = LoopRunner()
    c = SimController(make_spec(), run_fn=runner)
    c.configure(CONFIG)
    c.runner = runner
    return c


def test_per_loop_engage_and_old_calls(ctl):
    ctl.engage()
    assert ctl.status()["engaged"] and all(v["engaged"] for v in ctl.status()["loops"].values())
    ctl.disengage("lock0.fast")
    ctl.acquire(1e-3, ["v.other"])
    tog = ctl.runner.specs[-1].toggles
    assert tog["pid0f.feedback"] is False and tog["pid0s.feedback"] is True and tog["exc0.enabled"] is False
    assert not ctl.status()["engaged"]
    ctl.engage("lock0.fast")
    ctl.disengage()
    ctl.acquire(1e-3, ["v.other"])
    assert not any(v for k, v in ctl.runner.specs[-1].toggles.items() if k.endswith(".feedback"))
    with pytest.raises(KeyError):
        ctl.engage("nope")
    rec = ctl.log_record()
    assert rec["loops_engaged"] == {"lock0.fast": False, "lock0.slow": False}
    assert {"op": "disengage", "loop": "lock0.fast"} in rec["ops"]


def test_acquire_injections_on(ctl):
    ctl.set_params({"exc0.amplitude": 0.1})
    ctl.acquire(1e-3, ["v.err"], injections_on=True)
    assert ctl.runner.specs[-1].toggles["exc0.enabled"] is True
    ctl.acquire(1e-3, ["v.err"], injections_on=["err0"])
    assert ctl.runner.specs[-1].toggles["exc0.enabled"] is True
    ctl.acquire(1e-3, ["v.err"])                                       # default: off again
    assert ctl.runner.specs[-1].toggles["exc0.enabled"] is False
    with pytest.raises(KeyError):
        ctl.acquire(1e-3, ["v.err"], injections_on=["nope"])


def test_measure_loop_gain_stepped_and_band(ctl):
    ctl.engage()
    freqs = [500.0, 2.0e3, 8.0e3]
    res = ctl.measure_loop_gain("lock0.fast", "stepped_sine", freqs, amplitude=0.05, settle_s=2e-3, dwell_s=12e-3,
                                t_start_s=1e-3)
    spec = ctl.runner.specs[-1]
    p = leaf(spec, "exc0").params
    assert p["waveform"] == "stepped_sine" and p["amplitude"] == 0.05 and spec.toggles["exc0.enabled"]
    assert set(spec.measure["record"]["slots"]) == {"ctrl.exc0.x", "v.ref", "v.err"}
    assert spec.run.duration_s == pytest.approx(1e-3 + 3 * 12e-3 + 2 / FS)
    zi = np.exp(-2j * np.pi * res["f_hz"] / FS)
    np.testing.assert_allclose(res["l_hat"], K * zi / (1 - zi), rtol=1e-9)
    assert res["f_u_hz"] == pytest.approx(K * FS / (2 * math.pi), rel=0.05)
    assert {"f_hz", "l_hat", "coherence", "s_hat", "t_hat", "phase_margin_deg", "gain_margin_db"} <= set(res)
    res = ctl.measure_loop_gain("lock0.fast", "stepped_sine", band_hz=[300.0, 3.0e3], n_freqs=4, amplitude=0.05,
                                settle_s=2e-3, dwell_s=12e-3)
    assert res["f_hz"].size == 4
    ops = [o["op"] for o in ctl.log_record()["ops"]]
    assert ops.count("measure_loop_gain") == 2


def test_measure_loop_gain_errors(ctl):
    with pytest.raises(KeyError):
        ctl.measure_loop_gain("lock0.slow", "stepped_sine", [1e3], amplitude=0.1, settle_s=1e-3, dwell_s=5e-3)
    with pytest.raises(TypeError, match="settle_s"):
        ctl.measure_loop_gain("lock0.fast", "stepped_sine", [1e3], amplitude=0.1, dwell_s=5e-3)
    with pytest.raises(TypeError, match="unknown setting"):
        ctl.measure_loop_gain("lock0.fast", "stepped_sine", [1e3], amplitude=0.1, settle_s=1e-3, dwell_s=5e-3, foo=1)
    bad = dict(CONFIG, excitations={"err0": {**CONFIG["excitations"]["err0"], "reference": None}})
    with pytest.raises(ValueError, match="reference"):
        SimController(make_spec(), run_fn=LoopRunner()).configure(bad)


def test_base_interface_default_measure_loop_gain():
    class Old(ControllerInterface):  # a backend written before D85-7 (old signatures) still instantiates
        def configure(self, config): ...
        def set_params(self, params): return {}
        def get_params(self): return {}
        def engage(self): ...
        def disengage(self): ...
        def reset_integrators(self): ...
        def ramp(self, channel, *, amplitude, period_s, offset, enable): ...
        def acquire(self, duration_s, signals): return {}
        def measure_transfer_function(self, injection_point, freqs_hz): return np.zeros(0)
        def status(self) -> Mapping[str, Any]: return {}
        def log_record(self): return {}

    with pytest.raises(NotImplementedError):
        Old().measure_loop_gain("l", "stepped_sine", [1.0], amplitude=1.0)


def test_measure_loop_gain_refine_coarse_margin(ctl):
    """Phase 5 review I-3: a margin flagged coarse is re-swept with refine_n points around it and merged."""
    ctl.engage()
    kw = dict(amplitude=0.05, settle_s=2e-3, dwell_s=12e-3, t_start_s=1e-3)
    res = ctl.measure_loop_gain("lock0.fast", "stepped_sine", [500.0, 2.0e3, 8.0e3], margin_dphase_max_deg=1.0, **kw)
    assert res["f_u_coarse"] and "refined" not in res
    n0 = len(ctl.log_record()["ops"])
    res = ctl.measure_loop_gain("lock0.fast", "stepped_sine", [500.0, 2.0e3, 8.0e3], margin_dphase_max_deg=1.0,
                                refine_n=5, **kw)
    assert len(ctl.log_record()["ops"]) == n0 + 2                     # coarse run + one refinement run
    f = res["refined"]["f_hz"]
    assert f.size > 3 and np.all(np.diff(f) > 0) and f.min() == 500.0 and f.max() == 8.0e3
    zi = np.exp(-2j * np.pi * f / FS)
    np.testing.assert_allclose(res["refined"]["l_hat"], K * zi / (1 - zi), rtol=1e-9)
    assert res["f_u_hz"] == pytest.approx(K * FS / (2 * math.pi), rel=0.02)
    assert res["coarse_margins"]["f_u_coarse"]
