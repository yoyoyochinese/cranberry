"""SimController and the abstract interface (PLAN §2.10, T16 §4.6). The runner is replaced by a stand-in;
the real-runner test uses WP-F1's ``runner.run``."""

import json
import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import pytest

from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core import registry
from cranberry.core.block import BaseBlock, BlockSpec
from cranberry.core.pipeline import SeqSpec, StackSpec
from cranberry.hardware.interface import ControllerInterface
from cranberry.hardware.mapping import DeviceRangeWarning
from cranberry.hardware.sim import SimController, demodulate
from cranberry.runner.run import RunResult


class _Stub(BaseBlock):
    def __init__(self, kind, params: Mapping[str, Any], name, runtime, toggles):
        super().__init__(kind, name, params, runtime_params=runtime, toggle_names=toggles)

    def __call__(self, state, params):
        return state


@registry.register_block("test_hw_pid")
def make_test_hw_pid(params, *, name=None):
    return _Stub("test_hw_pid", params, name, ("kp", "ki_per_s", "kd_s"), ("feedback",))


@registry.register_block("test_hw_ramp")
def make_test_hw_ramp(params, *, name=None):
    return _Stub("test_hw_ramp", params, name, ("amplitude_v", "period_s", "offset_v"), ("enabled",))


@registry.register_block("test_hw_inj")
def make_test_hw_inj(params, *, name=None):
    return _Stub("test_hw_inj", params, name, ("freq_hz", "amplitude"), ("enabled",))


@registry.register_block("test_hw_plant")
def make_test_hw_plant(params, *, name=None):
    return _Stub("test_hw_plant", params, name, ("gain",), ())


F_SIM = 1.0e6
FC = 2.0e3  # corner of the stand-in plant seen by the injection


def make_spec():
    pid = BlockSpec("test_hw_pid", "lock1.pid", {"kp": -2.41e-4, "ki_per_s": -703.0, "kd_s": -2.216e-8, "toggles": {"feedback": False}})
    ramp = BlockSpec("test_hw_ramp", "lock1.ramp", {"amplitude_v": 0.0, "period_s": 0.1, "offset_v": 0.0, "toggles": {"enabled": False}})
    inj = BlockSpec("test_hw_inj", "inj1", {"freq_hz": 1.0, "amplitude": 0.0, "toggles": {"enabled": False}})
    plant = BlockSpec("test_hw_plant", "plant", {"gain": 1.0})
    pipeline = SeqSpec((inj, StackSpec((pid,), (plant,)), ramp))
    return RunSpec(run=RunSettings(1e-3, F_SIM, 1024, 7), pipeline=pipeline, preset="unit_test")


def leaf_params(spec, name):
    def walk(s):
        if isinstance(s, BlockSpec):
            yield s
        else:
            for c in getattr(s, "children", ()) or (*s.stages, *s.actuators):
                yield from walk(c)

    return next(b for b in walk(spec.pipeline) if b.name == name).params


class FakeRunner:
    """Records specs; returns ref = A cos(2πft), resp = H·ref with H = 1/(1 + i f/FC)."""

    def __init__(self):
        self.specs = []

    def __call__(self, spec):
        self.specs.append(spec)
        n = spec.run.n_samples
        t = np.arange(n) / spec.run.f_sim_hz
        p = leaf_params(spec, "inj1")
        on = spec.toggles.get("inj1.enabled", False)
        f, a = p["freq_hz"], p["amplitude"] if on else 0.0
        h = 1.0 / (1.0 + 1j * f / FC)
        series = {
            "v.inj": a * np.cos(2 * np.pi * f * t),
            "v.resp": a * abs(h) * np.cos(2 * np.pi * f * t + np.angle(h)),
            "ctrl.lock1.u": np.full(n, leaf_params(spec, "lock1.pid")["kp"]),
            "lock.lock1.state": np.ones(n, dtype=np.int64),
        }
        return RunResult(spec=spec, run_dir=None, t_s=t, series=series, events=[{"kind": "test"}])


CONFIG = {
    "loops": {"lock1": {"pid_block": "lock1.pid", "device": "dsc1"}},
    "ramps": {"out1": {"block": "lock1.ramp", "device": "dsc1"}},
    "injections": {
        "plant_in": {"block": "inj1", "reference": "v.inj", "response": "v.resp", "amplitude": 0.01, "settle_s": 1e-4, "periods": 4}
    },
}


@pytest.fixture
def ctl():
    runner = FakeRunner()
    c = SimController(make_spec(), run_fn=runner)
    c.configure(CONFIG)
    c.runner = runner
    return c


def test_interface_is_abstract():
    with pytest.raises(TypeError):
        ControllerInterface()
    assert issubclass(SimController, ControllerInterface)


def test_requires_configure():
    c = SimController(make_spec(), run_fn=FakeRunner())
    with pytest.raises(RuntimeError):
        c.engage()


def test_configure_validates_contract():
    c = SimController(make_spec(), run_fn=FakeRunner())
    with pytest.raises(KeyError):
        c.configure({"loops": {"l": {"pid_block": "nope", "device": None}}})
    with pytest.raises(ValueError, match="feedback|runtime"):
        c.configure({"loops": {"l": {"pid_block": "lock1.ramp", "device": None}}})
    with pytest.raises(ValueError):
        c.configure({"loops": {}, "bogus": {}})


def test_get_params_lists_runtime_params(ctl):
    p = ctl.get_params()
    assert p["lock1.pid.ki_per_s"] == -703.0 and "plant.gain" in p and "inj1.freq_hz" in p


def test_set_params_unknown_raises(ctl):
    with pytest.raises(KeyError):
        ctl.set_params({"lock1.pid.nope": 1.0})


def test_set_params_clips_like_dsc1(ctl):
    with pytest.warns(DeviceRangeWarning):
        applied = ctl.set_params({"lock1.pid.kp": -1.0e6, "plant.gain": 3.0})
    assert applied["lock1.pid.kp"] == -1.0e5  # |p| ≤ 1e5 (DSC1 p.6)
    assert applied["plant.gain"] == 3.0
    assert ctl.get_params()["lock1.pid.kp"] == -1.0e5
    assert ctl.status()["warnings"]


def test_set_params_in_range_unchanged(ctl):
    applied = ctl.set_params({"lock1.pid.ki_per_s": -500.0})
    assert math.isclose(applied["lock1.pid.ki_per_s"], -500.0, rel_tol=1e-12)


def test_acquire_applies_params_and_engage(ctl):
    ctl.set_params({"lock1.pid.kp": -1e-3, "run.seed": 11})
    ctl.engage()
    data = ctl.acquire(2e-4, ["ctrl.lock1.u"])
    spec = ctl.runner.specs[-1]
    assert spec.toggles["lock1.pid.feedback"] is True and spec.toggles["lock1.ramp.enabled"] is False
    assert spec.run.duration_s == 2e-4 and spec.run.seed == 11
    assert leaf_params(spec, "lock1.pid")["kp"] == -1e-3
    assert set(data) == {"ctrl.lock1.u", "t_s"} and data["t_s"].size == 200
    ctl.disengage()
    ctl.acquire(1e-4, [])
    assert ctl.runner.specs[-1].toggles["lock1.pid.feedback"] is False


def test_acquire_unknown_signal(ctl):
    with pytest.raises(KeyError, match="not recorded"):
        ctl.acquire(1e-4, ["v.nope"])


def test_ramp_clipped_to_device(ctl):
    with pytest.warns(DeviceRangeWarning):
        ctl.ramp("out1", amplitude=1.0, period_s=1.0, offset=0.0, enable=True)  # DSC1 ramp ≤ 400 ms (p.7)
    ctl.acquire(1e-4, [])
    spec = ctl.runner.specs[-1]
    assert leaf_params(spec, "lock1.ramp")["period_s"] == 0.4
    assert spec.toggles["lock1.ramp.enabled"] is True
    with pytest.raises(KeyError):
        ctl.ramp("nope", amplitude=1.0, period_s=0.1, offset=0.0, enable=True)


def test_measure_transfer_function(ctl):
    f = np.array([500.0, 2.0e3, 1.0e4])
    h = ctl.measure_transfer_function("plant_in", f)
    np.testing.assert_allclose(h, 1.0 / (1.0 + 1j * f / FC), rtol=2e-3, atol=2e-3)
    # the injection is switched off again afterwards
    ctl.acquire(1e-4, [])
    assert ctl.runner.specs[-1].toggles["inj1.enabled"] is False


def test_demodulate_phase_convention():
    t = np.arange(1000) / 1.0e5
    x = 2.0 * np.cos(2 * np.pi * 1.0e3 * t + 0.3)
    assert np.isclose(demodulate(x, t, 1.0e3), 2.0 * np.exp(0.3j))


def test_status_and_log_record_are_json(ctl):
    ctl.set_params({"lock1.pid.kp": -1e-3})
    ctl.engage()
    ctl.reset_integrators()
    ctl.acquire(1e-4, ["ctrl.lock1.u"])
    st = ctl.status()
    assert st["engaged"] and st["last_run"]["locks"] == {"lock.lock1.state": 1}
    rec = ctl.log_record()
    json.dumps(st)
    json.dumps(rec)
    assert rec["devices"]["dsc1"]["meta"]["device"] == "dsc1"
    assert [o["op"] for o in rec["ops"]][:4] == ["configure", "set_params", "engage", "reset_integrators"]


@pytest.mark.integration
def test_sim_controller_with_real_runner():
    c = SimController(make_spec())
    c.configure({"loops": {"lock1": {"pid_block": "lock1.pid", "device": "dsc1"}}})
    c.engage()
    data = c.acquire(1e-4, [])
    assert "t_s" in data
