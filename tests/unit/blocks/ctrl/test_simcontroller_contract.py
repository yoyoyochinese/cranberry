"""WP-F5 contract: ``pid`` / ``ramp`` / ``sine_injection`` work behind ``SimController`` (configure,
set_params with DSC1 clipping, get_params, engage, ramp, measure_transfer_function)."""

import numpy as np
import pytest

from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.pipeline import spec_from_obj
from cranberry.hardware.mapping import DeviceRangeWarning
from cranberry.hardware.sim import SimController

FS = 1.0e6


def spec():
    obj = [
        {"block": "sine_injection", "name": "inj", "params": {"in": None, "out": "v.x", "unit": "V", "freq_hz": 1e3,
                                                              "amplitude": 0.0, "toggles": {"enabled": False}}},
        {"block": "pid", "name": "pid", "params": {
            "in": "v.x", "out": "ctrl.pid.u", "device": None, "discretisation": "tustin", "anti_windup": "none",
            "kp": 0.2, "ki_per_s": 3e3, "kd_s": 0.0, "omega_is_rad_per_s": 0.0, "omega_ds_rad_per_s": 0.0,
            "ki2_per_s2": 0.0, "setpoint_v": 0.0, "offset_v": 0.0, "u_min_v": -1e3, "u_max_v": 1e3,
            "toggles": {"feedback": False, "p_term": True, "i_term": True, "d_term": True}}},
        {"block": "ramp", "name": "ramp", "params": {"out": "v.ramp", "waveform": "sawtooth", "amplitude_v": 1.0,
                                                     "period_s": 0.02, "offset_v": 0.0, "toggles": {"enabled": False}}},
    ]
    return RunSpec(run=RunSettings(duration_s=1e-3, f_sim_hz=FS, chunk_samples=1 << 14, seed=1),
                   pipeline=spec_from_obj(obj), measure={"record": {"slots": ["v.x", "ctrl.pid.u", "v.ramp"],
                                                                    "decim": 1}})


@pytest.mark.numba
def test_simcontroller_drives_b4_blocks():
    ctl = SimController(spec())
    ctl.configure({"loops": {"l1": {"pid_block": "pid", "device": "dsc1"}},
                   "ramps": {"r": {"block": "ramp", "device": "dsc1"}},
                   "injections": {"in": {"block": "inj", "reference": "v.x", "response": "ctrl.pid.u",
                                         "amplitude": 0.01, "settle_s": 2e-3, "periods": 5}}})
    assert ctl.get_params()["pid.kp"] == 0.2
    with pytest.warns(DeviceRangeWarning):
        applied = ctl.set_params({"pid.kp": 1e-7})  # below DSC1's |p| ≥ 1e-5 → disabled (0)
    assert applied["pid.kp"] == 0.0
    ctl.set_params({"pid.kp": 0.2, "pid.ki_per_s": 3e3})
    ctl.engage()
    f = np.array([500.0, 2e3])
    h = ctl.measure_transfer_function("in", f)
    np.testing.assert_allclose(h, 0.2 + 3e3 / (2j * np.pi * f), rtol=0.01)
    with pytest.warns(DeviceRangeWarning):
        ctl.ramp("r", amplitude=1.0, period_s=0.5, offset=0.0, enable=True)  # DSC1 clips the period to 0.4 s
    data = ctl.acquire(1e-3, ["v.ramp"])
    assert data["v.ramp"][0] == pytest.approx(-1.0)
    assert np.all(np.diff(data["v.ramp"]) > 0)
    ctl.disengage()
    data = ctl.acquire(1e-4, ["ctrl.pid.u"])
    np.testing.assert_array_equal(data["ctrl.pid.u"], 0.0)  # feedback off → offset only
