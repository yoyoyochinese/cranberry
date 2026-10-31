"""Canonical PID ↔ device settings (PLAN §2.10, T16 §4.3; DSC1 manual p.6, p.18–20; D53; Moku p.65)."""

import math
import warnings

import numpy as np
import pytest

from cranberry.core.block import MissingParameterError
from cranberry.hardware.devices import device_values, load_device
from cranberry.hardware.mapping import (
    CANONICAL_KEYS,
    DeviceRangeWarning,
    PidCoefficients,
    apply_device_limits,
    map_canonical,
    map_to_canonical,
)

DSC1 = device_values(load_device("dsc1"))
MOKU = device_values(load_device("moku_delta"))
FPGA = device_values(load_device("fpga_generic"))
FPGA["pid"]["update_hz"] = 10.0e6  # null in the preset (D34); a preset sets it


def close(a: PidCoefficients, b: PidCoefficients, rel=1e-12):
    assert a.double_integrator == b.double_integrator
    for k in ("kp", "ki_per_s", "kd_s", "omega_is_rad_per_s", "omega_ds_rad_per_s", "ki2_per_s2"):
        x, y = getattr(a, k), getattr(b, k)
        if x is None or y is None:
            assert x == y, k
        else:
            assert math.isclose(x, y, rel_tol=rel, abs_tol=0.0), (k, x, y)


# ------------------------------------------------------------------ DSC1


def test_d53_user_settings_map_to_documented_canonical_gains():
    """D53: p = 0.000241, i = 0.001406, d = 0.01108 on the device; units per p.20 (I = K_i·h, D = K_d/h,
    h = 2 µs) and e = S − v (p.18) ⇒ K_p = −2.41e-4, K_i = −703 /s, K_d = −2.216e-8 s."""
    settings = {"p": 0.000241, "i": 0.001406, "d": 0.01108, "p_enabled": True, "i_enabled": True, "d_enabled": True}
    c, msgs = map_to_canonical(settings, "dsc1", DSC1)
    assert not msgs
    close(c, PidCoefficients(-2.41e-4, -703.0, -2.216e-8, None, None, False))
    back, msgs = map_canonical(c, "dsc1", DSC1)
    assert not msgs
    for k in ("p", "i", "d"):
        assert math.isclose(back[k], settings[k], rel_tol=1e-12)
        assert back[f"{k}_enabled"]


def test_dsc1_sign_flip_from_error_definition():
    s, _ = map_canonical(PidCoefficients(1.0, 0.0, 0.0, None, None, False), "dsc1", DSC1)
    assert s["p"] == -1.0 and not s["i_enabled"] and not s["d_enabled"]


def test_dsc1_per_second_alternative_reading():
    dev = {**DSC1, "pid": {**DSC1["pid"], "ki_entry": "per_second", "kd_entry": "seconds"}}
    s, _ = map_canonical(PidCoefficients(-1e-3, -2.0, -3e-3, None, None, False), "dsc1", dev)
    assert (s["p"], s["i"], s["d"]) == (1e-3, 2.0, 3e-3)


@pytest.mark.parametrize(
    "c",
    [
        PidCoefficients(-2.41e-4, -703.0, -2.216e-8, None, None, False),
        PidCoefficients(0.5, 1.0e3, 0.0, None, None, False),
        PidCoefficients(0.0, -10.0, -1.0e-7, None, None, False),
    ],
)
def test_dsc1_round_trip(c):
    s, msgs = map_canonical(c, "dsc1", DSC1)
    assert not msgs
    close(map_to_canonical(s, "dsc1", DSC1)[0], c)


def test_dsc1_clipping_warns_and_returns_clipped():
    c = PidCoefficients(-1.0e6, -3.0e-6 / 2e-6, 0.0, None, None, False)  # |p| > 1e5; |i| = 3e-6 < 1e-5
    with pytest.warns(DeviceRangeWarning):
        s, msgs = map_canonical(c, "dsc1", DSC1)
    assert s["p"] == 1.0e5 and len(msgs) == 2
    assert s["i"] == 0.0 and not s["i_enabled"]  # 3e-6 rounds to "disabled" (nearer to 0 than to 1e-5)


def test_dsc1_drops_unsupported_terms_with_warning():
    c = PidCoefficients(1.0, 10.0, 0.0, 5.0, None, True, 1.0)
    with pytest.warns(DeviceRangeWarning):
        back, msgs = apply_device_limits(c, "dsc1", DSC1)
    assert len(msgs) == 2
    close(back, PidCoefficients(1.0, 10.0, 0.0, None, None, False))


def test_dsc1_mixed_signs_rejected():
    with pytest.raises(ValueError, match="sign"):
        map_canonical(PidCoefficients(1.0, -1.0, 0.0, None, None, False), "dsc1", DSC1)


def test_dsc1_disabled_term_is_zero():
    s = {"p": 1.0, "i": 5.0, "d": 1.0, "p_enabled": True, "i_enabled": False, "d_enabled": True}
    c, _ = map_to_canonical(s, "dsc1", DSC1)
    assert c.ki_per_s == 0.0


# ------------------------------------------------------------------ Moku:Delta


@pytest.mark.parametrize(
    "c",
    [
        PidCoefficients(0.5, 2.0e3, 1.0e-6, None, None, False),
        PidCoefficients(-2.0, -1.0e4, -2.0e-7, 1.0e2, 1.5e7, True, -1.0e6),
        PidCoefficients(0.0, 50.0, 0.0, None, None, False),  # P off: reference gain 1
    ],
)
def test_moku_round_trip(c):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeviceRangeWarning)
        s, msgs = map_canonical(c, "moku_delta", MOKU)
    assert not msgs
    close(map_to_canonical(s, "moku_delta", MOKU)[0], c)


def test_moku_crossover_definition():
    s, _ = map_canonical(PidCoefficients(2.0, 2.0 * 2 * math.pi * 1.0e3, 2.0 / (2 * math.pi * 1.0e5), None, None, False), "moku_delta", MOKU)
    assert math.isclose(s["kp_db"], 20 * math.log10(2.0))
    assert math.isclose(s["f_i_hz"], 1.0e3) and math.isclose(s["f_d_hz"], 1.0e5)
    s, _ = map_canonical(PidCoefficients(0.0, 2 * math.pi * 10.0, 0.0, None, None, False), "moku_delta", MOKU)
    assert not s["p_enabled"] and math.isclose(s["f_i_hz"], 10.0)


def test_moku_clipping():
    c = PidCoefficients(1.0e4, 1.0e4 * 2 * math.pi * 1.0e6, 0.0, None, None, False)  # +80 dB, f_I 1 MHz
    with pytest.warns(DeviceRangeWarning):
        s, msgs = map_canonical(c, "moku_delta", MOKU)
    assert math.isclose(s["kp_db"], 60.0) and s["f_i_hz"] == 3.125e5 and len(msgs) == 2


def test_moku_mixed_signs_rejected():
    with pytest.raises(ValueError):
        map_canonical(PidCoefficients(1.0, 1.0, -1.0e-6, None, None, False), "moku_delta", MOKU)


# ------------------------------------------------------------------ generic FPGA (Tustin)


@pytest.mark.parametrize(
    "c",
    [
        PidCoefficients(0.3, 1.0e4, 1.0e-7, 50.0, 1.0e6, True, 1.0e6),
        PidCoefficients(-0.3, -1.0e4, 0.0, None, None, False),
    ],
)
def test_fpga_round_trip(c):
    s, _ = map_canonical(c, "fpga_generic", FPGA)
    close(map_to_canonical(s, "fpga_generic", FPGA)[0], c, rel=1e-9)


def test_fpga_tustin_matches_continuous_at_low_frequency():
    c = PidCoefficients(0.3, 1.0e4, 1.0e-7, 50.0, 1.0e6, True, 1.0e6)
    s, _ = map_canonical(c, "fpga_generic", FPGA)
    f = 1.0e3
    z1 = np.exp(-2j * np.pi * f / s["fs_hz"])  # z⁻¹ at s = +i2πf

    def h(sec):
        return np.polyval(sec["b"][::-1], z1) / np.polyval(sec["a"][::-1], z1)

    total = s["kp"] + h(s["i"]) + h(s["i2"]) + h(s["d"])
    assert abs(total / c.transfer(f) - 1) < 1e-4


def test_fpga_needs_update_rate():
    with pytest.raises(MissingParameterError):
        map_canonical(PidCoefficients(1.0, 0.0, 0.0, None, None, False), "fpga_generic")


def test_fpga_pure_derivative_warns():
    with pytest.warns(DeviceRangeWarning, match="z = -1"):
        map_canonical(PidCoefficients(0.0, 0.0, 1e-6, None, None, False), "fpga_generic", FPGA)


# ------------------------------------------------------------------ canonical form


def test_params_round_trip_and_keys():
    c = PidCoefficients(1.0, 2.0, 3.0, None, 4.0, True, 5.0)
    p = c.as_params()
    assert set(p) <= set(CANONICAL_KEYS) and "omega_is_rad_per_s" not in p
    assert PidCoefficients.from_params(p) == c


def test_transfer_function_terms():
    c = PidCoefficients(2.0, 3.0, 0.0, None, None, False)
    s = 1j * 2 * np.pi * 7.0
    assert np.isclose(c.transfer(7.0), 2.0 + 3.0 / s)


def test_validation():
    with pytest.raises(ValueError):
        PidCoefficients(1.0, 0.0, 0.0, None, None, True)
    with pytest.raises(ValueError):
        PidCoefficients(1.0, 0.0, 0.0, -1.0, None, False)


def test_unknown_device():
    with pytest.raises(KeyError):
        map_canonical(PidCoefficients(1.0, 0.0, 0.0, None, None, False), "nope", {})
