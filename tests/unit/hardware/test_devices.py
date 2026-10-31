"""Device presets devices/*.yaml (PLAN §2.10, §7.3 'every value has a source'; D52, D53)."""

import math
import re

import pytest

from cranberry.hardware.devices import (
    DeviceFormatError,
    available_devices,
    device_values,
    load_device,
    sweep_of,
)

PRESETS = ("dsc1", "fpga_generic", "moku_delta")
# a citation: manual page, D-number, T-note, PLAN / conv / prompt section, or an explicit UNVERIFIED
_CITATION = re.compile(r"p\.\s?\d+|\bD\d{2}\b|\bT1\d\b|PLAN §|conv §|prompt §|UNVERIFIED|not a device fact")


def _leaves(node, path=""):
    if isinstance(node, dict) and "value" in node:
        yield path, node
    elif isinstance(node, dict):
        for k, v in node.items():
            yield from _leaves(v, f"{path}.{k}" if path else k)


def test_three_presets_available():
    assert set(PRESETS) <= set(available_devices())


@pytest.mark.parametrize("name", PRESETS)
def test_every_value_cites_a_source(name):
    tree = load_device(name)
    leaves = [lf for k, sec in tree.items() if k != "meta" for lf in _leaves(sec, k)]
    assert leaves
    for path, leaf in leaves:
        assert _CITATION.search(leaf["source"]), f"{name}:{path} source lacks a citation: {leaf['source']!r}"


def test_dsc1_values():
    tree = load_device("dsc1")
    d = device_values(tree)
    assert d["adc"]["bits"] == 16 and d["dac"]["bits"] == 12
    assert d["dac"]["range_v"] == [-4.0, 4.0] and d["adc"]["range_v"] == [-4.0, 4.0]
    assert d["adc"]["fs_hz"] == 1.0e6 and d["adc"]["avg_samples"] == 2
    assert d["dac"]["output_impedance_ohm"] == 220.0 and d["adc"]["input_impedance_ohm"] == 1.0e6
    assert d["pid"]["h_s"] == 2.0e-6 and d["pid"]["update_hz"] == 5.0e5
    assert d["pid"]["error_sign"] == "S - v"
    assert d["pid"]["ki_entry"] == "per_sample" and d["pid"]["kd_entry"] == "per_sample"  # p.20, D53
    assert "p.20" in tree["pid"]["ki_entry"]["source"] and "D53" in tree["pid"]["ki_entry"]["source"]
    assert d["pid"]["gain_abs_range"] == [1.0e-5, 1.0e5]
    assert d["loop"]["latency_s"] == 5.7e-6
    assert sweep_of(tree, "loop.latency_s") == [5.7e-6, 8.3e-6]  # D52 / PLAN §8.4
    assert d["pid"]["anti_windup"] is None and "UNVERIFIED" in tree["pid"]["anti_windup"]["source"]


def test_moku_delta_values():
    tree = load_device("moku_delta")
    d = device_values(tree)
    assert d["loop"]["latency_s"] <= 0.21e-6
    assert math.isclose(30.0 / (360.0 * d["loop"]["max_bw_hz"]), 0.208e-6, rel_tol=0.01)  # p.65, T16 §4.2
    assert d["pid"]["output_limit_v"] == [-0.5, 0.5]  # PID Controller, p.65
    assert 10.0 in d["dac"]["range_vpp_options"]  # hardware output range, p.7 / p.79
    assert "10vpp_multi_instrument" in d["pid"]["output_range_options"]
    src = tree["pid"]["output_range_options"]["source"]
    assert "p.79" in src and "UNVERIFIED" in src
    assert d["pid"]["kp_abs_range"] == [1.0e-3, 1.0e3]  # ±60 dB
    assert d["pid"]["f_i_hz_range"] == [3.125, 3.125e5] and d["pid"]["f_d_hz_range"] == [31.25, 3.125e6]


def test_fpga_generic_leaves_choices_to_presets():
    tree = load_device("fpga_generic")
    d = device_values(tree)
    assert d["pid"]["discretisation"] == "tustin" and d["pid"]["update_hz"] is None
    assert d["loop"]["latency_s"] is None and sweep_of(tree, "loop.latency_s") == [1.0e-7, 2.0e-6]


def _write(tmp_path, body):
    (tmp_path / "x.yaml").write_text("meta: {device: x}\n" + body, encoding="utf-8")
    return tmp_path


def test_value_without_source_rejected(tmp_path):
    with pytest.raises(DeviceFormatError, match="source"):
        load_device("x", _write(tmp_path, "adc: {bits: {value: 16}}\n"))


def test_bare_scalar_rejected(tmp_path):
    with pytest.raises(DeviceFormatError):
        load_device("x", _write(tmp_path, "adc: {bits: 16}\n"))


def test_unsigned_exponent_caught(tmp_path):
    with pytest.raises(DeviceFormatError, match="string"):
        load_device("x", _write(tmp_path, "adc: {fs_hz: {value: 1.0e6, source: 'p.1'}}\n"))


def test_meta_name_must_match(tmp_path):
    (tmp_path / "y.yaml").write_text("meta: {device: x}\nadc: {bits: {value: 1, source: p.1}}\n", encoding="utf-8")
    with pytest.raises(DeviceFormatError, match="meta.device"):
        load_device("y", tmp_path)


def test_missing_preset():
    with pytest.raises(FileNotFoundError):
        load_device("no_such_device")
