"""Phase 6 fix of reviews/phase6_wp6g_ab.md §5 (G-3..G-6, B-7): every ordered pair of the overlays {controller/*,
gate/*, lab/*_future, probe/detect_*, probe/freq_*} on F_liu20_nominal either builds with the right physics, or
config.load raises a ConfigError that names BOTH overlays.

"Right physics" (checked on the first point of every axis):
* the last controller overlay defines every controller choice, and the ADCs (also the ones a gate overlay restates),
  the fine PIDs (law, update rate), the DAC, chk.lock_stability and chk.gated_stability all follow it (G-3, D91, B-7);
* one stretcher mode list feeds pzt0 / pzt1, chk.lock_stability and chk.gated_stability (G-4);
* probe line (freq_*) and heterodyne beat (detect_heterodyne) are separate keys: probe offset = line + beat, the
  reference on the line (G-5);
* the controller-pole / extra-pole lists of both checks agree.
"""

from __future__ import annotations

import itertools
import warnings

import pytest

from cranberry import config
from cranberry.config import ConfigError, input_tree
from cranberry.config.spec import spec_to_tree
from cranberry.core.pipeline import spec_to_obj
from cranberry.runner.run import offline_block_specs
from cranberry.runner.sweep import axis_values, expand, find_markers, path_str

pytestmark = pytest.mark.integration
NOMINAL = "presets/future/F_liu20_nominal.yaml"
OVERLAYS = ["controller/dsc1", "controller/moku_delta", "controller/fpga_generic",
            "gate/borealis_window", "gate/chopper_wdm", "gate/intermittent_slow",
            "lab/acoustic_future", "lab/multimode915b_future",
            "probe/detect_heterodyne", "probe/detect_dither", "probe/detect_after_demux",
            "probe/freq_80MHz", "probe/freq_7GHz", "probe/freq_100GHz", "probe/freq_1THz"]
CTRL = {  # what each controller overlay must leave in the spec
    "controller/dsc1": {"ctrl_device": "dsc1", "pid_law": "dsc1", "chk_law": "dsc1", "pid_rate_div": 2,
                        "adc_range_v": 4.0, "adc_avg_samples": 2, "ctrl_dac_bits": 12, "out_range_v": 4.0,
                        "ctrl_latency_s": 5.7e-6, "ctrl_in_pole_hz": 1.2e5, "ctrl_out_pole_hz": 1.0e5, "zoh_delay_s": 0.0},
    "controller/moku_delta": {"ctrl_device": "moku_delta", "pid_law": "tustin", "chk_law": "continuous",
                              "pid_rate_div": 1, "adc_range_v": 0.5, "adc_avg_samples": 1, "ctrl_dac_bits": 14,
                              "out_range_v": 0.5, "ctrl_in_pole_hz": float("inf"), "zoh_delay_s": 0.5e-6},
    "controller/fpga_generic": {"ctrl_device": "fpga_generic", "pid_law": "tustin", "chk_law": "continuous",
                                "pid_rate_div": 1, "adc_range_v": 0.5, "adc_avg_samples": 1, "out_range_v": 0.5,
                                "ctrl_in_pole_hz": float("inf"), "zoh_delay_s": 0.5e-6},
}
LINE = {"probe/freq_80MHz": 8e7, "probe/freq_7GHz": 7e9, "probe/freq_100GHz": 1e11, "probe/freq_1THz": 1e12}
_CACHE: dict = {}


def blocks(node, out):
    if isinstance(node, dict):
        if "block" in node and "name" in node:
            out[node["name"]] = node
        for v in node.values():
            blocks(v, out)
    elif isinstance(node, list):
        for v in node:
            blocks(v, out)
    return out


def first_point(spec):
    if not config.sweep_markers(spec):
        return spec
    tree = input_tree(spec) or spec_to_tree(spec)
    axes = {path_str(tree, q): [axis_values(m, path_str(tree, q))[0]] for q, m in find_markers(tree)}
    for k in list(axes):  # D91 / B-7: the 5.7 µs + poles DSC1 model
        if k.endswith("dsc1_latency_s"):
            axes[k] = [5.7e-6]
    (pt,) = expand(spec, axes)
    return pt.spec


def load(overlays):
    key = tuple(overlays)
    if key not in _CACHE:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                _CACHE[key] = first_point(config.load({"extends": NOMINAL}, overlays=list(overlays)))
            except ConfigError as e:
                _CACHE[key] = e
    return _CACHE[key]


def physics(spec, overlays):
    p = spec.offline["params"]
    b = blocks(spec_to_obj(spec.pipeline), {})
    chk = {x.name: x.params for x in offline_block_specs(spec)}
    ctrl = [o for o in overlays if o.startswith("controller/")]
    if ctrl:
        for k, v in CTRL[ctrl[-1]].items():
            assert p[k] == pytest.approx(v), (overlays, k)
    for k in "01":
        adc = b[f"adc{k}"]["params"]
        assert adc["device"] == p["adc_device"] == p["ctrl_device"]
        assert adc["range_v"] == [-p["adc_range_v"], p["adc_range_v"]] and adc["avg_samples"] == p["adc_avg_samples"]
        pid = b[f"pid{k}f"]
        assert pid["params"]["discretisation"] == p["pid_law"] and pid.get("rate_div", 1) == p["pid_rate_div"]
        assert b[f"dac{k}f"].get("rate_div", 1) == p["pid_rate_div"] and b[f"dac{k}f"]["params"]["bits"] == p["ctrl_dac_bits"]
        assert b[f"pzt{k}"]["params"]["modes"] == p["pzt_modes"]
    lsc = chk["chk.lock_stability"]
    assert lsc["discretisation"] == p["chk_law"] and lsc["modes"] == p["pzt_modes"]
    assert lsc["delay_s"] == pytest.approx(p["ctrl_latency_s"] + p["zoh_delay_s"] + p["lock_extra_delay_s"])
    poles = [p["ctrl_in_pole_hz"], p["ctrl_out_pole_hz"], p["lock_extra_pole_hz"]]
    assert lsc["poles_hz"] == poles
    if "chk.gated_stability" in chk:
        g = chk["chk.gated_stability"]
        assert g["pid_discretisation"] == p["pid_law"] and g["pid_rate_div"] == p["pid_rate_div"]
        assert g["modes"] == p["pzt_modes"] and g["poles_hz"] == poles
        assert g["delay_s"] == pytest.approx(p["ctrl_latency_s"] - p["ctrl_hold_delay_s"] + p["zoh_delay_s"]
                                             + p["lock_extra_delay_s"])
    if "lab/multimode915b_future" in overlays:
        assert len(p["pzt_modes"]) == 4 and all("r" in m for m in p["pzt_modes"])
    lines = [LINE[o] for o in overlays if o in LINE]
    assert p["probe_line_offset_hz"] == (lines[-1] if lines else 0.0)
    shift = 8e7 if "probe/detect_heterodyne" in overlays else 0.0
    assert p["liu_f_shift_hz"] == shift and p["probe_freq_offset_hz"] == pytest.approx(p["probe_line_offset_hz"] + shift)
    beams = {x["id"]: x for x in spec.topology["beams"]}
    assert beams["pr1"]["freq_offset_hz"] == pytest.approx(p["probe_freq_offset_hz"])
    if shift:
        assert beams["ref1"]["freq_offset_hz"] == pytest.approx(p["probe_line_offset_hz"])


#: pairs that must raise (both orders): two overlays that replace liu0 / adc0's input differently, or two gate schemes
def must_raise(a, b):
    det = {"probe/detect_heterodyne", "probe/detect_dither"}
    pair = {a, b}
    if pair <= det or all(o.startswith("gate/") for o in pair):
        return True
    if pair & det and (any(o.startswith("gate/") for o in pair) or "probe/detect_after_demux" in pair):
        return True
    return False


PAIRS = list(itertools.permutations(OVERLAYS, 2))


@pytest.mark.parametrize("pair", PAIRS, ids=["+".join(p) for p in PAIRS])
def test_pair_builds_or_raises_naming_both(pair):
    res = load(pair)
    if isinstance(res, ConfigError):
        msg = str(res)
        assert all(o in msg for o in pair), msg
        assert must_raise(*pair), f"unexpected conflict: {msg}"
    else:
        assert not must_raise(*pair), f"{pair} should raise"
        physics(res, pair)


@pytest.mark.parametrize("a, b", [("controller/dsc1", "gate/borealis_window"), ("controller/moku_delta", "gate/chopper_wdm"),
                                  ("lab/multimode915b_future", "gate/borealis_window"),
                                  ("probe/freq_1THz", "probe/detect_heterodyne"), ("controller/dsc1", "lab/multimode915b_future")])
def test_order_independent(a, b):
    """G-3/G-4/G-5: the two orders give the same spec (apart from the applied_overlays record)."""
    t1, t2 = (spec_to_tree(load(x)) for x in ((a, b), (b, a)))
    for t in (t1, t2):
        t.pop("applied_overlays", None)
    assert t1 == t2


def test_controller_pair_last_wins_entirely():
    t1, t2 = spec_to_tree(load(("controller/dsc1", "controller/moku_delta"))), spec_to_tree(load(("controller/moku_delta",)))
    for t in (t1, t2):
        t.pop("applied_overlays", None)
        t["run"].pop("seed")  # the per-point seed depends on the axes
        t["offline"]["params"].pop("dsc1_latency_s", None)
        t["offline"]["params"].pop("pid_update_hz", None)
    assert t1 == t2
