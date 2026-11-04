"""Phase 6 review (WP6-G D90-4/6/8, WP6-AB D89): overlay compositions on F_liu20_nominal — which fail loudly and which
silently change the physics. For every composition: does config.load + expand succeed, and what do the ADC / DAC /
PID / checks end up with.

Run: uv run python reviews/scripts/p6_wp6g_compose.py
"""

from __future__ import annotations

import warnings

from cranberry import config
from cranberry.core.pipeline import spec_to_obj
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import offline_block_specs
from cranberry.config import input_tree
from cranberry.config.spec import spec_to_tree
from cranberry.runner.sweep import axis_values, expand, find_markers, path_str

NOMINAL = "presets/future/F_liu20_nominal.yaml"
PIN = {"offline.params.gate_period_s": [1e-4], "offline.params.gate_lock_fraction": [0.65]}


def leaves(node, out):
    if isinstance(node, dict):
        if "block" in node and "name" in node:
            out[node["name"]] = node
        for v in node.values():
            leaves(v, out)
    elif isinstance(node, list):
        for v in node:
            leaves(v, out)
    return out


def first(spec):
    if not config.sweep_markers(spec):
        return spec
    tree = input_tree(spec) or spec_to_tree(spec)
    axes = {path_str(tree, q): [axis_values(m, path_str(tree, q))[0]] for q, m in find_markers(tree)}
    axes.update(PIN)  # every axis at its first value (a full expand of the multimode grid is ~10^7 points)
    pts = expand(spec, axes)
    return pts[0].spec


def show(overlays):
    print(f"--- {overlays}")
    try:
        spec = first(config.load({"extends": NOMINAL}, overlays=list(overlays)))
    except Exception as e:  # noqa: BLE001
        print(f"    FAILS LOUDLY: {type(e).__name__}: {str(e).splitlines()[0][:200]}")
        return
    blocks = leaves(spec_to_obj(spec.pipeline), {})
    p = spec.offline["params"]
    adc = blocks.get("adc1", {}).get("params", {})
    dac = blocks.get("dac1f", {}).get("params", {})
    pid = blocks.get("pid1f", {}).get("params", {})
    liu = blocks.get("liu1", {}).get("params", {})
    chk = {b.name: b.params for b in offline_block_specs(spec)}
    g = chk.get("chk.gated_stability", {})
    lsc = chk.get("chk.lock_stability", {})
    print(f"    out_range_v {p.get('out_range_v')}, latency {p.get('latency_s')}, probe_freq_offset_hz "
          f"{p.get('probe_freq_offset_hz')}")
    print(f"    adc1: device {adc.get('device')}, range {adc.get('range_v')}, bits {adc.get('bits')} | dac1f: device "
          f"{dac.get('device')}, range {dac.get('range_v')} | pid1f: device {pid.get('device')}, disc "
          f"{pid.get('discretisation')} | liu1 mode {liu.get('mode')}")
    print(f"    chk.lock_stability: disc {lsc.get('discretisation', 'continuous')}, n_modes {len(lsc.get('modes', []))}"
          f" | chk.gated_stability: {'absent' if not g else ('pid ' + g['pid_discretisation'] + ', n_modes ' + str(len(g['modes'])))}")


def main():
    warnings.simplefilter("ignore")
    ensure_loaded()
    for ov in [
        ("gate/borealis_window",),
        ("controller/dsc1",),
        ("controller/dsc1", "gate/borealis_window"),          # D90-4: the documented silent case
        ("gate/borealis_window", "controller/dsc1"),          # D90-4: should fail loudly
        ("controller/moku_delta", "gate/borealis_window"),
        ("gate/borealis_window", "controller/moku_delta"),
        ("probe/detect_heterodyne", "gate/borealis_window"),
        ("gate/borealis_window", "probe/detect_heterodyne"),
        ("probe/detect_heterodyne", "probe/freq_1THz"),       # D90-6
        ("probe/freq_1THz", "probe/detect_heterodyne"),
        ("lab/multimode915b_future", "gate/borealis_window"),  # modes of the gated check
        ("gate/borealis_window", "lab/multimode915b_future"),
    ]:
        show(ov)


if __name__ == "__main__":
    main()
