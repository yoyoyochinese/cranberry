"""Phase 6 axis 1: which overlay combinations of the four bases build, and which sweep axes each leaves open.
Usage: uv run python reviews/scripts/p6_axis1_variants.py
"""
import sys
import warnings

from cranberry import config
from cranberry.runner.sweep import axis_values

warnings.simplefilter("ignore")
P = "presets/bench/pulsed_wg_1port.yaml"
COMBOS = [
    (P, []),
    (P, ["probe/detect_heterodyne"]),
    (P, ["probe/sep_wdm"]),
    (P, ["probe/sep_wdm", "probe/freq_100GHz"]),
    (P, ["probe/sep_wdm", "probe/freq_1THz"]),
    (P, ["probe/sep_wdm", "probe/freq_1THz", "probe/actuator_pathlength"]),
    (P, ["probe/sep_wdm", "probe/freq_1THz", "probe/detect_heterodyne"]),
    (P, ["gate/chopper_wdm"]),
    (P, ["probe/sep_wdm", "probe/freq_1THz", "gate/chopper_wdm"]),
    (P, ["probe/two_tone"]),
    (P, ["probe/two_tone_copol"]),
    (P, ["controller/dsc1"]),
    (P, ["controller/moku_delta"]),
    (P, ["lab/acoustic_future", "lab/multimode915b_future"]),
    (P, ["probe/detect_heterodyne", "controller/moku_delta"]),
    (P, ["gate/chopper_wdm", "controller/dsc1"]),
    (P, ["probe/two_tone", "controller/dsc1"]),
    ("presets/bench/refA_fiber_bench.yaml", ["bench/refA_laser_dx1", "bench/refA_fiber_acoustic", "bench/refA_sq_lo_offset"]),
    ("presets/bench/refA_fiber_bench.yaml", ["bench/dsc1_refA"]),
    ("presets/bench/refA_fiber_bench.yaml", ["bench/moku_delta_refA"]),
    ("presets/bench/refC_fiber_bench.yaml", ["bench/dsc1_refC"]),
    ("presets/bench/refC_fiber_bench.yaml", ["bench/moku_delta_refC"]),
    ("presets/bench/refD_lab.yaml", ["ident/wire_refD", "ident/in_lock_stepped"]),
]
only = sys.argv[1:]
for base, ovs in COMBOS:
    if only and not any(o in " ".join(ovs) for o in only):
        continue
    try:
        s = config.load(base, overlays=ovs)
        ax = {}
        for k, v in config.sweep_markers(s).items():
            try:
                ax[k] = axis_values(v, k)
            except Exception:  # noqa: BLE001 - expressions are not axes
                pass
        print("OK ", base.split("/")[-1], ovs, ax)
    except Exception as e:  # noqa: BLE001 - report the failing combination
        print("ERR", base.split("/")[-1], ovs, f"{type(e).__name__}: {e}"[:300])
