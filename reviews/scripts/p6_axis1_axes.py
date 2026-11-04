"""Phase 6 axis 1: list the open sweep axes (!sweep / !sweep_log / !choice) of a base + overlays, and the run settings.
Usage: uv run python reviews/scripts/p6_axis1_axes.py <preset> [overlay ...]
"""
import sys
import warnings

from cranberry import config
from cranberry.runner.sweep import axis_values

warnings.simplefilter("ignore")
base, *ovs = sys.argv[1:]
s = config.load(base, overlays=ovs)
for k, v in config.sweep_markers(s).items():
    try:
        vals = axis_values(v, k)
    except Exception:  # noqa: BLE001 - print the marker as is
        vals = repr(v)
    print(k, vals)
print("duration", s.run.duration_s, "f_sim", s.run.f_sim_hz, "seed", s.run.seed, s.run.seed_mode)
