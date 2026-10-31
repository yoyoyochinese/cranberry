"""M6 review: the ≈ 58 mrad mean of metric (ii) in every future base preset (D74-3's "σ = 101–154 mrad" is the rms
including it). Which source sets it? F_liu20 nominal, 30 ms, sources switched on one group at a time.
Usage: uv run python reviews/scripts/m6_metric_ii_mean.py
"""
import importlib.util
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.run import run  # noqa: E402

ALL = ["laser.fm_white", "env.mech.ports.mech", "env.mech.lo.mech", "env.mech.dist.mech", "env.mech.fs_lock.mech",
       "env.mech.fs_sq.mech", "dmx0.drive_phase_noise", "dmx1.drive_phase_noise", "dmx0.bias_drift", "dmx1.bias_drift",
       "env.thermal.F0.drift", "env.thermal.F1.drift", "env.thermal.LO.drift", "env.thermal.xtal.drift",
       "liu0.shot_noise", "liu1.shot_noise", "liu0.electronic_noise", "liu1.electronic_noise", "opo.noise.length_noise"]
cases = {"all off": [], "laser only": ["laser.fm_white"], "thermal F0/F1 only": ["env.thermal.F0.drift", "env.thermal.F1.drift"],
         "mech ports only": ["env.mech.ports.mech"], "fs_sq only": ["env.mech.fs_sq.mech"],
         "fs_lock only": ["env.mech.fs_lock.mech"], "dist only": ["env.mech.dist.mech"],
         "demux only": ["dmx0.drive_phase_noise", "dmx1.drive_phase_noise", "dmx0.bias_drift", "dmx1.bias_drift"],
         "opo only": ["opo.noise.length_noise"], "xtal only": ["env.thermal.xtal.drift"], "all on": ALL}
import sys  # noqa: E402
if len(sys.argv) > 1:
    cases = {k: v for k, v in cases.items() if k in sys.argv[1:]}
for lab, on in cases.items():
    tg = {t: (t in on) for t in ALL}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.03}, toggles=tg))
    st = r.metrics.metric_ii[1].stats
    e0, e1 = r.series["lock.L0.err"], r.series["lock.L1.err"]
    print(f"{lab:20s}: Theta mean {1e3*st.mean:+8.2f} mrad, std {1e3*st.std:7.2f} mrad | lock err mean L0 "
          f"{1e3*e0[len(e0)//2:].mean():+7.2f}, L1 {1e3*e1[len(e1)//2:].mean():+7.2f} mrad | "
          f"pzt0 {r.series['act.pzt0.dpsi'][-1]:+.4f}, pzt1 {r.series['act.pzt1.dpsi'][-1]:+.4f} rad, "
          f"tec0 {r.series['act.tec0.dT'][-1]:+.2e}, tec1 {r.series['act.tec1.dT'][-1]:+.2e} K")
