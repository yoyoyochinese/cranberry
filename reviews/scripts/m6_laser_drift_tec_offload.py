"""M6 review: slow laser-frequency drift + TEC offload + orthogonal (fast-axis) probe in PM fibre -> metric (ii).

The Liu lock of port k holds phi_pr,k - phi_ref,k. A laser frequency offset dnu puts 2 pi dnu (T_b + tau_k - T_ref)
~ 2 pi dnu tau_k into that error (probe crosses the port fibre once); the PZT takes it and the slow stage offloads it to
the spool TEC (F_common stack). The TEC changes the fibre temperature, which shifts the fast axis (probe) by
(k_T + |K_T|) L dT but the slow axis (squeezed) by k_T L dT: the squeezed light is left with
    |Theta_k| ~ (|K_T| / (k_T + |K_T|)) * 2 pi dnu tau_k        (~ 4.5 % of the offloaded phase; sign: phi = -2 pi int dnu,
    conv §4, so a rising nu gives a positive Theta here).
Run: F_liu20 nominal, every noise off, laser drift nu_dot (ramp), fast TEC (as test_bad_presets), orthogonal probe
vs probe/pol_same.  Expect the slope of Theta_k ~ -(|K_T|/(k_T+|K_T|)) 2 pi nu_dot tau_k once the offload has caught up.
Usage: uv run python reviews/scripts/m6_laser_drift_tec_offload.py
"""
import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.run import run  # noqa: E402

QUIET = {"laser.fm_white": False, "env.mech.ports.mech": False, "env.mech.lo.mech": False, "env.mech.dist.mech": False,
         "env.mech.fs_lock.mech": False, "env.mech.fs_sq.mech": False, "dmx0.drive_phase_noise": False,
         "dmx1.drive_phase_noise": False, "dmx0.bias_drift": False, "dmx1.bias_drift": False,
         "opo.noise.length_noise": False, "env.thermal.F0.drift": False, "env.thermal.F1.drift": False,
         "env.thermal.LO.drift": False, "env.thermal.xtal.drift": False, "laser.drift": True}
P = {"latency_s": 1.0e-5, "tec_rate_div": 1000, "tau_th_spool_s": 1.0, "tec_tau_th_s": 1.0, "f_12_hz": 0.1}
NU_DOT = 3.0e4    # Hz/s: 300 kHz after 10 s


def go(overlays=()):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(fs.load(fs.nominal("F_liu20"), params=P, run={"f_sim_hz": 1.0e5, "duration_s": 10.0}, overlays=overlays,
                           toggles=QUIET, pipeline=[{"name": "laser", "params": {"nu_dot_hz_per_s": NU_DOT}}]))


orth, same = go(), go(("probe/pol_same",))
p = fs.params(orth)
tau_k = p["tau_k_s"]
kt, KT = p["k_T_rad_per_m_K"], abs(p["K_T_pm_rad_per_m_K"])
expect = (KT / (kt + KT)) * 2 * math.pi * NU_DOT * tau_k
for lab, r in (("orthogonal (fast-axis) probe", orth), ("same-polarization probe", same)):
    th = r.series["metric.mii.port1"]
    s = fs.slope(r.t_s, th, 6.0)
    print(f"{lab:30s}: dTheta/dt = {s:+.4e} rad/s, Theta(10 s) = {th[-1]:+.4f} rad, "
          f"TEC dT(10 s) = {r.series['act.tec1.dT'][-1]:+.4e} K, PZT = {r.series['act.pzt1.dpsi'][-1]:+.3f} rad")
print(f"prediction (|K_T|/(k_T+|K_T|)) 2 pi nu_dot tau_k = {expect:+.4e} rad/s;  "
      f"=> 10 mrad needs dnu <= {0.01 / ((KT / (kt + KT)) * 2 * math.pi * tau_k):.3g} Hz of offloaded laser drift")
s_o, s_s = fs.slope(orth.t_s, orth.series["metric.mii.port1"], 6.0), fs.slope(same.t_s, same.series["metric.mii.port1"], 6.0)
assert abs(s_o - expect) < 0.3 * abs(expect), (s_o, expect)
assert abs(s_s) < 0.05 * abs(s_o)
print("self-checks passed")
