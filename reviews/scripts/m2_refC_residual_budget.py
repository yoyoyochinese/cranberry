"""Phase 5 M2 Ref C review (D67-4): non-dither residual of refC_nominal, its attribution, the axes that move
it, and the FREE-RUNNING differential phase of the model compared with Ref C Fig. 2 (right: 30 MHz, 1 km,
not locked; the trace stays near the squeezing level for the whole 1 s with ms-long excursions to
+2.5…+7 dB, i.e. |δ| ≲ 0.1–0.2 rad most of the time, bursts of 0.4–0.6 rad).

Run: uv run python reviews/scripts/m2_refC_residual_budget.py   (≈ 2–4 min)
"""
import importlib.util
import math
from pathlib import Path

import numpy as np

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refC_support", ROOT / "tests" / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

M = math.pi * 0.1 / 4.0
F30 = 3e7


def summary(res):
    nl = res.metrics.noise_levels["i.port0"]
    tot = rs.at(nl, "sigma_total_rad", F30)
    loop = rs.at(nl, "std_loop_rad", F30)
    return tot, math.sqrt(max(tot**2 - M**2 / 2, 0.0)), loop


def show(label, **kw):
    kw.setdefault("run", {"duration_s": 0.3})
    tot, rest, loop = summary(run(rs.load(**kw)))
    print(f"  {label:<44s} total {tot*1e3:7.2f} mrad  non-dither {rest*1e3:7.2f}  loop std {loop*1e3:7.2f}")
    return rest


print("1) closed loop, refC_nominal (analytic dither), 0.3 s")
show("nominal")
for name, tog in (("radiometer only", {k: v for k, v in rs.NOISE_OFF.items() if k != "qnl.radiometer_noise"}),
                  ("conduit mech only", {k: v for k, v in rs.NOISE_OFF.items() if k != "env.mech.conduit.mech"}),
                  ("spool mech only", {k: v for k, v in rs.NOISE_OFF.items() if k != "env.mech.spool.mech"}),
                  ("all noise off (dither + quantisation)", dict(rs.NOISE_OFF))):
    show(name, toggles=tog)
print("2) single axes (others nominal)")
for lab, p in (("f_u1 = 3 kHz", {"f_u1_hz": 3e3}), ("f_u1 = 1 kHz", {"f_u1_hz": 1e3}),
               ("h_ref_conduit = 40 (axis top)", {"h_ref_conduit_rad2_hz": 40.0}),
               ("mech_f_c = 10 kHz (axis top)", {"mech_f_c_hz": 1e4}),
               ("mech_gamma = 2", {"mech_gamma": 2.0}),
               ("rho_conduit = 0.99", {"rho_conduit": 0.99}),
               ("dL = 2 m, h0 = 1e4 Hz^2/Hz", {"dL_m": 2.0, "h0_hz2_per_hz": 1e4}),
               ("f_u1 = 1 kHz, mech_f_c = 10 kHz", {"f_u1_hz": 1e3, "mech_f_c_hz": 1e4})):
    show(lab, params=p)

print("3) free-running differential phase (all three PIDs feedback off), 1.2 s, compare Ref C Fig. 2 right")
OPEN = {"pid1.feedback": False, "pid2.feedback": False, "pid3.feedback": False}
REC = {"record": {"slots": ["lock.lockC.err"], "decim": 40, "antialias": "boxcar"}}
for lab, p in (("nominal fibre noise", {}), ("h_ref /10", {"h_ref_conduit_rad2_hz": 0.4, "h_ref_spool_rad2_hz": 0.004})):
    spec = rs.load(params=p, run={"duration_s": 1.2, "f_sim_hz": 1.0e6}, toggles=OPEN, measure=REC)
    res = run(spec)
    err = np.asarray(res.series["lock.lockC.err"], dtype=float)
    fr = spec.run.f_sim_hz / 40
    d = np.unwrap(2 * err) / 2
    d = d[int(0.2 * fr):]                      # last 1.0 s
    t = np.arange(d.size) / fr
    lin = np.polyval(np.polyfit(t, d, 1), t)
    spec_d = np.fft.rfft((d - lin) * np.hanning(d.size))
    f = np.fft.rfftfreq(d.size, 1 / fr)
    # rms in 1–300 Hz (≈ what a VBW-300 Hz zero-span trace can show within 1 s)
    band = (f >= 1) & (f <= 300)
    psd = 2 * np.abs(spec_d) ** 2 / (fr * np.sum(np.hanning(d.size) ** 2))
    rms_band = math.sqrt(np.sum(psd[band]) * (f[1] - f[0]))
    print(f"  {lab:<22s} peak-to-peak over 1 s {np.ptp(d):6.2f} rad, linear drift {(lin[-1]-lin[0]):+6.2f} rad, "
          f"rms after detrend {np.std(d - lin):5.2f} rad, rms 1–300 Hz {rms_band:5.2f} rad")
    # fraction of time within ±0.2 rad of the best squeezing angle (Fig. 2 right: most of the time)
    print(f"  {'':<22s} fraction of the 1 s with |δ − median| < 0.2 rad: {np.mean(np.abs(d - np.median(d)) < 0.2):.2f}")
