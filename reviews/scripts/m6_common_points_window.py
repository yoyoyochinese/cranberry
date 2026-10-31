"""M6 review: metric (ii) component 2 (common points PRS_k, REF_k) vs the record length.

F_common: env.mech.dist on PRsrc, PRS0/1 (2 m), REF, REF0/1 (1 m), h_ref = 0.04 rad^2 Hz at 1 km, gamma = 1,
f_lo = 0.1 Hz, f_c = 1 kHz (T08 eq. 2.7 shape). The Liu lock writes -(phi_pr,k - phi_ref,k) onto the squeezed path, so
PRS_k and REF_k (not PRsrc / REF, shared) enter Theta_k one to one. Their f^-2 spectrum puts almost all the variance
below 1/T_rec for the 30 ms runs of D74-3 (variance reported with the mean removed): the decomposition sees
~3.6 mrad, the 30 ms run shows a ~52 mrad MEAN. Analytic band variance int_{1/T}^inf S df vs T, plus a 3 s run
(f_sim 100 kHz, like the slow drift tests) to confirm the growth. No thermal driver exists on these segments at all
(drivers {strain} only): their k_T = 44 rad/(m K) drift is not modelled (2 m: 88 rad/K).
Usage: uv run python reviews/scripts/m6_common_points_window.py
"""
import importlib.util
import math
import warnings
from pathlib import Path

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.run import run  # noqa: E402

H, LREF, FLO, FC = 0.04, 1000.0, 0.1, 1000.0


def S(f, L):
    return (L / LREF) * H / f**2 / ((1 + (FLO / f) ** 2) * (1 + (f / FC) ** 2))


def band(L, f1):
    return quad(lambda f: S(f, L), f1, 1e3, limit=400)[0] + quad(lambda f: S(f, L), 1e3, 5e5, limit=400)[0]


print("analytic Var Theta_k from PRS0, PRS1 (2 m) + REF0, REF1 (1 m), variance above 1/T:")
for T in (0.03, 0.3, 3.0, 30.0, 1e9):
    v = 2 * band(2.0, 1 / T) + 2 * band(1.0, 1 / T)
    print(f"  T = {T:8.3g} s: {v:.3e} rad^2 = {1e3*math.sqrt(v):6.1f} mrad")

ALL = ["laser.fm_white", "env.mech.ports.mech", "env.mech.lo.mech", "env.mech.dist.mech", "env.mech.fs_lock.mech",
       "env.mech.fs_sq.mech", "dmx0.drive_phase_noise", "dmx1.drive_phase_noise", "dmx0.bias_drift", "dmx1.bias_drift",
       "env.thermal.F0.drift", "env.thermal.F1.drift", "env.thermal.LO.drift", "env.thermal.xtal.drift",
       "liu0.shot_noise", "liu1.shot_noise", "liu0.electronic_noise", "liu1.electronic_noise", "opo.noise.length_noise",
       "dac0f.quantization", "dac1f.quantization", "adc0.quantization", "adc1.quantization"]
tg = {t: t == "env.mech.dist.mech" for t in ALL}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for T, fsim in ((0.03, 1e6), (3.0, 1e5)):
        par = {"latency_s": 1.0e-5, "tec_rate_div": 1000} if fsim < 1e6 else {}
        r = run(fs.load(fs.nominal("F_liu20"), params=par, run={"duration_s": T + 0.005, "f_sim_hz": fsim}, toggles=tg))
        st = r.metrics.metric_ii[1].stats
        print(f"sim T = {T:5.2f} s (f_sim {fsim:.0e}): std {1e3*st.std:6.2f} mrad, mean {1e3*st.mean:+7.2f} mrad, "
              f"rms {1e3*math.sqrt(st.var + st.mean**2):6.2f} mrad")
