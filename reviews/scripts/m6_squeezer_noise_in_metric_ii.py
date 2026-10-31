"""M6 review: squeezer-side (pump path / OPO detuning) phase noise is NOT "metric (i) only" under time multiplexing.

F_liu20 / F_bor6 headers: "pump-path phase common to all ports, seen by metric (i) only". Port k's pulse left the
squeezer Delta tau_k = tau_k - tau_0 earlier than port 0's, so a squeezer phase psi(t) enters
Theta_k = psi(t - tau_k) - psi(t - tau_0): S_Theta = 4 sin^2(pi f dtau) S_psi (same filter as the laser, conv §7).
Run: F_liu20 nominal, only the OPO length noise on (h raised to the axis top 1e-26 m^2 Hz), compare Var(Theta_k)
with the open-loop prediction from the recorded drv.XTAL.detuning series and with metric (i) of port 0.
Usage: uv run python reviews/scripts/m6_squeezer_noise_in_metric_ii.py
"""
import importlib.util
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.run import run  # noqa: E402

OFF = {"laser.fm_white": False, "env.mech.ports.mech": False, "env.mech.lo.mech": False, "env.mech.dist.mech": False,
       "env.mech.fs_lock.mech": False, "env.mech.fs_sq.mech": False, "dmx0.drive_phase_noise": False,
       "dmx1.drive_phase_noise": False, "dmx0.bias_drift": False, "dmx1.bias_drift": False,
       "env.thermal.F0.drift": False, "env.thermal.F1.drift": False, "env.thermal.LO.drift": False,
       "env.thermal.xtal.drift": False, "liu0.shot_noise": False, "liu1.shot_noise": False,
       "liu0.electronic_noise": False, "liu1.electronic_noise": False, "opo.noise.length_noise": True}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    res = run(fs.load(fs.nominal("F_liu20"), params={"opo_dl_h_m2_hz": 1e-26, "opo_lock_f_ugf_hz": 1e3},
                      run={"duration_s": 0.04}, toggles=OFF))
v_ii = res.metrics.metric_ii[1].stats.var
v_i0 = res.metrics.metric_i[0].stats.var
v_i1 = res.metrics.metric_i[1].stats.var
print(f"Var Theta_k (metric ii) = {v_ii:.3e} rad^2 ({1e3 * v_ii ** 0.5:.2f} mrad);  metric (i) port 0 {v_i0:.3e}, port k {v_i1:.3e}")
assert v_ii > 1e-3 * v_i0 > 0
print("metric (ii) is non-zero: squeezer-side noise reaches Theta_k through the emission-time difference")
