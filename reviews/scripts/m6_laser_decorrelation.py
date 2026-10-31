"""M6 review of D74-4: laser decorrelation in metric (ii) of the future presets.

For each base preset (nominal point, port k = 31, only laser.fm_white on):
  * time-loop variance of Theta_k (locks on / locks off) and the ledger variance above f_N,
  * analytic open-loop references:  V_inf = 2 pi^2 h0 dtau  (conv §7, delayed self-interference),
    V_band(f<=f_N) = int_0^fN 4 sin^2(pi f dtau) h0/f^2 df, V_hf = int_fN^fmax (same),
  * analytic in-loop fraction 2 dtau f_u (white FM, |H|^2 flat below 1/dtau),
so the total is checked against V_inf to ~1 % (not 15 %) and the lock's in-band effect is isolated.
Usage: uv run python reviews/scripts/m6_laser_decorrelation.py [bases...] [--fu]
"""
import importlib.util
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry.runner.run import run  # noqa: E402

OFF = {"laser.fm_white": True, "env.mech.ports.mech": False, "env.mech.lo.mech": False, "env.mech.dist.mech": False,
       "env.mech.fs_lock.mech": False, "env.mech.fs_sq.mech": False, "dmx0.drive_phase_noise": False,
       "dmx1.drive_phase_noise": False, "dmx0.bias_drift": False, "dmx1.bias_drift": False,
       "env.thermal.F0.drift": False, "env.thermal.F1.drift": False, "env.thermal.LO.drift": False,
       "env.thermal.xtal.drift": False, "liu0.shot_noise": False, "liu1.shot_noise": False,
       "liu0.electronic_noise": False, "liu1.electronic_noise": False}
LOCK_OFF = {"pid0f.feedback": False, "pid1f.feedback": False, "pid0s.feedback": False, "pid1s.feedback": False}


def band(h0, dtau, lo, hi):
    f = lambda x: 4 * math.sin(math.pi * x * dtau) ** 2 * h0 / x**2  # noqa: E731
    edges = np.unique(np.concatenate([[lo], np.arange(math.ceil(lo * dtau), math.floor(hi * dtau) + 1)[:2000] / dtau, [hi]]))
    edges = edges[(edges >= lo) & (edges <= hi)]
    s = sum(quad(f, a, b, limit=200)[0] for a, b in zip(edges[:-1], edges[1:]))
    if hi * dtau > 2001:   # tail beyond the explicit periods: mean of sin^2 = 1/2
        s += 2 * h0 * (1 / edges[-1] - 1 / hi)
    return s


def one(base, extra_params=None, extra_toggles=None, seed=None):
    p = {"h0_hz2_per_hz": 64.0, "port_k": 31, **(extra_params or {})}
    r = {"duration_s": 0.04}
    if seed is not None:
        r["seed"] = seed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(fs.load(fs.nominal(base), params=p, run=r, toggles={**OFF, **(extra_toggles or {})}))
    pp = fs.params(res)
    return res, pp


def main(bases):
    print(f"{'base':8s} {'dtau[us]':>8s} {'f_u':>6s} {'V_inf':>9s} {'loop_on':>9s} {'loop_off':>9s} {'ledger':>9s} "
          f"{'an_band':>9s} {'an_hf':>9s} {'on/inf':>7s} {'off/inf':>7s} {'2dtau fu':>8s} {'sigma[mrad]':>11s}")
    for base in bases:
        res, pp = one(base)
        res_off, _ = one(base, extra_toggles=LOCK_OFF)
        h0 = 64.0
        dtau = pp["tau_k_s"] - pp["tau_0_s"]
        fN = res.spec.run.f_sim_hz / 2
        vinf = 2 * math.pi**2 * h0 * dtau
        von = res.metrics.metric_ii[1].stats.var
        voff = res_off.metrics.metric_ii[1].stats.var
        hf = res.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
        an_band = band(h0, dtau, 1.0 / 0.035, fN)
        an_hf = band(h0, dtau, fN, pp["f_max_metric_hz"])
        print(f"{base:8s} {dtau*1e6:8.3f} {pp['f_u_hz']:6.0f} {vinf:9.3e} {von:9.3e} {voff:9.3e} {hf:9.3e} {an_band:9.3e} "
              f"{an_hf:9.3e} {(von+hf)/vinf:7.3f} {(voff+hf)/vinf:7.3f} {2*dtau*pp['f_u_hz']:8.4f} "
              f"{1e3*math.sqrt(von+hf):11.1f}")
    if not FU_SWEEP:
        return
    # lock bandwidth dependence (F_liu20): in-band residual (1/(1+L)) -> larger f_u removes more (D77-1)
    print("\nF_liu20, f_u sweep (locks on), (loop+ledger)/V_inf and the 2 dtau f_u prediction of the removed share;"
          "\nwith the 915B-like stretcher resonance (20 kHz, Q 10) and with the resonance toggled off:")
    for tg, lab in (({}, "resonance on "), ({"pzt0.resonance": False, "pzt1.resonance": False}, "resonance off")):
        for fu in (1e3, 1e4, 3e4):
            res, pp = one("F_liu20", {"f_u_hz": fu}, tg)
            dtau = pp["tau_k_s"] - pp["tau_0_s"]
            vinf = 2 * math.pi**2 * 64.0 * dtau
            tot = res.metrics.metric_ii[1].stats.var + res.metrics.noise_levels["ii.port1"]["hf_variance_rad2"]
            print(f"  {lab} f_u = {fu:7.0f} Hz: ratio {tot/vinf:9.4f}   1 - 2 dtau f_u = {1-2*dtau*fu:.4f}   "
                  f"valid={getattr(res, 'valid', None)}")


FU_SWEEP = "--fu" in sys.argv

if __name__ == "__main__":
    main([a for a in sys.argv[1:] if not a.startswith("--")] or ["F_liu20", "F_wg100", "F_wg6", "F_bor6"])
