"""Phase 6 WP6-AB (D85-1): refD_nominal + lab/acoustic + lab/multimode915b against the Ref D measurement.

CALIBRATED ON THE REF D MEASUREMENT AFTER THE BLIND COMPARISON (D81, D85-1): the measured numbers below are copied
from predictions/refD/comparison/COMPARISON.md (§2 |S|², §3 (A) in-loop and free-running phase PSDs); the frozen
prediction is not touched. The D53 gains and every refD_nominal value are unchanged; only the overlay axes are set.

1. Analytic: |S|² = |1/(1 + L)|² of chk.lockD_stability (dsc1 law, parallel-sum modes) at refD_nominal's plant gain,
   and at the plant gains of the measured setpoints S = 0, ±0.5, ±1.5, +2 V (K_d ∝ √(V0² − S²), V0 of the model).
2. Time loop (short runs through the runner, same seed locked / free-running): |S|² = PSD(locked err) / PSD(free
   err) of the true phase error lock.lockD.err, the in-loop phase PSD at 300 Hz and 1 kHz (from v.pdb / K_d, the
   measured quantity, and from lock.lockD.err), the free-running phase PSD; for S = 0, 0.5, 1.5 V.

Run: uv run python reviews/scripts/p6_ab_check.py [--duration 2.0] [--seeds 2] [--plot out.png]
"""

from __future__ import annotations

import argparse
import dataclasses
import warnings

import numpy as np
from scipy.signal import welch

from cranberry import config
from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

NOMINAL = "presets/reference/refD_nominal.yaml"
OVERLAYS = ["lab/acoustic", "lab/multimode915b"]
# the point (overlay axes only): acoustic at the calibration (scale 1, γ = 1), the |S|² fit of the modes
POINT = {"acoustic_scale": 1.0, "acoustic_gamma": 1.0, "f_res_hz": 2.52e4, "q_res": 25.0,
         "mm2_f_hz": 5.4e3, "mm2_q": 15.0, "mm2_r": 0.11, "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07,
         "mm4_f_hz": 4.03e4, "mm4_q": 100.0, "mm4_r": 0.15}
# COMPARISON.md §2 (|S|², S = 0 / ±0.5 / ±1.5 / +2 V) and §3 (A)
MEAS_S2 = {"5.3k": [13.6, 9.2, 7.0, 3.8], "24.9k": [7.0, 15.4, 13.6, 9.1], "10k": [3.0] * 4, "40k": [4.0] * 4,
           "300": [0.020, 0.013, 0.028, 0.032], "1k": [0.16, 0.20, 0.28, 0.43], "3k": [0.62, 0.62, 0.80, 1.0]}
MEAS_INLOOP = {300.0: 7e-8, 1000.0: 1.1e-7}
MEAS_FREE = {100.0: 2e-5, 300.0: 3.6e-6, 1000.0: 6.7e-7, 3000.0: 9e-8}
BANDS = {"5.3k": (5.0e3, 5.8e3), "10k": (9.3e3, 1.07e4), "24.9k": (2.2e4, 2.75e4), "40k": (3.7e4, 4.3e4)}
SETPOINTS = [0.0, 0.5, 1.5, 2.0]
F_REC = 5.0e5


def spec_at(setpoint_v=0.0, seed=None, duration_s=None, feedback=True):
    tree = {"extends": NOMINAL, "overlays": OVERLAYS,
            "measure": {"record": {"slots": ["lock.lockD.err", "v.pdb"], "decim": 2, "antialias": "fir"}}}
    r = {}
    if duration_s:
        r["duration_s"] = duration_s
    if seed is not None:
        r["seed"] = seed
    if r:
        tree["run"] = r
    (pt,) = expand(config.load(tree), {f"offline.params.{k}": [v] for k, v in {**POINT, "setpoint_v": setpoint_v}.items()})
    spec = pt.spec
    if not feedback:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, "dsc1.pid.feedback": False})
    return spec


def check_params(spec):
    (b,) = [x for x in offline_block_specs(spec) if x.name == "chk.lockD_stability"]
    return dict(b.params)


def analytic():
    p = check_params(spec_at())
    f = np.geomspace(100.0, 6.0e4, 30000)
    ff = np.geomspace(1e-2, 1e7, 400001)
    stable, mx, f180 = nyquist_stable(ff, open_loop(ff, p), n_integrators(p))
    print(f"analytic: modes {p['modes']}\n  stable {stable}, max|L| at -180° {mx:.3f} at {f180:.0f} Hz (gain margin "
          f"{-20 * np.log10(mx):.1f} dB)")
    v0 = 2.86  # model K_d at refD_nominal (COMPARISON §0)
    out = {}
    for s_v in SETPOINTS:
        q = dict(p, k0=p["k0"] * np.sqrt(v0**2 - s_v**2) / v0)
        s2 = np.abs(1.0 / (1.0 + open_loop(f, q))) ** 2
        row = {}
        for k, (lo, hi) in BANDS.items():
            m = (f > lo) & (f < hi)
            row[k] = (f[m][np.argmax(s2[m])], s2[m].max())
        for k, x in (("300", 300.0), ("1k", 1e3), ("3k", 3e3)):
            row[k] = (x, s2[np.argmin(abs(f - x))])
        out[s_v] = row
    print("  |S|² model (measured) per setpoint S [V]:")
    for k in ("300", "1k", "3k", "5.3k", "10k", "24.9k", "40k"):
        cells = "  ".join(f"{out[s][k][1]:6.3g} ({MEAS_S2[k][i]:.3g})@{out[s][k][0] / 1e3:.2f}k"
                          for i, s in enumerate(SETPOINTS))
        print(f"    {k:>6}: {cells}")
    return f, {s: np.abs(1.0 / (1.0 + open_loop(f, dict(p, k0=p["k0"] * np.sqrt(v0**2 - s**2) / v0)))) ** 2
               for s in SETPOINTS}


def psd(x, nper=1 << 15):
    return welch(x, F_REC, nperseg=nper, detrend="linear")


def time_loop(duration, seeds, setpoints):
    res = {}
    for s_v in setpoints:
        acc = {"free": [], "lock": [], "pdb": []}
        for k in range(seeds):
            seed = 20261031 + k
            lk = run(spec_at(s_v, seed, duration, True))
            fr = run(spec_at(s_v, seed, duration, False))
            m = lk.t_s >= 0.5
            f, pl = psd(np.asarray(lk.series["lock.lockD.err"])[m])
            _, pf = psd(np.asarray(fr.series["lock.lockD.err"])[m])
            _, pv = psd(np.asarray(lk.series["v.pdb"])[m])
            # fringe slope at the lock point from the locked run (v.pdb against the true error; a 1.5 s free run
            # does not sweep a whole fringe in the model, so its max − min would under-state V0)
            e, v = (np.asarray(lk.series[k])[m] for k in ("lock.lockD.err", "v.pdb"))
            kd = abs(np.cov(v, e)[0, 1] / np.var(e))
            acc["lock"].append(pl)
            acc["free"].append(pf)
            acc["pdb"].append(pv / kd**2)
            print(f"  S = {s_v} V seed {seed}: valid {lk.valid}, events {[e['code'] for e in lk.events]}, "
                  f"K_d {kd:.3f} V/rad, err std {np.std(np.asarray(lk.series['lock.lockD.err'])[m]):.4f} rad")
        res[s_v] = (f, *(np.mean(acc[k], axis=0) for k in ("lock", "free", "pdb")))
    print("time loop (mean over seeds):")
    for s_v, (f, pl, pf, pv) in res.items():
        s2 = pl / pf

        def near(x, a, rel=0.1):
            m = (f > x / (1 + rel)) & (f < x * (1 + rel))
            return float(np.median(a[m]))

        def peak(lo, hi):
            m = (f > lo) & (f < hi)
            return f[m][np.argmax(s2[m])], s2[m].max()

        print(f"  S = {s_v} V: |S|² " + ", ".join(f"{k} {peak(*b)[1]:.3g}@{peak(*b)[0] / 1e3:.2f}k"
                                                 for k, b in BANDS.items())
              + f"; 300 Hz {near(300, s2):.3g}, 1 kHz {near(1e3, s2):.3g}, 3 kHz {near(3e3, s2):.3g}")
        print("    in-loop phase PSD [rad²/Hz] (v.pdb/K_d | true err) vs measured: "
              + ", ".join(f"{x:g} Hz {near(x, pv):.3g} | {near(x, pl):.3g} vs {MEAS_INLOOP[x]:.2g} "
                          f"(×{near(x, pv) / MEAS_INLOOP[x]:.2f})" for x in MEAS_INLOOP))
        print("    free-running phase PSD vs measured: "
              + ", ".join(f"{x:g} Hz {near(x, pf):.3g} vs {MEAS_FREE[x]:.2g}" for x in MEAS_FREE))
    return res


def plot(path, fa, s2a, res):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed (uv sync --extra plot): no plot")
        return
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for s_v, s2 in s2a.items():
        ax[0].loglog(fa, s2, lw=1, label=f"analytic S = {s_v} V")
    for s_v, (f, pl, pf, _) in res.items():
        ax[0].loglog(f[1:], (pl / pf)[1:], lw=0.6, alpha=0.7, label=f"time loop S = {s_v} V")
    fm = [300, 1e3, 3e3, 5.3e3, 1e4, 2.49e4, 4e4]
    ax[0].loglog(fm, [MEAS_S2[k][0] for k in ("300", "1k", "3k", "5.3k", "10k", "24.9k", "40k")], "ko",
                 label="measured S = 0 (COMPARISON §2)")
    ax[0].set(xlim=(100, 1e5), ylim=(1e-3, 100), xlabel="f [Hz]", ylabel="|S|²", title="loop sensitivity")
    ax[0].legend(fontsize=7)
    if res:
        f, pl, pf, pv = res[min(res)]
        ax[1].loglog(f[1:], pf[1:], lw=0.8, label="free-running (model)")
        ax[1].loglog(f[1:], pv[1:], lw=0.8, label="in-loop v.pdb/K_d (model)")
        ax[1].loglog(list(MEAS_FREE), list(MEAS_FREE.values()), "ks", label="free-running measured")
        ax[1].loglog(list(MEAS_INLOOP), list(MEAS_INLOOP.values()), "k^", label="in-loop measured")
    ax[1].set(xlim=(30, 1e5), xlabel="f [Hz]", ylabel="S_φ [rad²/Hz]", title="phase PSD, S = 0")
    ax[1].legend(fontsize=7)
    fig.suptitle("refD_nominal + lab/acoustic + lab/multimode915b (calibrated after the blind comparison, D81, D85-1)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"plot: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=2.0)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--setpoints", default="0,0.5,1.5")
    ap.add_argument("--no-time-loop", action="store_true")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--set", nargs="*", default=[], help="override POINT entries, e.g. q_res=30 mm4_r=0.1")
    a = ap.parse_args()
    POINT.update({k: float(v) for k, v in (x.split("=") for x in a.set)})
    warnings.simplefilter("ignore")
    ensure_loaded()
    fa, s2a = analytic()
    res = {} if a.no_time_loop else time_loop(a.duration, a.seeds, [float(x) for x in a.setpoints.split(",")])
    if a.plot:
        plot(a.plot, fa, s2a, res)


if __name__ == "__main__":
    main()
