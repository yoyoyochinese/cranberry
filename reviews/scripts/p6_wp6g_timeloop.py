"""Phase 6 review (WP6-G, D90): `gated_loop_stability_check` against the REAL time loop of F_liu20_nominal +
gate/borealis_window (not the toy loop of the unit test).

1. Unstable point (bad configuration, hold_mode hold: the held error keeps the integrator ramping): T = 1 ms, lock
   fraction 0.1, f_u = 1 kHz. The per-period growth of the time loop's true lock error lock.L1.err against ρ.
2. Stable point (gate/borealis_si, zero mode): the masked (gate.gw.meas) variance of lock.L1.err in a locked run against
   the check's var_meas with `disturbance` = the Welch PSD of lock.L1.err from a free-running run (pid1f feedback off,
   same seed). In-loop noises the check does not model (ADC / DAC quantisation, detector shot / electronic noise) are
   switched off in both runs, and the comparison is also made with them on.

Run: uv run python reviews/scripts/p6_wp6g_timeloop.py [--duration 0.2]
"""

from __future__ import annotations

import argparse
import dataclasses
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy.signal import welch

from cranberry import config
from cranberry.blocks.chk.gated_loop_stability import make_gated_loop_stability_check
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

NOMINAL = "presets/future/F_liu20_nominal.yaml"
REC = ["lock.L1.err", "gate.gw.meas", "gate.gw.lock"]
QUIET = {"adc1.quantization": False, "dac1f.quantization": False, "liu1.shot_noise": False,
         "liu1.electronic_noise": False, "adc0.quantization": False, "dac0f.quantization": False,
         "liu0.shot_noise": False, "liu0.electronic_noise": False}


def load(overlays, duration, params=None, toggles=None, seed=None):
    tree = {"extends": NOMINAL,
            "run": {"duration_s": duration, **({"seed": seed} if seed else {})},
            "offline": {"params": {"t_settle_s": 0.0}},
            "measure": {"record": {"slots": REC, "decim": 1, "antialias": "none"}}}
    spec = config.load(tree, overlays=list(overlays))
    if config.sweep_markers(spec) or params:  # the overlay axes are pinned AFTER the overlays
        (pt,) = expand(spec, {f"offline.params.{k}": [v] for k, v in (params or {}).items()})
        spec = pt.spec
    if toggles:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **toggles})
    return spec


def chk_params(spec):
    (b,) = [x for x in offline_block_specs(spec) if x.name == "chk.gated_stability"]
    return dict(b.params)


def hold_overlay(tmp):
    hold = Path(tmp) / "hold.yaml"
    hold.write_text("offline:\n  blocks:\n    - name: chk.gated_stability\n      params: {hold_mode: hold}\n"
                    "pipeline:\n  - name: adc1\n    seq:\n      - name: gh1\n        params: {mode: hold}\n"
                    "  - name: adc0\n    seq:\n      - name: gh0\n        params: {mode: hold}\n")
    return str(hold)


def unstable_point(hold):
    spec = load(["gate/borealis_window", hold], 12e-3,
                params={"gate_period_s": 1e-3, "gate_lock_fraction": 0.1, "f_u_hz": 1e3}, toggles=QUIET)
    p = chk_params(spec)
    r = make_gated_loop_stability_check(p, name="chk").result
    res = run(spec, mode="compiled")
    e = np.asarray(res.series["lock.L1.err"], float)
    n = int(round(p["period_s"] * p["f_sim_hz"]))
    rms = np.sqrt(np.mean(e[: e.size // n * n].reshape(-1, n) ** 2, axis=1))
    print(f"(1) hold mode T 1 ms, D 0.1, f_u 1 kHz: check ρ = {r['rho']:.4f} (mult {np.round(r['multipliers'][:3], 4)})")
    print("    time loop per-period rms of lock.L1.err:", np.array2string(rms, precision=4))
    print("    ratios:", np.array2string(rms[1:] / rms[:-1], precision=3), f"(valid {res.valid})")


def stable_point(duration, quiet, label="borealis_si", ov=("gate/borealis_window", "gate/borealis_si"), params=None):
    tg = QUIET if quiet else {}
    lk = run(load(ov, duration, params, toggles=tg), mode="compiled")
    fr = run(load(ov, duration, params, toggles={**tg, "pid1f.feedback": False, "pid0f.feedback": False,
                                                 "pid1s.feedback": False, "pid0s.feedback": False}), mode="compiled")
    spec = load(ov, duration, params)
    p = chk_params(spec)
    fs = p["f_sim_hz"]
    t = np.asarray(lk.t_s)
    keep = t >= 0.02  # after the start-up
    ef = np.asarray(fr.series["lock.L1.err"], float)[keep]
    f, s = welch(ef, fs, nperseg=min(ef.size // 4, 1 << 17), detrend="linear")
    f, s = f[1:], s[1:]
    r = make_gated_loop_stability_check({**p, "disturbance": {"f_hz": f.tolist(), "psd_rad2_per_hz": s.tolist()}},
                                        name="chk").result
    el = np.asarray(lk.series["lock.L1.err"], float)[keep]
    meas = np.asarray(lk.series["gate.gw.meas"], float)[keep] > 0.5
    el = el - el[meas].mean()
    lab = "noise / quantisation off" if quiet else "all on"
    print(f"(2) {label} ({lab}): time loop var(meas) {np.var(el[meas]):.4e}, var(lock win) "
          f"{np.var(el[~meas]):.4e} rad² | check var_meas {r['var_meas_rad2']:.4e}, var_lock {r['var_lock_rad2']:.4e}, "
          f"conv {r['var_conversion_rad2']:.3e}, hold drift {r['var_hold_drift_rad2']:.3e}, ungated {r['var_ungated_rad2']:.3e}"
          f" | metric_ii[1] std {lk.metrics.metric_ii[1].stats.std:.4e} rad (n {lk.metrics.metric_ii[1].stats.n})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=0.2)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    ensure_loaded()
    with tempfile.TemporaryDirectory() as tmp:
        hold = hold_overlay(tmp)
        unstable_point(hold)
        stable_point(a.duration, True)
        stable_point(a.duration, False)
        stable_point(a.duration, True, "zero T 1 ms D 0.1", ("gate/borealis_window",),
                     {"gate_period_s": 1e-3, "gate_lock_fraction": 0.1})
        stable_point(a.duration, True, "hold T 300 us D 0.3", ("gate/borealis_window", hold),
                     {"gate_period_s": 3e-4, "gate_lock_fraction": 0.3})


if __name__ == "__main__":
    main()
