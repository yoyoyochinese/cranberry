"""Phase 6 axis 2, Procedure 5: constrained BoTorch (runner.optimize.optimize) per fibre length.

  uv run python reviews/scripts/p6_axis2_optimize.py --case liu --base F_liu20 --k 31 [--n-init 6 --n-iter 14]

Objective: metric (ii) of port k with its static mean removed, σ = √(std_loop² + ledger) [rad] (the quantity the
decomposition sums; a static per-port offset is calibrated out at every t_cal and is reported separately), mean over
3 seeds (common random numbers), 3 s runs.
Constraints (worst seed; feasible when ≤ 0):
  lock      runner.optimize.lock_maintained (no lock loss in any supervisor, run valid: no loop_unstable CHK).
  bump      max |S(f)|² of lock.L1.err (in-loop / free-running PSD, same seed; free run = every feedback toggle off)
            over 100 Hz – 100 kHz ≤ 4, i.e. |S| ≤ 2 (6 dB): a modulus margin ≥ 0.5 guarantees gain margin ≥ 6 dB and
            phase margin ≥ 29° whatever the plant's phase (robust against the 915B mode set, which the Ref D fit does not
            fix uniquely, D96-8, and against a ±1 dB Moku gain calibration error, D93-3 / D97-5), and caps the
            out-of-loop noise amplification near the bump at 6 dB.
Variables (bounds; f_u upper bound = the P7-identified stability limit f_u·10^(GM/20), runs/p6_axis2/ident.json; the 6 dB robustness is the bump constraint):
  liu      f_u_hz [1e2, f_u,max] log, f_12_hz (PZT → TEC actuator split) [3e-3, 3] log
  p2       + gate_period_s [1e-5, 1e-2] log, gate_lock_fraction [0.1, 0.9]
  p3       + gate_f_ch_hz [1e3, 3e4] log (the chopper / demodulation rate), lpf ratio gate_lpf_hz/f_ch [0.05, 0.45]
Writes runs/p6_axis2/opt_<case>_<base>_k<k>.json (every evaluation, the best feasible point).
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import math
import warnings
from pathlib import Path

import numpy as np

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

from cranberry.runner.explore import sensitivity_s2, servo_bump  # noqa: E402
from cranberry.runner.optimize import lock_maintained, optimize  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

LABEL = {"liu": "base", "p2": "P2_borealis_si", "p3": "P3_wdm_chop6k"}
SEEDS = (1, 2, 3)
BUMP_LIMIT = 4.0
FEEDBACK_OFF = {t: False for t in ("pid0f.feedback", "pid1f.feedback", "pid0s.feedback", "pid1s.feedback")}


def objective(r) -> float:
    q = r.summary()
    return math.sqrt(q["ii.port1.std_loop_rad"] ** 2 + q["ii.port1.hf_variance_rad2"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=sorted(LABEL), default="liu")
    ap.add_argument("--base", default="F_liu20")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--n-init", type=int, default=6)
    ap.add_argument("--n-iter", type=int, default=14)
    ap.add_argument("--duration", type=float, default=3.0)
    a = ap.parse_args()
    warnings.simplefilter("ignore")
    ident = json.loads((lib.OUT / "ident.json").read_text())
    fu_max = ident.get(f"{a.base}/fpga", {}).get("f_u_stab_hz", 1.0e4)   # identified stability limit (GM → 0 dB)
    if not math.isfinite(fu_max):
        fu_max = 1.0e4
    bounds = {"f_u_hz": (1.0e2, fu_max, "log"), "f_12_hz": (3.0e-3, 3.0, "log")}
    if a.case == "p2":
        bounds |= {"gate_period_s": (1.0e-5, 1.0e-2, "log"), "gate_lock_fraction": (0.1, 0.9)}
    if a.case == "p3":
        bounds |= {"gate_f_ch_hz": (1.0e3, 3.0e4, "log"), "lpf_ratio": (0.05, 0.45)}
    label = LABEL[a.case]
    free = {}

    def pins(x):
        e = {k: v for k, v in x.items() if k != "lpf_ratio"}
        if a.case == "p2":
            e |= {"win_period_s": x["gate_period_s"], "win_duty_lock": x["gate_lock_fraction"]}
        if a.case == "p3":
            e["gate_lpf_hz"] = x["lpf_ratio"] * x["gate_f_ch_hz"]
        return e

    def evaluate(x):
        objs, c_lock, c_bump = [], -math.inf, -math.inf
        for s in SEEDS:
            spec = lib.build(a.base, label, port_k=a.k, duration_s=a.duration, seed=s, extra=pins(x))
            if s not in free:  # the free-running error does not depend on the controller parameters
                free[s] = run(dataclasses.replace(spec, toggles={**spec.toggles, **FEEDBACK_OFF}))
            r = run(spec)
            objs.append(objective(r))
            c_lock = max(c_lock, lock_maintained(r))
            fs_ = 1.0 / float(np.diff(r.t_s[:2])[0])
            f, s2 = sensitivity_s2(r.series["lock.L1.err"], free[s].series["lock.L1.err"], fs_, nperseg=4096)
            c_bump = max(c_bump, servo_bump(f, s2, 1.0e2, 1.0e5)[0] - BUMP_LIMIT)
        return float(np.mean(objs)), {"lock": c_lock, "bump": c_bump}

    def cb(ev):
        print(f"{'F' if ev.feasible else '-'} obj {1e3 * ev.objective:8.2f} mrad  c {ev.constraints}  x "
              + " ".join(f"{k}={v:.4g}" for k, v in ev.x.items()) + (f"  ERR {ev.error}" if ev.error else ""), flush=True)

    res = optimize(evaluate, bounds, n_init=a.n_init, n_iter=a.n_iter, seed=0, callback=cb)
    # reference: the nominal point with the same objective / constraints (the gate overlays' own points: Borealis SI
    # 100 µs / 0.65, P3 6 kHz with the LPF at 0.3 f_ch)
    nominal = evaluate({"liu": {}, "p2": {"gate_period_s": 1.0e-4, "gate_lock_fraction": 0.65},
                        "p3": {"gate_f_ch_hz": 6.0e3, "lpf_ratio": 0.3}}[a.case])
    out = {"case": a.case, "base": a.base, "k": a.k, "bounds": {k: list(v) for k, v in bounds.items()},
           "bump_limit": BUMP_LIMIT, "seeds": SEEDS, "duration_s": a.duration,
           "nominal": {"objective": nominal[0], "constraints": nominal[1]},
           "history": [dataclasses.asdict(e) for e in res.history],
           "best": dataclasses.asdict(res.best) if res.best else None}
    (lib.OUT / f"opt_{a.case}_{a.base}_k{a.k}.json").write_text(json.dumps(out, indent=1, default=float))
    print("best", out["best"], "nominal", out["nominal"], flush=True)


if __name__ == "__main__":
    main()
