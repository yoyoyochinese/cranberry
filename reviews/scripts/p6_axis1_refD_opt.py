"""Phase 6 axis 1 (d): BoTorch tuning of the DSC1 PID gains on refD_lab (the calibrated bench, D89/D91-2), under the
915B side-mode non-uniqueness (D96-8: the 5.4 kHz mode at Q 15 / r 0.11 and at Q 5 / r 0.34 both match the measured
|S|²).

Variables: the entered DSC1 values dsc1_p, dsc1_i, dsc1_d (DSC1 p.20 units, D53; log bounds around the D53 point).
Objective: minimise log10 of the true in-loop phase error rms lock.lockD.err (t ≥ T0), mean over models × seeds.
Constraints (worst over models × seeds): lock maintained (optimize.lock_maintained: no loss, valid — valid includes
chk.lockD_stability's loop_unstable), servo bump max |S|² ≤ S2_MAX in [F_LO, F_HI] (|S|² = PSD locked / PSD free, same
seed, explore.sensitivity_s2) and of chk.lockD_stability's |1/(1 + L)|², the larger (on log10(|S|²/S2_MAX)).
S2_MAX = 4 (6 dB): the classical robustness bound M_s ≤ 2, which guarantees a gain margin ≥ 6 dB and a phase margin
≥ 29° whatever the plant details; the measured bench sits at 13.6 (11.3 dB) at 5.3 kHz (D81-4) with a 3.3 dB gain
margin (D89-4), i.e. one r step from instability (D96-8).

Usage:
  uv run python reviews/scripts/p6_axis1_refD_opt.py opt <q15|q5|robust> [n_init n_iter]
  uv run python reviews/scripts/p6_axis1_refD_opt.py validate [n_seeds]
  uv run python reviews/scripts/p6_axis1_refD_opt.py analytic
"""
import dataclasses
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

MODELS = {"q15": {"offline.params.mm2_q": 15.0, "offline.params.mm2_r": 0.11},
          "q5": {"offline.params.mm2_q": 5.0, "offline.params.mm2_r": 0.34}}
D53 = {"offline.params.dsc1_p": 2.41e-4, "offline.params.dsc1_i": 1.406e-3, "offline.params.dsc1_d": 1.108e-2}
BOUNDS = {"offline.params.dsc1_p": (1.0e-5, 3.0e-3, "log"), "offline.params.dsc1_i": (1.0e-4, 1.0e-2, "log"),
          "offline.params.dsc1_d": (1.0e-4, 5.0e-2, "log")}
S2_MAX, F_LO, F_HI = 4.0, 1.0e3, 4.5e4
DUR, T0, NPERSEG = 0.2, 0.02, 2048
SEEDS_OPT = (1, 2)
OUT = L.RUNS


def spec(model, gains, seed, feedback=True, duration_s=DUR):
    s = L.point(**L.d_point(seed=seed, duration_s=duration_s, axes={**MODELS[model], **gains}))
    if not feedback:
        s = dataclasses.replace(s, toggles={**dict(s.toggles), "dsc1.pid.feedback": False})
    return s


def go(s):
    from cranberry.runner.run import run

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(s)


_FREE = {}


def free(model, seed):
    if (model, seed) not in _FREE:
        r = go(spec(model, D53, seed, feedback=False))
        _FREE[(model, seed)] = np.asarray(r.series["lock.lockD.err"], dtype=float)
    return _FREE[(model, seed)]


def measure(model, gains, seed):
    """{err_rms, s2_max, f_s2_max, lock c, valid, n_slip, range} of one locked run, plus (f, |S|²)."""
    from cranberry.runner.explore import run_quantities, sensitivity_s2, servo_bump
    from cranberry.runner.optimize import lock_maintained

    r = go(spec(model, gains, seed))
    t = np.asarray(r.t_s)
    e = np.asarray(r.series["lock.lockD.err"], dtype=float)
    fs = 1.0 / float(t[1] - t[0])
    f, s2 = sensitivity_s2(e, free(model, seed), fs, nperseg=NPERSEG)
    try:
        bump, fb = servo_bump(f, s2, F_LO, F_HI)
    except ValueError:
        bump, fb = float("inf"), float("nan")
    q = run_quantities(r)
    return {"err_rms": float(np.std(e[t >= T0])), "s2_max": bump, "f_s2_max": fb, "lock_c": lock_maintained(r),
            "valid": bool(r.valid), "n_slip": q.get("lock.dsc1.sup.n_slip", np.nan),
            "range": q.get("range.stretcher915b", np.nan), "events": [str(ev.get("code", ev)) for ev in r.events][:4]}, (f, s2)


def analytic(model, gains):
    """chk.lockD_stability model: stable, gain margin [dB], max |S|² in [F_LO, F_HI] and its frequency, f_u, (f, |S|²)."""
    from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop
    from cranberry.runner.run import offline_block_specs

    s = spec(model, gains, 1)
    (b,) = [x for x in offline_block_specs(s) if x.name == "chk.lockD_stability"]
    p = dict(b.params)
    ff = np.geomspace(1e-2, 1e7, 200001)
    stable, mx, f180 = nyquist_stable(ff, open_loop(ff, p), n_integrators(p))
    f = np.geomspace(100.0, 6.0e4, 6000)
    lg = open_loop(f, p)
    s2 = np.abs(1.0 / (1.0 + lg)) ** 2
    m = (f >= F_LO) & (f <= F_HI)
    i = int(np.argmax(np.where(m, s2, -np.inf)))
    iu = int(np.argmin(np.abs(np.abs(lg) - 1.0)))
    return {"stable": bool(stable), "gm_db": float(-20 * np.log10(mx)) if mx > 0 else float("inf"), "f180": float(f180),
            "s2_max": float(s2[i]), "f_s2_max": float(f[i]), "f_u": float(f[iu])}, (f, s2)


def opt(which, n_init=8, n_iter=22):
    from cranberry.runner.optimize import optimize

    models = ["q15", "q5"] if which == "robust" else [which]
    log = []

    def evaluate(x):
        rows = [measure(m, dict(x), s)[0] for m in models for s in SEEDS_OPT]
        ana = [analytic(m, dict(x))[0] for m in models]
        log.append({"x": dict(x), "rows": rows, "analytic": ana})
        # log scales (unstable points reach 1e6 in |S|² and rad in the error: keep the GP standardisation sane);
        # the bump is the larger of the time-loop and the chk.lockD_stability |S|² (conservative)
        bump = max([r["s2_max"] for r in rows] + [a["s2_max"] for a in ana])
        return float(np.log10(np.mean([r["err_rms"] for r in rows]))), {
            "lock": max(r["lock_c"] for r in rows), "bump": float(np.log10(bump / S2_MAX))}

    def cb(ev):
        print(f"  {len(log):3d} obj {10 ** ev.objective * 1e3:8.3f} mrad  c {ev.constraints}  feasible {ev.feasible}  "
              f"x {[f'{v:.3g}' for v in ev.x.values()]}  {ev.error or ''}", flush=True)

    d53 = evaluate(D53)
    print(f"D53 gains ({which}): obj {10 ** d53[0] * 1e3:.3f} mrad, c {d53[1]}", flush=True)
    res = optimize(evaluate, BOUNDS, n_init=n_init, n_iter=n_iter, minimize=True, seed=0, callback=cb)
    best = res.best
    out = {"which": which, "d53": {"objective": d53[0], "constraints": d53[1]},
           "best": None if best is None else {"x": best.x, "objective": best.objective, "constraints": best.constraints},
           "history": [{"x": e.x, "objective": e.objective, "constraints": e.constraints, "feasible": e.feasible,
                        "error": e.error} for e in res.history], "log": log}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"refD_opt_{which}.json").write_text(json.dumps(out, default=str))
    print("best", out["best"])


def validate(n_seeds=10):
    import matplotlib

    matplotlib.use("Agg")
    from cranberry.viz.plots import plot_sensitivity

    gainsets = {"D53": D53}
    for w in ("q15", "q5", "robust"):
        p = OUT / f"refD_opt_{w}.json"
        if p.exists() and json.loads(p.read_text())["best"]:
            gainsets[f"opt_{w}"] = json.loads(p.read_text())["best"]["x"]
    summ, curves = {}, {}
    for gname, g in gainsets.items():
        for m in MODELS:
            rows, s2s = [], []
            for s in range(1, n_seeds + 1):
                r, (f, s2) = measure(m, g, s)
                rows.append(r)
                s2s.append(s2)
            ana, (fa, s2a) = analytic(m, g)
            arr = lambda k: np.array([r[k] for r in rows], dtype=float)  # noqa: E731
            summ[f"{gname}|{m}"] = {"gains": g, "analytic": ana, "n": n_seeds,
                                    **{k: {"mean": float(arr(k).mean()), "std": float(arr(k).std(ddof=1)),
                                           "min": float(arr(k).min()), "max": float(arr(k).max())}
                                       for k in ("err_rms", "s2_max", "f_s2_max", "n_slip", "range")},
                                    "lock_ok": int(sum(r["lock_c"] <= 0 for r in rows)),
                                    "events": sorted({e for r in rows for e in r["events"]})}
            curves[f"{gname}|{m}"] = (f, np.median(s2s, axis=0), fa, s2a)
            e = summ[f"{gname}|{m}"]
            print(f"{gname:11s} {m:4s} err {L.fmt(e['err_rms'], 1e3)} mrad  |S|²max {L.fmt(e['s2_max'])} @ "
                  f"{e['f_s2_max']['mean'] / 1e3:.2f} kHz  lock {e['lock_ok']}/{n_seeds}  range {e['range']['max']:.3f}  "
                  f"analytic: stable {ana['stable']} GM {ana['gm_db']:.1f} dB |S|² {ana['s2_max']:.2f} @ "
                  f"{ana['f_s2_max'] / 1e3:.2f} kHz f_u {ana['f_u'] / 1e3:.2f} kHz", flush=True)
    L.FIG.mkdir(parents=True, exist_ok=True)
    (L.FIG / "refD_opt_validate.json").write_text(json.dumps(summ, indent=1, default=str))
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), layout="constrained")
    for ax, m in zip(axs, MODELS):
        plot_sensitivity({f"{g} (time loop, median of {n_seeds})": (c[0], c[1]) for g0, c in curves.items()
                          for g in [g0.split("|")[0]] if g0.endswith("|" + m)}, limit=S2_MAX, ax=ax,
                         title=f"refD_lab, 5.4 kHz side mode {m} ({'Q 15, r 0.11' if m == 'q15' else 'Q 5, r 0.34'})")
        for g0, c in curves.items():
            if g0.endswith("|" + m):
                ax.loglog(c[2], c[3], lw=0.6, ls="--", color="0.4")
        ax.set_xlim(300, 5e4)
        ax.set_ylim(1e-3, 50)
    fig.savefig(L.FIG / "d_refD_opt_sensitivity.png", dpi=130)
    print("saved", L.FIG / "d_refD_opt_sensitivity.png")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "opt":
        opt(sys.argv[2], *(int(a) for a in sys.argv[3:5]))
    elif cmd == "validate":
        validate(int(sys.argv[2]) if len(sys.argv) > 2 else 10)
    elif cmd == "analytic":
        for m in MODELS:
            print(m, analytic(m, D53)[0])
