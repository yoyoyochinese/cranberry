"""Phase 6 axis 2: tables and figures from runs/p6_axis2/*.json (written by p6_axis2_run / _ablate / _optimize /
_ident / _convergence).

  uv run python reviews/scripts/p6_axis2_report.py

Outputs (committed): docs/figures/phase6/axis2/*.png, docs/figures/phase6/axis2/summary.csv (one row per job:
mean ± sample std over the seeds of every headline quantity), runs/p6_axis2/summary.md (the same as a markdown table,
for the notes). Metric (ii) here is the decomposition total √(Σ components) [mrad]: Var(Θ_k) with the static mean
removed plus the ledger rows; the static mean (a per-port calibration offset) is the separate column `mean_ii`.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

FIG = lib.ROOT / "docs" / "figures" / "phase6" / "axis2"
MR = (lambda v: 1e3 * np.sqrt(np.clip(v, 0, None)))  # rad² → mrad
CLABEL = {"locked_path_residual": "(1) locked path", "common_points": "(2) common points",
          "laser_decorrelation_oob": "(3) laser decorrelation", "squeezed_only": "(4) squeezed only"}


def headline(rows):
    """{column: stats} of one job (lib.stats over the seeds)."""
    h = {"ii_total": lib.stats(rows, "ii.port1.total_rad2", MR)}
    for c in lib.COMP:
        h[c] = lib.stats(rows, f"ii.port1.{c}_rad2", MR)
    for key, col, f in (("ii.port1.mean_loop_rad", "mean_ii", lambda v: 1e3 * v),
                        ("i.port1.sigma_total_rad", "sigma_i", lambda v: 1e3 * v),
                        ("ii.port1.v_minus_db", "Vm_ii_db", None), ("eff.port1.eta_gbs", "eta_gbs", None),
                        ("off.stray1.photons_per_mode", "stray_mode", None),
                        ("stray1.photons_per_mode", "stray_mode", None),
                        ("stray1.clicks_per_gate", "stray_clicks", None),
                        ("off.winstray1.photons_per_mode_mean", "win_mode_mean", None),
                        ("off.winstray1.clicks_per_gate_mean", "win_clicks_mean", None),
                        ("off.winstray1.photons_per_mode_max", "win_mode_max", None),
                        ("off.frefstray1.photons_per_mode_mean", "fref_mode", None),
                        ("off.frefstray1.clicks_per_gate_mean", "fref_clicks", None),
                        ("lock.sup1.n_slip", "slips", None), ("lock.sup1.n_loss", "losses", None),
                        ("lock.sup1.lock_duration_s", "lock_s", None), ("range.pzt1", "range_pzt1", None),
                        ("range.tec1", "range_tec1", None), ("valid", "valid", None), ("L_k_m", "L_k_m", None),
                        ("ii.port1.linearity_ratio", "linearity", None), ("ii.port1.decomposition_ok", "decomp_ok", None)):
        if any(key in r for r in rows):
            h[col] = lib.stats(rows, key, f)
    return h


def table(res):
    out = []
    for tag, d in res.items():
        if not tag.startswith(("S_", "L_")):
            continue
        h = headline(d["rows"])
        out.append({"tag": tag, "base": d["base"], "label": d["label"], "k": d["port_k"], "dur": d["duration_s"],
                    "n": len(d["rows"]), "events": " ".join(d["events"]), **{k: v for k, v in h.items()}})
    return out


def fmt(s, p=1):
    if not isinstance(s, dict) or not math.isfinite(s["mean"]):
        return "–"
    return f"{s['mean']:.{p}f} ± {s['std']:.{p}f}"


def write_tables(tab):
    FIG.mkdir(parents=True, exist_ok=True)
    cols = ["ii_total", *lib.COMP, "mean_ii", "sigma_i", "eta_gbs", "stray_mode", "stray_clicks", "win_mode_mean",
            "win_mode_max", "fref_mode", "slips", "losses", "lock_s", "range_pzt1", "range_tec1", "valid", "L_k_m",
            "linearity", "decomp_ok"]
    with open(FIG / "summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag", "base", "label", "k", "duration_s", "n_seeds", "chk_events"]
                   + [f"{c}_{s}" for c in cols for s in ("mean", "std")])
        for r in tab:
            w.writerow([r["tag"], r["base"], r["label"], r["k"], r["dur"], r["n"], r["events"]]
                       + [f"{r[c][s]:.6g}" if c in r else "" for c in cols for s in ("mean", "std")])
    lines = ["| base | config | k | L_k [m] | T [s] | (ii) total [mrad] | (1) | (2) | (3) | (4) | static mean [mrad] | "
             "(i) [mrad] | η_gbs | stray /mode | slips | range PZT/TEC | CHK |", "|" + "---|" * 17]
    for r in sorted(tab, key=lambda r: (r["dur"], r["base"], r["label"], r["k"])):
        lines.append(f"| {r['base']} | {r['label']}{'' if r['tag'].endswith(f'k{r['k']}') else ' ' + r['tag'].split(f'_k{r['k']}')[-1]} "
                     f"| {r['k']} | {r['L_k_m']['mean']:.0f} | {r['dur']:g} | {fmt(r['ii_total'])} | "
                     + " | ".join(fmt(r[c]) for c in lib.COMP)
                     + f" | {fmt(r.get('mean_ii'))} | {fmt(r.get('sigma_i'))} | {fmt(r.get('eta_gbs'), 3)} | "
                     f"{r['stray_mode']['mean']:.2g} | {r['slips']['mean']:.1f} | "
                     f"{r['range_pzt1']['mean']:.2f}/{r['range_tec1']['mean']:.2g} | {r['events'] or '–'} |"
                     if all(k in r for k in ('stray_mode', 'slips', 'range_pzt1', 'range_tec1')) else
                     f"| {r['base']} | {r['label']} | {r['k']} | – | {r['dur']:g} | {fmt(r['ii_total'])} | (incomplete row) |")
    (lib.OUT / "summary.md").write_text("\n".join(lines) + "\n")


def fig_components_vs_length(tab):
    import matplotlib.pyplot as plt

    labels = ("base", "P2_borealis_si", "P3_wdm_chop6k", "P4_wdm_R0.3", "P8_two_tone", "P8_intermittent")
    for base in ("F_liu20", "F_wg6"):
        fig, axes = plt.subplots(1, 5, figsize=(17, 3.8), layout="constrained", sharex=True)
        for j, lab in enumerate(labels):
            rr = sorted([r for r in tab if r["base"] == base and r["label"] == lab and r["dur"] < 10
                         and r["tag"].endswith(f"k{r['k']}")], key=lambda r: r["L_k_m"]["mean"])
            if not rr:
                continue
            L = [r["L_k_m"]["mean"] for r in rr]
            for i, c in enumerate(("ii_total", *lib.COMP)):
                axes[i].errorbar(L, [r[c]["mean"] for r in rr], yerr=[r[c]["std"] for r in rr], marker="o", ms=3,
                                 capsize=2, lw=1, label=lab, color=f"C{j}")
        for ax, c in zip(axes, ("ii_total", *lib.COMP)):
            ax.set_title("metric (ii) total" if c == "ii_total" else CLABEL[c], fontsize=9)
            if c in ("ii_total", "locked_path_residual"):
                ax.set_yscale("log")
            ax.set_xscale("log")
            ax.set_xlabel("L_k [m]")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("mrad rms (10 seeds, ±1σ)")
        axes[-1].legend(fontsize=7)
        fig.suptitle(f"{base}: metric (ii) and its four components vs port fibre length (1 s runs; P8 intermittent 3 s)", fontsize=10)
        fig.savefig(FIG / f"components_vs_length_{base}.png", dpi=130)
        plt.close(fig)


def fig_stray_vs_gate(tab, res):
    import matplotlib.pyplot as plt

    rr = [r for r in tab if r["label"] == "P2_borealis_si" and "_T" in r["tag"]]
    if not rr:
        return
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8), layout="constrained")
    for d, mk in ((0.35, "o"), (0.65, "s")):
        sel = sorted([r for r in rr if r["tag"].endswith(f"_d{d:g}")], key=lambda r: float(r["tag"].split("_T")[1].split("_")[0]))
        T = [float(r["tag"].split("_T")[1].split("_")[0]) for r in sel]
        ax[0].loglog(T, [r["win_mode_mean"]["mean"] for r in sel], marker=mk, label=f"window mean, lock fraction {d}")
        ax[0].loglog(T, [r["win_mode_max"]["mean"] for r in sel], marker=mk, ls="--", label=f"worst slot, {d}")
        ax[1].errorbar(T, [r["locked_path_residual"]["mean"] for r in sel], [r["locked_path_residual"]["std"] for r in sel],
                       marker=mk, capsize=2, label=f"lock fraction {d}")
        ax[2].errorbar(T, [r["ii_total"]["mean"] for r in sel], [r["ii_total"]["std"] for r in sel], marker=mk, capsize=2,
                       label=f"lock fraction {d}")
    b = next((r for r in tab if r["tag"] == "S_F_liu20_base_k31_T3"), None)   # same 3 s run length
    if b:
        ax[0].axhline(b["stray_mode"]["mean"], color="k", lw=0.8, label="ungated Liu probe (CW on, 3 s)")
        ax[1].axhline(b["locked_path_residual"]["mean"], color="k", lw=0.8, label="ungated")
        ax[2].axhline(b["ii_total"]["mean"], color="k", lw=0.8, label="ungated")
    ax[0].axhline(5.8e-4, color="r", lw=0.8, ls=":", label="T17 limit 5.8e-4")
    for a_, t in zip(ax, ("probe stray photons per mode in the window", "(1) locked-path residual [mrad]",
                          "metric (ii) total [mrad]")):
        a_.set_xscale("log")
        a_.set_xlabel("window period T [s]")
        a_.set_title(t, fontsize=9)
        a_.grid(alpha=0.3)
        a_.legend(fontsize=6)
    fig.suptitle("P2 (Borealis windows) on F_liu20 port 31: stray light and phase error vs window period", fontsize=10)
    fig.savefig(FIG / "p2_stray_vs_gate_period.png", dpi=130)
    plt.close(fig)


def fig_sweeps(tab):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for key, a_, xl in (("_fch", ax[0], "P3 chopper rate f_ch [Hz]"), ("P8_intermittent", ax[1], "P8 hold cycle T_c [s]")):
        if key == "_fch":
            sel = sorted([r for r in tab if key in r["tag"]], key=lambda r: float(r["tag"].split("_fch")[1]))
            x = [float(r["tag"].split("_fch")[1]) for r in sel]
        else:
            sel = [r for r in tab if r["base"] == "F_liu20" and r["label"] == key and r["k"] == 31 and r["dur"] < 10]
            x = [float(r["tag"].split("_T")[1]) if "_T" in r["tag"] else 0.1 for r in sel]
            o = np.argsort(x)
            sel, x = [sel[i] for i in o], [x[i] for i in o]
        for c in ("ii_total", "locked_path_residual"):
            a_.errorbar(x, [r[c]["mean"] for r in sel], [r[c]["std"] for r in sel], marker="o", capsize=2,
                        label="metric (ii) total" if c == "ii_total" else CLABEL[c])
        a_.set_xscale("log")
        a_.set_yscale("log")
        a_.set_xlabel(xl)
        a_.set_ylabel("mrad rms")
        a_.grid(alpha=0.3)
        a_.legend(fontsize=7)
    fig.suptitle("F_liu20 port 31 (3 s, 10 seeds): gate-rate windows of P3 and P8", fontsize=10)
    fig.savefig(FIG / "p3_p8_gate_rate.png", dpi=130)
    plt.close(fig)


def fig_ranking(tab):
    import matplotlib.pyplot as plt

    for dur, name in ((1.0, "screen_k31"), (30.0, "long_k31")):
        for base in ("F_liu20", "F_wg6", "F_bor6", "F_wg100"):
            rr = [r for r in tab if r["base"] == base and r["k"] == 31 and r["dur"] == dur and r["tag"].endswith("k31")
                  and "ii_total" in r and math.isfinite(r["ii_total"]["mean"])]
            b = next((r for r in rr if r["label"] == "base"), None)
            if not b or len(rr) < 3:
                continue
            rr.sort(key=lambda r: -r["ii_total"]["mean"])
            fig, ax = plt.subplots(1, 3, figsize=(15, 0.3 * len(rr) + 1.6), layout="constrained", sharey=True)
            y = np.arange(len(rr))
            for c, mk in (("ii_total", "o"), *zip(lib.COMP, ("s", "^", "v", "D"))):
                ax[0].errorbar([r[c]["mean"] for r in rr], y + (0 if c == "ii_total" else 0.0),
                               xerr=[r[c]["std"] for r in rr], fmt=mk, ms=4 if c == "ii_total" else 3, capsize=2,
                               color="k" if c == "ii_total" else None, label="total" if c == "ii_total" else CLABEL[c])
            ax[0].axvline(b["ii_total"]["mean"], color="k", lw=0.8, ls="--")
            ax[0].set_xscale("log")
            ax[0].set_yticks(y, [r["label"] for r in rr], fontsize=7)
            ax[0].set_xlabel("metric (ii) and components [mrad rms] (±1σ over seeds; dashed: Liu baseline)")
            ax[0].legend(fontsize=6)
            ax[0].grid(alpha=0.3)
            st = [r["win_mode_mean"]["mean"] if "win_mode_mean" in r else r.get("stray_mode", {"mean": np.nan})["mean"]
                  for r in rr]
            ax[1].barh(y, st, color="C7")
            ax[1].axvline(5.8e-4, color="r", lw=0.8, ls=":")
            ax[1].set_xscale("log")
            ax[1].set_xlabel("probe stray photons / mode (window mean if gated; red: T17 limit)")
            ax[2].barh(y, [r["eta_gbs"]["mean"] for r in rr], color="C2")
            ax[2].set_xlim(0.3, 0.75)
            ax[2].set_xlabel("η_k at the GBS input")
            fig.suptitle(f"{base} port 31, {dur:g} s runs: schemes ranked by metric (ii)", fontsize=10)
            fig.savefig(FIG / f"ranking_{name}_{base}.png", dpi=130)
            plt.close(fig)


def fig_botorch(res):
    import matplotlib.pyplot as plt

    opts = {k: v for k, v in res.items() if k.startswith("opt_")}
    if not opts:
        return
    fig, ax = plt.subplots(1, 1, figsize=(7, 4), layout="constrained")
    for i, (k, d) in enumerate(sorted(opts.items())):
        h = d["history"]
        y = np.array([e["objective"] * 1e3 if e["feasible"] else np.nan for e in h])
        best = np.fmin.accumulate(np.where(np.isfinite(y), y, np.inf))
        ax.plot(np.arange(1, len(h) + 1), np.where(np.isfinite(best), best, np.nan), color=f"C{i}", label=k[4:])
        ax.plot(np.arange(1, len(h) + 1), y, "o", ms=3, color=f"C{i}", alpha=0.5)
        ax.plot([len(h) + 0.5], [d["nominal"]["objective"] * 1e3], "x", color=f"C{i}", ms=8)
    ax.set_xlabel("evaluation")
    ax.set_ylabel("metric (ii), mean removed [mrad] (3 seeds)")
    ax.set_title("BoTorch traces (dots: feasible evaluations; line: best so far; ×: nominal point)", fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.savefig(FIG / "botorch_traces.png", dpi=130)
    plt.close(fig)


def fig_ablation(res):
    import matplotlib.pyplot as plt

    abl = {k: v for k, v in res.items() if k.startswith("ablate_")}
    if not abl:
        return
    fig, axes = plt.subplots(1, len(abl), figsize=(4.2 * len(abl), 6), layout="constrained")
    axes = np.atleast_1d(axes)
    for ax, (k, d) in zip(axes, sorted(abl.items())):
        imp = {t: v["ii.port1.sigma_total_rad"] * 1e3 for t, v in d["impact"].items()
               if math.isfinite(v["ii.port1.sigma_total_rad"])}
        top = sorted(imp, key=lambda t: -abs(imp[t]))[:15]
        ax.barh(range(len(top)), [imp[t] for t in top], color=["C3" if imp[t] > 0 else "C0" for t in top])
        ax.set_yticks(range(len(top)), top, fontsize=6)
        ax.invert_yaxis()
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_xlabel("Δ metric (ii) σ_total when switched off [mrad]")
        ax.set_title(f"{d['base']} (baseline {1e3 * d['baseline']['ii.port1.sigma_total_rad']:.0f} mrad)", fontsize=9)
        ax.grid(alpha=0.3)
    fig.suptitle("Ablation (one toggle off at a time, port 31, 3 s, seed 1)", fontsize=10)
    fig.savefig(FIG / "ablation.png", dpi=130)
    plt.close(fig)


def main():
    import matplotlib

    matplotlib.use("Agg")
    res = lib.load_results()
    tab = table(res)
    write_tables(tab)
    for f in (fig_components_vs_length, fig_ranking):
        f(tab)
    fig_stray_vs_gate(tab, res)
    fig_sweeps(tab)
    fig_botorch(res)
    fig_ablation(res)
    print(f"{len(tab)} jobs tabulated → {FIG}")


if __name__ == "__main__":
    main()
