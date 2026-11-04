"""Phase 6 axis 1: figures and tables of notes/findings_axis1.md from the sweep summaries (viz.plots).

Usage: uv run python reviews/scripts/p6_axis1_plots.py [A] [C] [P] [Pdec]
Reads docs/figures/phase6/axis1/sweep_<case>_summary.json (and runs/p6_axis1/sweep_<case>.json for the PSDs),
writes PNGs to docs/figures/phase6/axis1/ and prints the tables (Markdown).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cranberry.viz.plots import plot_metric_ii_components, plot_psd, plot_vs_length  # noqa: E402

LEN = [1.0, 3.1622776601683795, 10.0, 31.622776601683793, 100.0]


def summ(case):
    return json.loads((L.FIG / f"sweep_{case}_summary.json").read_text())


def series(S, fmt_label, key, scale=1.0, lengths=LEN):
    m, s = [], []
    for Lm in lengths:
        e = S.get(fmt_label.format(L=f"{Lm:g}"), {}).get(key)
        m.append(np.nan if e is None else e["mean"] * scale)
        s.append(np.nan if e is None else e["std"] * scale)
    return np.array(m), np.array(s)


def psd_rows(case, labels, slot):
    rows = L.load(f"sweep_{case}")
    out = {}
    for r in rows:
        if r.get("label") in labels and r.get("psd") and slot in r["psd"]:
            f, p = r["psd"][slot]
            out[labels[r["label"]]] = (np.array(f), np.array(p))
    return out


def cell(e, scale=1e3, nd=2):
    return "—" if e is None else f"{e['mean'] * scale:.{nd}f} ± {e['std'] * scale:.{nd}f}"


def fig_A():
    S = summ("A")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.3), layout="constrained")
    for ax, route in zip(axs[:2], ("fiber_matched", "free_space")):
        ser = {f"{lz} / {ac}": series(S, f"A|{route}|{lz}|{ac}|fpga|{{L}}", "i.port0.sigma_total_rad", 1e3)
               for lz in ("koh", "dx1") for ac in ("quiet", "lab")}
        plot_vs_length(LEN, ser, ylabel="σ_total of metric (i) [mrad]", logy=True, ax=ax,
                       title=f"(a) refA_fiber_bench, LO {route}, fpga_generic")
    ser = {f"{c} / {lz}-{ac}": series(S, f"A|fiber_matched|{lz}|{ac}|{c}|{{L}}", "i.port0.sigma_total_rad", 1e3)
           for c in ("fpga", "moku", "dsc1") for lz, ac in (("koh", "quiet"), ("dx1", "lab"))}
    plot_vs_length(LEN, ser, ylabel="σ_total of metric (i) [mrad]", logy=True, ax=axs[2],
                   title="(a) controllers (LO fiber_matched)")
    fig.savefig(L.FIG / "a_sigma_vs_length.png", dpi=130)
    # squeezing V′ [dB]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    ser = {f"{route} / {lz}-{ac}": series(S, f"A|{route}|{lz}|{ac}|fpga|{{L}}", "i.port0.r_minus_db")
           for route in ("fiber_matched", "free_space") for lz, ac in (("koh", "quiet"), ("dx1", "lab"))}
    plot_vs_length(LEN, ser, ylabel="R′− (squeezing with jitter) [dB]", ax=ax, title="(a) squeezing at the BHD")
    fig.savefig(L.FIG / "a_squeezing_vs_length.png", dpi=130)
    # sq-LO offset
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    offs = ("as_built", "zero", "fc_opt", "fc_opt_m10ns", "fc_opt_p10ns")
    ser = {o: series(S, f"Aoff|fiber_matched|dx1|lab|{o}|{{L}}", "i.port0.sigma_total_rad", 1e3, [1.0, 10.0, 100.0])
           for o in offs}
    plot_vs_length([1.0, 10.0, 100.0], ser, ylabel="σ_total [mrad]", ax=ax,
                   title="(a) DX1 + lab acoustics: squeezed − LO delay (bench/refA_sq_lo_offset)")
    fig.savefig(L.FIG / "a_sq_lo_offset.png", dpi=130)
    # PSD of metric (i)
    lab = {f"A|fiber_matched|{lz}|{ac}|fpga|{Lm:g}": f"{lz}-{ac} {Lm:g} m" for lz, ac in (("koh", "quiet"), ("dx1", "lab"))
           for Lm in (1.0, 100.0)}
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    plot_psd(psd_rows("A", lab, "metric.mi.port0"), ylabel="PSD of metric (i) [rad²/Hz]", ax=ax,
             title="(a) time-loop part of metric (i) (seed 1, t ≥ 15 ms)")
    fig.savefig(L.FIG / "a_psd_metric_i.png", dpi=130)
    # table
    print("\n### (a) σ_total [mrad] / R′− [dB] / η (fpga, mean ± std over seeds)\n")
    print("| route | laser / noise | " + " | ".join(f"{Lm:g} m" for Lm in LEN) + " |")
    print("|---|---|" + "---|" * len(LEN))
    for route in ("fiber_matched", "free_space"):
        for lz in ("koh", "dx1"):
            for ac in ("quiet", "lab"):
                row = [S.get(f"A|{route}|{lz}|{ac}|fpga|{Lm:g}", {}) for Lm in LEN]
                print(f"| {route} | {lz} / {ac} | " + " | ".join(
                    f"{cell(r.get('i.port0.sigma_total_rad'))} / {r['i.port0.r_minus_db']['mean']:.2f} dB"
                    for r in row) + " |")
    print("\n### (a) controllers (fiber_matched), σ_total [mrad], max range use EOM2 / stretcher\n")
    for c in ("fpga", "moku", "dsc1"):
        for lz, ac in (("koh", "quiet"), ("dx1", "lab")):
            row = [S.get(f"A|fiber_matched|{lz}|{ac}|{c}|{Lm:g}", {}) for Lm in LEN]
            print(f"| {c} | {lz}-{ac} | " + " | ".join(
                f"{cell(r.get('i.port0.sigma_total_rad'))} (EOM {r.get('range.eom2', {}).get('max', float('nan')):.2f},"
                f" fs {r.get('range.fs2', {}).get('max', float('nan')):.2f}, valid {r.get('valid', {}).get('min', 0):.0f})"
                for r in row) + " |")
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    ser = {"Δν_FC 3.3 MHz, fc_opt": series(S, "Aoff|fiber_matched|dx1|lab|fc_opt|{L}", "i.port0.sigma_total_rad", 1e3,
                                          [1.0, 10.0, 100.0])}
    for v in ("fc33MHz_as_built", "fc33MHz_zero", "fc33MHz_fc_opt", "nofc_as_built", "nofc_zero"):
        ser[v] = series(S, f"Afc|fiber_matched|dx1|lab|{v}|{{L}}", "i.port0.sigma_total_rad", 1e3, [1.0, 10.0, 100.0])
    plot_vs_length([1.0, 10.0, 100.0], ser, ylabel="σ_total [mrad]", logy=True, ax=ax,
                   title="(a) DX1 + lab: LO filter cavity width and delay offset")
    fig.savefig(L.FIG / "a_dx1_remedies.png", dpi=130)
    print("\n### (a) DX1 remedies (fiber_matched, lab), σ_total [mrad] / R′− [dB] / valid\n")
    for v in ("fc33MHz_as_built", "fc33MHz_zero", "fc33MHz_fc_opt", "nofc_as_built", "nofc_zero"):
        row = [S.get(f"Afc|fiber_matched|dx1|lab|{v}|{Lm:g}", {}) for Lm in (1.0, 10.0, 100.0)]
        print(f"| {v} | " + " | ".join(f"{cell(r.get('i.port0.sigma_total_rad'))} / {r['i.port0.r_minus_db']['mean']:.2f} dB"
                                       f" / {r['valid']['min']:.0f}" for r in row) + " |")
    print("\n### (a) squeezed − LO offset (DX1 + lab), σ_total [mrad]\n")
    for o in offs:
        row = [S.get(f"Aoff|fiber_matched|dx1|lab|{o}|{Lm:g}", {}) for Lm in (1.0, 10.0, 100.0)]
        print(f"| {o} | " + " | ".join(f"{cell(r.get('i.port0.sigma_total_rad'))} / "
                                       f"{r['i.port0.r_minus_db']['mean']:.2f} dB" for r in row) + " |")


def fig_C():
    S = summ("C")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), layout="constrained")
    for ax, key, yl in zip(axs, ("i.port0.std_rad", "i.port0.sigma_hf_rad"),
                           ("time-loop std of metric (i) [mrad]", "σ_hf (ledger) [mrad]")):
        ser = {f"{c} / h0 {h}": series(S, f"C|{c}|{h}|{{L}}", key, 1e3) for c in ("fpga", "moku", "dsc1")
               for h in ("nom", "dx1")}
        plot_vs_length(LEN, ser, ylabel=yl, logy=key.endswith("std_rad"), ax=ax, title="(b) refC_fiber_bench")
    fig.savefig(L.FIG / "b_sigma_vs_length.png", dpi=130)
    lab = {f"C|{c}|nom|{Lm:g}": f"{c} {Lm:g} m" for c in ("fpga", "dsc1") for Lm in (1.0, 100.0)}
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    plot_psd(psd_rows("C", lab, "metric.mi.port0"), ylabel="PSD of metric (i) [rad²/Hz]", ax=ax,
             title="(b) metric (i) time-loop PSD (seed 1)")
    fig.savefig(L.FIG / "b_psd_metric_i.png", dpi=130)
    print("\n### (b) time-loop std of metric (i) [mrad] (σ_hf, η, slips)\n")
    print("| controller | h0 | " + " | ".join(f"{Lm:g} m" for Lm in LEN) + " |")
    print("|---|---|" + "---|" * len(LEN))
    for c in ("fpga", "moku", "dsc1"):
        for h in ("nom", "dx1"):
            row = [S.get(f"C|{c}|{h}|{Lm:g}", {}) for Lm in LEN]
            print(f"| {c} | {h} | " + " | ".join(
                f"{cell(r.get('i.port0.std_rad'), nd=3)} (hf {r['i.port0.sigma_hf_rad']['mean'] * 1e3:.2f}, "
                f"η {r['i.port0.eta']['mean']:.3f}, slip {r.get('lock.sup.n_slip', {}).get('max', float('nan')):.0f})"
                for r in row) + " |")


P_LOCKS = ("liu", "het", "wdm100G", "wdm1T", "wdm1T_path", "P3", "P8", "P8copol")


def fig_P():
    S = summ("P")
    fig, axs = plt.subplots(2, 3, figsize=(17, 8.5), layout="constrained")
    for ax, (c, lab) in zip(axs.flat, (("fpga", "quiet"), ("fpga", "lab"), ("dsc1", "lab"), ("moku", "lab"),
                                       ("dsc1", "quiet"), ("moku", "quiet"))):
        ser = {lk: series(S, f"P|{lk}|{c}|{lab}|{{L}}", "ii.port1.std_loop_rad", 1e3) for lk in P_LOCKS}
        plot_vs_length(LEN, ser, ylabel="time-loop std of metric (ii) [mrad]", logy=True, ax=ax,
                       title=f"(c) pulsed_wg_1port, {c}, {lab}")
    fig.savefig(L.FIG / "c_metric_ii_vs_length.png", dpi=130)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), layout="constrained")
    for ax, (c, lab) in zip(axs, (("fpga", "lab"), ("dsc1", "lab"))):
        ser = {lk: series(S, f"P|{lk}|{c}|{lab}|{{L}}", "i.port1.std_loop_rad", 1e3) for lk in P_LOCKS}
        plot_vs_length(LEN, ser, ylabel="time-loop std of metric (i), port 1 [mrad]", logy=True, ax=ax,
                       title=f"(c) metric (i) at the bench port, {c}, {lab}")
    fig.savefig(L.FIG / "c_metric_i_vs_length.png", dpi=130)
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    ser = {lk: series(S, f"P|{lk}|fpga|lab|{{L}}", "ii.port1.mean_loop_rad", 1e3) for lk in P_LOCKS}
    plot_vs_length(LEN, ser, ylabel="static mean of metric (ii) [mrad]", ax=ax,
                   title="(c) static bias 2πΔf·Δτ (off the port grid), fpga, lab")
    fig.savefig(L.FIG / "c_metric_ii_static_bias.png", dpi=130)
    lab = {f"P|{lk}|fpga|lab|100": f"{lk} 100 m" for lk in ("liu", "het", "P3", "P8")}
    fig, ax = plt.subplots(figsize=(6.4, 4.2), layout="constrained")
    plot_psd(psd_rows("P", lab, "metric.mii.port1"), ylabel="PSD of metric (ii) [rad²/Hz]", ax=ax,
             title="(c) metric (ii) time-loop PSD, fpga + lab noise (seed 1)")
    fig.savefig(L.FIG / "c_psd_metric_ii.png", dpi=130)
    print("\n### (c) metric (ii) time-loop std [mrad] (static mean in brackets; port 1 vs the 2 m port 0)\n")
    print("| lock | ctrl | noise | " + " | ".join(f"{Lm:g} m" for Lm in LEN) + " | slips (max) | valid (min) |")
    print("|---|---|---|" + "---|" * len(LEN) + "---|---|")
    for lk in P_LOCKS:
        for c in ("fpga", "dsc1", "moku"):
            for lab in ("quiet", "lab"):
                row = [S.get(f"P|{lk}|{c}|{lab}|{Lm:g}", {}) for Lm in LEN]
                slips = max(r.get("lock.sup1.n_slip", {}).get("max", np.nan) for r in row)
                valid = min(r.get("valid", {}).get("min", np.nan) for r in row)
                errs = [r.get("_errors") for r in row if r.get("_errors")]
                print(f"| {lk} | {c} | {lab} | " + " | ".join((cell(r.get("ii.port1.std_loop_rad"), nd=1) + " (" + cell(r.get("ii.port1.mean_loop_rad"), nd=0) + ")") for r in row)
                      + f" | {slips:.0f} | {valid:.0f} |" + (f" ERR {errs[0][0][:80]}" if errs else ""))


def fig_Pdec():
    S = summ("Pdec")
    comps = ["locked_path_residual", "common_points", "laser_decorrelation_oob", "squeezed_only"]
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.3), layout="constrained")
    for ax, Lm in zip(axs, ("10", "100")):
        var, spr = {}, {}
        for lk in ("liu", "het", "wdm1T", "P3", "P8"):
            for c in ("fpga", "dsc1"):
                e = S.get(f"Pdec|{lk}|{c}|lab|{Lm}", {})
                var[f"{lk}/{c}"] = {k: e.get(f"ii.port1.{k}_rad2", {}).get("mean", 0.0) for k in comps}
                spr[f"{lk}/{c}"] = {k: e.get(f"ii.port1.{k}_rad2", {}).get("std", 0.0) for k in comps}
        plot_metric_ii_components(var, spread_rad2=spr, ax=ax, title=f"(c) metric (ii) components, lab noise, {Lm} m")
    fig.savefig(L.FIG / "c_components.png", dpi=130)
    print("\n### (c) metric (ii) components [mrad rms] (lab noise)\n")
    print("| config | L | " + " | ".join(comps) + " | total | linearity |")
    print("|---|---|" + "---|" * (len(comps) + 2))
    for lk in ("liu", "het", "wdm1T", "P3", "P8"):
        for c in ("fpga", "dsc1"):
            for Lm in ("10", "100"):
                e = S.get(f"Pdec|{lk}|{c}|lab|{Lm}", {})
                vals = [e.get(f"ii.port1.{k}_rad2") for k in comps + ["total"]]
                print(f"| {lk}/{c} | {Lm} m | " + " | ".join(
                    "—" if v is None else f"{np.sqrt(max(v['mean'], 0)) * 1e3:.1f}" for v in vals)
                      + f" | {e.get('ii.port1.linearity_ratio', {}).get('mean', float('nan')):.2f} |")


if __name__ == "__main__":
    for c in sys.argv[1:] or ("A", "C", "P", "Pdec"):
        {"A": fig_A, "C": fig_C, "P": fig_P, "Pdec": fig_Pdec}[c]()
