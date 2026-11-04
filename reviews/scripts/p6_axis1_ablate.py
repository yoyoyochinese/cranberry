"""Phase 6 axis 1, §8 step 1–2: ablation at each base's headline point, then the combination of the effects that
mattered (explore.combine_effects logic, run in parallel here).

Each toggle (runner.ablate.list_toggles, every toggle currently on) is switched off alone with the SAME seed (common
random numbers, G11), so the change isolates the effect. "Mattered": |Δq| > 5 % of the baseline for the case's
quantity. The mattering toggles are then rerun with 3 more seeds (spread check), and one run with ONLY them on
(every other candidate off) gives the combined value (non-additivity = combined − baseline).

Usage: uv run python reviews/scripts/p6_axis1_ablate.py <case> [...]   (cases: A_lab, A_nom, C, P, P_lab, D)
Writes runs/p6_axis1/ablate_<case>.json and docs/figures/phase6/axis1/ablate_<case>.json (summary).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

CASES = {
    # our lab case for (a): DX1 laser + lab acoustics, 10 m fibre, LO in a matched fibre
    "A_lab": (L.a_point(10.0, dx1=True, acoustic=True), "i.port0.sigma_total_rad"),
    "A_nom": (L.a_point(10.0), "i.port0.sigma_total_rad"),
    # DX1 + lab, LO filter cavity removed (Δν_FC 10 GHz) and τ_sq = τ_lo: what is left of the laser term
    "A_nofc": (L.a_point(10.0, dx1=True, acoustic=True, offset="zero",
                         pipeline=[{"name": "fc", "params": {"dnu_fc_hz": 1.0e10}}]), "i.port0.sigma_total_rad"),
    "C": (L.c_point(10.0), "i.port0.sigma_phase_rad"),   # time loop + HF ledger, excluding the analytic 1 MHz dither
    "P": (L.p_point(10.0), "ii.port1.sigma_total_rad"),
    "P_lab": (L.p_point(10.0, lab=True), "ii.port1.sigma_total_rad"),
    "D": (L.d_point(), "series.lock.lockD.err.std"),   # true phase error, t >= 20 ms (JOBX)
}
REL = 0.05
JOBX = {"D": {"record": ["lock.lockD.err"], "t_from": 0.02}}


def toggles_of(pt):
    from cranberry.runner.ablate import list_toggles

    return list_toggles(L.point(**pt))


def main(case):
    pt, q = CASES[case]
    tg = toggles_of(pt)
    print(f"{case}: {len(tg)} toggles, quantity {q}", flush=True)
    jx = JOBX.get(case, {})
    jobs = [{"label": "baseline", "point": pt, **jx}] + [{"label": t, "point": pt, "toggles": {t: False}, **jx} for t in tg]
    rows = L.pool(jobs, f"ablate_{case}")
    base = rows[0][q]
    impact = {r["label"]: r.get(q, float("nan")) - base for r in rows[1:]}
    errs = {r["label"]: r["error"] for r in rows[1:] if "error" in r}
    mattered = sorted([t for t, d in impact.items() if abs(d) > REL * abs(base)], key=lambda t: -abs(impact[t]))
    print(f"baseline {q} = {base:.6g}; mattered (>{REL:.0%}): {len(mattered)}", flush=True)
    # spread check: seeds 2..4 for the baseline and each mattering toggle
    jobs = []
    for s in (2, 3, 4):
        p2 = {**pt, "seed": s}
        jobs.append({"label": "baseline", "point": p2, **jx})
        jobs += [{"label": t, "point": p2, "toggles": {t: False}, **jx} for t in mattered]
    rest = [t for t in tg if t not in mattered]
    for s in (1, 2, 3, 4):
        jobs.append({"label": "combined", "point": {**pt, "seed": s}, "toggles": {t: False for t in rest}, **jx})
    rows2 = L.pool(jobs, f"ablate_{case}_seeds")
    allrows = rows + rows2
    st = {lab: L.stats(allrows, lab, [q])[q] for lab in ["baseline", "combined"] + mattered}
    summary = {"case": case, "quantity": q, "point": pt, "n_toggles": len(tg), "baseline_seed1": base,
               "impact_seed1": impact, "errors": errs, "mattered": mattered,
               "stats": {k: {kk: v[kk] for kk in ("mean", "std", "min", "max", "n")} for k, v in st.items()}}
    L.FIG.mkdir(parents=True, exist_ok=True)
    (L.FIG / f"ablate_{case}.json").write_text(json.dumps(summary, indent=1, default=str))
    print(f"\n{case}: {q}  baseline {L.fmt(st['baseline'], 1e3, 3)} mrad-or-unit×1e3, combined {L.fmt(st['combined'], 1e3, 3)}")
    for t in mattered:
        print(f"  {t:55s} seed1 Δ {impact[t] * 1e3:+9.3f}   4 seeds: {L.fmt(st[t], 1e3, 3)}"
              f"  differs={L.differs(st[t], st['baseline'])}")
    if errs:
        print("errors:", errs)


if __name__ == "__main__":
    for c in sys.argv[1:]:
        main(c)
