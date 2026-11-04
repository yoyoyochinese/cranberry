"""Phase 6 axis 2: the job lists of the campaign, run on a pool of at most 5 processes (the axis-1 campaign shares the
machine). Each job = one configuration × one port index × 10 seeds (common random numbers), with the four-component
decomposition of metric (ii) per seed; it writes runs/p6_axis2/<tag>.json and is skipped when that file exists.

  uv run python reviews/scripts/p6_axis2_run.py --phase screen [--procs 5] [--dry]

Phases:
  screen  3 s runs (TEC τ 10 s partly inside; component 2's band down to 1/t_cal is the env.dist.window ledger row):
          baseline length sweep (4 bases × k = 1, 4, 16, 31), lock presets, pair / lab / controller variants, the
          schemes P1–P8 alone at k = 31, scheme length sweep (F_liu20, F_wg6; k = 4, 16), gate sweeps (P2 window,
          P3 chopper rate, P8 hold period) on F_liu20 k = 31.
  long    30 s runs = t_cal_s (the window row is 0; the slow layers act in the run): the headline set on F_liu20 k = 31
          and the P5 family on F_wg6 k = 31 (D78: TEC and drift need seconds; D98: P5 at 30 s).
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEEDS = list(range(1, 11))


def _lib():
    sp = importlib.util.spec_from_file_location("p6_axis2_lib", HERE / "p6_axis2_lib.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


COMP_SEEDS = list(range(1, 6))   # seeds with the four-component reruns (the total uses all 10)


def job(base, label, k=31, dur=1.0, extra=None, suffix="", seeds=SEEDS, comp_seeds=COMP_SEEDS):
    tag = f"{'L' if dur >= 10 else 'S'}_{base}_{label}_k{k}{suffix}"
    return {"base": base, "label": label, "port_k": k, "duration_s": dur, "seeds": list(seeds), "decompose": True,
            "component_seeds": list(comp_seeds), "extra": extra or {}, "tag": tag}


def screen_jobs():
    """Screening: 1 s runs (D78 lower end; the 3 s / 10 s / 30 s check on F_liu20 base is the length artefact check),
    3 s for the window sweeps and the slow P8 cycle; 10 seeds for every quantity and the metric-(ii) total, the four
    components on seeds 1–5 (cost: 5 runs per decomposed seed)."""
    B, J = ("F_wg100", "F_wg6", "F_bor6", "F_liu20"), []
    two = ("F_liu20", "F_wg6")
    # P2 window sweep on F_liu20 k = 31 (3 s: ≥ 300 windows at 10 ms)
    for T in (1e-5, 1e-4, 1e-3, 1e-2):
        for d in (0.35, 0.65):
            J.append(job("F_liu20", "P2_borealis_si", dur=3.0, extra={"gate_period_s": T, "gate_lock_fraction": d,
                                                                      "win_period_s": T, "win_duty_lock": d},
                         suffix=f"_T{T:g}_d{d:g}"))
    # baseline length sweep and run-length check
    J += [job(b, "base", k) for b in B for k in (1, 4, 16, 31)]
    J += [job("F_liu20", "base", dur=3.0, suffix="_T3"), job("F_liu20", "base", dur=10.0, suffix="_T10",
                                                             comp_seeds=[1, 2, 3])]
    # schemes alone at k = 31 (the slow layers of P5 and its fusions run in the long phase only: in a short run the
    # declared window row of component 2 would be counted although the slow layer suppresses that band)
    J += [job(b, s) for s in ("P2_borealis_si", "P3_wdm_chop6k", "P4_wdm_R0.3", "P8_two_tone") for b in B]
    J += [job(b, s, dur=3.0) for s in ("P8_intermittent", "F_P8int+P3") for b in two]
    J += [job(b, s) for s in ("P1_guardslot", "P6_steps") for b in two]
    # P1 in calibration frames: the time loop is the baseline's; one seed for the offline stray table only
    J += [{**job(b, "P1_calframe", seeds=[1]), "decompose": False} for b in B]
    J += [job(b, s) for s in ("het_ref", "P8_two_tone_copol", "P8_two_tone_detune0.3", "P3_wdm_chop6k_R1.5",
                              "P4_wdm_R1.5", "P4_wdm_Rm0.5", "P4_wdm_implicit_bad") for b in two]
    # baseline variants
    J += [job(b, s) for s in ("L-LiuPD-A", "L-LiuPD-B", "L-A-seed", "L-LOdeliv") for b in two]
    J += [job(b, s) for s in ("pair_dk1", "pair_dk16", "pair_dk16_rho99", "lab_acoustic_g0", "lab_acoustic_g1",
                              "lab_mm915b", "lab_both_g0", "ctrl_moku", "ctrl_dsc1") for b in two]
    # scheme length sweep (F_liu20)
    J += [job("F_liu20", s, k) for s in ("P2_borealis_si", "P3_wdm_chop6k", "P4_wdm_R0.3", "P8_two_tone")
          for k in (4, 16)]
    J += [job("F_liu20", "P8_intermittent", k, dur=3.0) for k in (4, 16)]
    for f in (1e3, 3e3, 1e4, 3e4):
        J.append(job("F_liu20", "P3_wdm_chop6k", extra={"gate_f_ch_hz": f}, suffix=f"_fch{f:g}"))
    for T in (0.01, 1.0):
        J.append(job("F_liu20", "P8_intermittent", dur=3.0, extra={"gate_period_s": T, "win_period_s": T},
                     suffix=f"_T{T:g}"))
    # D12: the LO-delivery lock starts at the top of its fringe (error −π/2, zero slope) and pulls in over ~1 s, longer
    # than the metric settle time 5 ms: metric (i) evaluated after 1 s (3 s runs), with the baseline at the same window
    for b in two:
        for s in ("L-LOdeliv", "base"):
            J.append(job(b, s, dur=3.0, extra={"t_settle_s": 1.0}, suffix="_settle1", comp_seeds=[1, 2, 3]))
    return J


def long_jobs():
    """30 s = t_cal_s (the window row is 0, the slow layers act in the run; D78 / D98). 10 seeds for the totals and
    every quantity, the four components on seeds 1–3."""
    J = [job("F_liu20", s, dur=30.0, comp_seeds=[1, 2, 3])
         for s in ("base", "P5_dbs", "F_P3+P5", "P3_wdm_chop6k", "F_P2+P5", "P5_dbs_pumped_bad", "F_P1slow+Liu")]
    J += [job("F_wg6", s, dur=30.0, comp_seeds=[1, 2, 3]) for s in ("base", "P5_dbs")]
    # E41 (orthogonal-probe PM birefringence drift) sets component 1 at 30 s (49 mrad on F_liu20): the co-polarised
    # two-tone probe (P8) avoids it; its Rayleigh backscatter is ±675 MHz off the OPO mode (stray 1.5e-8 at 1 s)
    J += [job("F_liu20", "P8_two_tone_copol", dur=30.0, comp_seeds=[1, 2, 3])]
    return J


def _run(j):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    lib = _lib()
    try:
        out = lib.run_job(j)
        return j["tag"], f"ok {out['wall_s']:.0f} s"
    except Exception as e:  # noqa: BLE001 - a failed job is reported, the campaign goes on
        return j["tag"], f"ERROR {type(e).__name__}: {str(e)[:300]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("screen", "long"), required=True)
    ap.add_argument("--procs", type=int, default=5)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--only", default=None, help="substring filter on the tag")
    a = ap.parse_args()
    procs = min(a.procs, 5)  # hard cap: the axis-1 campaign runs at the same time
    lib = _lib()
    jobs = screen_jobs() if a.phase == "screen" else long_jobs()
    jobs = [j for j in jobs if (a.only is None or a.only in j["tag"]) and not (lib.OUT / f"{j['tag']}.json").exists()]
    print(f"{len(jobs)} jobs to run", flush=True)
    if a.dry:
        for j in jobs:
            print(j["tag"])
        return
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = [ex.submit(_run, j) for j in jobs]
        for f in as_completed(futs):
            tag, msg = f.result()
            print(f"{tag}: {msg}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
