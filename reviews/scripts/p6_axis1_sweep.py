"""Phase 6 axis 1, §8 step 3: main-axis sweeps of the bases (a) refA_fiber_bench, (b) refC_fiber_bench,
(c) pulsed_wg_1port, with N_SEEDS seeds per configuration (common seeds 1..N across configurations).

Usage: uv run python reviews/scripts/p6_axis1_sweep.py <A|C|P|Pdec> [n_seeds]
Writes runs/p6_axis1/sweep_<case>.json (raw rows) and docs/figures/phase6/axis1/sweep_<case>_summary.json
(mean / std / min / max per configuration and quantity).
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

LENGTHS = [1.0, 3.1622776601683795, 10.0, 31.622776601683793, 100.0]

P_LOCKS = {
    "liu": [],
    "het": ["probe/detect_heterodyne"],
    "wdm100G": ["probe/sep_wdm", "probe/freq_100GHz"],
    "wdm1T": ["probe/sep_wdm", "probe/freq_1THz"],
    "wdm1T_path": ["probe/sep_wdm", "probe/freq_1THz", "probe/actuator_pathlength"],
    "P3": ["probe/sep_wdm", "probe/freq_1THz", "gate/chopper_wdm"],
    "P8": ["probe/two_tone"],
    "P8copol": ["probe/two_tone", "probe/two_tone_copol"],
}
P_AXES = {
    "wdm1T_path": {"offline.params.bias_case": "R_+0.3"},
    "P3": {"offline.params.gate_f_ch_hz": 6.0e3, "offline.params.gate_jitter_s": 0.0},   # P3's 6 kHz (gate/chopper_p3)
    "P8": {"offline.params.tt_dnu_target_hz": 7.0e9, "offline.params.tt_on_grid": True, "offline.params.tt_r": 0.1,
           "offline.params.tt_dL0_detune_rad": 0.0},
}
P_AXES["P8copol"] = P_AXES["P8"]


def configs_A():
    out = {}
    for route in ("fiber_matched", "free_space"):
        for dx1 in (False, True):
            for ac in (False, True):
                for Lm in LENGTHS:
                    out[f"A|{route}|{'dx1' if dx1 else 'koh'}|{'lab' if ac else 'quiet'}|fpga|{Lm:g}"] = \
                        L.a_point(Lm, route, dx1=dx1, acoustic=ac)
    for ctrl, cn in (("bench/moku_delta_refA", "moku"), ("bench/dsc1_refA", "dsc1")):
        for dx1, ac in ((False, False), (True, True)):
            for Lm in LENGTHS:
                out[f"A|fiber_matched|{'dx1' if dx1 else 'koh'}|{'lab' if ac else 'quiet'}|{cn}|{Lm:g}"] = \
                    L.a_point(Lm, dx1=dx1, acoustic=ac, ctrl=ctrl)
    for off in ("as_built", "zero", "fc_opt", "fc_opt_m10ns", "fc_opt_p10ns"):
        for Lm in (1.0, 10.0, 100.0):
            out[f"Aoff|fiber_matched|dx1|lab|{off}|{Lm:g}"] = L.a_point(Lm, dx1=True, acoustic=True, offset=off)
    # DX1 remedies: the LO filter cavity's laser low-pass is the DX1 limit (ablation A_lab): a 10x wider cavity, or none
    # (Δν_FC = 10 GHz ≈ no cavity: time loop AND the koheras ledger pair read fc.dnu_fc_hz; the toggle
    # fc.laser_lowpass alone would leave the ledger's FC term in; LO excess noise at the sideband is NOT modelled)
    wide = [{"name": "fc", "params": {"dnu_fc_hz": 3.3e7}}]
    none = [{"name": "fc", "params": {"dnu_fc_hz": 1.0e10}}]
    for Lm in (1.0, 10.0, 100.0):
        for off in ("as_built", "zero", "fc_opt"):
            out[f"Afc|fiber_matched|dx1|lab|fc33MHz_{off}|{Lm:g}"] = L.a_point(Lm, dx1=True, acoustic=True, offset=off,
                                                                              pipeline=wide)
        for off in ("as_built", "zero"):
            out[f"Afc|fiber_matched|dx1|lab|nofc_{off}|{Lm:g}"] = L.a_point(Lm, dx1=True, acoustic=True, offset=off,
                                                                           pipeline=none)
    return out


def configs_C():
    out = {}
    for ctrl, cn in ((None, "fpga"), ("bench/moku_delta_refC", "moku"), ("bench/dsc1_refC", "dsc1")):
        for h0, hn in ((None, "nom"), (2.546e4, "dx1")):
            for Lm in LENGTHS:
                out[f"C|{cn}|{hn}|{Lm:g}"] = L.c_point(Lm, ctrl=ctrl,
                                                      axes={"offline.params.h0_hz2_per_hz": h0} if h0 else None)
    return out


def configs_P():
    out = {}
    for lk, ovs in P_LOCKS.items():
        for ctrl, cn in ((None, "fpga"), ("controller/dsc1", "dsc1"), ("controller/moku_delta", "moku")):
            for lab in (False, True):
                for Lm in LENGTHS:
                    out[f"P|{lk}|{cn}|{'lab' if lab else 'quiet'}|{Lm:g}"] = \
                        L.p_point(Lm, ovs, ctrl=ctrl, lab=lab, axes=P_AXES.get(lk))
    return out


def configs_Pdec():
    """The four metric-(ii) components (G13) for the headline lock variants, lab noise on, fpga and DSC1."""
    out = {}
    for lk in ("liu", "het", "wdm1T", "P3", "P8"):
        for ctrl, cn in ((None, "fpga"), ("controller/dsc1", "dsc1")):
            for Lm in (10.0, 100.0):
                out[f"Pdec|{lk}|{cn}|lab|{Lm:g}"] = L.p_point(Lm, P_LOCKS[lk], ctrl=ctrl, lab=True, axes=P_AXES.get(lk))
    return out


def configs_P5():
    """P5 two-layer lock (lock/dual_band_slow at its D98 nominal: seed after the squeezer, calibration frames) against
    the Liu baseline, 2 s runs (slow layer f_c = 2 Hz), lab noise, fpga. PRELIMINARY: seed_psa / cal_slot are not yet
    Phase 5 reviewed (D98-6)."""
    out = {}
    for Lm in (10.0, 100.0):
        out[f"P5|liu|fpga|lab|{Lm:g}"] = L.p_point(Lm, [], lab=True, duration_s=2.0)
        out[f"P5|dbs|fpga|lab|{Lm:g}"] = L.p_point(Lm, ["lock/dual_band_slow", "lock/points/dual_band_slow_nominal.yaml"],
                                                   lab=True, duration_s=2.0)
    return out


RECORD = {"P5": [],"A": ["metric.mi.port0", "v.dac1", "v.dac2"], "Aoff": ["metric.mi.port0"], "Afc": ["metric.mi.port0"],
          "C": ["metric.mi.port0"], "P": ["metric.mi.port1", "metric.mii.port1"], "Pdec": []}
PSD = {"P5": [], "A": ["metric.mi.port0"], "Aoff": [], "Afc": [], "C": ["metric.mi.port0"], "P": ["metric.mii.port1", "metric.mi.port1"],
       "Pdec": []}


KEEP = re.compile(r"^((i|ii)\.port\d\.(sigma_total_rad|sigma_phase_rad|std_rad|std_loop_rad|sigma_hf_rad|sigma_dither_rad|mean_loop_rad|r_minus_db|r_plus_db|eta)"
                  r"|ii\.port\d\.\w+_rad2|ii\.port\d\.linearity_ratio|eff\..*|lock\.[\w.]+\.(n_slip|n_loss|lock_duration_s)"
                  r"|range\..*|valid|validity\..*settle_frac|series\..*\.(std|absmax)|stray\..*|wall_s|n_events)$")


def main(case, n_seeds):
    cfgs = {"A": configs_A, "C": configs_C, "P": configs_P, "Pdec": configs_Pdec, "P5": configs_P5}[case]()
    jobs = []
    for s in range(1, n_seeds + 1):
        for lab, pt in cfgs.items():
            kind = lab.split("|")[0]
            jobs.append({"label": lab, "point": {**pt, "seed": s}, "record": RECORD[kind],
                         "t_from": L.T_SETTLE_A if kind.startswith("A") else 0.0,
                         "psd": PSD[kind] if s == 1 else None, "decompose": kind == "Pdec"})
    print(f"{case}: {len(cfgs)} configurations × {n_seeds} seeds = {len(jobs)} runs", flush=True)
    rows = L.load(f"sweep_{case}") if RESUM else L.pool(jobs, f"sweep_{case}")
    keys = sorted({k for r in rows for k, v in r.items() if isinstance(v, (int, float)) and KEEP.match(k)})
    summ = {}
    for lab in cfgs:
        st = L.stats(rows, lab, keys)
        summ[lab] = {k: {kk: v[kk] for kk in ("mean", "std", "min", "max", "n")} for k, v in st.items()}
        summ[lab]["_errors"] = [r["error"] for r in rows if r.get("label") == lab and "error" in r][:3]
        summ[lab]["_events"] = sorted({e for r in rows if r.get("label") == lab for e in r.get("events", [])})
    L.FIG.mkdir(parents=True, exist_ok=True)
    (L.FIG / f"sweep_{case}_summary.json").write_text(json.dumps(summ, indent=0))
    print("done", case)


RESUM = "--resum" in sys.argv   # rebuild the summary from runs/p6_axis1/sweep_<case>.json without running

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--resum"]
    main(args[0], int(args[1]) if len(args) > 1 else 10)
