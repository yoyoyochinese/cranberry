"""Phase 6 axis 2 (GBS integration) campaign: named configurations of the future presets and one job runner.

A configuration is (base preset, overlays in order, pins). Pins fix every sweep axis the overlays leave open, so a
configuration is one point (`build` refuses a spec with an axis left open, naming it): no axis is silently taken at
its first value. Pins go through `runner.sweep.expand` (the spec's unresolved tree), so dependent `!expr` values
follow (D64). Port index `port_k` (fibre length L_k = L0 + k(N_b + n_guard)c/(f_rep n_g)) is a pin like any other.
A pin "toggle:<block>.<toggle>" sets a RunSpec toggle instead.

Used by p6_axis2_run.py (jobs, ≤ 5-process pool) and p6_axis2_report.py (tables, figures). Raw JSON in runs/p6_axis2/.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import math
import time
import warnings
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

OUT = ROOT / "runs" / "p6_axis2"
BASES = ("F_wg100", "F_wg6", "F_bor6", "F_liu20")
OPO = ("F_bor6", "F_liu20")
COMP = ("locked_path_residual", "common_points", "laser_decorrelation_oob", "squeezed_only")

# ---- representative points of the overlay axes (each a documented centre / paper value, not an optimum, D34) -------
# D89-4 / D96-8: the Ref D-calibrated 915B mode set (one example consistent with the measurement, not unique).
MM915B = {"pzt_f_res_hz": 2.52e4, "pzt_q": 25.0, "mm2_f_hz": 5.4e3, "mm2_q": 15.0, "mm2_r": 0.11,
          "mm3_f_hz": 1.0e4, "mm3_q": 30.0, "mm3_r": 0.07, "mm4_f_hz": 4.03e4, "mm4_q": 100.0, "mm4_r": 0.15,
          "dac_output_noise_v": 0.0}
ACOUSTIC = {"acoustic_scale": 1.0, "acoustic_hf_slope": 2.0, "acoustic_rho_arms": 0.0}   # D89-3: arms uncorrelated
# P2 tail table (stray/window_tail): AOM 100 ns rise, −50 dB off (F_common eps_aom), no discrete reflection, mean.
WINTAIL = {"win_edge_s": 1.0e-7, "win_eps_off": 1.0e-5, "win_limit_statistic": "mean", "win_refl_r": 0.0}
# P1 forward reference (ref/aurora_calframe): 1 µW peak, 1 ns pulses, −50 dB carving, 50 dB slow gate (one fibre AOM),
# 1 ns edges, SNSPD 43 ns recovery, 1e-3 afterpulse at 100 ns (all UNVERIFIED axis values, D88-4).
AURORA = {"fref_p_peak_w": 1.0e-6, "fref_t_pulse_s": 1.0e-9, "fref_eps_mod": 1.0e-5, "fref_eps_gate": 1.0e-5,
          "fref_edge_tau_s": 1.0e-9, "det_p_afterpulse": 1.0e-3, "det_tau_afterpulse_s": 1.0e-7, "det_recovery_s": 4.3e-8}
CHOPPER_P3 = {"gate_f_ch_hz": 6.0e3, "gate_jitter_s": 0.0}   # gate/chopper_p3 (P3's best rate)
INTERMITTENT = {"gate_period_s": 0.1, "gate_lock_fraction": 0.2, "gate_settle_s": 1.0e-5}   # [AE] UNVERIFIED (Luo: no values)
P6 = {"prog_dwell_s": 1.0e-2, "pipeline.squeezer.seq.prog.params.step_rad": 0.7853981633974483}   # π/4 steps, 10 ms
# P8 (probe/two_tone): Luo's 675 MHz on the probe-axis grid, r = 0.1 (middle of 0–0.2), ΔL_0 trimmed (D94-5, D97-3).
TWO_TONE = {"tt_dnu_target_hz": 6.75e8, "tt_on_grid": True, "tt_r": 0.1, "tt_dL0_detune_rad": 0.0}
DSC1 = {"dsc1_latency_s": 5.7e-6}   # D96-2: τ_L 5.7 µs with the 120/100 kHz poles
WIRE = "lock/dual_band_slow_gate_wire"


def _p(d: dict[str, Any]) -> dict[str, Any]:
    """Short names → offline.params paths (a key with a dot is a full sweep path)."""
    return {(k if "." in k else f"offline.params.{k}"): v for k, v in d.items()}


def configs(base: str) -> dict[str, tuple[tuple[str, ...], dict[str, Any]]]:
    """{label: (overlays, pins)} of every axis-2 configuration on ``base`` (labels shared across bases)."""
    opo = base in OPO
    refmode = (f"ref/opo_mode_{base.split('_')[1]}",) if opo else ()
    straymode = (f"stray/opo_mode_{base.split('_')[1]}",) if opo else ()
    wdm = ("probe/sep_wdm", "probe/freq_1THz", "probe/actuator_pathlength")
    chop = ("gate/chopper_wdm", "gate/chopper_p3")
    bor = ("stray/window_tail", *straymode, "gate/borealis_window", "gate/borealis_si")
    inter = ("stray/window_tail", *straymode, "gate/intermittent_slow")
    dbs = ("lock/dual_band_slow", "lock/points/dual_band_slow_nominal")
    borp = {**WINTAIL, "win_period_s": 1.0e-4, "win_duty_lock": 0.65}
    interp = {**WINTAIL, **INTERMITTENT, "win_period_s": 0.1, "win_duty_lock": 0.2}
    return {
        # ---- baseline (Procedure 2) ----
        "base": ((), {}),
        "L-LiuPD-A": (("presets/locks/L-LiuPD-A.yaml",), {}),
        "L-LiuPD-B": (("presets/locks/L-LiuPD-B.yaml",), {}),
        "L-A-seed": (("presets/locks/L-A-seed.yaml",), {}),
        "L-LOdeliv": (("presets/locks/L-LOdeliv.yaml",), {}),
        "pair_dk1": (("ports/pair",), {"pair_dk": 1, "rho_ports": 0.0}),
        "pair_dk16": (("ports/pair",), {"pair_dk": 16, "rho_ports": 0.0}),
        "pair_dk16_rho99": (("ports/pair",), {"pair_dk": 16, "rho_ports": 0.99}),
        "lab_acoustic_g0": (("lab/acoustic_future",), {**ACOUSTIC, "acoustic_gamma": 0.0}),
        "lab_acoustic_g1": (("lab/acoustic_future",), {**ACOUSTIC, "acoustic_gamma": 1.0}),   # EXTRAPOLATED > 3 kHz
        "lab_mm915b": (("lab/multimode915b_future",), MM915B),
        "lab_both_g0": (("lab/acoustic_future", "lab/multimode915b_future"), {**ACOUSTIC, "acoustic_gamma": 0.0, **MM915B}),
        "ctrl_moku": (("controller/moku_delta",), {}),
        "ctrl_dsc1": (("controller/dsc1",), DSC1),
        # ---- P1 (forward reference pulses in calibration frames; offline stray table; same time loop as base) ----
        "P1_calframe": (("ref/aurora_calframe", *refmode), AURORA),
        "P1_guardslot": (("ref/aurora_guardslot", *refmode), {k: v for k, v in AURORA.items() if k != "fref_eps_gate"}),
        # ---- P2 (Borealis windows on the Liu probe) ----
        "P2_borealis_si": (bor, borp),
        # ---- P3 (WDM 1 THz probe + chopper anti-sync + LPF, path-length actuator, R nominal) ----
        "P3_wdm_chop6k": ((*wdm, *chop), {**CHOPPER_P3, "bias_case": "R_+0.3"}),
        "P3_wdm_chop6k_R1.5": ((*wdm, *chop), {**CHOPPER_P3, "bias_case": "R_+1.5"}),
        # ---- P4 (different-wavelength co-propagating reference, no gate) ----
        "P4_wdm_R0.3": (wdm, {"bias_case": "R_+0.3"}),
        "P4_wdm_R1.5": (wdm, {"bias_case": "R_+1.5"}),
        "P4_wdm_Rm0.5": (wdm, {"bias_case": "R_-0.5"}),
        "P4_wdm_implicit_bad": (wdm, {"bias_case": "implicit"}),
        # ---- P5 (two-layer: Liu fast + in-band seed counting slow; gated pump nominal, D98-7) ----
        "P5_dbs": (dbs, {}),
        "P5_dbs_pumped_bad": (dbs, {"dbs_cal_pump": "pumped", "dbs_mu_per_gate": 4.2}),
        "P5_dbs_guard": ((*dbs, "lock/points/dbs_cal_guard"), {}),
        # ---- P6 (program port phases through the heterodyne RF phase, settle mask) ----
        "P6_steps": (("probe/detect_heterodyne", "program/p6_steps"), P6),
        "het_ref": (("probe/detect_heterodyne",), {}),
        # ---- P8 (two-tone probe; intermittent slow lock) ----
        "P8_two_tone": (("probe/two_tone",), TWO_TONE),
        "P8_two_tone_copol": (("probe/two_tone", "probe/two_tone_copol"), TWO_TONE),
        "P8_two_tone_detune0.3": (("probe/two_tone",), {**TWO_TONE, "tt_dL0_detune_rad": 0.3}),
        "P8_intermittent": (inter, interp),
        # ---- fusions (a gate first, then the slow layer, its point and the wire file, D96-9) ----
        "F_P3+P5": ((*wdm, *chop, *dbs, WIRE), {**CHOPPER_P3, "bias_case": "R_+0.3"}),
        "F_P4+P5": ((*wdm, *dbs), {"bias_case": "R_+0.3"}),
        "F_P2+P5": ((*bor, *dbs, WIRE), borp),
        # P1 as the slow layer: bright classical reference pulses in the calibration frames, sensed without counting
        # noise or click saturation (an upper bound of a classical-PD reference; the slow-layer loop is P5's)
        "F_P1slow+Liu": (("ref/aurora_calframe", *refmode, *dbs),
                         {**AURORA, "toggle:pcf1.counting_noise": False, "toggle:pcf1.saturation": False}),
        # P8's intermittent hold on the P3/P4 WDM sensor (P8-a's second sensor is not a configuration, D85-15)
        "F_P8int+P3": ((*wdm, *inter), {**interp, "bias_case": "R_+0.3"}),
    }


def build(base: str, label: str, *, port_k: int = 31, duration_s: float = 3.0, seed: int = 1,
          extra: dict[str, Any] | None = None):
    """One single-run spec of configuration ``label`` on ``base`` (module docstring)."""
    from cranberry.runner.sweep import expand

    overlays, pins = configs(base)[label]
    allp = {**pins, "port_k": port_k, **(extra or {})}
    toggles = {k[7:]: bool(v) for k, v in allp.items() if k.startswith("toggle:")}
    spec = fs.load(fs.nominal(base), run={"duration_s": float(duration_s), "seed": int(seed), "seed_mode": "common"},
                   overlays=overlays)
    axes = {k: [v] for k, v in _p({k: v for k, v in allp.items() if not k.startswith("toggle:")}).items()}
    pts = expand(spec, axes, limit=2)
    if len(pts) != 1:
        raise ValueError(f"{base}/{label}: axes left open (pin them): {_open_axes(spec, axes)}")
    out = pts[0].spec
    if toggles:
        out = dataclasses.replace(out, toggles={**out.toggles, **toggles})
    return out


def _open_axes(spec, axes) -> list[str]:
    from cranberry.config import input_tree
    from cranberry.config.resolve import find_markers, path_str

    tree = input_tree(spec)
    return sorted(n for n in (path_str(tree, p) for p, _ in find_markers(tree)) if n not in axes)


def quantities(res, spec) -> dict[str, float]:
    """run_quantities + L_k + offline stray-table summaries (window_stray_light / forward_ref_stray_table check())."""
    from cranberry.runner.explore import run_quantities

    q = run_quantities(res)
    q["L_k_m"] = float(spec.offline["params"]["L_k_m"])
    q["L_0_m"] = float(spec.offline["params"].get("L_0_m", math.nan))
    q["n_chk"] = float(len(res.events))
    for b in res.offline.blocks:
        r = getattr(b, "result", None)
        if isinstance(r, dict):
            for key, v in r.items():
                if isinstance(v, (int, float)):
                    q[f"off.{b.name}.{key}"] = float(v)
    return q


def decompose_with_result(spec, components: bool = True):
    """runner.decompose.decompose_metric_ii, returning also the all-on RunResult (the spec itself with the validity
    slots recorded), so the comparison quantities need no extra run. Same groups, validity, ledger and 5 % rule.
    ``components`` False: the all-on run only (total = time loop + ledger, components NaN, status "total_only")."""
    from cranberry.blocks.meas.metrics import METRIC_II_COMPONENTS, MetricIIComponents
    from cranberry.runner.decompose import _var_with_ledger, validity_violations
    from cranberry.runner.explore import with_recorded
    from cranberry.runner.run import INVALIDATING_CHK, run

    cfg = spec.measure["decompose"]
    groups, validity, mask = cfg["groups"], cfg["validity"], cfg.get("mask_slot")
    owner = {t: c for c, ts in groups.items() for t in ts}
    spec = with_recorded(spec, [v["slot"] for v in validity] + ([mask] if mask and validity else []))

    def variant(on):
        return dataclasses.replace(spec, toggles={**dict(spec.toggles),
                                                  **{t: (on is None or owner[t] == on) for t in owner}})

    runs = {None: run(variant(None))}
    if not components:
        return runs[None], {k: MetricIIComponents(dict.fromkeys(METRIC_II_COMPONENTS, math.nan),
                                                  _var_with_ledger(runs[None], k), math.nan, "total_only", None)
                            for k in runs[None].metrics.metric_ii}
    for c in METRIC_II_COMPONENTS:
        if groups[c]:
            runs[c] = run(variant(c))
    reasons = []
    for key, r in runs.items():
        bad = validity_violations(r, validity, mask if validity else None)
        if bad:
            reasons.append(f"{key or 'all-on'}: " + "; ".join(bad))
        if any(e.get("code") in INVALIDATING_CHK for e in r.events):
            reasons.append(f"{key or 'all-on'}: unstable")
    out = {}
    for k in runs[None].metrics.metric_ii:
        total = _var_with_ledger(runs[None], k)
        var = {c: (_var_with_ledger(runs[c], k) if c in runs else 0.0) for c in METRIC_II_COMPONENTS}
        ratio = sum(var.values()) / total if total > 0 else 1.0
        why = list(reasons) + ([f"linearity {ratio:.4f}"] if abs(ratio - 1) > 0.05 else [])
        out[k] = MetricIIComponents(var, total, ratio, "decomposition_invalid" if why else "ok", "; ".join(why) or None)
    return runs[None], out


def run_job(job: dict[str, Any]) -> dict[str, Any]:
    """job: {base, label, port_k, duration_s, seeds, decompose, extra, tag}. Writes runs/p6_axis2/<tag>.json."""
    from cranberry.runner.explore import actuator_ranges, with_recorded
    from cranberry.runner.run import run

    warnings.simplefilter("ignore")
    t0 = time.time()
    rows, events = [], set()
    for s in job["seeds"]:
        spec = build(job["base"], job["label"], port_k=job.get("port_k", 31), duration_s=job["duration_s"], seed=s,
                     extra=job.get("extra"))
        spec = with_recorded(spec, [slot for slot, _, _ in actuator_ranges(spec).values()])
        if job.get("decompose", True):
            res, comps = decompose_with_result(spec, components=s in job.get("component_seeds", job["seeds"]))
        else:
            res, comps = run(spec), {}
        events |= {e["code"] for e in res.events}
        row = quantities(res, spec)
        try:  # the probe's CW stray light (stray1, both D46 units)
            st = fs.stray_light(res)
            row["stray1.photons_per_mode"] = float(st["photons_per_mode"])
            row["stray1.clicks_per_gate"] = float(st["clicks_per_gate"])
        except (StopIteration, KeyError):
            pass
        for k, c in comps.items():
            for name, v in c.variances_rad2.items():
                row[f"ii.port{k}.{name}_rad2"] = float(v)
            row[f"ii.port{k}.total_rad2"] = float(c.total_var_rad2)
            row[f"ii.port{k}.linearity_ratio"] = float(c.linearity_ratio)
            if c.status != "total_only":
                row[f"ii.port{k}.decomposition_ok"] = float(c.status == "ok")
            if c.reason:
                row[f"ii.port{k}.reason"] = c.reason
        rows.append(row)
    out = {**job, "rows": rows, "events": sorted(events), "wall_s": time.time() - t0}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{job['tag']}.json").write_text(json.dumps(out, indent=0, default=float))
    return out


def stats(rows: list[dict[str, Any]], key: str, f=None) -> dict[str, float]:
    """mean / sample std / min / max of rows[key] (optionally mapped by f, e.g. rad² → mrad)."""
    import numpy as np

    v = np.array([r.get(key, math.nan) for r in rows], dtype=float)
    if f is not None:
        v = f(v)
    return {"mean": float(np.nanmean(v)), "std": float(np.nanstd(v, ddof=1)) if v.size > 1 else 0.0,
            "min": float(np.nanmin(v)), "max": float(np.nanmax(v)), "n": int(v.size)}


def load_results(pattern: str = "*") -> dict[str, dict[str, Any]]:
    return {p.stem: json.loads(p.read_text()) for p in sorted(OUT.glob(f"{pattern}.json"))}
