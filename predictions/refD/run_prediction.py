"""Ref D blind prediction (end of Phase 5; PLAN §6.5, DECISIONS D04, D61-4, D65-3, D66-5, D68-5, D75).

Built ONLY from ``presets/reference/refD.yaml`` / ``refD_nominal.yaml`` exactly as committed (no parameter is
changed or tuned; ``blind/`` is never read). Every run goes through the experiment runner
(``cranberry.runner.run.run``) and the sweep expander (``cranberry.runner.sweep.expand``), so each point is the
preset with its own ``offline.params`` axis values substituted and every dependent ``!expr`` re-evaluated (D64).

Run (from the repository root)::

    uv run python predictions/refD/run_prediction.py --workers 5            # all groups, cached / resumable
    uv run python predictions/refD/run_prediction.py --aggregate-only       # rebuild results/ from the cache

Groups (``--groups``):

* ``nominal``  refD_nominal, NOMINAL_SEEDS, DAC quantisation on and off (D65-3), locked;
* ``unlocked`` refD_nominal with the DSC1 PID feedback off (free-running), first 5 nominal seeds;
* ``fsim2``    refD_nominal at f_sim = 2 MHz (DSC1 device rates kept: ADC 1 MHz, PID / DAC 500 kHz, D64),
               FSIM2_SEEDS, quantisation on (convergence check, PLAN §6.4);
* ``oat``      one-at-a-time: every sweep axis of refD over the preset's own grid values, others at
               refD_nominal; quantisation on and off;
* ``lhs``      Latin hypercube (scipy.stats.qmc.LatinHypercube, seed LHS_SEED, N = LHS_N) over all axes'
               full ranges (``!sweep`` linear, ``!sweep_log`` logarithmic, ``!choice`` equiprobable);
               quantisation on and off;
* ``resgrid``  the preset's 4 × 4 (f_res, Q) grid, others nominal, quantisation on (time-loop check of the
               analytic loop_stability_check verdict, D75);
* ``stabgrid`` (offline only, no time loop) the loop_stability_check open loop of the preset on a fine
               (f_res, Q) grid for the DSC1 latency variants and the extreme plant gains (D75 region).

Per run (after T_SETTLE_S): Welch PSDs (Hann, linear detrend per segment, 50 % overlap, NPERSEG at the record
rate 500 kHz → 0.24 Hz resolution, ≈ 27 averages over the preset's 60 s) of the true phase error
``lock.lockD.err``, the PDB450C output ``v.pdb`` (DSC1 input) and the DSC1 output ``v.dsc1.out`` (tee to the
914 and the oscilloscope), recorded with decim 2 and ``antialias: fir`` (D68-5) and divided by the recorder
FIR's |H(f)|²; log-binned (100 bins/decade); automatic peak detection (no prior frequency: running-median
baseline over ±0.5 decade in log f, a bin is a peak if it is a local maximum whose excess over the baseline
is ≥ max(3 dB, 6 σ of the Welch estimate in that bin)); band rms; scope-like time features; lock bookkeeping
(lock.dsc1.sup.*), the CV saturation fraction and the loop_stability_check verdict of the point.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REFD = "presets/reference/refD.yaml"
NOMINAL = "presets/reference/refD_nominal.yaml"

SLOTS = ["lock.lockD.err", "v.pdb", "v.dsc1.out", "v.dsc1.cv"]
SPEC_SIGNALS = ["lock.lockD.err", "v.pdb", "v.dsc1.out"]
F_REC_HZ = 5.0e5  # record rate: decim = f_sim / F_REC_HZ (2 at 1 MHz, 4 at 2 MHz)
T_SETTLE_S = 0.5  # lock acquisition + recorder FIR start-up discarded
NPERSEG = 2**21  # 4.19 s at 500 kHz → Δf = 0.238 Hz
F_REPORT = (1.0, 2.0e5)  # reported band
BIN_EDGES = 10.0 ** np.arange(np.log10(0.2), np.log10(2.5e5) + 1e-9, 0.01)  # 100 bins / decade
ERR_BANDS = [(1.0, 10.0), (10.0, 100.0), (100.0, 1.0e3), (1.0e3, 1.0e4), (1.0e4, 1.0e5), (1.0e5, 2.0e5),
             (1.0, 2.0e5)]
PEAK_MIN_DB, PEAK_NSIGMA, PEAK_BASE_DEC = 3.0, 6.0, 0.5
WELCH_C = 2.05  # Hann bin-to-bin and 50 % overlap correlation factor of the Welch variance (VT-14 docstring)
CV_LIMIT_V = 4.0  # D35
ROOT_SEED = 20261031  # refD run.seed
NOMINAL_SEEDS = [ROOT_SEED + k for k in range(8)]
UNLOCKED_SEEDS = NOMINAL_SEEDS[:5]
FSIM2_SEEDS = NOMINAL_SEEDS[:4]
LHS_SEED = 20261101
LHS_N = 128


# ----------------------------------------------------------------------------------------------- axes
def axes_of_preset():
    """[(name, kind, lo, hi, options)] of every offline.params sweep axis of refD (preset order) and the
    refD_nominal value of each."""
    import yaml

    from cranberry import config
    from cranberry.config.yaml_io import Choice, Sweep, SweepLog

    spec = config.load(REFD)
    axes = []
    for path, m in config.sweep_markers(spec).items():
        if not path.startswith("offline.params.") or not isinstance(m, (Sweep, SweepLog, Choice)):
            continue
        name = path.split(".", 2)[2]
        if isinstance(m, Choice):
            axes.append({"name": name, "kind": "choice", "options": list(m.options)})
        else:
            axes.append({"name": name, "kind": "sweep" if isinstance(m, Sweep) else "sweep_log",
                         "lo": float(m.lo), "hi": float(m.hi), "n": int(m.n)})
    nominal = yaml.safe_load((ROOT / NOMINAL).read_text())["offline"]["params"]
    for a in axes:
        a["nominal"] = nominal[a["name"]]
    return axes


def grid_values(a):
    from cranberry.runner.sweep import axis_values
    from cranberry.config.yaml_io import Choice, Sweep, SweepLog

    if a["kind"] == "choice":
        return list(a["options"])
    cls = Sweep if a["kind"] == "sweep" else SweepLog
    return axis_values(cls(a["lo"], a["hi"], a["n"]))


def lhs_design(axes, n, seed):
    from scipy.stats import qmc

    u = qmc.LatinHypercube(d=len(axes), seed=seed).random(n)
    pts = []
    for row in u:
        v = {}
        for a, x in zip(axes, row):
            if a["kind"] == "choice":
                v[a["name"]] = a["options"][min(int(x * len(a["options"])), len(a["options"]) - 1)]
            elif a["kind"] == "sweep":
                v[a["name"]] = a["lo"] + x * (a["hi"] - a["lo"])
            else:
                v[a["name"]] = a["lo"] * (a["hi"] / a["lo"]) ** x
        pts.append(v)
    return u, pts


def jobs_for(groups, axes):
    nom = {a["name"]: a["nominal"] for a in axes}
    jobs = []
    if "nominal" in groups:
        for s in NOMINAL_SEEDS:
            for q in (True, False):
                jobs.append({"id": f"nominal_s{s}_q{int(q)}", "group": "nominal", "base": "nominal", "seed": s,
                             "quant": q, "feedback": True, "f_sim": 1.0e6, "values": nom})
    if "unlocked" in groups:
        for s in UNLOCKED_SEEDS:
            jobs.append({"id": f"unlocked_s{s}", "group": "unlocked", "base": "nominal", "seed": s, "quant": True,
                         "feedback": False, "f_sim": 1.0e6, "values": nom})
    if "fsim2" in groups:
        for s in FSIM2_SEEDS:
            jobs.append({"id": f"fsim2_s{s}", "group": "fsim2", "base": "nominal", "seed": s, "quant": True,
                         "feedback": True, "f_sim": 2.0e6, "values": nom})
    if "oat" in groups:
        lower = [a["name"].lower() for a in axes]
        for k, a in enumerate(axes):
            # case-insensitive file systems: k_T_rad_per_m_K (thermal) and K_T_rad_per_m_K (birefringence)
            # would share a cache file name, so such axes get their axis index in the id
            tag = a["name"] + (f"_ax{k:02d}" if lower.count(a["name"].lower()) > 1 else "")
            for i, v in enumerate(grid_values(a)):
                if v == a["nominal"] or (not isinstance(v, str) and not isinstance(a["nominal"], str)
                                         and math.isclose(float(v), float(a["nominal"]), rel_tol=1e-12, abs_tol=0)):
                    continue
                for q in (True, False):
                    jobs.append({"id": f"oat_{tag}_{i}_q{int(q)}", "group": "oat", "axis": a["name"],
                                 "base": "refD", "seed": None, "quant": q, "feedback": True, "f_sim": 1.0e6,
                                 "values": {**nom, a["name"]: v}})
    if "lhs" in groups:
        _, pts = lhs_design(axes, LHS_N, LHS_SEED)
        for k, v in enumerate(pts):
            for q in (True, False):
                jobs.append({"id": f"lhs_{k:03d}_q{int(q)}", "group": "lhs", "lhs_index": k, "base": "refD",
                             "seed": None, "quant": q, "feedback": True, "f_sim": 1.0e6, "values": v})
    if "resgrid" in groups:
        fa = next(a for a in axes if a["name"] == "f_res_hz")
        qa = next(a for a in axes if a["name"] == "q_res")
        for i, fr in enumerate(grid_values(fa)):
            for j, qq in enumerate(grid_values(qa)):
                jobs.append({"id": f"resgrid_{i}{j}", "group": "resgrid", "base": "refD", "seed": None, "quant": True,
                             "feedback": True, "f_sim": 1.0e6, "values": {**nom, "f_res_hz": fr, "q_res": qq}})
    return jobs


# ----------------------------------------------------------------------------------------------- spec
def build_spec(job, duration_s=None):
    """The single-run spec of ``job`` (preset values only; record / toggles / f_sim are run settings)."""
    from cranberry import config
    from cranberry.runner.sweep import expand

    f_sim = float(job["f_sim"])
    decim = int(round(f_sim / F_REC_HZ))
    tree = {"extends": NOMINAL if job["base"] == "nominal" else REFD,
            "measure": {"record": {"slots": SLOTS, "decim": decim, "antialias": "fir"}}}
    run = {}
    if duration_s is not None:
        run["duration_s"] = float(duration_s)
    if f_sim != 1.0e6:
        k = int(round(f_sim / 1.0e6))
        run["f_sim_hz"] = f_sim
        # DSC1 device rates stay physical (measure.convergence.keep_rate of refD; D64)
        tree["pipeline"] = [{"name": "dsc1.adc", "rate_div": k}, {"name": "dsc1.pid", "rate_div": 2 * k},
                            {"name": "dsc1.sup", "rate_div": 2 * k}, {"name": "dsc1.dac", "rate_div": 2 * k}]
    if job["base"] == "nominal":
        run["seed"] = int(job["seed"])
        if run:
            tree["run"] = run
        spec = config.load(tree)
    else:
        if run:
            tree["run"] = run
        base = config.load(tree)
        (pt,) = expand(base, {f"offline.params.{k}": [v] for k, v in job["values"].items()})
        spec = pt.spec
    tog = dict(spec.toggles)
    if not job["quant"]:
        tog["dsc1.dac.quantization"] = False
    if not job["feedback"]:
        tog["dsc1.pid.feedback"] = False
    return dataclasses.replace(spec, toggles=tog)


# ----------------------------------------------------------------------------------------------- analysis
def welch(x, fs, nperseg):
    """One-sided PSD (Hann, linear detrend per segment, 50 % overlap); (f, P, n_segments)."""
    from scipy.signal import get_window

    w = get_window("hann", nperseg)
    step = nperseg // 2
    nseg = (x.size - nperseg) // step + 1
    t = np.arange(nperseg, dtype=float) - (nperseg - 1) / 2.0
    tt = float(t @ t)
    acc = np.zeros(nperseg // 2 + 1)
    for k in range(nseg):
        s = x[k * step:k * step + nperseg]
        s = s - s.mean() - t * (float(t @ s) / tt)
        acc += np.abs(np.fft.rfft(s * w)) ** 2
    p = acc / nseg * 2.0 / (fs * float(w @ w))
    p[0] /= 2.0
    p[-1] /= 2.0
    return np.fft.rfftfreq(nperseg, 1.0 / fs), p, nseg


def log_bin(f, p):
    idx = np.digitize(f, BIN_EDGES) - 1
    ok = (idx >= 0) & (idx < BIN_EDGES.size - 1) & (f > 0)
    nb = BIN_EDGES.size - 1
    cnt = np.bincount(idx[ok], minlength=nb).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        pb = np.bincount(idx[ok], weights=p[ok], minlength=nb) / cnt
        fb = np.bincount(idx[ok], weights=f[ok], minlength=nb) / cnt
    return fb, pb, cnt


def _local_line(lf, lp, x):
    """Baseline at log f = x: straight line in (log f, dB) fitted over |log f − x| ≤ PEAK_BASE_DEC, refitted
    twice without the points more than 3 dB above it (robust to the peak itself, and unbiased on a sloped
    spectrum whose bins are unevenly spread in log f — a running median is not: at low f the raw bins are
    linear in f, so the median sits at higher f and a red spectrum shows spurious 'peaks')."""
    m = np.abs(lf - x) <= PEAK_BASE_DEC
    xs, ys = lf[m] - x, lp[m]
    keep = np.ones(xs.size, bool)
    for _ in range(3):
        if keep.sum() < 3:
            break
        a, b = np.polyfit(xs[keep], ys[keep], 1)
        keep = ys - (a * xs + b) <= 3.0
    return float(b)


def detect_peaks(fb, pb, cnt, nseg):
    """Automatic peaks of a log-binned PSD (module docstring). Returns a list of dicts."""
    from scipy.signal import find_peaks

    ok = np.isfinite(pb) & (pb > 0) & (cnt > 0)
    f, p, c = fb[ok], pb[ok], cnt[ok]
    lf, lp = np.log10(f), 10.0 * np.log10(p)
    base = np.array([_local_line(lf, lp, x) for x in lf])
    exc = lp - base
    thr = np.maximum(PEAK_MIN_DB, 10.0 * np.log10(1.0 + PEAK_NSIGMA * np.sqrt(WELCH_C / (nseg * c))))
    idx, _ = find_peaks(exc)
    blin = 10.0 ** (base / 10.0)
    e = p - blin
    out = []
    for i in idx:
        if not (F_REPORT[0] <= f[i] <= F_REPORT[1]) or exc[i] < thr[i]:
            continue
        half = e[i] / 2.0
        lo = i
        while lo > 0 and e[lo] > half:
            lo -= 1
        hi = i
        while hi < f.size - 1 and e[hi] > half:
            hi += 1
        fwhm = math.nan
        if e[lo] <= half and e[hi] <= half:
            fl = 10 ** np.interp(half, [e[lo], e[lo + 1]], [lf[lo], lf[lo + 1]]) if lo < i else f[i]
            fh = 10 ** np.interp(half, [e[hi], e[hi - 1]], [lf[hi], lf[hi - 1]]) if hi > i else f[i]
            fwhm = float(fh - fl)
        if not math.isfinite(fwhm):
            continue  # not a peak: the excess does not fall to half on both sides (edge of a slope)
        # integrated excess power over the contiguous region where the PSD is above the baseline
        a, b = i, i
        while a > 0 and e[a - 1] > 0:
            a -= 1
        while b < f.size - 1 and e[b + 1] > 0:
            b += 1
        seg = slice(a, b + 1)
        power = float(np.trapezoid(np.clip(e[seg], 0, None), f[seg])) if b > a else float(e[i] * f[i] * 0.023)
        out.append({"f_hz": float(f[i]), "psd": float(p[i]), "excess_db": float(exc[i]), "threshold_db": float(thr[i]),
                    "fwhm_hz": fwhm, "bin_width_rel": 0.023, "power": power})
    return out


def band_rms(f, p, lo, hi):
    m = (f >= lo) & (f <= hi)
    return float(math.sqrt(max(np.trapezoid(p[m], f[m]), 0.0))) if m.sum() > 1 else math.nan


def scope_features(x, fs):
    """Oscilloscope-like numbers of a recorded voltage (record rate 500 kHz, FIR-limited to 250 kHz)."""
    out = {"mean": float(x.mean()), "std": float(x.std()), "ptp_full": float(np.ptp(x))}
    for w_s, tag in ((1e-3, "1ms"), (1e-2, "10ms"), (1e-1, "100ms")):
        n = int(round(w_s * fs))
        k = x.size // n
        pp = np.ptp(x[: k * n].reshape(k, n), axis=1)
        out[f"ptp_{tag}_median"] = float(np.median(pp))
        out[f"ptp_{tag}_p95"] = float(np.percentile(pp, 95))
        out[f"ptp_{tag}_max"] = float(pp.max())
    return out


def analyse(res, spec):
    from scipy.signal import freqz

    from cranberry.runner.run import antialias_taps

    f_sim = spec.run.f_sim_hz
    decim = int(round(f_sim / F_REC_HZ))
    t = res.t_s
    m = t >= T_SETTLE_S
    out = {"scalars": {}, "peaks": {}, "spectra": {}}
    sc = out["scalars"]
    taps = antialias_taps("fir", decim)
    for name in SPEC_SIGNALS:
        x = np.asarray(res.series[name], float)[m]
        f, p, nseg = welch(x, F_REC_HZ, NPERSEG)
        h2 = np.abs(freqz(taps, worN=f, fs=f_sim)[1]) ** 2
        p = p / np.maximum(h2, 1e-6)
        fb, pb, cnt = log_bin(f, p)
        out["spectra"][name] = pb
        out["fb"] = fb
        out["cnt"] = cnt
        out["nseg"] = nseg
        out["peaks"][name] = detect_peaks(fb, pb, cnt, nseg)
        for lo, hi in ERR_BANDS:
            sc[f"{name}.rms_{lo:g}_{hi:g}"] = band_rms(f, p, lo, hi)
        sel = (fb >= 10.0) & (fb <= F_REPORT[1]) & np.isfinite(pb)
        sc[f"{name}.f_max_fS_hz"] = float(fb[sel][np.argmax((fb * pb)[sel])])
        if name != "lock.lockD.err":
            for k, v in scope_features(x, F_REC_HZ).items():
                sc[f"{name}.{k}"] = v
        del x
    err = np.asarray(res.series["lock.lockD.err"], float)[m]
    mean = float(err.mean())
    branch = int(round(mean / math.pi))
    sc["err.mean"] = mean
    sc["err.branch"] = branch  # 0: lock at +π/2 (setpoint), ±1: the other fringe slope (−π/2)
    sc["err.offset"] = mean - branch * math.pi  # lock-point offset from the nearest stable point
    sc["err.std"] = float(err.std())
    sc["err.max_abs_dev"] = float(np.max(np.abs(err - mean)))
    cv = np.asarray(res.series["v.dsc1.cv"], float)
    sc["cv.frac_saturated"] = float(np.mean(np.abs(cv[m]) >= CV_LIMIT_V - 1e-9))
    sc["cv.frac_saturated_all"] = float(np.mean(np.abs(cv) >= CV_LIMIT_V - 1e-9))
    early = (t >= T_SETTLE_S) & (t < T_SETTLE_S + 0.5)
    sc["cv.u0"] = float(cv[early].mean())
    sc["cv.mean"] = float(cv[m].mean())
    sc["cv.end"] = float(cv[-int(F_REC_HZ):].mean())
    st = res.final_state
    sc["sup.n_loss"] = int(st["lock.dsc1.sup.n_loss"])
    sc["sup.t_loss"] = float(st["lock.dsc1.sup.t_loss"])
    sc["sup.n_slip"] = int(st["lock.dsc1.sup.n_slip"])
    chk = next(b for b in res.offline.blocks if getattr(b, "kind", "") == "loop_stability_check")
    out["stability"] = {k: (v if not isinstance(v, float) or math.isfinite(v) else None)
                        for k, v in chk.result.items() if k != "modes"}
    out["valid"] = bool(res.valid)
    out["events"] = [e.get("code") for e in res.events]
    return out


# ----------------------------------------------------------------------------------------------- worker
def run_job(job, cache_dir, duration_s=None):
    os.chdir(ROOT)
    warnings.simplefilter("ignore")
    from cranberry.config import spec_hash
    from cranberry.runner.run import run

    t0 = time.time()
    spec = build_spec(job, duration_s)
    res = run(spec, keep_state=True)
    a = analyse(res, spec)
    rec = {k: v for k, v in job.items()}
    rec.update({"seed_used": int(spec.run.seed), "spec_hash": spec_hash(spec), "duration_s": spec.run.duration_s,
                "f_sim_hz": spec.run.f_sim_hz, "mode": res.mode, "wall_s": time.time() - t0,
                "scalars": a["scalars"], "peaks": a["peaks"], "stability": a["stability"], "valid": a["valid"],
                "events": a["events"], "nseg": a["nseg"], "toggles_override": {
                    k: v for k, v in spec.toggles.items() if k in ("dsc1.dac.quantization", "dsc1.pid.feedback")}})
    np.savez_compressed(Path(cache_dir) / f"{job['id']}.npz", fb=a["fb"], cnt=a["cnt"],
                        **{n: a["spectra"][n] for n in SPEC_SIGNALS})
    (Path(cache_dir) / f"{job['id']}.json").write_text(json.dumps(rec, default=_json_default))
    return job["id"], rec["wall_s"]


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# ----------------------------------------------------------------------------------------------- stability grid
def stability_grid(n_f=41, n_q=41):
    """loop_stability_check verdict of refD on a fine (f_res, Q) grid (D75), for the DSC1 latency variants and
    the smallest / nominal / largest plant gain k0 reachable inside the preset's axis ranges. Uses the
    check's own open-loop model and Nyquist test (cranberry.blocks.chk.loop_stability) with the chk params
    of a resolved refD point, so the formula is the preset's."""
    os.chdir(ROOT)
    warnings.simplefilter("ignore")
    from cranberry import config
    from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable, open_loop
    from cranberry.runner.run import offline_block_specs
    from cranberry.runner.sweep import expand

    axes = axes_of_preset()
    nom = {a["name"]: a["nominal"] for a in axes}
    gains = {
        "k0_min": {"il_915b_db": 2.0, "coupler_excess_db": 0.3, "responsivity_a_per_w": 1.0, "r_914_ohm": 1.0e4,
                   "kappa1": 0.485, "kappa2": 0.485},
        "k0_nominal": {},
        "k0_max": {"il_915b_db": 0.0, "coupler_excess_db": 0.0, "responsivity_a_per_w": 1.03, "r_914_ohm": 1.0e7,
                   "kappa1": 0.5, "kappa2": 0.5},
    }
    variants = {"5.7us_poles": {"latency_s": 5.7e-6, "in_out_poles": True},
                "8.3us_nopoles": {"latency_s": 8.3e-6, "in_out_poles": False},
                "5.7us_nopoles": {"latency_s": 5.7e-6, "in_out_poles": False},
                "8.3us_poles": {"latency_s": 8.3e-6, "in_out_poles": True}}
    fa = next(a for a in axes if a["name"] == "f_res_hz")
    qa = next(a for a in axes if a["name"] == "q_res")
    f_res = np.geomspace(fa["lo"], fa["hi"], n_f)
    q_res = np.geomspace(qa["lo"], qa["hi"], n_q)
    base = config.load(REFD)
    out = {"f_res_hz": f_res.tolist(), "q_res": q_res.tolist(), "cases": {}}
    for gname, gv in gains.items():
        for vname, vv in variants.items():
            vals = {**nom, **gv, **vv}
            (pt,) = expand(base, {f"offline.params.{k}": [v] for k, v in vals.items()})
            (cs,) = [b for b in offline_block_specs(pt.spec) if b.kind == "loop_stability_check"]
            p = dict(cs.params)
            lo, hi, n = p["f_grid_hz"]
            fgrid = np.geomspace(float(lo), float(hi), int(n))
            stable = np.zeros((n_f, n_q), bool)
            gm = np.full((n_f, n_q), np.nan)
            f180 = np.full((n_f, n_q), np.nan)
            nint = n_integrators(p)
            for i, fr in enumerate(f_res):
                for j, qq in enumerate(q_res):
                    p["modes"] = [{"f_hz": float(fr), "q": float(qq)}]
                    lg = open_loop(fgrid, p)
                    ok, mx, fx = nyquist_stable(fgrid, lg, nint)
                    stable[i, j], gm[i, j], f180[i, j] = ok, mx, fx
            out["cases"][f"{gname}/{vname}"] = {"k0": float(p["k0"]), "stable": stable.tolist(),
                                                "max_l_at_180": gm.tolist(), "f_180_hz": f180.tolist()}
    return out


# ----------------------------------------------------------------------------------------------- aggregation
def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def quantiles(x):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    if x.size == 0:
        return {"n": 0}
    return {"n": int(x.size), "median": float(np.median(x)), "p05": float(np.percentile(x, 5)),
            "p95": float(np.percentile(x, 95)), "min": float(x.min()), "max": float(x.max()),
            "mean": float(x.mean()), "std": float(x.std(ddof=1)) if x.size > 1 else 0.0}


def k_act_g(v):
    """915B gain [rad/V] (dopl_full 16 µm / λ0 / 10 V) and the DAC load divider of point values ``v``."""
    k_act = 2 * math.pi * 1.6e-5 / (1.5488e-6 * 10.0)
    r = 1.0 / (1.0 / 1.0e6 + 1.0 / float(v["r_914_ohm"]))
    return k_act, r / (r + 220.0)


def t_loss_extrapolated(v, u0):
    """Lock-loss time [s] of the thermal ramp of the preset's model (segment_thermal: ΔT = r_T (t − τ(1 −
    e^{−t/τ})), differential phase k_T·ΔL·ΔT, compensated by the 915B until the CV reaches ∓4 V; VT-21 form):
    the CV moves from its early value u0 towards −4 V (loop_sign +1) or +4 V (−1)."""
    from scipy.optimize import brentq

    k_act, g = k_act_g(v)
    rate = float(v["k_T_rad_per_m_K"]) * float(v["dL_m"]) * float(v["r_T_k_per_s"])  # rad/s at full ramp
    if rate <= 0:
        return math.inf
    reach_v = (u0 + CV_LIMIT_V) if float(v["loop_sign"]) > 0 else (CV_LIMIT_V - u0)
    reach = reach_v * k_act * g
    if reach <= 0:
        return 0.0
    tau = float(v["tau_th_s"])
    fn = lambda t: rate * (t - tau * (1 - math.exp(-t / tau))) - reach  # noqa: E731
    hi = 10.0
    while fn(hi) < 0:
        hi *= 2
    return float(brentq(fn, 1e-6, hi))


LOCK_STD_MAX_RAD, LOCK_OFFSET_MAX_RAD, LOCK_SAT_MAX = 0.3, 0.3, 0.01


def time_loop_locked(sc):
    """The time loop holds a lock point: true-phase std < 0.3 rad, |offset from the nearest stable point|
    < 0.3 rad and the CV at its limit < 1 % of the time (an oscillating / slipping loop has std ≳ 1 rad)."""
    return (sc["err.std"] < LOCK_STD_MAX_RAD and abs(sc["err.offset"]) < LOCK_OFFSET_MAX_RAD
            and sc["cv.frac_saturated"] < LOCK_SAT_MAX)


def aggregate(cache_dir, out_dir, axes, stab, wall_total):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    recs = {}
    for pth in sorted(Path(cache_dir).glob("*.json")):
        if pth.stem != "stability_grid":
            recs[pth.stem] = json.loads(pth.read_text())
    ids = sorted(recs)
    fb = None
    arrays = {n: [] for n in SPEC_SIGNALS}
    for i in ids:
        z = np.load(Path(cache_dir) / f"{i}.npz")
        fb = z["fb"] if fb is None else fb
        for n in SPEC_SIGNALS:
            arrays[n].append(z[n])
            # peaks re-detected from the cached spectrum with the current rule (the worker's copy is replaced)
            recs[i]["peaks"][n] = detect_peaks(z["fb"], z[n], z["cnt"], recs[i]["nseg"])
    np.savez_compressed(out_dir / "spectra.npz", run_ids=np.array(ids), f_bin_hz=fb, bin_edges_hz=BIN_EDGES,
                        **{f"psd.{n}": np.array(arrays[n]) for n in SPEC_SIGNALS})
    # derived per-run numbers
    for i in ids:
        r = recs[i]
        sc = r["scalars"]
        v = r["values"]
        sc["t_loss_extrap_s"] = t_loss_extrapolated(v, sc["cv.u0"]) if r["feedback"] else None
        tau = 1.468 * float(v["dL_m"]) / 299792458.0
        h0 = float(v["h0_hz2_per_hz"])
        sc["err.sigma_oob_rad"] = math.sqrt(max(2 * math.pi**2 * h0 * tau - 4 * math.pi**2 * h0 * tau**2
                                                * r["f_sim_hz"] / 2, 0.0))
        r["stable"] = bool(r["stability"]["stable"])
        # time-loop outcome, independent of the analytic verdict (D75 asks for the verdict; both are kept)
        r["locked"] = bool(r["feedback"] and time_loop_locked(sc))
        # lock.dsc1.sup.n_slip counts changes of round(err / 2π) around the +π/2 setpoint: a loop that locks
        # on the other slope (err ≈ −π, branch −1) toggles that boundary with every noise excursion, so the
        # counter is meaningful only for branch 0 (recorded as is, flagged here)
        sc["sup.n_slip_meaningful"] = bool(sc["err.branch"] == 0)
    with open(out_dir / "runs.json", "w") as fh:
        json.dump({"runs": [recs[i] for i in ids]}, fh, default=_json_default, indent=0)
    summary = summarise(recs, ids, fb, arrays, axes)
    summary["stability_grid_summary"] = summarise_stab(stab)
    summary["wall_time_s_total_jobs"] = wall_total
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, default=_json_default, indent=1)
    with open(out_dir / "stability_grid.json", "w") as fh:
        json.dump(stab, fh, default=_json_default)
    return summary


SCALAR_KEYS = None  # all


def _sel(recs, ids, **kw):
    out = []
    for i in ids:
        r = recs[i]
        if all((r.get(k) == v) for k, v in kw.items()):
            out.append(i)
    return out


def scalar_bands(recs, sel):
    keys = sorted({k for i in sel for k in recs[i]["scalars"]})
    return {k: quantiles([recs[i]["scalars"].get(k) if isinstance(recs[i]["scalars"].get(k), (int, float))
                          else None for i in sel]) for k in keys}


def spectral_bands(arrays, ids, sel, fb):
    pos = [ids.index(i) for i in sel]
    out = {}
    for n in SPEC_SIGNALS:
        a = np.array(arrays[n])[pos] if pos else np.zeros((0, fb.size))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[n] = {"median": np.nanmedian(a, axis=0).tolist() if pos else [],
                      "p05": np.nanpercentile(a, 5, axis=0).tolist() if pos else [],
                      "p95": np.nanpercentile(a, 95, axis=0).tolist() if pos else []}
    return out


def peak_summary(recs, sel):
    """Dominant peak (largest integrated excess power) per run and signal, and all peaks pooled by decade."""
    out = {}
    for n in SPEC_SIGNALS:
        dom, pooled = [], []
        for i in sel:
            pk = recs[i]["peaks"][n]
            pooled += pk
            if pk:
                dom.append(max(pk, key=lambda d: d["power"]))
        out[n] = {"n_runs": len(sel), "n_runs_with_peak": len(dom),
                  "dominant_f_hz": quantiles([d["f_hz"] for d in dom]),
                  "dominant_excess_db": quantiles([d["excess_db"] for d in dom]),
                  "dominant_psd": quantiles([d["psd"] for d in dom]),
                  "dominant_fwhm_hz": quantiles([d["fwhm_hz"] for d in dom]),
                  "dominant_rms": quantiles([math.sqrt(d["power"]) for d in dom]),
                  "pooled_by_decade": {f"1e{k}": sum(1 for d in pooled if 10**k <= d["f_hz"] < 10**(k + 1))
                                       for k in range(0, 6)}}
    return out


def summarise(recs, ids, fb, arrays, axes):
    s = {"f_bin_hz": fb.tolist(), "groups": {}}
    defs = {
        "nominal_q1": dict(group="nominal", quant=True), "nominal_q0": dict(group="nominal", quant=False),
        "unlocked": dict(group="unlocked"), "fsim2": dict(group="fsim2"), "resgrid": dict(group="resgrid"),
    }
    for g, kw in defs.items():
        sel = _sel(recs, ids, **kw)
        if sel:
            s["groups"][g] = {"ids": sel, "scalars": scalar_bands(recs, sel), "peaks": peak_summary(recs, sel),
                              "spectra": spectral_bands(arrays, ids, sel, fb),
                              "peaks_per_run": {i: recs[i]["peaks"] for i in sel}}
    for g in ("lhs", "oat"):
        for q in (True, False):
            allsel = _sel(recs, ids, group=g, quant=q)
            if not allsel:
                continue
            cls = {"stable": lambda r: r["stable"], "unstable": lambda r: not r["stable"], "all": lambda r: True,
                   "stable_locked": lambda r: r["stable"] and r["locked"],
                   "stable_notlocked": lambda r: r["stable"] and not r["locked"],
                   "unstable_locked": lambda r: (not r["stable"]) and r["locked"]}
            for st, fn in cls.items():
                sel = [i for i in allsel if fn(recs[i])]
                key = f"{g}_q{int(q)}_{st}"
                s["groups"][key] = {"ids": sel, "scalars": scalar_bands(recs, sel), "peaks": peak_summary(recs, sel)}
                if g == "lhs":
                    s["groups"][key]["spectra"] = spectral_bands(arrays, ids, sel, fb)
            if g == "lhs":
                # T15 pairing of the DSC1 latency model (5.7 µs + poles / 8.3 µs pure delay), stable only
                sel = [i for i in allsel if recs[i]["stability"]["stable"] and
                       (recs[i]["values"]["latency_s"] < 7e-6) == bool(recs[i]["values"]["in_out_poles"])]
                s["groups"][f"lhs_q{int(q)}_stable_t15paired"] = {"ids": sel, "scalars": scalar_bands(recs, sel),
                                                                 "peaks": peak_summary(recs, sel)}
    # OAT table: per axis value the key observables (quant on / off)
    keys = ["err.std", "err.offset", "err.branch", "lock.lockD.err.rms_1_200000", "v.pdb.std", "v.dsc1.out.std",
            "v.dsc1.out.ptp_10ms_median", "cv.frac_saturated", "sup.n_loss", "sup.t_loss", "t_loss_extrap_s"]
    oat = {}
    for i in _sel(recs, ids, group="oat"):
        r = recs[i]
        ax = r["axis"]
        val = r["values"][ax]
        d = oat.setdefault(ax, {}).setdefault(str(val), {})
        d[f"q{int(r['quant'])}"] = {**{k: r["scalars"].get(k) for k in keys}, "stable": r["stability"]["stable"],
                                    "max_l_at_180": r["stability"]["max_l_at_180"],
                                    "dsc1_out_dominant_peak_hz": (max(r["peaks"]["v.dsc1.out"], key=lambda d: d["power"])["f_hz"]
                                                                  if r["peaks"]["v.dsc1.out"] else None)}
    s["oat_table"] = oat
    # quantisation ablation: paired differences (same point, same seed)
    abl = {}
    for g in ("nominal", "oat", "lhs"):
        pairs = []
        for i in _sel(recs, ids, group=g, quant=True):
            j = i[:-1] + "0"
            if j in recs:
                pairs.append((i, j))
        if not pairs:
            continue
        def d2(k):
            out = []
            for i, j in pairs:
                a, b = recs[i]["scalars"].get(k), recs[j]["scalars"].get(k)
                if a is not None and b is not None:
                    out.append(math.copysign(math.sqrt(abs(a * a - b * b)), a * a - b * b))
            return out
        stable_pairs = [(i, j) for i, j in pairs if recs[i]["stability"]["stable"]]
        abl[g] = {"n_pairs": len(pairs), "n_stable_pairs": len(stable_pairs),
                  "quadrature_increment_err_std_rad": quantiles(d2("err.std")),
                  "quadrature_increment_err_std_rad_stable": quantiles(
                      [math.copysign(math.sqrt(abs(recs[i]["scalars"]["err.std"]**2 - recs[j]["scalars"]["err.std"]**2)),
                                     recs[i]["scalars"]["err.std"]**2 - recs[j]["scalars"]["err.std"]**2)
                       for i, j in stable_pairs]),
                  "ratio_err_std_on_off_stable": quantiles([recs[i]["scalars"]["err.std"] / recs[j]["scalars"]["err.std"]
                                                            for i, j in stable_pairs]),
                  "ratio_dsc1_out_std_on_off_stable": quantiles(
                      [recs[i]["scalars"]["v.dsc1.out.std"] / recs[j]["scalars"]["v.dsc1.out.std"] for i, j in stable_pairs]),
                  "ratio_pdb_std_on_off_stable": quantiles(
                      [recs[i]["scalars"]["v.pdb.std"] / recs[j]["scalars"]["v.pdb.std"] for i, j in stable_pairs])}
        if g == "lhs":
            abl[g]["per_point"] = [{"lhs_index": recs[i]["lhs_index"], "stable": recs[i]["stability"]["stable"],
                                    "dL_m": recs[i]["values"]["dL_m"], "f_res_hz": recs[i]["values"]["f_res_hz"],
                                    "q_res": recs[i]["values"]["q_res"],
                                    "err_std_on": recs[i]["scalars"]["err.std"], "err_std_off": recs[j]["scalars"]["err.std"]}
                                   for i, j in pairs]
    s["quantisation_ablation"] = abl
    # LHS stability classification in (f_res, Q)
    lhs = _sel(recs, ids, group="lhs", quant=True)
    s["lhs_stability"] = {"n": len(lhs), "n_stable": sum(recs[i]["stability"]["stable"] for i in lhs),
                          "points": [{"lhs_index": recs[i]["lhs_index"], "f_res_hz": recs[i]["values"]["f_res_hz"],
                                      "q_res": recs[i]["values"]["q_res"], "latency_s": recs[i]["values"]["latency_s"],
                                      "in_out_poles": recs[i]["values"]["in_out_poles"],
                                      "stable": recs[i]["stability"]["stable"],
                                      "max_l_at_180": recs[i]["stability"]["max_l_at_180"],
                                      "time_loop_locked": recs[i]["locked"], "err_std": recs[i]["scalars"]["err.std"]}
                                     for i in lhs]}
    # analytic verdict vs time-loop outcome (every feedback-on run)
    cont = {}
    for g in ("nominal", "oat", "lhs", "resgrid", "fsim2"):
        sel = [i for i in _sel(recs, ids, group=g) if recs[i]["feedback"]]
        if sel:
            cont[g] = {f"{'stable' if a else 'unstable'}_{'locked' if b else 'notlocked'}":
                       sum(1 for i in sel if recs[i]["stable"] == a and recs[i]["locked"] == b)
                       for a in (True, False) for b in (True, False)}
            cont[g]["disagreeing_runs"] = [i for i in sel if recs[i]["stable"] != recs[i]["locked"]]
    s["verdict_vs_time_loop"] = cont
    # convergence: f_sim x2 and seeds (nominal, quantisation on)
    conv = {}
    nq = _sel(recs, ids, group="nominal", quant=True)
    f2 = _sel(recs, ids, group="fsim2")
    ckeys = ["err.std", "err.offset", "lock.lockD.err.rms_1_200000", "lock.lockD.err.rms_1000_10000",
             "lock.lockD.err.rms_10000_100000", "lock.lockD.err.rms_100000_200000", "v.pdb.std",
             "v.pdb.rms_10000_100000", "v.dsc1.out.std", "v.dsc1.out.ptp_10ms_median", "v.dsc1.out.rms_10000_100000"]
    for k in ckeys:
        a = [recs[i]["scalars"][k] for i in nq]
        b = [recs[i]["scalars"][k] for i in f2]
        if a and b:
            ma, sa = float(np.mean(a)), float(np.std(a, ddof=1))
            mb, sb = float(np.mean(b)), float(np.std(b, ddof=1)) if len(b) > 1 else 0.0
            conv[k] = {"f1MHz_mean": ma, "f1MHz_seed_std": sa, "f1MHz_sem": sa / math.sqrt(len(a)),
                       "f2MHz_mean": mb, "f2MHz_seed_std": sb, "rel_change": (mb - ma) / ma if ma else None,
                       "change_in_sem_units": (mb - ma) / math.hypot(sa / math.sqrt(len(a)), sb / math.sqrt(len(b)))
                       if (sa or sb) else None}
    # spectral f_sim x2 comparison: ratio of seed-mean PSDs per decade band
    pos1 = [ids.index(i) for i in nq]
    pos2 = [ids.index(i) for i in f2]
    spec_conv = {}
    if pos1 and pos2:
        for n in SPEC_SIGNALS:
            a1 = np.nanmean(np.array(arrays[n])[pos1], axis=0)
            a2 = np.nanmean(np.array(arrays[n])[pos2], axis=0)
            spec_conv[n] = {f"{lo:g}-{hi:g}": float(np.nanmean(a2[(fb >= lo) & (fb < hi)]) /
                                                    np.nanmean(a1[(fb >= lo) & (fb < hi)]))
                            for lo, hi in [(1, 10), (10, 100), (100, 1e3), (1e3, 1e4), (1e4, 5e4), (5e4, 1.25e5),
                                           (1.25e5, 2e5)]}
    s["convergence"] = {"fsim2_vs_fsim1_scalars": conv, "fsim2_vs_fsim1_psd_ratio_by_band": spec_conv,
                        "n_seeds_1MHz": len(nq), "n_seeds_2MHz": len(f2)}
    # lock-loss over the documented drift range (nominal u0, model formula)
    u0 = float(np.mean([recs[i]["scalars"]["cv.u0"] for i in nq])) if nq else -0.25
    nom = {a["name"]: a["nominal"] for a in axes}
    table = []
    for dl in (5.0, 9.0, 15.0):
        for rt in (1e-5, 1e-4, 1e-3):
            for tau in (1.0, 300.0, 3000.0):
                for kt in (40.0, 44.0, 48.0):
                    table.append({"dL_m": dl, "r_T_k_per_s": rt, "tau_th_s": tau, "k_T": kt,
                                  "t_loss_s": t_loss_extrapolated({**nom, "dL_m": dl, "r_T_k_per_s": rt,
                                                                   "tau_th_s": tau, "k_T_rad_per_m_K": kt}, u0)})
    s["lock_loss_table"] = {"u0_v": u0, "rows": table}
    # extrapolation vs simulated loss (points that lost lock within the run while stable)
    chk = []
    for i in ids:
        r = recs[i]
        if r["feedback"] and r["scalars"]["sup.n_loss"] > 0 and r["stability"]["stable"]:
            chk.append({"id": i, "t_loss_sim": r["scalars"]["sup.t_loss"], "t_loss_extrap": r["scalars"]["t_loss_extrap_s"]})
    s["lock_loss_extrap_check"] = chk
    return s


def summarise_stab(stab):
    out = {}
    f, q = np.array(stab["f_res_hz"]), np.array(stab["q_res"])
    for name, c in stab["cases"].items():
        st = np.array(c["stable"])
        # for each f_res: largest stable Q (None if none), and smallest f_res stable at every Q
        qmax = [float(q[st[i]].max()) if st[i].any() else None for i in range(f.size)]
        all_q = [float(f[i]) for i in range(f.size) if st[i].all()]
        out[name] = {"k0": c["k0"], "frac_stable": float(st.mean()), "q_max_stable_per_f_res": qmax,
                     "f_res_min_stable_all_q": min(all_q) if all_q else None,
                     "f_res_min_stable_any_q": float(f[st.any(axis=1)].min()) if st.any() else None}
    return out


# ----------------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--groups", default="nominal,unlocked,fsim2,oat,lhs,resgrid,stabgrid")
    ap.add_argument("--cache", type=Path, default=HERE / "_cache")
    ap.add_argument("--out", type=Path, default=HERE / "results")
    ap.add_argument("--duration", type=float, default=None, help="smoke tests only (default: the preset's 60 s)")
    ap.add_argument("--limit", type=int, default=None, help="smoke tests only: first N jobs")
    ap.add_argument("--aggregate-only", action="store_true")
    a = ap.parse_args(argv)
    os.chdir(ROOT)
    warnings.simplefilter("ignore")
    groups = set(a.groups.split(","))
    a.cache.mkdir(parents=True, exist_ok=True)
    axes = axes_of_preset()
    jobs = jobs_for(groups, axes)
    if a.limit:
        jobs = jobs[: a.limit]
    t_start = time.time()
    stab_path = a.cache / "stability_grid.json"
    if not a.aggregate_only:
        todo = [j for j in jobs if not (a.cache / f"{j['id']}.json").exists()]
        print(f"{len(jobs)} jobs, {len(todo)} to run, {a.workers} workers", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(run_job, j, str(a.cache), a.duration): j["id"] for j in todo}
            if "stabgrid" in groups and not stab_path.exists():
                futs[ex.submit(stability_grid)] = "stabgrid"
            for k, fu in enumerate(as_completed(futs), 1):
                name = futs[fu]
                r = fu.result()
                if name == "stabgrid":
                    stab_path.write_text(json.dumps(r))
                    print(f"[{k}/{len(futs)}] stabgrid done", flush=True)
                else:
                    print(f"[{k}/{len(futs)}] {r[0]} {r[1]:.0f} s (elapsed {time.time() - t_start:.0f} s)", flush=True)
    stab = json.loads(stab_path.read_text()) if stab_path.exists() else {"f_res_hz": [], "q_res": [], "cases": {}}
    walls = [json.loads(p.read_text())["wall_s"] for p in a.cache.glob("*_*.json") if p.name != "stability_grid.json"]
    a.out.mkdir(parents=True, exist_ok=True)
    u, pts = lhs_design(axes, LHS_N, LHS_SEED)
    design = {"axes": axes, "n_axes": len(axes), "lhs": {"n": LHS_N, "seed": LHS_SEED, "method":
              "scipy.stats.qmc.LatinHypercube (centered=False); sweep: lo + u(hi−lo), sweep_log: lo(hi/lo)^u, "
              "choice: options[floor(u·n)]", "unit_cube": u.tolist(), "points": pts},
              "nominal_seeds": NOMINAL_SEEDS, "unlocked_seeds": UNLOCKED_SEEDS, "fsim2_seeds": FSIM2_SEEDS,
              "sweep_seed_rule": "run.seed_mode per_point: cranberry.runner.sweep.point_seed(20261031, coordinates)",
              "jobs": [{k: v for k, v in j.items() if k != "values"} for j in jobs],
              "analysis": {"T_settle_s": T_SETTLE_S, "f_record_hz": F_REC_HZ, "record_antialias": "fir",
                           "welch_nperseg": NPERSEG, "welch": "hann, linear detrend per segment, 50 % overlap",
                           "fir_correction": "PSD / |H_fir(f)|^2 (cranberry.runner.run.antialias_taps)",
                           "bins_per_decade": 100, "peak_rule": {"min_excess_db": PEAK_MIN_DB, "n_sigma": PEAK_NSIGMA,
                                                                  "baseline_halfwidth_decades": PEAK_BASE_DEC},
                           "report_band_hz": F_REPORT},
              "provenance": {"git_commit": git("rev-parse", "HEAD"), "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
                             "git_dirty_outside_predictions": bool(git("status", "--porcelain", "--", "src", "presets",
                                                                       "devices", "tests")),
                             "preset_sha256": {REFD: sha256(ROOT / REFD), NOMINAL: sha256(ROOT / NOMINAL),
                                               "devices/dsc1.yaml": sha256(ROOT / "devices/dsc1.yaml"),
                                               "devices/pdb450c.yaml": sha256(ROOT / "devices/pdb450c.yaml")},
                             "python": sys.version, "platform": platform.platform(),
                             "blind_read": False, "parameters_tuned": False}}
    with open(a.out / "design.json", "w") as fh:
        json.dump(design, fh, default=_json_default, indent=1)
    summary = aggregate(a.cache, a.out, axes, stab, float(np.sum(walls)) if walls else 0.0)
    print(json.dumps({k: summary["groups"][k]["scalars"]["err.std"] for k in summary["groups"]
                      if "err.std" in summary["groups"][k]["scalars"]}, indent=0, default=_json_default)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
