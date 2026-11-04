"""Shared helpers of the Phase 6 axis-1 scripts (reviews/scripts/p6_axis1_*.py, notes/findings_axis1.md).

* :func:`point` — one runnable spec: base preset + overlays, run settings, every listed sweep axis pinned
  (``runner.sweep.expand``; dependent ``!expr`` values follow). Raises when an axis is left open.
* :func:`job` / :func:`pool` — run many (label, spec-args) jobs in at most ``N_PROC`` processes (CPU budget of the
  task: ≤ 6), each returning ``explore.run_quantities`` plus the requested series statistics.
* :func:`stats` — mean / sample std / min / max per quantity over seeds (the ``explore.compare_seeds`` reduction, done
  here because the runs are spread over processes).
Raw results go to ``runs/p6_axis1/`` (gitignored); summaries are printed and saved as small JSON next to the figures.
"""
from __future__ import annotations

import json
import math
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "p6_axis1"
FIG = ROOT / "docs" / "figures" / "phase6" / "axis1"
N_PROC = int(os.environ.get("P6_NPROC", "6"))


def point(base: str, overlays=(), axes: dict | None = None, *, duration_s: float | None = None, seed: int = 1,
          extra: dict | None = None, measure_esa_off: bool = True):
    """Runnable spec (module docstring). ``extra``: a dict merged on top as a last overlay (e.g. pipeline params)."""
    from cranberry import config
    from cranberry.runner.sweep import expand

    doc: dict[str, Any] = {"extends": str(ROOT / base) if not base.startswith("/") else base,
                           "run": {"seed": int(seed), "seed_mode": "common"}}
    if duration_s is not None:
        doc["run"]["duration_s"] = float(duration_s)
    if measure_esa_off:
        doc["measure"] = {"esa": None}
    ovs = list(overlays)
    if extra:
        RUNS.mkdir(parents=True, exist_ok=True)
        import hashlib

        import yaml

        txt = yaml.safe_dump(extra, sort_keys=True)
        p = RUNS / f"extra_{hashlib.sha1(txt.encode()).hexdigest()[:12]}.yaml"
        p.write_text(txt)
        ovs.append(str(p))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = config.load(doc, overlays=ovs)
        pts = expand(spec, {k: [v] for k, v in (axes or {}).items()})
    if len(pts) != 1:
        left = sorted(config.sweep_markers(spec))
        raise ValueError(f"{base} {ovs}: {len(pts)} points, open axes {left}")
    return pts[0].spec


def _series_stats(res, slots, t_from):
    out = {}
    t = np.asarray(res.t_s)
    for s in slots:
        if s in res.series:
            x = np.asarray(res.series[s], dtype=float)[t >= t_from]
            if x.size:
                out[f"series.{s}.std"] = float(np.std(x))
                out[f"series.{s}.mean"] = float(np.mean(x))
                out[f"series.{s}.absmax"] = float(np.max(np.abs(x)))
    return out


def job(args: dict) -> dict:
    """Worker: build ``point(**args['point'])``, run it, return quantities (+ series stats, wall time, events)."""
    from cranberry.runner.explore import run_quantities, with_recorded, actuator_ranges
    from cranberry.runner.run import run

    t0 = time.perf_counter()
    try:
        pt = dict(args["point"])
        tog = {**pt.pop("_toggles", {}), **(args.get("toggles") or {})}
        spec = point(**pt)
        if tog:
            args = {**args, "toggles": tog}
            import dataclasses

            spec = dataclasses.replace(spec, toggles={**dict(spec.toggles), **args["toggles"]})
        spec = with_recorded(spec, [s for s, _, _ in actuator_ranges(spec).values()] + list(args.get("record", [])))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = run(spec)
        q = run_quantities(res)
        q.update(_series_stats(res, args.get("record", []), args.get("t_from", 0.0)))
        for k in [k[: -len(".std_rad")] for k in q if k.endswith(".std_rad")]:
            if f"{k}.sigma_total_rad" not in q and f"{k}.sigma_hf_rad" in q:   # cw_sideband metrics (refC): no V±
                q[f"{k}.sigma_phase_rad"] = math.hypot(q[f"{k}.std_rad"], q[f"{k}.sigma_hf_rad"])   # excl. the dither
        q["n_events"] = float(len(res.events))
        q["events"] = sorted({str(getattr(e, "code", e))[:60] for e in res.events})[:8]
        if args.get("psd"):
            from cranberry.core.fft import psd_onesided

            psd = {}
            t = np.asarray(res.t_s)
            fs = 1.0 / float(t[1] - t[0])
            for s in args["psd"]:
                if s not in res.series:
                    continue
                x = np.asarray(res.series[s], dtype=float)[t >= args.get("t_from", 0.0)]
                f, p = psd_onesided(x - x.mean(), fs, nperseg=args.get("nperseg"))
                psd[s] = (f.tolist(), p.tolist())
            q["psd"] = psd
        if args.get("decompose"):
            from cranberry.runner.explore import _decomposed

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                q.update(_decomposed(spec, "auto"))
        if args.get("measurements"):
            q["measurements"] = _jsonable(res.measurements)
    except Exception as e:  # noqa: BLE001 - a failed point is data (D75)
        q = {"error": f"{type(e).__name__}: {e}"[:400]}
    q["wall_s"] = time.perf_counter() - t0
    q["label"] = args.get("label")
    q["seed"] = args["point"].get("seed")
    return q


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if np.iscomplexobj(x):
            return {"re": x.real.tolist(), "im": x.imag.tolist()}
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, complex):
        return {"re": x.real, "im": x.imag}
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    return repr(x)


def pool(jobs: list[dict], name: str, n_proc: int = N_PROC) -> list[dict]:
    """Run ``jobs`` in parallel, save the raw list to runs/p6_axis1/<name>.json and return it."""
    RUNS.mkdir(parents=True, exist_ok=True)
    out = []
    with ProcessPoolExecutor(max_workers=min(n_proc, max(len(jobs), 1))) as ex:
        for i, r in enumerate(ex.map(job, jobs)):
            out.append(r)
            err = r.get("error")
            print(f"  [{i + 1}/{len(jobs)}] {r['label']} seed {r['seed']}: {r['wall_s']:.0f} s"
                  + (f" ERROR {err}" if err else ""), flush=True)
    (RUNS / f"{name}.json").write_text(json.dumps(out))
    return out


def load(name: str) -> list[dict]:
    return json.loads((RUNS / f"{name}.json").read_text())


def stats(rows: list[dict], label: str, keys) -> dict[str, dict[str, float]]:
    """{key: {mean, std, min, max, n, values}} of the rows with this label (errors skipped)."""
    sel = [r for r in rows if r.get("label") == label and "error" not in r]
    out = {}
    for k in keys:
        v = np.array([r.get(k, math.nan) for r in sel], dtype=float)
        if v.size:
            out[k] = {"mean": float(np.nanmean(v)), "std": float(np.nanstd(v, ddof=1)) if v.size > 1 else 0.0,
                      "min": float(np.nanmin(v)), "max": float(np.nanmax(v)), "n": int(v.size), "values": v.tolist()}
    return out


def differs(a, b, k=1.0):
    return abs(a["mean"] - b["mean"]) > k * (a["std"] + b["std"])


def fmt(s, scale=1.0, nd=2):
    return f"{s['mean'] * scale:.{nd}f} ± {s['std'] * scale:.{nd}f}"


# ------------------------------------------------------------------------------------------------
# the four bases at their headline points (notes/findings_axis1.md §0)

A_BASE = "presets/bench/refA_fiber_bench.yaml"
C_BASE = "presets/bench/refC_fiber_bench.yaml"
P_BASE = "presets/bench/pulsed_wg_1port.yaml"
D_BASE = "presets/bench/refD_lab.yaml"
T_SETTLE_A = 0.015   # D97-1: bench statistics from t >= 15 ms
A_EXTRA = {"pipeline": [{"name": "mi", "params": {"t_settle_s": T_SETTLE_A}}]}
# lab/acoustic(-type) overlays pinned at the D89 calibration (scale 1, gamma 1, f^-2 slope, rho 0, no DAC dither)
LAB_ACOUSTIC = {"offline.params.acoustic_scale": 1.0, "offline.params.acoustic_gamma": 1.0}
LAB_FUTURE = {**LAB_ACOUSTIC, "offline.params.acoustic_hf_slope": 2.0, "offline.params.acoustic_rho_arms": 0.0,
              "offline.params.dac_output_noise_v": 0.0,
              "offline.params.pzt_f_res_hz": 2.52e4, "offline.params.pzt_q": 25.0,
              "offline.params.mm2_f_hz": 5.4e3, "offline.params.mm2_q": 15.0, "offline.params.mm2_r": 0.11,
              "offline.params.mm3_f_hz": 1.0e4, "offline.params.mm3_q": 30.0, "offline.params.mm3_r": 0.07,
              "offline.params.mm4_f_hz": 4.03e4, "offline.params.mm4_q": 100.0, "offline.params.mm4_r": 0.15}
D_PIN = {"offline.params.acoustic_hf_slope": 2.0, "offline.params.acoustic_rho_arms": 0.0,
         "offline.params.dac_output_noise_v": 0.0}


def a_point(L, route="fiber_matched", *, dx1=False, acoustic=False, ctrl=None, offset=None, eta_out=0.98, seed=1,
            duration_s=0.05, extra_overlays=(), axes=None, pipeline=(), toggles=None):
    ovs, ax = [], {"offline.params.L_sq_fiber_m": L, "offline.params.lo_route": route,
                   "offline.params.eta_fiber_out": eta_out}
    if ctrl:
        ovs.append(ctrl)
    if dx1:
        ovs.append("bench/refA_laser_dx1")
    if acoustic:
        ovs.append("bench/refA_fiber_acoustic")
        ax.update(LAB_ACOUSTIC)
    if offset:
        ovs.append("bench/refA_sq_lo_offset")
        ax["offline.params.sq_lo_offset"] = offset
    ovs += list(extra_overlays)
    ax.update(axes or {})
    out = dict(base=A_BASE, overlays=ovs, axes=ax, duration_s=duration_s, seed=seed,
               extra={"pipeline": A_EXTRA["pipeline"] + list(pipeline)})
    if toggles:
        out["_toggles"] = dict(toggles)
    return out


def c_point(L, *, ctrl=None, seed=1, duration_s=0.1, axes=None):
    return dict(base=C_BASE, overlays=[ctrl] if ctrl else [], axes={"offline.params.L_fiber_m": L, **(axes or {})},
                duration_s=duration_s, seed=seed)


def p_point(L, lock=(), *, ctrl=None, lab=False, seed=1, duration_s=0.05, axes=None):
    ovs = list(lock)
    ax = {"offline.params.L_port_m": L}
    if ctrl:
        ovs.append(ctrl)
        if ctrl == "controller/dsc1":
            ax["offline.params.dsc1_latency_s"] = 5.7e-6
    if lab:
        ovs += ["lab/acoustic_future", "lab/multimode915b_future"]
        ax.update(LAB_FUTURE)
    ax.update(axes or {})
    return dict(base=P_BASE, overlays=ovs, axes=ax, duration_s=duration_s, seed=seed)


def d_point(*, seed=1, duration_s=0.2, axes=None, overlays=()):
    return dict(base=D_BASE, overlays=list(overlays), axes={**D_PIN, **(axes or {})}, duration_s=duration_s, seed=seed)
