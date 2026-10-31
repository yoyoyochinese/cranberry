"""Shared helpers of the future-preset tests (WP-I5): spec variants of the F_* presets, the manifest, CHK collection,
stray-light / GAWBS readouts and the 32-port scale estimate. Not a test module (loaded by path: ``tests`` is not a
package under ``--import-mode=importlib``).

Spec variants go through ``config.load`` with an inline tree that ``extends`` a preset, so every ``!expr`` (port
lengths, loop gains, τ_L − 1/f_sim, coefficients) is re-evaluated (D64).
"""

from __future__ import annotations

import copy
import dataclasses
import math
import re
import time
import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

from cranberry import config
from cranberry.config.spec import RunSpec
from cranberry.config.yaml_io import load_yaml

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "presets" / "future" / "manifest.yaml"
H_PLANCK = 6.62607015e-34
C_LIGHT = 299792458.0


def manifest() -> dict[str, Any]:
    return load_yaml(MANIFEST)


def nominal(base: str) -> str:
    return f"presets/future/{base}_nominal.yaml"


def combinations() -> list[tuple[str, tuple[str, ...]]]:
    """Every (preset file, overlays) the manifest lists: each base alone, each overlay alone on every applicable base,
    the multi-overlay combinations and the bad presets."""
    m = manifest()
    out: list[tuple[str, tuple[str, ...]]] = []
    for b in m["base"]:
        out.append((nominal(b), ()))
        extra = m["overlays_opo"] if b in m["opo"] else m["overlays_waveguide"]
        out += [(nominal(b), (o,)) for o in [*m["overlays_all"], *extra]]
    out += [(nominal(c["base"]), tuple(c["overlays"])) for c in m["combinations"]]
    out += [(f"presets/bad/{b}.yaml", ()) for b in m["bad"]]
    return out


def load(preset: str, params: dict[str, Any] | None = None, run: dict[str, Any] | None = None,
         toggles: dict[str, bool] | None = None, overlays: tuple[str, ...] = (), pipeline: list | None = None,
         measure: dict[str, Any] | None = None) -> RunSpec:
    """``preset`` with ``offline.params`` / ``run`` / pipeline (keyed by name) / ``measure`` overrides and overlays."""
    tree: dict[str, Any] = {"extends": str(ROOT / preset)}
    if params:
        tree["offline"] = {"params": dict(params)}
    if run:
        tree["run"] = dict(run)
    if pipeline:
        tree["pipeline"] = list(pipeline)
    if measure:
        tree["measure"] = dict(measure)
    spec = config.load(tree, overlays=list(overlays))
    if toggles:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **toggles})
    return spec


def first_point(spec: RunSpec) -> RunSpec:
    """The spec itself, or its first sweep point when overlays left axes open."""
    if not config.sweep_markers(spec):
        return spec
    from cranberry.runner.sweep import expand

    return expand(spec)[0].spec


def build_check(preset: str, overlays: tuple[str, ...], params: dict[str, Any] | None = None):
    """Build and validate ``preset`` + ``overlays`` through the full runner (offline stage, topology, pipeline checks,
    metrics) for three samples of the pure-Python loop (no compilation)."""
    from cranberry.runner.run import run

    settle0 = [{"name": "mi", "params": {"t_settle_s": 0.0}}, {"name": "mii", "params": {"t_settle_s": 0.0}}]
    spec = first_point(load(preset, params=params, run={"duration_s": 3.0e-6}, overlays=overlays, pipeline=settle0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(spec, mode="reference")


def chk_codes(result) -> list[str]:
    """CHK event codes of a run (G16-j): ``RunResult.events`` — the construction-time events of the offline blocks
    and the ``check()`` of the blocks that check on demand (demux_design, sbs_check, stray_light_table with a
    ``limit``), which the runner calls itself (Phase 5 M6-12)."""
    return [e["code"] for e in result.events]


def block(result, name: str):
    return next(b for b in result.offline.blocks if b.name == name)


def params(result) -> dict[str, Any]:
    return dict(result.spec.offline["params"])


def stray_light(result, name: str = "stray1") -> dict[str, float]:
    """Stray photons of the probe in the squeezed direction at the run's P̄ (B10 ``stray_photons``; D46 units)."""
    from cranberry.blocks.fib.scattering import stray_photons

    b = block(result, name)
    p = params(result)
    return stray_photons(b.params, b._toggles, float(p["p_pr_avg_w"]))


def gawbs_photons_per_mode(result, power_w: float, name: str = "gawbs1") -> float:
    """Forward-GAWBS photons per squeezed mode from a co-propagating field of ``power_w`` (T17 §6:
    N_G = Σ_m (P σ_m²/(2hν0))(|ũ(f_m)|² + |ũ(−f_m)|²)), with the source's optimal mode."""
    from cranberry.blocks.source import temporal_modes as tm

    t = block(result, name).tables()["gawbs_peaks"]
    f_m, s2 = t.values[:, 0], t.values[:, 1]
    so = result.offline.source_outputs["src"]
    if so.basis == "temporal":
        u, _ = tm.optimal_lo_mode(so)
        fg, sg = tm.mode_spectrum(u, so.grid, 16)
        both = np.interp(f_m, fg, sg, left=0.0, right=0.0) + np.interp(-f_m, fg, sg, left=0.0, right=0.0)
    else:  # single-mode stand-in of a ps pulse: |ũ(f)|² ≈ τ_eff ≈ τ_mode, flat over the GAWBS lines (≪ 1/τ_p)
        both = 2.0 * float(params(result)["tau_mode_s"]) * np.ones_like(f_m)
    hnu = H_PLANCK * C_LIGHT / float(params(result)["lambda0_m"])
    return float(np.sum(power_w * s2 / (2 * hnu) * both))


def slope(t: np.ndarray, x: np.ndarray, t_from: float) -> float:
    """Least-squares slope of x(t) for t ≥ t_from."""
    m = t >= t_from
    return float(np.polyfit(t[m], x[m], 1)[0])


# ------------------------------------------------------------------------------------------- closed-loop path model


def liu_lock_ratio(p: dict[str, Any], f_u_hz: float, f_sim_hz: float) -> float:
    """R(f_u) = locked / unlocked white-FM variance of Θ_k (time loop ≤ f_N plus the same ledger part above it) from
    the closed-loop path model of the Phase 5 M6 review §2.2 (``reviews/scripts/m6_liu_lock_long_fiber.py``; the
    M4 review §2.1 form with the actuator upstream, T_px = 0): per port j
    θ_j = e^{−iω(T_a + τ_j)} − T(ω)·[e^{−iω(T_b + τ_j)} − e^{−iωT_ref}]·e^{−iωτ_j}, T = L/(1 + L),
    L = (2π f_u/s)·e^{−sτ_L} (no stretcher resonance), τ_L = latency + 0.5/f_sim. Delays from the F_common
    geometry (offline params ``p``): T_a squeezer → demux, T_b probe source → PBS_k, T_ref reference → det_k,
    τ_j = D_j + F_j (the common section)."""
    c = C_LIGHT
    t_a = (1.0003 * (p["L_pump_m"] + p["L_up_m"]) + p["n_g_xtal"] * p["L_xtal_m"]) / c
    t_b = p["n_g"] * (p["L_prsrc_m"] + p["L_prs_m"]) / c
    t_ref = p["n_g"] * (p["L_ref_m"] + p["L_refk_m"]) / c
    tau_k = p["tau_k_s"] + 1.0003 * p["L_d_m"] / c
    tau_0 = p["tau_0_s"] + 1.0003 * p["L_d_m"] / c
    f_n = 0.5 * f_sim_hz
    f = np.geomspace(1.0, f_n, 400001)
    w = 2 * np.pi * f
    lg = (2 * np.pi * f_u_hz / (1j * w)) * np.exp(-1j * w * (p["latency_s"] + 0.5 / f_sim_hz))
    t = lg / (1 + lg)

    def th(tau):
        return (np.exp(-1j * w * (t_a + tau))
                - t * (np.exp(-1j * w * (t_b + tau)) - np.exp(-1j * w * t_ref)) * np.exp(-1j * w * tau))

    s_phi = 1.0 / f**2                                            # h0 cancels
    free = np.trapezoid(np.abs(np.exp(-1j * w * tau_k) - np.exp(-1j * w * tau_0)) ** 2 * s_phi, f)
    locked = np.trapezoid(np.abs(th(tau_k) - th(tau_0)) ** 2 * s_phi, f)
    total = 2 * np.pi**2 * (tau_k - tau_0)                       # 2π² Δτ (per unit h0): band + ledger
    return (total - free + locked) / total


# ------------------------------------------------------------------------------------------- 32-port scale estimate

_PORT1 = re.compile(r"\b(DMX|D|F|G|PRS|REF|LO|sqz|pr|ref|lo|det|L|liu|adc|pid|sup|dac|dl|pzt|tec|stage|act|lock|delay)1(?![0-9])")
_PORT1_SEGS = ("DMX1", "D1", "F1", "G1", "PRS1", "REF1", "LO1")


def _to_port2(x: Any) -> Any:
    """Copy of a config subtree with every port-1 name (segments, beams, blocks, slots) renamed to port 2; expressions
    are kept (port 2 duplicates port k: same lengths, gains and noise)."""
    if isinstance(x, str):
        return _PORT1.sub(r"\g<1>2", x)
    if isinstance(x, Mapping):
        return {_to_port2(k): _to_port2(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_to_port2(v) for v in x]
    return copy.deepcopy(x)


def _add_segment_copies(params: dict[str, Any]) -> None:
    """Noise / thermal block params: add the port-2 copy of every port-1 segment in ``segments`` and in per-segment
    mappings (tau_th_s, spatial_coherence)."""
    segs = params.get("segments")
    if not isinstance(segs, list) or not any(x in _PORT1_SEGS for x in segs):
        return
    params["segments"] = segs + [_to_port2(x) for x in segs if x in _PORT1_SEGS]
    for v in params.values():
        if isinstance(v, dict):
            v.update({_to_port2(k): copy.deepcopy(val) for k, val in list(v.items()) if k in _PORT1_SEGS})


def three_port_tree(preset: str) -> dict[str, Any]:
    """The unresolved tree of ``preset`` with a THIRD port (index 2, address port_k − 1) that duplicates port k: its
    segments, beams, lock, port entry, coefficients, the Liu detector / ADC / lock stack / delay line, a demux
    switch port and its loss element, and its segments added to the noise blocks (Phase 5 M6-09)."""
    from cranberry.config.yaml_io import Expr

    tree = copy.deepcopy(dict(config.load(ROOT / preset).unresolved))
    topo = tree["topology"]
    topo["segments"] += [_to_port2(x) for x in topo["segments"] if x["id"] in _PORT1_SEGS]
    topo["beams"] += [_to_port2(x) for x in topo["beams"] if x["id"] in ("sqz1", "pr1", "ref1", "lo1")]
    topo["observers"] += [_to_port2(x) for x in topo["observers"] if x["id"] == "det1"]
    topo["locks"] += [_to_port2(x) for x in topo["locks"] if x["id"] == "L1"]
    topo["ports"] += [_to_port2(x) for x in topo["ports"] if x["index"] == 1]
    topo["ports"][-1]["index"] = 2
    topo["coef"] += [_to_port2(x) for x in topo["coef"] if x["segment"] in _PORT1_SEGS]
    k2 = Expr("offline.params.port_k - 1")
    for b in tree["offline"]["blocks"]:
        if b.get("block") == "demux_design":
            b["params"]["port_segments"]["DMX2"] = k2
    pipe = tree["pipeline"]
    pipe += [_to_port2(x) for x in pipe if x.get("name") in ("liu1", "adc1", "lock1.stack", "delay1")]

    def walk(node):
        if isinstance(node, list):
            for x in node:
                walk(x)
        elif isinstance(node, Mapping):
            if node.get("block") == "demux_switch":
                node["params"]["ports"].append({**_to_port2(node["params"]["ports"][1]), "port": k2})
            elif "block" in node and isinstance(node.get("params"), dict):
                _add_segment_copies(node["params"])
            for key in ("seq", "stages", "actuators"):
                walk(node.get(key))
            walk(node.get("stack"))
    walk(pipe)
    return tree


def scale_estimate(preset: str = "F_liu20", durations_s: tuple[float, float] = (0.05, 0.3), n_ports_target: int = 32,
                   repeats: int = 3) -> dict[str, Any]:
    """Memory / time per simulated port from the MARGINAL cost of a real third port (Phase 5 M6-09; PLAN §1.2 M6: an
    estimate, not a 32-port run). The 2-port nominal configuration (reference port 0 + port k) and the same with a
    third port (:func:`three_port_tree`) each run compiled at two lengths, interleaved ``repeats`` times (minimum
    wall per length: shared machine), which separates build / compile from the per-sample loop time. Per port:
    Δ(µs/sample) and Δ(state) from 2 → 3 ports; 32 ports: 2-port value + 30 × Δ (the M6 review measured
    1.85 µs/sample per port → ≈ 59 µs/sample). Port 0 is cheaper than port k (n_sub = 2), so apportioning the
    2-port time by leaves / state share (the old estimate) underestimated a port."""
    from cranberry.runner.run import run

    trees = {2: {"extends": str(ROOT / nominal(preset))}, 3: three_port_tree(nominal(preset))}
    walls = {(n, d): [] for n in trees for d in durations_s}
    res: dict[int, Any] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for n, tr in trees.items():
            run(config.load({**tr, "run": {**tr.get("run", {}), "duration_s": 0.006}}))   # compile cache warm-up
        for _ in range(repeats):
            for n, tr in trees.items():
                for d in durations_s:
                    spec = config.load({**tr, "run": {**tr.get("run", {}), "duration_s": d}})
                    t0 = time.perf_counter()
                    res[n] = run(spec)
                    walls[(n, d)].append(time.perf_counter() - t0)
    f_sim = float(res[2].spec.run.f_sim_hz)
    ns = [int(round(d * f_sim)) for d in durations_s]
    us = {n: (min(walls[(n, durations_s[1])]) - min(walls[(n, durations_s[0])])) / (ns[1] - ns[0]) * 1e6 for n in trees}
    elems = {n: sum(r.final_state.layout.slot(s).size for s in r.final_state.layout.slot_names()) for n, r in res.items()}
    inc = us[3] - us[2]
    out = dict(f_sim_hz=f_sim, us_per_sample_2=us[2], us_per_sample_3=us[3], us_per_sample_per_port=inc,
               state_elems_2=elems[2], port_elems=elems[3] - elems[2], target_ports=n_ports_target,
               metric_ii_ports=sorted(res[3].metrics.metric_ii))
    out["target_us_per_sample"] = us[2] + (n_ports_target - 2) * inc
    out["target_state_bytes"] = 8 * (elems[2] + (n_ports_target - 2) * out["port_elems"])
    out["target_wall_s_per_sim_s"] = out["target_us_per_sample"] * 1e-6 * f_sim
    return out


def fmt_scale(est: dict[str, Any]) -> str:
    return (f"2 ports {est['us_per_sample_2']:.2f} µs/sample, 3 ports {est['us_per_sample_3']:.2f} µs/sample → "
            f"{est['us_per_sample_per_port']:.2f} µs/sample and {est['port_elems']} state elements per extra port. "
            f"{est['target_ports']} ports: ≈ {est['target_state_bytes'] / 1e6:.2f} MB state, "
            f"{est['target_us_per_sample']:.1f} µs/sample → {est['target_wall_s_per_sim_s']:.1f} s wall per simulated s "
            f"at f_sim = {est['f_sim_hz'] / 1e6:g} MHz")


def math_db(x: float) -> float:
    return 10.0 * math.log10(x)
