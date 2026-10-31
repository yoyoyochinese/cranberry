"""Shared helpers of the future-preset tests (WP-I5): spec variants of the F_* presets, the manifest, CHK collection,
stray-light / GAWBS readouts and the 32-port scale estimate. Not a test module (loaded by path: ``tests`` is not a
package under ``--import-mode=importlib``).

Spec variants go through ``config.load`` with an inline tree that ``extends`` a preset, so every ``!expr`` (port
lengths, loop gains, τ_L − 1/f_sim, coefficients) is re-evaluated (D64).
"""

from __future__ import annotations

import dataclasses
import math
import time
import warnings
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
    """CHK event codes of a run (G16-j): construction-time events of the offline blocks plus the ``check()`` of the
    blocks that check on demand (demux_design, sbs_check). The runner does not collect them (proposal, WP-I5)."""
    codes: list[str] = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for b in result.offline.blocks:
            if hasattr(b, "check"):
                b.check()
            log = getattr(b, "chk_log", None)
            if log is not None:
                codes += [e.code for e in log.events]
    return codes


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


# ------------------------------------------------------------------------------------------- 32-port scale estimate


def scale_estimate(preset: str = "F_liu20", durations_s: tuple[float, float] = (0.02, 0.1),
                   n_ports_target: int = 32) -> dict[str, Any]:
    """Memory / time per simulated port and the extrapolation to ``n_ports_target`` ports (PLAN §1.2 M6: an
    estimate, not a 32-port run). Two compiled runs of the 2-port configuration (reference port 0 + port k) of
    different lengths separate the build / compile time from the per-sample loop time; the port-k share of the state
    (slots of the port's segments, beams, loop blocks, with their delay rings) apportions the loop time to a port
    (cost ∝ state touched per sample). The shared part (laser, squeezer, offline stage) is paid once.
    """
    from cranberry.runner.run import run

    walls, res = [], None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run(load(nominal(preset), run={"duration_s": 0.006}))    # compile cache warm-up (past the 5 ms settle)
        for d in durations_s:
            t0 = time.perf_counter()
            res = run(load(nominal(preset), run={"duration_s": d}))
            walls.append(time.perf_counter() - t0)
    f_sim = float(res.spec.run.f_sim_hz)
    n = [int(round(d * f_sim)) for d in durations_s]
    us = (walls[1] - walls[0]) / (n[1] - n[0]) * 1e6
    lay = res.final_state.layout
    sizes = {name: lay.slot(name).size for name in lay.slot_names()}
    total = sum(sizes.values())
    port = sum(v for name, v in sizes.items() if _port_owned(name, "1"))
    port0 = sum(v for name, v in sizes.items() if _port_owned(name, "0"))
    from cranberry.config.spec import iter_block_specs

    leaves = [b.name for _, b in iter_block_specs(res.spec.pipeline)]
    port_leaves = [x for x in leaves if _port_owned(x, "1")]
    # time share of a port: the larger of its share of the leaf blocks and of the state (conservative; the path
    # blocks' per-port work, Σ beams × subsegments, is inside the shared leaves and scales the same way)
    frac = max(len(port_leaves) / len(leaves), port / total)
    out = dict(state_elems=total, state_bytes=8 * total, port_k_elems=port, port_0_elems=port0, f_sim_hz=f_sim,
               leaves=len(leaves), port_leaves=len(port_leaves),
               build_s=walls[0] - n[0] * us * 1e-6, us_per_sample=us, us_per_sample_per_port=us * frac,
               target_ports=n_ports_target)
    shared_elems = total - port - port0
    out["target_state_bytes"] = 8 * (shared_elems + n_ports_target * port)
    out["target_us_per_sample"] = us * (1 - 2 * frac) + n_ports_target * us * frac
    out["target_wall_s_per_sim_s"] = out["target_us_per_sample"] * 1e-6 * f_sim
    return out


def _port_owned(slot: str, idx: str) -> bool:
    """Slots of port index ``idx`` (its segments, beams, loop blocks and delay line)."""
    keys = {f"F{idx}", f"D{idx}", f"G{idx}", f"DMX{idx}", f"PRS{idx}", f"REF{idx}", f"LO{idx}", f"sqz{idx}",
            f"pr{idx}", f"ref{idx}", f"lo{idx}", f"L{idx}", f"liu{idx}", f"adc{idx}", f"pid{idx}f", f"pid{idx}s",
            f"sup{idx}", f"dac{idx}f", f"dl{idx}f", f"dl{idx}s", f"pzt{idx}", f"tec{idx}", f"dmx{idx}", f"delay{idx}",
            f"port{idx}"}
    return any(p in keys for p in slot.split("."))


def fmt_scale(est: dict[str, Any]) -> str:
    return (f"2-port run: state {est['state_elems']} elements ({est['state_bytes'] / 1e3:.0f} kB), port k "
            f"{est['port_k_elems']}, port 0 {est['port_0_elems']}; {est['port_leaves']} of {est['leaves']} leaf blocks per port; "
            f"loop {est['us_per_sample']:.2f} µs/sample "
            f"({est['us_per_sample_per_port']:.2f} µs/sample/port), build+compile {est['build_s']:.1f} s. "
            f"{est['target_ports']} ports: ≈ {est['target_state_bytes'] / 1e6:.2f} MB state, "
            f"{est['target_us_per_sample']:.1f} µs/sample → {est['target_wall_s_per_sim_s']:.1f} s wall per simulated s "
            f"at f_sim = {est['f_sim_hz'] / 1e6:g} MHz")


def math_db(x: float) -> float:
    return 10.0 * math.log10(x)
