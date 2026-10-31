"""``convergence(spec, what)``: numerical convergence checks (PLAN §6.4, prompt §8). Owner: WP-F4.

``what`` ∈ {"f_sim", "chunk", "seed", "n_sub"}: f_sim ×2/×4 (change ≤ 1/3 of tolerance), chunk length
(bit-identical), seed ensemble (10–30 seeds; mean and spread), subsegment count ×2 (≤ 2 % PSD change).

Checks (WP-F4). Quantities are ``RunResult.summary()`` keys (metric moments, R′∓ [dB], recorded-series
mean / std / rms), with their tolerances::

    measure:
      convergence: {quantities: {<summary key>: <tolerance>, ...}, seeds: [<int>, ...], keep_rate: [<block>, ...]}

* ``chunk``: the spec and the same spec with another ``chunk_samples`` (``chunk_samples`` argument;
  None → ⌊chunk/3⌋ + 1) give bit-identical recorded series and final state (G11; a difference is a bug).
* ``f_sim``: f_sim × each of ``f_sim_factors`` (2, 4) through ``runner.sweep.expand``, same root seed;
  passed when |q(k·f_sim) − q(f_sim)| ≤ tolerance/3 for every quantity. ``run.f_sim_hz`` is set in
  the spec's unresolved tree and the tree is resolved again (WP-F2c, D63-9), so expressions depending
  on it, e.g. ``tau_s = τ_L − 1/f_sim`` of D55 rule 1, are re-evaluated (a spec without an unresolved
  tree warns, ``config.FrozenExpressionsWarning``). Blocks named in ``keep_rate`` (argument or
  ``measure.convergence.keep_rate``) keep their physical update rate: their ``rate_div`` is multiplied
  by the factor (``core.rates.scale_rate_divs``, as ``with_f_sim``; the factors must then be
  integers). A ``rate_div`` written as an ``!expr`` of ``run.f_sim_hz`` follows by itself and must not
  be listed.
* ``seed``: one run per seed of ``seeds``; mean and standard deviation of every quantity (a conclusion
  must exceed the spread, prompt §8); passed when 10 ≤ len(seeds) ≤ 30.
* ``n_sub``: every topology segment's ``n_sub`` doubled (subsegment indices remapped: point k → 2k,
  subs [k] → [2k, 2k+1], path ends and observers to the matching end of the finer grid; in
  ``RunSpec.topology`` and in explicit ``topology`` params); passed when every quantity changes by
  ≤ 2 % (relative; absolute ≤ tolerance/3 when the reference is 0). Blocks whose own params depend on
  N_sub (per-subsegment noise of shape (N_sub,)) must be adapted by their preset; a build failure is
  reported as ``error``.
"""

from __future__ import annotations

import copy
import dataclasses
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["CONVERGENCE_CHECKS", "ConvergenceReport", "convergence", "double_n_sub"]

CONVERGENCE_CHECKS: tuple[str, ...] = ("f_sim", "chunk", "seed", "n_sub")


@dataclass
class ConvergenceReport:
    """Per check: {check: {"passed": bool, "values": ..., "criterion": str}}."""

    checks: dict[str, dict[str, Any]] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        """True when every check ran and passed."""
        return all(c.get("passed") is True for c in self.checks.values())


def _identical(a: RunResult, b: RunResult) -> tuple[bool, float]:
    same = set(a.series) == set(b.series)
    diff = 0.0
    for k in a.series:
        if k in b.series:
            x, y = np.asarray(a.series[k]), np.asarray(b.series[k])
            if x.shape != y.shape or x.tobytes() != y.tobytes():
                same = False
                if x.shape == y.shape and x.size:
                    diff = max(diff, float(np.max(np.abs(x - y))))
    sa, sb = a.final_state, b.final_state
    if sa is not None and sb is not None:
        # every slot except the recorder's (sized by the chunk) and the DrawBuffer fill state ``rng.*``
        # (refilled per piece, so its contents depend on the chunk pattern even when every drawn number is
        # identical — WP-I1/I2 finding)
        skip = ("rec.", "rng.")
        names = [n for n in sa.layout.slot_names() if not n.startswith(skip)]
        if names != [n for n in sb.layout.slot_names() if not n.startswith(skip)]:
            return False, diff
        for n in names:
            if np.asarray(sa[n]).tobytes() != np.asarray(sb[n]).tobytes():
                same = False
    return same, diff


def _remap_point(k: int) -> int:
    return 2 * int(k)


def _double_topology(topo: Mapping[str, Any]) -> dict[str, Any]:
    t = copy.deepcopy(dict(topo))
    n_of: dict[str, int] = {}
    for s in t.get("segments", []):
        n = int(s["n_sub"])
        n_of[s["id"]] = n
        s["n_sub"] = 2 * n
        if s.get("actuator_sub") is not None:
            s["actuator_sub"] = _remap_point(s["actuator_sub"])
        for drv, pl in list((s.get("drivers") or {}).items()):
            if isinstance(pl, Mapping) and "point" in pl:
                s["drivers"][drv] = {"point": _remap_point(pl["point"])}
            elif isinstance(pl, Mapping) and "subs" in pl:
                s["drivers"][drv] = {"subs": [j for k in pl["subs"] for j in (2 * int(k), 2 * int(k) + 1)]}
    end_dir: dict[str, int] = {}
    for b in t.get("beams", []):
        for st in b.get("path", []):
            d = int(st["direction"])
            if "sub_in" in st:
                st["sub_in"] = 2 * int(st["sub_in"]) + (0 if d == 1 else 1)
            if "sub_out" in st:
                st["sub_out"] = 2 * int(st["sub_out"]) + (1 if d == 1 else 0)
        if b.get("path"):
            end_dir.setdefault(b["observer"], int(b["path"][-1]["direction"]))
    for o in t.get("observers", []):
        if o.get("sub") is not None:
            o["sub"] = 2 * int(o["sub"]) + (1 if end_dir.get(o["id"], -1) == 1 else 0)
    return t


def _double_params(obj: Any) -> Any:
    if isinstance(obj, list):
        return [_double_params(x) for x in obj]
    if isinstance(obj, dict):
        out = {k: _double_params(v) for k, v in obj.items()}
        if "block" in out and isinstance(out.get("params"), dict) and isinstance(out["params"].get("topology"), Mapping):
            out["params"] = {**out["params"], "topology": _double_topology(out["params"]["topology"])}
        return out
    return obj


def double_n_sub(spec: RunSpec) -> RunSpec:
    """``spec`` with every segment's ``n_sub`` doubled and subsegment indices remapped (module docstring)."""
    from cranberry.config.spec import spec_from_tree, spec_to_tree

    tree = spec_to_tree(spec)
    if tree.get("topology"):
        tree["topology"] = _double_topology(tree["topology"])
    tree["pipeline"] = _double_params(tree["pipeline"])
    if tree.get("offline", {}).get("blocks"):
        tree["offline"]["blocks"] = _double_params(tree["offline"]["blocks"])
    return spec_from_tree(tree)


def _compare(base: Mapping[str, float], other: Mapping[str, float], quantities: Mapping[str, float]) -> dict[str, Any]:
    missing = [q for q in quantities if q not in base or q not in other]
    if missing:
        raise KeyError(f"convergence quantities {missing} not in the run summary (known: {sorted(base)})")
    return {q: {"reference": base[q], "value": other[q], "change": other[q] - base[q], "tolerance": tol}
            for q, tol in quantities.items()}


def convergence(
    spec: RunSpec,
    what: Sequence[str] = CONVERGENCE_CHECKS,
    *,
    quantities: Mapping[str, float] | None = None,
    seeds: Sequence[int] | None = None,
    f_sim_factors: Sequence[float] = (2, 4),
    chunk_samples: int | None = None,
    mode: str = "auto",
    keep_rate: Sequence[str] | None = None,
) -> ConvergenceReport:
    """Run the requested convergence checks and return the report (module docstring).

    ``quantities`` ({summary key: tolerance}), ``seeds`` and ``keep_rate`` (block names) default to
    ``spec.measure["convergence"]``.
    """
    from cranberry.core.rates import scale_rate_divs
    from cranberry.runner.run import _as_runspec, run
    from cranberry.runner.sweep import expand

    rs = _as_runspec(spec)
    unknown = [w for w in what if w not in CONVERGENCE_CHECKS]
    if unknown:
        raise ValueError(f"unknown convergence check(s) {unknown}; known {CONVERGENCE_CHECKS}")
    cfg = (rs.measure or {}).get("convergence") or {}
    quantities = dict(cfg.get("quantities") or {}) if quantities is None else dict(quantities)
    seeds = cfg.get("seeds") if seeds is None else seeds
    keep_rate = list(cfg.get("keep_rate") or ()) if keep_rate is None else list(keep_rate)
    if any(w in ("f_sim", "n_sub") for w in what) and not quantities:
        raise KeyError("convergence f_sim / n_sub need quantities {summary key: tolerance} "
                       "(argument or measure.convergence.quantities)")
    if "seed" in what and not seeds:
        raise KeyError("convergence 'seed' needs seeds (argument or measure.convergence.seeds)")
    report = ConvergenceReport()
    base = run(rs, mode=mode)
    ref = base.summary()

    if "chunk" in what:
        c0 = rs.run.chunk_samples
        c1 = int(chunk_samples) if chunk_samples is not None else c0 // 3 + 1
        other = run(dataclasses.replace(rs, run=dataclasses.replace(rs.run, chunk_samples=c1)), mode=mode)
        same, diff = _identical(base, other)
        report.checks["chunk"] = {"passed": same, "values": {"chunk_samples": [c0, c1], "max_abs_diff": diff},
                                  "criterion": "bit-identical series and final state (G11)"}

    if "f_sim" in what:
        f0 = rs.run.f_sim_hz
        common = dataclasses.replace(rs, run=dataclasses.replace(rs.run, seed_mode="common"))
        pts = expand(common, {"run.f_sim_hz": [f0 * float(k) for k in f_sim_factors]})
        vals, ok = {}, True
        for k, p in zip(f_sim_factors, pts):
            ps = p.spec
            if keep_rate:
                if isinstance(k, bool) or int(k) != k:
                    raise ValueError(f"convergence f_sim with keep_rate needs integer factors, got {k!r}")
                ps = dataclasses.replace(ps, pipeline=scale_rate_divs(ps.pipeline, int(k), keep_rate))
            cmp = _compare(ref, run(ps, mode=mode).summary(), quantities)
            ok &= all(abs(c["change"]) <= c["tolerance"] / 3 for c in cmp.values())
            vals[f"x{k:g}"] = cmp
        report.checks["f_sim"] = {"passed": bool(ok), "values": vals,
                                  "criterion": "|q(k·f_sim) − q(f_sim)| <= tolerance/3 (PLAN §6.4)"}

    if "seed" in what:
        rows = [run(dataclasses.replace(rs, run=dataclasses.replace(rs.run, seed=int(s))), mode=mode).summary()
                for s in seeds]
        keys = list(quantities) if quantities else sorted(set.intersection(*(set(r) for r in rows)))
        vals = {q: {"values": [r[q] for r in rows], "mean": float(np.mean([r[q] for r in rows])),
                    "std": float(np.std([r[q] for r in rows], ddof=1)) if len(rows) > 1 else 0.0} for q in keys}
        report.checks["seed"] = {"passed": 10 <= len(seeds) <= 30, "values": {"seeds": list(seeds), **vals},
                                 "criterion": "10–30 seeds; conclusions must exceed the spread (prompt §8)"}

    if "n_sub" in what:
        if not rs.topology or not rs.topology.get("segments"):
            report.checks["n_sub"] = {"passed": None, "values": {}, "criterion": "skipped: no topology segments"}
        else:
            try:
                cmp = _compare(ref, run(double_n_sub(rs), mode=mode).summary(), quantities)
            except (ValueError, KeyError, TypeError) as e:
                report.checks["n_sub"] = {"passed": False, "values": {}, "error": f"{type(e).__name__}: {e}",
                                          "criterion": "≤ 2 % change with N_sub × 2 (T10 test 3)"}
            else:
                ok = all(abs(c["change"]) <= (0.02 * abs(c["reference"]) if c["reference"] else c["tolerance"] / 3)
                         for c in cmp.values())
                report.checks["n_sub"] = {"passed": bool(ok), "values": cmp,
                                          "criterion": "≤ 2 % change with N_sub × 2 (T10 test 3)"}
    return report
