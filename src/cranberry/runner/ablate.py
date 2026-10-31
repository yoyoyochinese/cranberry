"""``ablate(spec, toggles)``: switch toggles off one at a time (PLAN §2.8, prompt §8). Owner: WP-F4 (D54).

The baseline is the spec as given; each ablated run differs by one ``RunSpec.toggles`` override
``"<block name>.<toggle>": False`` (pipeline or ``offline.blocks`` block, ``runner.run``). The seed and
therefore every other noise realisation are unchanged (G11: streams keyed by block name), so the
difference isolates the toggle. :func:`ablation_impact` tabulates the change of chosen scalar results
(``RunResult.summary()``: metric moments, R′∓ [dB], η, recorded-series statistics).
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["ablate", "ablation_impact", "list_toggles"]


def list_toggles(spec: RunSpec | Any, *, only_on: bool = True) -> list[str]:
    """``"<block>.<toggle>"`` keys of every pipeline and ``offline.blocks`` block (after the spec's
    overrides); ``only_on``: just those currently True (the candidates for ablation)."""
    from cranberry.config.spec import apply_toggles, iter_block_specs
    from cranberry.core.registry import ensure_loaded, make_block
    from cranberry.runner.run import _as_runspec, inject_topology, offline_block_specs

    rs = _as_runspec(spec)
    ensure_loaded()
    off = offline_block_specs(rs)
    off_names = {b.name for b in off}
    pipe_t = {k: v for k, v in rs.toggles.items() if k.rpartition(".")[0] not in off_names}
    off_t = {k: v for k, v in rs.toggles.items() if k not in pipe_t}
    pipe = apply_toggles(inject_topology(rs.pipeline, rs.topology), pipe_t)
    specs = [b for _, b in iter_block_specs(pipe)] + list(off)
    out = []
    for b in specs:
        for t, on in make_block(b).toggles.items():
            key = f"{b.name}.{t}"
            if b.name in off_names and key in off_t:
                on = off_t[key]
            if on or not only_on:
                out.append(key)
    return out


def ablate(spec: RunSpec, toggles: Sequence[str], *, mode: str = "auto") -> dict[str, RunResult]:
    """Return {"baseline": result, "<block>.<toggle>": result with that toggle off, ...}.
    Toggle names are ``RunSpec.toggles`` keys (``"<block name>.<toggle>"``; :func:`list_toggles`)."""
    from cranberry.runner.run import _as_runspec, run

    rs = _as_runspec(spec)
    keys = list(toggles)
    bad = [t for t in keys if not isinstance(t, str) or "." not in t or t == "baseline"]
    if bad or len(set(keys)) != len(keys):
        raise ValueError(f"ablate: toggles must be distinct '<block>.<toggle>' names, got {keys!r}")
    out = {"baseline": run(rs, mode=mode)}
    for t in keys:
        out[t] = run(dataclasses.replace(rs, toggles={**dict(rs.toggles), t: False}), mode=mode)
    return out


def ablation_impact(
    results: Mapping[str, RunResult], quantities: Iterable[str] | None = None
) -> dict[str, dict[str, dict[str, float]]]:
    """{toggle: {quantity: {"baseline", "ablated", "delta"}}} from :func:`ablate` results.

    ``quantities``: ``RunResult.summary()`` keys (e.g. ``"i.port0.r_minus_db"``, ``"ii.port1.std_rad"``);
    None → every key of the baseline summary. Keys missing in a run are skipped.
    """
    if "baseline" not in results:
        raise KeyError("ablation_impact needs the 'baseline' result of ablate()")
    base = results["baseline"].summary()
    names = list(base) if quantities is None else list(quantities)
    table: dict[str, dict[str, dict[str, float]]] = {}
    for key, res in results.items():
        if key == "baseline":
            continue
        s = res.summary()
        table[key] = {q: {"baseline": base[q], "ablated": s[q], "delta": s[q] - base[q]}
                      for q in names if q in base and q in s}
    return table
