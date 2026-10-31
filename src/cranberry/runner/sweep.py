"""``sweep(spec, axes)``: expand sweep markers / explicit axes into runs (PLAN §2.8, D16).

Owner: WP-F2 (D54; seed mode and records WP-F2b).

Axes:

* every ``!sweep [lo, hi, n]`` (``numpy.linspace``), ``!sweep_log [lo, hi, n]`` (``numpy.geomspace``)
  and ``!choice [...]`` left in the spec is an axis, including one inside an ``!expr``'s own
  ``where`` or a ``!fn``'s ``args``; a range marker without ``n`` cannot be put on a grid and is an
  error unless ``axes`` gives its values (D34: no default point count);
* ``axes = {dotted path: values}`` adds axes or replaces a marker's values; paths are those printed
  by ``config.resolve.path_str`` (``topology.segments.arm2.where.dL``, ``pipeline.lock1.pid.params.kp``,
  ``run.f_sim_hz``; ``[i]`` for unnamed list elements). A dB parameter is addressed by its name as
  written (``pipeline.x.params.il_db``, values in dB) or by its linear name (``...il_lin``, linear
  values; written back in dB so that expressions reading the dB name follow, equal to the axis value
  up to float rounding); each point is converted (``config.db``).

The points are the cartesian product of the axes in tree order (explicit axes not on a marker
follow, in the given order).

Re-evaluation (WP-F2c, D63-9): each point starts from the spec's *unresolved* tree
(``config.input_tree``: the input document with ``extends`` / overlays merged, expressions still
symbolic, ``run:`` from ``spec.run``), gets the axis values, and is then resolved, so every ``!expr`` /
``!fn`` depending on an axis (e.g. ``tau_s: !expr "1e-6 - 1/run.f_sim_hz"``, arm lengths from a
swept ΔL) is evaluated at that point, and dB → linear (D56-1) happens after. An axis may name a value
that an expression computes (the axis value replaces the expression). A spec without a usable tree
(built with ``RunSpec(...)`` / ``spec_from_tree``, or edited outside ``run:`` with
``dataclasses.replace``) is expanded from its resolved values: expressions already evaluated stay
frozen, with a :class:`~cranberry.config.FrozenExpressionsWarning` when there is an axis. A spec
loaded from a resolved document (``dump_yaml`` output, a record's ``config.resolved.yaml``) holds
values only, so its former expressions are frozen too (no warning: nothing in it is symbolic); sweep
the input document (``config.read_input``, ``config.input.yaml``) instead. Each point carries its own
unresolved tree, so it can be expanded / rescaled again.

Seeds (``run.seed_mode``, D34 option; required when there is an axis, no default):

* ``per_point``: a **deterministic seed derived from the root seed and the point's coordinates**
  (:func:`point_seed`): independent of the grid's size and order, reproducible, different between
  points;
* ``common``: every point uses the root seed (common random numbers, so that points differ only by
  the parameters).

A seed axis (``axes={"run.seed": [...]}``) is used verbatim instead (then ``seed_mode`` is not
needed). A spec without any axis gives one point with the ``per_point`` derivation unless
``seed_mode: common``.
"""

from __future__ import annotations

import copy
import itertools
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.config import ConfigError, FrozenExpressionsWarning, input_tree, resolve_input, sweep_markers
from cranberry.config.db import linear_name, linear_value_to_db
from cranberry.config.resolve import (
    _convert_pi_frac,
    find_markers,
    get_path,
    normalize,
    parse_path,
    path_str,
    set_path,
)
from cranberry.config.spec import SEED_MODES, spec_to_tree
from cranberry.config.yaml_io import Choice, Sweep, SweepLog, canonical_json
from cranberry.core.rng import child_seed_sequence

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["SweepPoint", "axis_values", "point_seed", "expand", "sweep"]

_SEED_PATH = "run.seed"


@dataclass(frozen=True)
class SweepPoint:
    """One point of a sweep: its coordinates {axis path: value} and the single-run spec."""

    values: Mapping[str, Any]
    spec: RunSpec


def axis_values(marker: Any, path: str = "") -> list[Any]:
    """Grid values of a sweep marker (Python scalars; ``!pi_frac`` options as floats)."""
    if isinstance(marker, Choice):
        return [_convert_pi_frac(v) for v in marker.options]
    if isinstance(marker, (Sweep, SweepLog)):
        if marker.n is None:
            tag = "!sweep" if isinstance(marker, Sweep) else "!sweep_log"
            raise ConfigError(
                f"{path}: {tag} [{marker.lo}, {marker.hi}] has no point count; write {tag} [lo, hi, n] "
                f"or pass axes={{'{path}': [...]}} (D34: no default grid)"
            )
        space = np.linspace if isinstance(marker, Sweep) else np.geomspace
        return [float(v) for v in space(marker.lo, marker.hi, marker.n)]
    raise TypeError(f"not a sweep marker: {marker!r}")


def point_seed(root_seed: int, values: Mapping[str, Any]) -> int:
    """Seed of the sweep point with coordinates ``values`` under ``root_seed``: a 63-bit integer from
    ``core.rng.child_seed_sequence(root, "sweep:" + canonical JSON of values)``."""
    key = "sweep:" + canonical_json({k: values[k] for k in sorted(values)})
    state = child_seed_sequence(int(root_seed), key).generate_state(1, np.uint64)[0]
    return int(state) >> 1


def _db_axis(tree: Any, name: str) -> tuple[tuple[Any, ...], str, str] | None:
    """(parent path, axis key, key to drop) for an explicit axis naming a dB parameter by the other
    name than the tree holds: the dB name while a resolved spec holds the linear one, or the linear
    name while an unresolved tree holds the dB one; None if ``name`` is not such a path."""
    parent, _, key = name.rpartition(".")
    if not parent:
        return None
    try:
        p = parse_path(tree, parent)
    except ConfigError:
        return None
    node = get_path(tree, p)
    if not isinstance(node, Mapping) or key in node:
        return None
    lin = linear_name(key)
    if lin is not None and lin in node:
        return p, key, lin
    db = [k for k in node if isinstance(k, str) and linear_name(k) == key]
    return (p, key, db[0]) if db else None


def expand(spec: RunSpec, axes: Mapping[str, Sequence[Any]] | None = None) -> list[SweepPoint]:
    """Expand the axes of ``spec`` (+ explicit ``axes``) into single-run specs (module docstring)."""
    source = input_tree(spec)
    tree = source if source is not None else spec_to_tree(spec)
    grid: dict[str, list[Any]] = {}
    paths: dict[str, tuple[Any, ...]] = {}
    db_axes: dict[str, tuple[tuple[Any, ...], str, str]] = {}
    explicit = {k: list(normalize(list(v), k)) for k, v in (axes or {}).items()}
    for k, v in explicit.items():
        if not v:
            raise ConfigError(f"axis {k!r} has no values")
    for p, marker in find_markers(tree):
        name = path_str(tree, p)
        paths[name] = p
        grid[name] = explicit[name] if name in explicit else axis_values(marker, name)
    for name, values in explicit.items():
        if name in grid:
            continue
        try:
            paths[name] = parse_path(tree, name)
        except ConfigError:
            hit = _db_axis(tree, name)
            if hit is None:
                raise
            db_axes[name] = hit
        grid[name] = values
    names = list(grid)
    seed_axis = _SEED_PATH in grid
    mode = spec.run.seed_mode
    if names and not seed_axis and mode is None:
        raise ConfigError(
            f"sweep over {names}: run.seed_mode is required ({' | '.join(SEED_MODES)}; D34: no default): "
            "per_point = a seed derived from the root seed and the point's coordinates, "
            "common = the root seed for every point (common random numbers)"
        )
    if names and source is None:
        why = ("carries no unresolved tree (built with RunSpec(...) / spec_from_tree)" if spec.unresolved is None
               else "was edited outside run: after loading (its unresolved tree no longer matches)")
        warnings.warn(
            f"sweep over {names}: the spec {why}, so expressions are frozen at their resolved values; "
            "load it with config.load / config.to_spec to re-evaluate them (D63-9)",
            FrozenExpressionsWarning, stacklevel=2,
        )
    points: list[SweepPoint] = []
    for combo in itertools.product(*(grid[n] for n in names)):
        values = dict(zip(names, combo))
        point_tree = copy.deepcopy(tree)
        for name, value in values.items():
            if name in db_axes:
                parent, key, other = db_axes[name]
                node = get_path(point_tree, parent)
                if linear_name(key) == other:  # dB axis, resolved tree holds the linear name
                    node.pop(other)
                    node[key] = value  # converted to linear by resolve() below
                else:  # linear axis, unresolved tree holds the dB name: expressions read the dB value
                    node[other] = linear_value_to_db(other, value, name)
            else:
                set_path(point_tree, paths[name], value)
        if not seed_axis:
            point_tree["run"]["seed"] = spec.run.seed if mode == "common" else point_seed(spec.run.seed, values)
        point = resolve_input(point_tree)
        left = sweep_markers(point)
        if left:
            raise ConfigError(f"sweep point {values}: unresolved {sorted(left)}")
        points.append(SweepPoint(values, point))
    return points


def sweep(
    spec: RunSpec,
    axes: Mapping[str, Sequence[Any]] | None = None,
    *,
    record_root: str | Path | None = None,
    config_input: Any = None,
    invocation: Mapping[str, Any] | None = None,
) -> list[RunResult]:
    """Run the grid of ``axes`` ({dotted param path: values}); None → expand the spec's own
    ``!sweep`` / ``!sweep_log`` / ``!choice`` markers. Each point is run by ``runner.run.run``
    (WP-F1); results are in point order (see :func:`expand` for the coordinates).

    ``record_root``: write one run record per point under it (``runner.record``); each record's
    ``invocation.json`` holds the point's coordinates and the seed mode, ``config.input.yaml`` the
    sweep's input document (``config_input``, e.g. ``config.read_input(path)``; None → ``spec.unresolved`` if
    present, else the sweep spec),
    plus the caller's ``invocation`` context (CLI argv, overrides).
    """
    from cranberry.runner.run import run

    pts = expand(spec, axes)
    if record_root is None:
        return [run(p.spec) for p in pts]
    # keep the expressions in config.input.yaml when the spec still carries its input tree (D64)
    inp = config_input if config_input is not None else (getattr(spec, "unresolved", None) or spec)
    base = {**dict(invocation or {}), "seed_mode": spec.run.seed_mode, "n_points": len(pts),
            "axes": {k: list(v) for k, v in (axes or {}).items()}}
    return [
        run(p.spec, record_root=record_root, config_input=inp, invocation={**base, "sweep_point": dict(p.values)})
        for p in pts
    ]
