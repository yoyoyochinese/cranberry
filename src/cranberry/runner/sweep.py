"""``sweep(spec, axes)``: expand sweep markers / explicit axes into runs (PLAN §2.8, D16).

Owner: WP-F2 (D54).

Axes:

* every ``!sweep [lo, hi, n]`` (``numpy.linspace``), ``!sweep_log [lo, hi, n]`` (``numpy.geomspace``)
  and ``!choice [...]`` left in the spec is an axis; a range marker without ``n`` cannot be put on a
  grid and is an error unless ``axes`` gives its values (D34: no default point count);
* ``axes = {dotted path: values}`` adds axes or replaces a marker's values; paths are those printed
  by ``config.resolve.path_str`` (``topology.segments.arm2.where.dL``, ``pipeline.lock1.pid.params.kp``,
  ``run.f_sim_hz``; ``[i]`` for unnamed list elements).

The points are the cartesian product of the axes in tree order (explicit axes not on a marker
follow, in the given order). Each point is resolved again (expressions depending on an axis are
evaluated) and gets a **deterministic seed derived from the root seed and the point's coordinates**
(:func:`point_seed`): independent of the grid's size and order, reproducible, and different
between points. A seed axis (``axes={"run.seed": [...]}``) is used verbatim instead.
"""

from __future__ import annotations

import itertools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.config import ConfigError, sweep_markers
from cranberry.config.resolve import _convert_pi_frac, find_markers, normalize, parse_path, path_str, resolve, set_path
from cranberry.config.spec import spec_from_tree, spec_to_tree
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


def expand(spec: RunSpec, axes: Mapping[str, Sequence[Any]] | None = None) -> list[SweepPoint]:
    """Expand the axes of ``spec`` (+ explicit ``axes``) into single-run specs (module docstring)."""
    tree = spec_to_tree(spec)
    grid: dict[str, list[Any]] = {}
    paths: dict[str, tuple[Any, ...]] = {}
    explicit = {k: list(normalize(list(v), k)) for k, v in (axes or {}).items()}
    for k, v in explicit.items():
        if not v:
            raise ConfigError(f"axis {k!r} has no values")
    for p, marker in find_markers(tree):
        name = path_str(tree, p)
        paths[name] = p
        grid[name] = explicit[name] if name in explicit else axis_values(marker, name)
    for name, values in explicit.items():
        if name not in grid:
            paths[name] = parse_path(tree, name)
            grid[name] = values
    names = list(grid)
    points: list[SweepPoint] = []
    for combo in itertools.product(*(grid[n] for n in names)):
        values = dict(zip(names, combo))
        point_tree = spec_to_tree(spec)
        for name, value in values.items():
            set_path(point_tree, paths[name], value)
        if _SEED_PATH not in values:
            point_tree["run"]["seed"] = point_seed(spec.run.seed, values)
        point = spec_from_tree(resolve(point_tree))
        left = sweep_markers(point)
        if left:
            raise ConfigError(f"sweep point {values}: unresolved {sorted(left)}")
        points.append(SweepPoint(values, point))
    return points


def sweep(spec: RunSpec, axes: Mapping[str, Sequence[Any]] | None = None) -> list[RunResult]:
    """Run the grid of ``axes`` ({dotted param path: values}); None → expand the spec's own
    ``!sweep`` / ``!sweep_log`` / ``!choice`` markers. Each point is run by ``runner.run.run``
    (WP-F1); results are in point order (see :func:`expand` for the coordinates)."""
    from cranberry.runner.run import run

    return [run(p.spec) for p in expand(spec, axes)]

