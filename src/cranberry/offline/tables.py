"""Offline lookup tables and their cache (G16-g, PLAN §2.2, §4.4). Owner: WP-F3.

A :class:`Table` is a named read-only array on a rectilinear grid with unit and provenance. Tables are
cached as ``runs/_cache/<key>.npz`` + ``<key>.json`` metadata; the key is the hash of the config
subtree that determines the table plus the code version (:func:`spec_hash`). The time loop only reads
tables (interpolation), it never builds them.

Cache contract (G16-g)
    * ``key = spec_hash(subtree, code_version)``. The subtree must contain everything the builder
      reads (including the table name if one subtree yields several tables). A change of any value
      in the subtree, or of ``code_version``, gives a new key, so a stale table is never reused
      (invalidation on hash change). Old entries stay on disk (sweeps reuse them); delete with
      :meth:`TableCache.invalidate` / :meth:`TableCache.clear`.
    * ``code_version`` convention: :func:`code_version` = ``"<cranberry version>|<builder>/<rev>"``;
      a builder bumps its ``rev`` whenever its numerics change.
    * The npz holds ``values`` and one ``grid__<axis>`` array per axis (``allow_pickle=False``); the
      JSON holds name, unit, axis order, value shape, provenance, ``spec_hash`` and a format version.
      A corrupt or mismatching entry (JSON ``spec_hash`` ≠ key, unreadable file, shape mismatch) is
      rebuilt and overwritten. Writes are atomic (temporary file + ``os.replace``).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["Table", "TableCache", "spec_hash", "code_version", "canonical_json", "CACHE_DIR", "CACHE_FORMAT"]

#: Cache directory relative to the repository root (gitignored with runs/).
CACHE_DIR = Path("runs") / "_cache"
#: Version of the on-disk layout (npz + JSON); bump on incompatible changes.
CACHE_FORMAT = 1


@dataclass(frozen=True)
class Table:
    """Read-only table.

    Attributes:
        name: table name (e.g. ``"vpm_bhd"``).
        grid: ordered {axis name with unit suffix, e.g. ``"f_hz"``: 1-D ascending grid}.
        values: array of shape ``tuple(len(g) for g in grid.values()) + value_shape``.
        unit: unit of ``values``.
        provenance: {"spec_hash": ..., "code_version": ..., "sources": [...], ...}.

    Arrays are copied and made read-only on construction. Axes must be strictly ascending; an axis
    of length 1 is allowed (a single-point table, e.g. the CW limit θ_off = 0) and then requires the
    coordinate to equal that point.
    """

    name: str
    grid: Mapping[str, np.ndarray]
    values: np.ndarray
    unit: str
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("table name must be non-empty")
        grid: dict[str, np.ndarray] = {}
        for axis, g in self.grid.items():
            a = np.array(g, dtype=float)
            if a.ndim != 1 or a.size < 1:
                raise ValueError(f"table {self.name!r}: axis {axis!r} must be a non-empty 1-D array")
            if a.size > 1 and not np.all(np.diff(a) > 0):
                raise ValueError(f"table {self.name!r}: axis {axis!r} must be strictly ascending")
            if not np.all(np.isfinite(a)):
                raise ValueError(f"table {self.name!r}: axis {axis!r} has non-finite values")
            a.flags.writeable = False
            grid[str(axis)] = a
        if not grid:
            raise ValueError(f"table {self.name!r}: needs at least one grid axis")
        values = np.array(self.values)
        lead = tuple(a.size for a in grid.values())
        if values.shape[: len(lead)] != lead:
            raise ValueError(f"table {self.name!r}: values shape {values.shape} does not start with grid shape {lead}")
        values.flags.writeable = False
        object.__setattr__(self, "grid", grid)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "provenance", dict(self.provenance))

    @property
    def axes(self) -> tuple[str, ...]:
        """Axis names in order."""
        return tuple(self.grid)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """Trailing shape of one table entry."""
        return self.values.shape[len(self.grid) :]

    def interp(self, **coords: Any) -> np.ndarray:
        """Multilinear interpolation at ``coords`` (one keyword per grid axis).

        Coordinates broadcast against each other; the result has shape
        ``broadcast(coords).shape + value_shape``. Points outside the grid raise ``ValueError``
        (no extrapolation; CHK ``grid_edge_monitor`` is the place for edge warnings).
        """
        missing = set(self.grid) - set(coords)
        extra = set(coords) - set(self.grid)
        if missing or extra:
            raise ValueError(f"table {self.name!r}: need coordinates {self.axes}, got {tuple(coords)}")
        pts = np.broadcast_arrays(*[np.asarray(coords[a], dtype=float) for a in self.axes])
        shape = pts[0].shape
        flat = [p.ravel() for p in pts]

        # per axis: lower index and weight of the upper neighbour
        idx, wts = [], []
        for axis, x in zip(self.axes, flat, strict=True):
            g = self.grid[axis]
            if g.size == 1:
                if not np.allclose(x, g[0], rtol=1e-12, atol=0.0):
                    raise ValueError(f"table {self.name!r}: axis {axis!r} has the single point {g[0]}, got {x}")
                idx.append(np.zeros(x.size, dtype=np.intp))
                wts.append(np.zeros(x.size))
                continue
            if np.any(x < g[0]) or np.any(x > g[-1]) or np.any(~np.isfinite(x)):
                raise ValueError(
                    f"table {self.name!r}: axis {axis!r} coordinate outside [{g[0]}, {g[-1]}] (no extrapolation)"
                )
            i = np.clip(np.searchsorted(g, x, side="right") - 1, 0, g.size - 2)
            idx.append(i)
            wts.append((x - g[i]) / (g[i + 1] - g[i]))

        out = np.zeros((flat[0].size,) + self.value_shape)
        n_ax = len(self.axes)
        for corner in range(2**n_ax):
            w = np.ones(flat[0].size)
            sel = []
            for d in range(n_ax):
                up = (corner >> d) & 1
                if self.grid[self.axes[d]].size == 1 and up:
                    w = None
                    break
                w = w * (wts[d] if up else 1.0 - wts[d])
                sel.append(idx[d] + up)
            if w is None:
                continue
            v = self.values[tuple(sel)]
            out += w.reshape((-1,) + (1,) * len(self.value_shape)) * v
        return out.reshape(shape + self.value_shape)


# --------------------------------------------------------------------------------------------------
# hashing
# --------------------------------------------------------------------------------------------------


def _to_jsonable(obj: Any, *, path: str = "$") -> Any:
    """Canonical JSON-compatible form used for hashing and metadata (deterministic, type-strict)."""
    if obj is None or isinstance(obj, bool | str):
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, int | np.integer):
        return int(obj)
    if isinstance(obj, float | np.floating):
        f = float(obj)
        if math.isnan(f) or math.isinf(f):
            return {"__float__": repr(f)}
        return f
    if isinstance(obj, complex | np.complexfloating):
        return {"__complex__": [_to_jsonable(obj.real, path=path), _to_jsonable(obj.imag, path=path)]}
    if isinstance(obj, np.ndarray):
        return {"__ndarray__": _to_jsonable(obj.tolist(), path=path), "dtype": str(obj.dtype), "shape": list(obj.shape)}
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, Mapping):
        out = {}
        for k, v in obj.items():
            if not isinstance(k, str):
                raise TypeError(f"{path}: mapping keys must be str for hashing, got {k!r}")
            out[k] = _to_jsonable(v, path=f"{path}.{k}")
        return out
    if isinstance(obj, list | tuple):
        return [_to_jsonable(v, path=f"{path}[{i}]") for i, v in enumerate(obj)]
    if isinstance(obj, set | frozenset):
        items = [_to_jsonable(v, path=path) for v in obj]
        return sorted(items, key=lambda x: json.dumps(x, sort_keys=True))
    if hasattr(obj, "to_json") and callable(obj.to_json):
        return _to_jsonable(obj.to_json(), path=path)
    raise TypeError(
        f"{path}: cannot hash value of type {type(obj).__name__}; give the config subtree in its textual / "
        "numeric form (callables as their !fn / !expr source, not as Python objects)"
    )


def canonical_json(obj: Any) -> str:
    """Deterministic JSON text of ``obj`` (sorted keys, no whitespace, exact float repr)."""
    return json.dumps(_to_jsonable(obj), sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def spec_hash(subtree: Any, code_version: str) -> str:
    """Stable hash of a resolved config subtree + code version (SHA-256 hex, 64 chars).

    Independent of mapping order and process (no Python ``hash``). Floats hash by their exact
    ``repr``; numpy scalars/arrays hash like the equal Python values plus dtype/shape for arrays.
    """
    if not isinstance(code_version, str) or not code_version:
        raise ValueError("code_version must be a non-empty str")
    text = canonical_json({"code_version": code_version, "subtree": subtree})
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def code_version(builder: str, rev: int | str) -> str:
    """Code-version string for :func:`spec_hash`: ``"<cranberry version>|<builder>/<rev>"``."""
    from cranberry import __version__

    if not builder:
        raise ValueError("builder name must be non-empty")
    return f"{__version__}|{builder}/{rev}"


# --------------------------------------------------------------------------------------------------
# cache
# --------------------------------------------------------------------------------------------------


class TableCache:
    """npz + JSON cache keyed by :func:`spec_hash`."""

    def __init__(self, root: Path = CACHE_DIR):
        self.root = Path(root)

    def paths(self, key: str) -> tuple[Path, Path]:
        """(npz path, JSON path) of ``key``."""
        _check_key(key)
        return self.root / f"{key}.npz", self.root / f"{key}.json"

    def has(self, key: str) -> bool:
        """True if a readable, matching entry for ``key`` exists."""
        return self.load(key) is not None

    def load(self, key: str) -> Table | None:
        """Cached table for ``key`` or None (missing, unreadable or mismatching entries → None)."""
        npz_path, json_path = self.paths(key)
        if not (npz_path.is_file() and json_path.is_file()):
            return None
        try:
            meta = json.loads(json_path.read_text(encoding="utf-8"))
            if meta.get("format") != CACHE_FORMAT or meta.get("spec_hash") != key:
                return None
            with np.load(npz_path, allow_pickle=False) as data:
                grid = {a: data[f"grid__{i}"] for i, a in enumerate(meta["axes"])}
                values = data["values"]
            if list(values.shape[len(grid) :]) != list(meta["value_shape"]):
                return None
            return Table(name=meta["name"], grid=grid, values=values, unit=meta["unit"], provenance=meta["provenance"])
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
            return None

    def store(self, key: str, table: Table) -> Table:
        """Write ``table`` under ``key`` (atomic); returns the table with ``spec_hash`` in provenance."""
        if not isinstance(table, Table):
            raise TypeError(f"table builder must return a Table, got {type(table).__name__}")
        prov = dict(table.provenance)
        if prov.setdefault("spec_hash", key) != key:
            raise ValueError(f"table {table.name!r}: provenance spec_hash {prov['spec_hash']!r} != cache key {key!r}")
        table = Table(name=table.name, grid=table.grid, values=table.values, unit=table.unit, provenance=prov)
        meta = {
            "format": CACHE_FORMAT,
            "spec_hash": key,
            "name": table.name,
            "unit": table.unit,
            "axes": list(table.axes),
            "value_shape": list(table.value_shape),
            "provenance": _to_jsonable(prov, path="provenance"),
        }
        npz_path, json_path = self.paths(key)
        self.root.mkdir(parents=True, exist_ok=True)
        arrays = {f"grid__{i}": table.grid[a] for i, a in enumerate(table.axes)}
        arrays["values"] = table.values
        _atomic_write(npz_path, lambda f: np.savez(f, **arrays))
        _atomic_write(json_path, lambda f: f.write(json.dumps(meta, indent=1, sort_keys=True).encode("utf-8")))
        return table

    def get_or_build(self, key: str, build: Callable[[], Table]) -> Table:
        """Return the cached table for ``key`` or build, store and return it (G16-g)."""
        cached = self.load(key)
        if cached is not None:
            return cached
        return self.store(key, build())

    def invalidate(self, key: str) -> bool:
        """Delete the entry for ``key``; True if something was removed."""
        removed = False
        for p in self.paths(key):
            if p.exists():
                p.unlink()
                removed = True
        return removed

    def clear(self) -> int:
        """Delete all cache entries under ``root``; returns the number of keys removed."""
        if not self.root.is_dir():
            return 0
        keys = {p.stem for p in self.root.iterdir() if p.suffix in (".npz", ".json") and _is_key(p.stem)}
        for k in keys:
            self.invalidate(k)
        return len(keys)


def _is_key(key: str) -> bool:
    return len(key) == 64 and all(c in "0123456789abcdef" for c in key)


def _check_key(key: str) -> None:
    if not isinstance(key, str) or not _is_key(key):
        raise ValueError(f"cache key must be a spec_hash (64 lowercase hex chars), got {key!r}")


def _atomic_write(path: Path, write: Callable[[Any], Any]) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as f:
            write(f)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
