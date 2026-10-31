"""Offline lookup tables and their cache (G16-g, PLAN §2.2, §4.4). Owner: WP-F3.

A :class:`Table` is a named read-only array on a rectilinear grid with unit and provenance. Tables are
cached as ``runs/_cache/<key>.npz`` + JSON metadata; the key is the hash of the config subtree that
determines the table plus the code version (:func:`spec_hash`). The time loop only reads tables
(interpolation), it never builds them.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["Table", "TableCache", "spec_hash", "CACHE_DIR"]

#: Cache directory relative to the repository root (gitignored with runs/).
CACHE_DIR = Path("runs") / "_cache"


@dataclass(frozen=True)
class Table:
    """Read-only table.

    Attributes:
        name: table name (e.g. ``"vpm_bhd"``).
        grid: ordered {axis name with unit suffix, e.g. ``"f_hz"``: 1-D ascending grid}.
        values: array of shape ``tuple(len(g) for g in grid.values()) + value_shape``.
        unit: unit of ``values``.
        provenance: {"spec_hash": ..., "code_version": ..., "sources": [...], ...}.
    """

    name: str
    grid: Mapping[str, np.ndarray]
    values: np.ndarray
    unit: str
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def interp(self, **coords: Any) -> np.ndarray:
        """Multilinear interpolation at ``coords`` (one keyword per grid axis). WP-F3."""
        raise NotImplementedError("offline.tables.Table.interp is implemented by WP-F3")


class TableCache:
    """npz + JSON cache keyed by :func:`spec_hash`."""

    def __init__(self, root: Path = CACHE_DIR):
        self.root = Path(root)

    def get_or_build(self, key: str, build: Callable[[], Table]) -> Table:
        """Return the cached table for ``key`` or build, store and return it. WP-F3."""
        raise NotImplementedError("offline.tables.TableCache.get_or_build is implemented by WP-F3")


def spec_hash(subtree: Any, code_version: str) -> str:
    """Stable hash of a resolved config subtree + code version (hex). WP-F3."""
    raise NotImplementedError("offline.tables.spec_hash is implemented by WP-F3")
