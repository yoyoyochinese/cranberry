"""Dependency resolution of callable / cross-referencing parameters (G12, prompt §5.5). Owner: WP-F2.

``!expr`` / ``!fn`` values may reference other parameters by dotted path. These references form a
directed graph that is resolved in topological order; a cycle is an error naming the cycle.
Time-varying callables (``time_varying: true``, e.g. L_ic(x, P, t)) are tabulated offline on a grid
and interpolated in the loop.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

__all__ = ["resolve"]


def resolve(tree: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve all callables / references in a merged config tree (sweep markers untouched). WP-F2."""
    raise NotImplementedError("config.resolve.resolve is implemented by WP-F2")
