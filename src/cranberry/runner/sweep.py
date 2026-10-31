"""``sweep(spec, axes)``: expand sweep markers / explicit axes into runs (PLAN §2.8, D16).

Owner: not assigned explicitly in PLAN §4.2; proposed WP-F1 (runner).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["sweep"]


def sweep(spec: RunSpec, axes: Mapping[str, Sequence[Any]] | None = None) -> list[RunResult]:
    """Run the grid of ``axes`` ({dotted param path: values}); None → expand the spec's own
    ``!sweep`` / ``!sweep_log`` / ``!choice`` markers. Proposed owner WP-F1."""
    raise NotImplementedError("runner.sweep is implemented by WP-F1 (proposed; PLAN §4.2 names only F1 for run)")
