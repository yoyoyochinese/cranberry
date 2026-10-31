"""``ablate(spec, toggles)``: switch toggles off one at a time (PLAN §2.8, prompt §8).

Owner: not assigned explicitly in PLAN §4.2; proposed WP-F1 (runner) or WP-F4.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["ablate"]


def ablate(spec: RunSpec, toggles: Sequence[str]) -> dict[str, RunResult]:
    """Return {"baseline": result, "<block>.<toggle>": result with that toggle off, ...}.
    Toggle names are ``RunSpec.toggles`` keys (``"<block name>.<toggle>"``)."""
    raise NotImplementedError("runner.ablate is implemented by WP-F1 (proposed; unassigned in PLAN §4.2)")
