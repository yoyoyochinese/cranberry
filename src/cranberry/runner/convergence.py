"""``convergence(spec, what)``: numerical convergence checks (PLAN §6.4, prompt §8). Owner: WP-F4.

``what`` ∈ {"f_sim", "chunk", "seed", "n_sub"}: f_sim ×2/×4 (change ≤ 1/3 of tolerance), chunk length
(bit-identical), seed ensemble (10–30 seeds; mean and spread), subsegment count ×2 (≤ 2 % PSD change).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec

__all__ = ["CONVERGENCE_CHECKS", "ConvergenceReport", "convergence"]

CONVERGENCE_CHECKS: tuple[str, ...] = ("f_sim", "chunk", "seed", "n_sub")


@dataclass
class ConvergenceReport:
    """Per check: {check: {"passed": bool, "values": ..., "criterion": str}}."""

    checks: dict[str, dict[str, Any]] = field(default_factory=dict)


def convergence(spec: RunSpec, what: Sequence[str] = CONVERGENCE_CHECKS) -> ConvergenceReport:
    """Run the requested convergence checks and return the report. WP-F4."""
    raise NotImplementedError("runner.convergence is implemented by WP-F4")
