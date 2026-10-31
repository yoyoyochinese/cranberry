"""Loss registry: the single place where transmissions η are combined (G06, PLAN §2.2). Owner: WP-F3.

Every block with a lossy element on a squeezed / LO / probe path returns :class:`LossElement` records
from ``Block.losses()``; blocks never multiply η into a state themselves (PLAN §7.3). The registry
then gives η per beam path and per metric target:

    squeezed → GBS_k ("gbs"), squeezed → BHD_k ("bhd"), LO_k → BHD_k ("lo"), probe ("probe")

Rules (G06): series losses multiply; each overlap kind (``mode_overlap_spatial``,
``mode_overlap_temporal``, ``lo_shape``) is counted once per path pair, a second registration is an
error (VT-F3); modulators / actuators / probe in-out couplers on the squeezed path must be registered
(prompt §5.5, D10). Pure-loss action on quadrature variance: V → η V + (1 − η)/2 (conv §2), applied
by ``offline.gaussian.pure_loss``. The registry is saved with provenance tags (``loss_registry.json``).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np

from cranberry.core.units import il_db_to_eta

__all__ = ["LOSS_KINDS", "LOSS_TARGETS", "LossElement", "LossRegistry"]

LOSS_KINDS: tuple[str, ...] = (
    "connector",
    "splice",
    "fiber",
    "component",
    "mode_overlap_spatial",
    "mode_overlap_temporal",
    "lo_shape",
    "electronic",
    "detector",
)
LOSS_TARGETS: tuple[str, ...] = ("gbs", "bhd", "lo", "probe")


@dataclass(frozen=True)
class LossElement:
    """One loss element.

    Attributes:
        id: unique element id (e.g. ``"arm1.915b.il"``).
        eta: power transmission in [0, 1] (use :meth:`from_il_db` for dB input).
        kind: one of :data:`LOSS_KINDS`.
        applies_to: subset of :data:`LOSS_TARGETS` the element counts for (GBS input vs BHD, D46).
        source: provenance tag ([Ref A], [manual p.N], [AE], UNVERIFIED, Dnn; PLAN §7.3).
        enabled: False keeps the element in the record but with η treated as 1 (toggle).
        segment: segment id the element sits on (path elements), or None.
        observer: observer id (detector / electronic elements not on a segment), or None.
    """

    id: str
    eta: float
    kind: str
    applies_to: frozenset[str]
    source: str
    enabled: bool
    segment: str | None = None
    observer: str | None = None

    def __post_init__(self) -> None:
        if not (0.0 <= float(self.eta) <= 1.0):
            raise ValueError(f"loss element {self.id!r}: eta must be in [0, 1], got {self.eta!r}")
        if self.kind not in LOSS_KINDS:
            raise ValueError(f"loss element {self.id!r}: kind {self.kind!r} not in {LOSS_KINDS}")
        object.__setattr__(self, "applies_to", frozenset(self.applies_to))
        bad = self.applies_to - set(LOSS_TARGETS)
        if bad or not self.applies_to:
            raise ValueError(f"loss element {self.id!r}: applies_to must be a non-empty subset of {LOSS_TARGETS}")

    @classmethod
    def from_il_db(cls, id: str, il_db: float, **kwargs: Any) -> LossElement:
        """Build from insertion loss IL [dB] (η = 10^(−IL/10))."""
        return cls(id=id, eta=float(il_db_to_eta(il_db)), **kwargs)


class LossRegistry:
    """Offline table of all loss elements; computes η per (beam path, metric target)."""

    def __init__(self, elements: Iterable[LossElement] = ()):
        self.elements: list[LossElement] = []
        for e in elements:
            self.add(e)

    def add(self, element: LossElement) -> None:
        """Register one element (duplicate ids and double-counted overlaps are errors, VT-F3)."""
        raise NotImplementedError("offline.loss_registry.LossRegistry.add is implemented by WP-F3")

    def eta(self, path: str, metric: str) -> float:
        """Total transmission of beam ``path`` (Topology beam id) for target ``metric`` ∈ LOSS_TARGETS."""
        raise NotImplementedError("offline.loss_registry.LossRegistry.eta is implemented by WP-F3")

    def eta_array(self, paths: Iterable[str], metric: str) -> np.ndarray:
        """Vectorised :meth:`eta`."""
        raise NotImplementedError("WP-F3")

    def to_json(self) -> dict[str, Any]:
        """JSON-serialisable record with provenance (``loss_registry.json``)."""
        raise NotImplementedError("WP-F3")
