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

Path model (WP-F3): a :class:`BeamPath` lists the segment ids a beam traverses and the observer ids it
reaches (Topology, WP-F6, provides them). An element belongs to a path through its ``segment`` or its
``observer``. GBS input vs BHD (D46, T12 V13) is separated by ``applies_to``: BHD-only elements (LO
overlap, QE, circuit noise) carry ``{"bhd"}`` and so never enter η(path, "gbs").
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from cranberry.core.units import eta_to_il_db, il_db_to_eta
from cranberry.offline import DoubleRegistrationError, MissingRegistrationError

__all__ = ["LOSS_KINDS", "LOSS_TARGETS", "OVERLAP_KINDS", "LossElement", "BeamPath", "LossRegistry"]

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
#: Mode-overlap kinds: counted at most once per path pair (beam path × metric target), G06.
OVERLAP_KINDS: tuple[str, ...] = ("mode_overlap_spatial", "mode_overlap_temporal", "lo_shape")


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

    @property
    def il_db(self) -> float:
        """Insertion loss IL = −10·log10(η) [dB] (inf for η = 0)."""
        return float("inf") if self.eta <= 0.0 else float(eta_to_il_db(self.eta))

    @property
    def location(self) -> tuple[str | None, str | None]:
        """(segment, observer) where the element sits."""
        return self.segment, self.observer

    def to_json(self) -> dict[str, Any]:
        """JSON-serialisable record of the element (with provenance tag)."""
        il = self.il_db
        return {
            "id": self.id,
            "eta": float(self.eta),
            "il_db": None if il == float("inf") else il,
            "kind": self.kind,
            "applies_to": sorted(self.applies_to),
            "source": self.source,
            "enabled": bool(self.enabled),
            "segment": self.segment,
            "observer": self.observer,
        }


@dataclass(frozen=True)
class BeamPath:
    """A beam path as the registry sees it: the segments it traverses and the observers it reaches.

    Attributes:
        id: beam-path id (Topology beam id, e.g. ``"sqz1"``, ``"lo1"``, ``"pr1"``).
        segments: segment ids along the path (order irrelevant for passive losses, T12 §1.2).
        observers: observer ids at its end (detectors ``"bhd1"``, GBS input ``"gbs1"``, ...).
    """

    id: str
    segments: tuple[str, ...]
    observers: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.id:
            raise ValueError("beam path id must be non-empty")
        object.__setattr__(self, "segments", tuple(self.segments))
        object.__setattr__(self, "observers", tuple(self.observers))
        for what, seq in (("segments", self.segments), ("observers", self.observers)):
            if len(set(seq)) != len(seq):
                raise ValueError(f"beam path {self.id!r}: duplicate entries in {what}")

    def contains(self, element: LossElement) -> bool:
        """True if ``element`` sits on this path (its segment or its observer belongs to it)."""
        return (element.segment is not None and element.segment in self.segments) or (
            element.observer is not None and element.observer in self.observers
        )

    def to_json(self) -> dict[str, Any]:
        return {"id": self.id, "segments": list(self.segments), "observers": list(self.observers)}


class LossRegistry:
    """Offline table of all loss elements; computes η per (beam path, metric target).

    An element counts for path p and metric m when it sits on p (:meth:`BeamPath.contains`) and
    m ∈ ``applies_to``; η(p, m) is the product of the enabled ones (disabled elements stay in the
    record with η treated as 1). Paths are declared with :meth:`define_path` (or the ``paths``
    constructor argument: {id: BeamPath | {"segments": [...], "observers": [...]}}).

    Build-time errors (VT-F3, :class:`~cranberry.offline.DoubleRegistrationError`):
        * a second element with the same ``id``;
        * a second element of the same overlap kind (:data:`OVERLAP_KINDS`) on the same path pair —
          "path pair" = (beam path, metric target), i.e. the beam and the reference mode it is
          overlapped with at that target (squeezed ↔ LO_k at ``bhd``, squeezed ↔ GBS mode at
          ``gbs``). Checked against every defined path, and path-independently for elements at the
          same segment / observer with intersecting ``applies_to``. All registrations count,
          enabled or not.
    Every element must sit somewhere (``segment`` or ``observer`` given), else ``ValueError``.
    """

    def __init__(
        self,
        elements: Iterable[LossElement] = (),
        paths: Mapping[str, BeamPath | Mapping[str, Any]] | None = None,
    ):
        self.elements: list[LossElement] = []
        self.paths: dict[str, BeamPath] = {}
        for pid, p in (paths or {}).items():
            if isinstance(p, BeamPath):
                if p.id != pid:
                    raise ValueError(f"paths[{pid!r}] has id {p.id!r}")
                self.define_path(pid, p.segments, p.observers)
            else:
                self.define_path(pid, p.get("segments", ()), p.get("observers", ()))
        for e in elements:
            self.add(e)

    # ---------------------------------------------------------------------------------------------
    # registration

    def add(self, element: LossElement) -> None:
        """Register one element (duplicate ids and double-counted overlaps are errors, VT-F3)."""
        if not isinstance(element, LossElement):
            raise TypeError(f"expected LossElement, got {type(element).__name__}")
        if element.segment is None and element.observer is None:
            raise ValueError(f"loss element {element.id!r}: give its segment or its observer")
        if any(e.id == element.id for e in self.elements):
            raise DoubleRegistrationError(f"loss element id {element.id!r} registered twice")
        if element.kind in OVERLAP_KINDS:
            for e in self.elements:
                if (
                    e.kind == element.kind
                    and e.applies_to & element.applies_to
                    and (
                        (e.segment is not None and e.segment == element.segment)
                        or (e.observer is not None and e.observer == element.observer)
                    )
                ):
                    raise DoubleRegistrationError(
                        f"{element.kind} counted twice at {element.location}: {e.id!r} and {element.id!r} (G06)"
                    )
        self.elements.append(element)
        try:
            for p in self.paths.values():
                self._check_overlaps(p)
        except DoubleRegistrationError:
            self.elements.pop()
            raise

    def add_all(self, elements: Iterable[LossElement]) -> None:
        """Register several elements (e.g. ``block.losses()``)."""
        for e in elements:
            self.add(e)

    def define_path(self, path: str, segments: Iterable[str], observers: Iterable[str] = ()) -> BeamPath:
        """Declare beam path ``path`` (redefinition is an error)."""
        if path in self.paths:
            raise DoubleRegistrationError(f"beam path {path!r} defined twice")
        bp = BeamPath(path, tuple(segments), tuple(observers))
        self._check_overlaps(bp)
        self.paths[path] = bp
        return bp

    def _check_overlaps(self, path: BeamPath) -> None:
        for metric in LOSS_TARGETS:
            for kind in OVERLAP_KINDS:
                hits = [e.id for e in self.elements if e.kind == kind and metric in e.applies_to and path.contains(e)]
                if len(hits) > 1:
                    raise DoubleRegistrationError(
                        f"{kind} counted {len(hits)}× on path pair ({path.id!r}, {metric!r}): {hits} (G06)"
                    )

    # ---------------------------------------------------------------------------------------------
    # queries

    def _path(self, path: str) -> BeamPath:
        try:
            return self.paths[path]
        except KeyError:
            raise KeyError(f"beam path {path!r} not defined (known: {sorted(self.paths)})") from None

    @staticmethod
    def _check_metric(metric: str) -> None:
        if metric not in LOSS_TARGETS:
            raise ValueError(f"metric target {metric!r} not in {LOSS_TARGETS}")

    def elements_on(self, path: str, metric: str, *, include_disabled: bool = False) -> list[LossElement]:
        """Elements that count for (``path``, ``metric``), in registration order."""
        self._check_metric(metric)
        bp = self._path(path)
        return [
            e for e in self.elements if metric in e.applies_to and bp.contains(e) and (include_disabled or e.enabled)
        ]

    def eta(self, path: str, metric: str) -> float:
        """Total transmission of beam ``path`` (Topology beam id) for target ``metric`` ∈ LOSS_TARGETS."""
        return float(np.prod([e.eta for e in self.elements_on(path, metric)]))

    def eta_array(self, paths: Iterable[str], metric: str) -> np.ndarray:
        """Vectorised :meth:`eta`."""
        return np.array([self.eta(p, metric) for p in paths], dtype=float)

    def loss(self, path: str, metric: str) -> float:
        """Total loss L = 1 − η(path, metric)."""
        return 1.0 - self.eta(path, metric)

    def breakdown(self, path: str, metric: str) -> list[dict[str, Any]]:
        """Per-element contributions (incl. disabled ones, flagged) for logging / budget tables."""
        return [e.to_json() for e in self.elements_on(path, metric, include_disabled=True)]

    def missing_registrations(self, path: str, required: Iterable[str]) -> list[str]:
        """Owners in ``required`` without a registered element on ``path``.

        An owner (block / component name, e.g. ``"arm1.915b"``) counts as registered if some element
        with ``id == owner`` or ``id`` starting with ``owner + "."`` sits on ``path`` (enabled or not).
        """
        bp = self._path(path)
        on_path = [e.id for e in self.elements if bp.contains(e)]
        return [o for o in required if not any(i == o or i.startswith(o + ".") for i in on_path)]

    def require_registered(self, path: str, required: Iterable[str]) -> None:
        """Raise :class:`~cranberry.offline.MissingRegistrationError` if a mandatory element is missing.

        Use for modulators / actuators / probe in-out couplers on the squeezed path (prompt §5.5,
        D10, T12 §8.2): they must have a ``LossElement`` on that path even if η = 1.
        """
        missing = self.missing_registrations(path, list(required))
        if missing:
            raise MissingRegistrationError(
                f"no loss element registered on path {path!r} for {missing} (squeezed-path elements must be "
                "registered, prompt §5.5 / D10 / G06)"
            )

    # ---------------------------------------------------------------------------------------------
    # record

    def to_json(self) -> dict[str, Any]:
        """JSON-serialisable record with provenance (``loss_registry.json``)."""
        totals: dict[str, dict[str, Any]] = {}
        for pid in self.paths:
            totals[pid] = {}
            for m in LOSS_TARGETS:
                els = self.elements_on(pid, m)
                if not els:
                    continue
                eta = self.eta(pid, m)
                totals[pid][m] = {
                    "eta": eta,
                    "loss": 1.0 - eta,
                    "il_db": None if eta <= 0.0 else float(eta_to_il_db(eta)),
                    "elements": [e.id for e in els],
                }
        return {
            "loss_kinds": list(LOSS_KINDS),
            "targets": list(LOSS_TARGETS),
            "elements": [e.to_json() for e in self.elements],
            "paths": {pid: p.to_json() for pid, p in self.paths.items()},
            "eta": totals,
        }

