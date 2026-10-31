"""Consistency warnings that never stop a run (CHK, G16-j, D34).

CHK blocks (``blocks/chk``) and other code report physically questionable configurations (SBS
threshold, repetition-rate tail limit, demux feasibility, GNLSE grid edge, ...) as :class:`ChkEvent`
records. They are emitted as Python warnings of category :class:`ChkWarning` and collected by the
runner into the run record; the computation continues (bad configurations are allowed on purpose,
D34).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any

__all__ = ["ChkWarning", "ChkEvent", "ChkLog"]


class ChkWarning(UserWarning):
    """Warning category for CHK events."""


@dataclass(frozen=True)
class ChkEvent:
    """One consistency warning.

    Attributes:
        code: stable identifier, e.g. ``"sbs_threshold"``, ``"rep_rate_tail"``.
        source: emitting block / module name.
        message: human-readable explanation with numbers.
        context: machine-readable values (JSON-serialisable).
    """

    code: str
    source: str
    message: str
    context: dict[str, Any] = field(default_factory=dict)


class ChkLog:
    """Collector of CHK events for one run (serialised into the run record by WP-F2)."""

    def __init__(self) -> None:
        self.events: list[ChkEvent] = []

    def warn(self, code: str, source: str, message: str, **context: Any) -> ChkEvent:
        """Record an event and emit it as a :class:`ChkWarning`; never raises."""
        ev = ChkEvent(code, source, message, dict(context))
        self.events.append(ev)
        warnings.warn(f"[{code}] {source}: {message}", ChkWarning, stacklevel=2)
        return ev

    def to_json(self) -> list[dict[str, Any]]:
        return [{"code": e.code, "source": e.source, "message": e.message, "context": e.context} for e in self.events]
