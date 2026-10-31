"""Segment × driver × beam coefficients (conv §9, D47, G16-a). Owner: WP-F6.

``coef[(segment, driver, beam, effect)]`` maps a driver value on a segment (ΔT [K], strain [1],
detuning [rad/s], tec [K], stretcher / eom [rad or V per preset]) to the phase [rad] or group delay [s]
that this beam acquires. Phase and group delay have separate coefficients (``effect``). Values may be a
constant, a callable of named coordinates, or a Table (interpolated), G16-a. Equal coefficients on a
shared segment make common-mode noise cancel in the computed sums; unequal ones do not (VT-23).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, NamedTuple, Union

from cranberry.offline.tables import Table

__all__ = ["DRIVERS", "EFFECTS", "CoefKey", "CoefValue", "evaluate_coef"]

#: Drivers written into ``drv.<segment>.<driver>`` (PLAN §2.3.1–§2.3.2).
DRIVERS: tuple[str, ...] = ("dT", "strain", "detuning", "tec", "stretcher", "eom")
#: What a coefficient produces.
EFFECTS: tuple[str, ...] = ("phase", "group_delay")


class CoefKey(NamedTuple):
    segment: str
    driver: str
    beam: str
    effect: str


#: constant | callable(**coords) | Table
CoefValue = Union[float, Callable[..., float], Table]


def evaluate_coef(value: CoefValue, **coords: Any) -> float:
    """Evaluate a coefficient at ``coords`` (ignored for constants). Implemented by WP-F6."""
    raise NotImplementedError("path.coef.evaluate_coef is implemented by WP-F6")
