"""Segment × driver × beam coefficients (conv §9, D47, G16-a). Owner: WP-F6.

``coef[(segment, driver, beam, effect)]`` maps a driver value on a segment (ΔT [K], strain [1],
detuning [rad/s], tec [K], stretcher / eom [rad or V per preset]) to the phase [rad] or group delay [s]
that this beam acquires. Phase and group delay have separate coefficients (``effect``). Values may be a
constant, a callable of named coordinates, or a Table (interpolated), G16-a. Equal coefficients on a
shared segment make common-mode noise cancel in the computed sums; unequal ones do not (VT-23).

Meaning of one coefficient (WP-F6): the response of the **whole segment** to a unit driver applied
uniformly over it (e.g. ``k_T·L`` [rad/K] for ``dT``). A driver acting on part of the segment enters
with the weight of its placement (``Segment.drivers``: ``uniform`` / ``per_sub`` → 1/N_sub per
subsegment, ``subs`` → 1/N_sub on the listed subsegments, ``point`` → 1 at one subsegment), so the same
coefficient serves a uniform temperature and a TEC on part of the segment (G08). A point actuator
(stretcher, EOM) has weight 1: its coefficient is its full response [rad per driver unit].

Coefficients are evaluated once at build time (``evaluate_coef``) at the coordinates given in the
topology spec (``coords`` of the coefficient entry); state-dependent responses (e.g. the OPO angle map
on (Δ_s, Δ_p, x_ss)) are blocks that write ``psi.*`` themselves (``opo_angle_map``, WP-B11).
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import Any, NamedTuple, Union

import numpy as np

from cranberry.offline.tables import Table

__all__ = ["DRIVERS", "EFFECTS", "CoefKey", "CoefValue", "evaluate_coef"]

#: Drivers written into ``drv.<segment>.<driver>`` (PLAN §2.3.1–§2.3.2).
DRIVERS: tuple[str, ...] = ("dT", "strain", "detuning", "tec", "stretcher", "eom", "dither", "mirror")
#: What a coefficient produces.
EFFECTS: tuple[str, ...] = ("phase", "group_delay")


class CoefKey(NamedTuple):
    segment: str
    driver: str
    beam: str
    effect: str


#: constant | callable(**coords) | Table
CoefValue = Union[float, Callable[..., float], Table]


def _call_with(fn: Callable[..., Any], coords: dict[str, Any]) -> Any:
    """Call ``fn`` with the subset of ``coords`` it accepts (all of them if it takes ``**kwargs``)."""
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):  # builtins without a signature
        return fn(**coords)
    params = sig.parameters.values()
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params):
        return fn(**coords)
    names = [p.name for p in params if p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)]
    required = [p.name for p in params if p.default is inspect.Parameter.empty
                and p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)]
    missing = [n for n in required if n not in coords]
    if missing:
        raise ValueError(f"coefficient callable needs coordinate(s) {missing}; given {sorted(coords)}")
    return fn(**{n: coords[n] for n in names if n in coords})


def evaluate_coef(value: CoefValue, **coords: Any) -> float:
    """Evaluate a coefficient at ``coords`` (ignored for constants).

    Accepted values (G16-a):

    * a real number (bool rejected);
    * a :class:`~cranberry.offline.tables.Table` with a scalar entry: multilinear interpolation at the
      coordinates of its axes (missing axes raise; extra coordinates are ignored);
    * a config ``!expr`` with ``args`` / ``!fn`` marker (``config.yaml_io.Expr`` / ``Fn``): turned into
      a function with ``config.resolve.as_callable`` and called with its arguments from ``coords``;
    * any Python callable: called with the coordinates it names (all of them if it takes ``**kw``).

    Returns a finite float; anything else raises ``ValueError`` / ``TypeError``.
    """
    from cranberry.config.yaml_io import Expr, Fn

    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"coefficient must be a number, callable or Table, got bool {value!r}")
    if isinstance(value, (int, float, np.integer, np.floating)):
        out = float(value)
    elif isinstance(value, Table):
        if value.value_shape != ():
            raise ValueError(f"coefficient table {value.name!r} must have scalar entries, got {value.value_shape}")
        missing = [a for a in value.axes if a not in coords]
        if missing:
            raise ValueError(f"coefficient table {value.name!r} needs coordinate(s) {missing}")
        out = float(value.interp(**{a: coords[a] for a in value.axes}))
    elif isinstance(value, Expr):
        from cranberry.config.resolve import as_callable

        args = tuple(value.args)
        missing = [a for a in args if a not in coords]
        if missing:
            raise ValueError(f"coefficient !expr {value.source!r} needs coordinate(s) {missing}")
        out = float(as_callable(value, args)(*[coords[a] for a in args]))
    elif isinstance(value, Fn):
        from cranberry.core.registry import get_fn

        out = float(_call_with(_bind_kwargs(get_fn(value.name), dict(value.args)), dict(coords)))
    elif callable(value):
        out = float(_call_with(value, dict(coords)))
    else:
        raise TypeError(f"coefficient must be a number, callable or Table, got {type(value).__name__}: {value!r}")
    if not np.isfinite(out):
        raise ValueError(f"coefficient evaluated to a non-finite value {out!r}")
    return out


def _bind_kwargs(fn: Callable[..., Any], fixed: dict[str, Any]) -> Callable[..., Any]:
    """``fn`` with keyword arguments ``fixed`` bound, keeping the signature of the remaining ones."""
    sig = inspect.signature(fn)
    rest = [p for p in sig.parameters.values() if p.name not in fixed]

    def bound(**kw: Any) -> Any:
        return fn(**fixed, **kw)

    bound.__signature__ = sig.replace(parameters=rest)  # type: ignore[attr-defined]
    return bound
