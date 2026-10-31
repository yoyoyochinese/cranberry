"""dB → linear conversion when a config is loaded (G16-i, ``core.units``). Owner: WP-F2b.

Code below the config layer works in linear units only. Configs may give logarithmic values under
names carrying a dB suffix (``core.units.is_db_name``); :func:`linearize` converts each such key and
renames it (longest suffix wins):

=================  ======================  ================================  ==============================
written as         becomes                 value                             e.g.
=================  ======================  ================================  ==============================
``<x>_db``         ``<x>_lin``             10^(v/10)   (power ratio)         ``il_db: 3`` → ``il_lin: 1.995``
``<x>_db_per_hz``  ``<x>_per_hz``          10^(v/10)   [1/Hz]                ``rin_db_per_hz: -140`` → 1e-14
``<x>_dbc_per_hz`` ``<x>_rad2_per_hz``     2·10^(v/10) one-sided S_φ (conv   ``pn_dbc_per_hz: -100`` → 2e-10
                                           §7: ℒ = S_φ/2)
=================  ======================  ================================  ==============================

Notes on the convention:

* ``_db`` is always a *power* ratio, 10·log10. An insertion loss ``il_db`` (≥ 0) therefore becomes
  the loss factor ``il_lin = 10^(IL/10) ≥ 1`` (transmission η = 1/il_lin); a noise level relative to
  shot noise ``sq_db`` becomes V/V_SN. Amplitude-ratio quantities must be written linear.
* Converted are numbers and (nested) lists of numbers; a callable ``!expr`` (with ``args``) is
  wrapped, ``10 ** ((<source>) / 10)``. Anything else under a dB name (string, bool, mapping, ``!fn``,
  null) is a :class:`ConfigError`, as is a dB key whose linear name is also present.
* Where: every mapping key in ``topology``, ``offline``, ``pipeline`` (block ``params``, incl. ``!fn``
  ``args``) and ``measure``; not in ``run``/``toggles``, and not the bindings of an ``!expr`` (they are
  names inside the expression).
* When: at the end of ``config.resolve.resolve`` and only for a fully resolved tree (no sweep axis
  or deferred ``!expr`` left). A sweep spec keeps the dB names as written, so its axes are dB axes
  (``...params.il_db``) and expressions see the dB values; ``runner.sweep`` converts each point.
  Expressions are always evaluated before conversion, on the values as written.

Lossless round trip: a resolved spec stores the linear names and values, so ``dump_yaml`` writes
``il_lin: 1.9952623149688795`` (exact float ``repr``), which reloads unchanged (it has no dB suffix):
``load(dump_yaml(s)) == s`` with no dB → linear → dB float error. The value as written in dB is
kept in the run record's ``config.input.yaml``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from cranberry.config.spec import ConfigError
from cranberry.config.yaml_io import Expr, Fn
from cranberry.core.units import db_to_lin, dbc_per_hz_to_s_phi, lin_to_db, s_phi_to_dbc_per_hz

__all__ = ["DB_RENAMES", "linear_name", "db_value_to_linear", "linear_value_to_db", "linearize", "LINEARIZED_SECTIONS"]

#: dB suffix → (linear suffix, scalar conversion, expression template). Longest suffix first.
DB_RENAMES: tuple[tuple[str, str, Any, str], ...] = (
    ("_dbc_per_hz", "_rad2_per_hz", lambda v: float(dbc_per_hz_to_s_phi(v)), "2 * 10 ** (({}) / 10)"),
    ("_db_per_hz", "_per_hz", lambda v: float(db_to_lin(v)), "10 ** (({}) / 10)"),
    ("_db", "_lin", lambda v: float(db_to_lin(v)), "10 ** (({}) / 10)"),
)

#: Top-level sections whose keys are linearised.
LINEARIZED_SECTIONS: tuple[str, ...] = ("topology", "offline", "pipeline", "measure")


def _rule(name: str) -> tuple[str, str, Any, str] | None:
    lname = name.lower()
    for rule in DB_RENAMES:
        if lname.endswith(rule[0]) and len(lname) > len(rule[0]):
            return rule
    return None


def linear_name(name: str) -> str | None:
    """The linear parameter name of a dB name (``il_db`` → ``il_lin``), or None if not a dB name."""
    rule = _rule(name)
    return None if rule is None else name[: -len(rule[0])] + rule[1]


def db_value_to_linear(name: str, value: Any, where: str = "") -> Any:
    """Linear value of ``value`` given under dB name ``name`` (table in the module docstring)."""
    rule = _rule(name)
    if rule is None:
        raise ConfigError(f"{where or name}: {name!r} is not a dB parameter name")
    _, _, conv, template = rule

    def one(v: Any, at: str) -> Any:
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return conv(v)
        if isinstance(v, list):
            return [one(x, f"{at}[{i}]") for i, x in enumerate(v)]
        if isinstance(v, Expr) and v.args:
            return Expr(template.format(v.source), dict(v.where), v.args)
        raise ConfigError(
            f"{at}: a dB parameter needs a number, a list of numbers or a callable !expr (converted to linear "
            f"on load, G16-i), got {v!r}; give the linear value under {linear_name(name)!r} instead"
        )

    return one(value, where or name)


#: Inverse conversions (linear value → value as written under the dB name), by dB suffix.
_TO_DB = {"_dbc_per_hz": lambda v: float(s_phi_to_dbc_per_hz(v)), "_db_per_hz": lambda v: float(lin_to_db(v)),
          "_db": lambda v: float(lin_to_db(v))}


def linear_value_to_db(name: str, value: Any, where: str = "") -> Any:
    """Inverse of :func:`db_value_to_linear` for numbers and (nested) lists of numbers: the value to
    write under dB name ``name`` so that it loads as the linear ``value`` (used by ``runner.sweep``
    for an axis named by the linear name of a parameter written in dB, WP-F2c)."""
    rule = _rule(name)
    if rule is None:
        raise ConfigError(f"{where or name}: {name!r} is not a dB parameter name")
    inv = _TO_DB[rule[0]]

    def one(v: Any, at: str) -> Any:
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return inv(v)
        if isinstance(v, list):
            return [one(x, f"{at}[{i}]") for i, x in enumerate(v)]
        raise ConfigError(f"{at}: a linear value for the dB parameter {name!r} must be a number or a list, got {v!r}")

    return one(value, where or name)


def _walk(node: Any, path: str) -> Any:
    if isinstance(node, Mapping):
        out: dict[str, Any] = {}
        renamed = {k: linear_name(k) for k in node if isinstance(k, str) and k != "where"}
        for k, v in node.items():
            at = f"{path}.{k}" if path else str(k)
            lin = renamed.get(k)
            if lin is None or k == "where":
                out[k] = v if k == "where" else _walk(v, at)
                continue
            if lin in node:
                raise ConfigError(f"{at}: both {k!r} (dB) and {lin!r} (linear) are given; keep one")
            out[lin] = db_value_to_linear(k, v, at)
        return out
    if isinstance(node, list):
        return [_walk(v, f"{path}[{i}]") for i, v in enumerate(node)]
    if isinstance(node, Fn):
        return Fn(node.name, _walk(node.args, f"{path}.args"))
    return node


def linearize(tree: Mapping[str, Any]) -> dict[str, Any]:
    """Return ``tree`` with every dB-named key of :data:`LINEARIZED_SECTIONS` converted and renamed."""
    return {k: _walk(v, k) if k in LINEARIZED_SECTIONS else v for k, v in tree.items()}
