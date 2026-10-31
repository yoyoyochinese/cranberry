"""YAML reading/writing with the project tags (G12, prompt §5.2). Owner: WP-F2.

Tags (names are the contract; presets written by WP-I* use exactly these):

    !sweep [a, b]         linear sweep axis            → :class:`Sweep`
    !sweep_log [a, b]     logarithmic sweep axis       → :class:`SweepLog`
    !choice [x, y, ...]   discrete sweep axis          → :class:`Choice`
    !expr "<expression>"  restricted expression (``config.expr``), with ``where:`` bindings → :class:`Expr`
    !fn {name: <registered>, args: {...}}   registered callable (``core.registry.register_fn``) → :class:`Fn`
    !pi_frac <x>          x·π (e.g. setpoints)         → :class:`PiFrac`

A single run requires every sweep marker to be resolved to one value (``runner.sweep`` expands them).
YAML is untrusted input: loading uses a safe loader; ``!expr`` is evaluated only by the restricted
evaluator, never by ``eval`` (G12).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["YAML_TAGS", "Sweep", "SweepLog", "Choice", "Expr", "Fn", "PiFrac", "load_yaml", "dump_yaml_obj"]

YAML_TAGS: tuple[str, ...] = ("!sweep", "!sweep_log", "!choice", "!expr", "!fn", "!pi_frac")


@dataclass(frozen=True)
class Sweep:
    lo: float
    hi: float


@dataclass(frozen=True)
class SweepLog:
    lo: float
    hi: float


@dataclass(frozen=True)
class Choice:
    options: tuple[Any, ...]


@dataclass(frozen=True)
class Expr:
    source: str
    where: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Fn:
    name: str
    args: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PiFrac:
    x: float


def load_yaml(path: str | Path) -> Any:
    """Parse a YAML file into plain data + tag marker objects (safe loader). WP-F2."""
    raise NotImplementedError("config.yaml_io.load_yaml is implemented by WP-F2")


def dump_yaml_obj(obj: Any) -> str:
    """Serialise plain data + tag markers back to YAML text (inverse of load_yaml). WP-F2."""
    raise NotImplementedError("config.yaml_io.dump_yaml_obj is implemented by WP-F2")
