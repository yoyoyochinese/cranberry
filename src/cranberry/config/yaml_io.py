"""YAML reading/writing with the project tags (G12, prompt §5.2). Owner: WP-F2.

Tags (names are the contract; presets written by WP-I* use exactly these):

    !sweep [a, b]         linear sweep axis            → :class:`Sweep`
    !sweep [a, b, n]      ... with n grid points (needed for grid expansion by ``runner.sweep``)
    !sweep_log [a, b]     logarithmic sweep axis       → :class:`SweepLog`  (also ``[a, b, n]``)
    !choice [x, y, ...]   discrete sweep axis          → :class:`Choice`
    !expr "<expression>"  restricted expression (``config.expr``), with ``where:`` bindings → :class:`Expr`
    !expr {source: "<expression>", where: {...}, args: [x]}   long form; ``args`` = callable arguments
    !fn {name: <registered>, args: {...}}   registered callable (``core.registry.register_fn``) → :class:`Fn`
    !pi_frac <x>          x·π (e.g. setpoints)         → :class:`PiFrac`

A single run requires every sweep marker to be resolved to one value (``runner.sweep`` expands them).
YAML is untrusted input: loading uses a ``yaml.SafeLoader`` subclass; ``!expr`` is evaluated only by
the restricted evaluator, never by ``eval`` (G12). Duplicate mapping keys are an error, and the
expanded size of a document (anchors/aliases) is bounded.

Differences from plain PyYAML ``safe_load`` (on purpose, so YAML and Python values agree):

* ``1e-8`` / ``1E+3`` (exponent without a dot) are floats, not strings (YAML 1.2 behaviour).
* Timestamps are not converted to ``datetime`` (they stay strings; the config tree is JSON-able).
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

__all__ = [
    "YAML_TAGS",
    "Sweep",
    "SweepLog",
    "Choice",
    "Expr",
    "Fn",
    "PiFrac",
    "SWEEP_MARKERS",
    "MAX_NODES",
    "load_yaml",
    "loads_yaml",
    "dump_yaml_obj",
    "canonical_json",
]

YAML_TAGS: tuple[str, ...] = ("!sweep", "!sweep_log", "!choice", "!expr", "!fn", "!pi_frac")

#: Upper bound on the number of nodes of a loaded document after alias expansion ("billion laughs").
MAX_NODES = 1_000_000


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _check_range(tag: str, lo: Any, hi: Any, n: Any) -> None:
    if not (_is_number(lo) and _is_number(hi)):
        raise ValueError(f"{tag}: bounds must be numbers, got {lo!r}, {hi!r}")
    if n is not None and (isinstance(n, bool) or not isinstance(n, int) or n < 1):
        raise ValueError(f"{tag}: number of points must be an int >= 1, got {n!r}")


@dataclass(frozen=True)
class Sweep:
    """Linear sweep axis ``[lo, hi]``; ``n`` grid points (None: range only, e.g. for optimisers)."""

    lo: float
    hi: float
    n: int | None = None

    def __post_init__(self) -> None:
        _check_range("!sweep", self.lo, self.hi, self.n)


@dataclass(frozen=True)
class SweepLog:
    """Logarithmic sweep axis ``[lo, hi]`` (same sign, non-zero); ``n`` grid points or None."""

    lo: float
    hi: float
    n: int | None = None

    def __post_init__(self) -> None:
        _check_range("!sweep_log", self.lo, self.hi, self.n)
        if self.lo == 0 or self.hi == 0 or (self.lo > 0) != (self.hi > 0):
            raise ValueError(f"!sweep_log: bounds must be non-zero and of the same sign, got {self.lo!r}, {self.hi!r}")


@dataclass(frozen=True)
class Choice:
    """Discrete sweep axis."""

    options: tuple[Any, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.options, tuple) or not self.options:
            raise ValueError("!choice needs a non-empty list of options")


@dataclass(frozen=True)
class Expr:
    """Restricted expression (``config.expr``).

    Attributes:
        source: the expression text.
        where: own bindings {name: value or dotted-path string}; usually the sibling ``where:`` key of
            the mapping holding the expression is used instead (PLAN §2.9 example).
        args: names of callable arguments (e.g. ``("x",)`` for Ref B's L(x)); empty → the expression
            is evaluated to a value at load time.
    """

    source: str
    where: dict[str, Any] = field(default_factory=dict)
    args: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("!expr needs a non-empty expression string")
        if not isinstance(self.where, dict) or not all(isinstance(k, str) for k in self.where):
            raise ValueError("!expr where: must be a mapping of names")
        if not isinstance(self.args, tuple) or not all(isinstance(a, str) and a.isidentifier() for a in self.args):
            raise ValueError(f"!expr args must be a list of identifiers, got {self.args!r}")


@dataclass(frozen=True)
class Fn:
    """Registered callable ``core.registry.get_fn(name)`` with keyword arguments ``args``."""

    name: str
    args: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("!fn needs a non-empty registered function name")
        if not isinstance(self.args, dict) or not all(isinstance(k, str) for k in self.args):
            raise ValueError("!fn args must be a mapping of keyword names")


@dataclass(frozen=True)
class PiFrac:
    """x·π (resolved to a float by ``config.resolve``)."""

    x: float

    def __post_init__(self) -> None:
        if not _is_number(self.x):
            raise ValueError(f"!pi_frac needs a number, got {self.x!r}")


#: Markers that define sweep axes (expanded by ``runner.sweep``).
SWEEP_MARKERS = (Sweep, SweepLog, Choice)


# --------------------------------------------------------------------------------------------- loader

_FLOAT_RE = re.compile(
    r"""^(?:[-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |[-+]?\.[0-9_]+(?:[eE][-+]?[0-9]+)?
    |[-+]?\.(?:inf|Inf|INF)
    |\.(?:nan|NaN|NAN))$""",
    re.X,
)
_FLOAT_TAG = "tag:yaml.org,2002:float"
_TIMESTAMP_TAG = "tag:yaml.org,2002:timestamp"


class _Loader(yaml.SafeLoader):
    """SafeLoader with YAML-1.2 floats, no timestamps, duplicate-key errors and the project tags."""


_Loader.yaml_implicit_resolvers = {
    ch: [(tag, rx) for tag, rx in resolvers if tag not in (_FLOAT_TAG, _TIMESTAMP_TAG)]
    for ch, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
_Loader.add_implicit_resolver(_FLOAT_TAG, _FLOAT_RE, list("-+0123456789."))


def _construct_mapping(loader: _Loader, node: yaml.MappingNode, deep: bool = False) -> dict[Any, Any]:
    loader.flatten_mapping(node)
    out: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            hash(key)
        except TypeError:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark, "found unhashable key", key_node.start_mark
            ) from None
        if key in out:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping", node.start_mark, f"found duplicate key {key!r}", key_node.start_mark
            )
        out[key] = loader.construct_object(value_node, deep=deep)
    return out


def _construct_yaml_map(loader: _Loader, node: yaml.MappingNode):
    data: dict[Any, Any] = {}
    yield data
    data.update(_construct_mapping(loader, node))


_Loader.add_constructor("tag:yaml.org,2002:map", _construct_yaml_map)


def _err(node: yaml.Node, msg: str) -> yaml.constructor.ConstructorError:
    return yaml.constructor.ConstructorError(None, None, msg, node.start_mark)


def _range_marker(cls: type, tag: str):
    def construct(loader: _Loader, node: yaml.Node):
        try:
            if isinstance(node, yaml.SequenceNode):
                items = loader.construct_sequence(node, deep=True)
                if len(items) not in (2, 3):
                    raise ValueError(f"{tag} takes [lo, hi] or [lo, hi, n], got {items!r}")
                return cls(*items)
            if isinstance(node, yaml.MappingNode):
                m = _construct_mapping(loader, node, deep=True)
                if set(m) - {"lo", "hi", "n"} or not {"lo", "hi"} <= set(m):
                    raise ValueError(f"{tag} mapping form is {{lo, hi, n?}}, got {sorted(m)}")
                return cls(m["lo"], m["hi"], m.get("n"))
            raise ValueError(f"{tag} needs a sequence [lo, hi(, n)]")
        except ValueError as e:
            raise _err(node, str(e)) from None

    return construct


def _construct_choice(loader: _Loader, node: yaml.Node) -> Choice:
    if not isinstance(node, yaml.SequenceNode):
        raise _err(node, "!choice needs a sequence of options")
    try:
        return Choice(tuple(loader.construct_sequence(node, deep=True)))
    except ValueError as e:
        raise _err(node, str(e)) from None


def _construct_expr(loader: _Loader, node: yaml.Node) -> Expr:
    try:
        if isinstance(node, yaml.ScalarNode):
            return Expr(loader.construct_scalar(node))
        if isinstance(node, yaml.MappingNode):
            m = _construct_mapping(loader, node, deep=True)
            if set(m) - {"source", "where", "args"} or "source" not in m:
                raise ValueError(f"!expr mapping form is {{source, where?, args?}}, got {sorted(m)}")
            return Expr(m["source"], dict(m.get("where") or {}), tuple(m.get("args") or ()))
        raise ValueError("!expr needs a string or {source, where?, args?}")
    except (ValueError, TypeError) as e:
        raise _err(node, str(e)) from None


def _construct_fn(loader: _Loader, node: yaml.Node) -> Fn:
    try:
        if isinstance(node, yaml.ScalarNode):
            return Fn(loader.construct_scalar(node))
        if isinstance(node, yaml.MappingNode):
            m = _construct_mapping(loader, node, deep=True)
            if set(m) - {"name", "args"} or "name" not in m:
                raise ValueError(f"!fn mapping form is {{name, args?}}, got {sorted(m)}")
            return Fn(m["name"], dict(m.get("args") or {}))
        raise ValueError("!fn needs a name or {name, args?}")
    except (ValueError, TypeError) as e:
        raise _err(node, str(e)) from None


def _construct_pi_frac(loader: _Loader, node: yaml.Node) -> PiFrac:
    if not isinstance(node, yaml.ScalarNode):
        raise _err(node, "!pi_frac needs a number")
    text = loader.construct_scalar(node)
    try:
        value = int(text) if re.fullmatch(r"[-+]?[0-9]+", text) else float(text)
    except ValueError:
        raise _err(node, f"!pi_frac needs a number, got {text!r}") from None
    return PiFrac(value)


_Loader.add_constructor("!sweep", _range_marker(Sweep, "!sweep"))
_Loader.add_constructor("!sweep_log", _range_marker(SweepLog, "!sweep_log"))
_Loader.add_constructor("!choice", _construct_choice)
_Loader.add_constructor("!expr", _construct_expr)
_Loader.add_constructor("!fn", _construct_fn)
_Loader.add_constructor("!pi_frac", _construct_pi_frac)


def _check_size(obj: Any, limit: int = MAX_NODES) -> None:
    """Count nodes with alias expansion; stop early once ``limit`` is exceeded."""
    count = 0
    stack = [obj]
    while stack:
        x = stack.pop()
        count += 1
        if count > limit:
            raise ValueError(f"YAML document too large after alias expansion (> {limit} nodes)")
        if isinstance(x, Mapping):
            stack.extend(x.values())
        elif isinstance(x, (list, tuple)):
            stack.extend(x)
        elif isinstance(x, Choice):
            stack.extend(x.options)
        elif isinstance(x, (Expr, Fn)):
            stack.extend((x.where if isinstance(x, Expr) else x.args).values())


def loads_yaml(text: str) -> Any:
    """Parse YAML text into plain data + tag marker objects (safe loader, one document)."""
    docs = list(yaml.load_all(text, Loader=_Loader))  # noqa: S506 - _Loader is a SafeLoader subclass
    if len(docs) > 1:
        raise ValueError("expected a single YAML document")
    obj = docs[0] if docs else None
    _check_size(obj)
    return obj


def load_yaml(path: str | Path) -> Any:
    """Parse a YAML file into plain data + tag marker objects (safe loader)."""
    return loads_yaml(Path(path).read_text(encoding="utf-8"))


# --------------------------------------------------------------------------------------------- dumper


class _Dumper(yaml.SafeDumper):
    """SafeDumper emitting the project tags; uses the loader's implicit resolvers so that strings
    such as ``"1e-8"`` are quoted (and read back as strings)."""

    def ignore_aliases(self, data: Any) -> bool:
        return True


_Dumper.yaml_implicit_resolvers = _Loader.yaml_implicit_resolvers


def _represent_float(dumper: _Dumper, data: float) -> yaml.Node:
    if math.isnan(data):
        text = ".nan"
    elif math.isinf(data):
        text = ".inf" if data > 0 else "-.inf"
    else:
        text = repr(float(data))
        if "e" in text and "." not in text:  # 1e-08 -> 1.0e-08: readable by YAML-1.1 parsers too
            text = text.replace("e", ".0e", 1)
    return dumper.represent_scalar(_FLOAT_TAG, text)


def _represent_range(tag: str):
    def represent(dumper: _Dumper, data: Sweep | SweepLog) -> yaml.Node:
        items = [data.lo, data.hi] + ([data.n] if data.n is not None else [])
        return dumper.represent_sequence(tag, items, flow_style=True)

    return represent


def _represent_expr(dumper: _Dumper, data: Expr) -> yaml.Node:
    if not data.where and not data.args:
        return dumper.represent_scalar("!expr", data.source, style='"')
    m: dict[str, Any] = {"source": data.source}
    if data.where:
        m["where"] = data.where
    if data.args:
        m["args"] = list(data.args)
    return dumper.represent_mapping("!expr", m, flow_style=True)


def _represent_fn(dumper: _Dumper, data: Fn) -> yaml.Node:
    return dumper.represent_mapping("!fn", {"name": data.name, "args": data.args}, flow_style=True)


def _represent_pi_frac(dumper: _Dumper, data: PiFrac) -> yaml.Node:
    return dumper.represent_scalar("!pi_frac", repr(data.x))


_Dumper.add_representer(float, _represent_float)
_Dumper.add_representer(tuple, lambda d, data: d.represent_list(list(data)))
_Dumper.add_representer(Sweep, _represent_range("!sweep"))
_Dumper.add_representer(SweepLog, _represent_range("!sweep_log"))
_Dumper.add_representer(Choice, lambda d, data: d.represent_sequence("!choice", list(data.options), flow_style=True))
_Dumper.add_representer(Expr, _represent_expr)
_Dumper.add_representer(Fn, _represent_fn)
_Dumper.add_representer(PiFrac, _represent_pi_frac)


def dump_yaml_obj(obj: Any) -> str:
    """Serialise plain data + tag markers back to YAML text (inverse of :func:`loads_yaml`)."""
    return yaml.dump(obj, Dumper=_Dumper, sort_keys=False, allow_unicode=True, default_flow_style=False, width=120)


# ------------------------------------------------------------------------------------- canonical form


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (Sweep, SweepLog)):
        return {f"!{type(obj).__name__.lower()}": [obj.lo, obj.hi, obj.n]}
    if isinstance(obj, Choice):
        return {"!choice": _jsonable(obj.options)}
    if isinstance(obj, Expr):
        return {"!expr": [obj.source, _jsonable(obj.where), list(obj.args)]}
    if isinstance(obj, Fn):
        return {"!fn": [obj.name, _jsonable(obj.args)]}
    if isinstance(obj, PiFrac):
        return {"!pi_frac": obj.x}
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    raise TypeError(f"not a config value: {type(obj).__name__} {obj!r}")


def canonical_json(obj: Any) -> str:
    """Deterministic JSON text of a config tree (sorted keys, markers encoded); used for hashes."""
    return json.dumps(_jsonable(obj), sort_keys=True, separators=(",", ":"), allow_nan=True)
