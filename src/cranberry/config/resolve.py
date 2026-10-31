"""Dependency resolution of callable / cross-referencing parameters (G12, prompt §5.5). Owner: WP-F2.

``!expr`` / ``!fn`` values may reference other parameters by dotted path. These references form a
directed graph that is resolved in topological order; a cycle is an error naming the cycle.
Time-varying callables (``time_varying: true``, e.g. L_ic(x, P, t)) are tabulated offline on a grid
and interpolated in the loop.

Name lookup for a name ``a.b.c`` used in an ``!expr`` (lexical scoping, innermost first):

1. the expression's own bindings (``!expr {source, where}``): the full dotted name, then ``a``;
2. its callable arguments (``args``) — these stay free and make the value a callable;
3. walking from the mapping that holds the expression up to the root, at each level:
   a mapping's ``where:`` bindings, then the mapping's own keys; in a list, the element whose
   ``id`` / ``name`` equals ``a`` (or the longest dotted prefix ``a.b``, block names contain dots);
   the rest of the name navigates into the found node (keys; ``id``/``name`` in lists);
4. the constants ``pi``, ``e``, ``inf``.

A ``where:`` binding is a literal, a marker (``!sweep`` ...: the expression is *deferred* until
``runner.sweep`` fixes the value), or a string = dotted reference resolved from the enclosing scope.
``where`` is therefore a reserved key; it is removed from a mapping once nothing under that mapping
still needs it (resolved configs carry values, not bindings). References must point at scalar values
or at other expressions, never at whole mappings/lists or ``!fn`` values.

An ``!expr``'s *own* ``where`` (long form) may also hold a sweep marker or a nested ``!expr``: the
marker is a sweep axis at ``<path of the expr>.where.<name>`` (found by :func:`find_markers`,
expanded by ``runner.sweep``) and a nested expression is resolved first (WP-F2b).

``!fn`` arguments (WP-F2b) may contain references: write the argument as an ``!expr`` (a bare
reference is ``!expr "opo.L0"``); it is resolved like any other expression, with the names looked up
from the mapping holding the ``!fn`` (and the sibling arguments), and the value replaces it in
``args``. Plain strings in ``args`` stay strings (registered functions may take string options)::

    loss: !fn {name: refB_loss_vs_x, args: {L0: !expr "opo.L0", k: !expr "2 * k_half"}}

Paths into these are ``<path>.where.<name>`` / ``<path>.args.<name>`` (:func:`path_str`,
:func:`parse_path`, :func:`get_path`, :func:`set_path` step into ``Expr.where`` / ``Fn.args``).

After resolution: expressions without ``args`` are replaced by their value; expressions with
``args`` stay :class:`~cranberry.config.yaml_io.Expr` markers whose ``where`` holds the resolved
constants (blocks turn them into functions with :func:`as_callable`); ``!pi_frac x`` becomes the
float x·π; sweep markers and expressions depending on them are left untouched. A fully resolved tree
(no sweep axis left) then has its dB-named parameters converted to linear (``config.db``, G16-i).
"""

from __future__ import annotations

import copy
import functools
import graphlib
import math
import re
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from cranberry.config.db import linearize
from cranberry.config.expr import CONSTANTS, ExprError, ExprFunction, compile_expr, free_names, safe_eval
from cranberry.config.spec import ConfigError
from cranberry.config.yaml_io import SWEEP_MARKERS, Choice, Expr, Fn, PiFrac, _check_size

__all__ = [
    "resolve",
    "TreePath",
    "find_markers",
    "get_path",
    "set_path",
    "parse_path",
    "path_str",
    "as_callable",
    "normalize",
]

#: A location in a config tree: keys (str) and list indices (int).
TreePath = tuple[Any, ...]

_DEFERRED = object()
_MARKERS = (*SWEEP_MARKERS, Expr, Fn, PiFrac)


# ------------------------------------------------------------------------------------ tree navigation


def _child(node: Any, p: Any) -> Any:
    """``node[p]``; also steps into an ``!expr``'s own ``where`` and a ``!fn``'s ``args`` (WP-F2b)."""
    if isinstance(node, Expr) and p == "where":
        return node.where
    if isinstance(node, Fn) and p == "args":
        return node.args
    return node[p]


def get_path(tree: Any, path: TreePath) -> Any:
    node = tree
    for p in path:
        node = _child(node, p)
    return node


def set_path(tree: Any, path: TreePath, value: Any) -> None:
    get_path(tree, path[:-1])[path[-1]] = value


def _item_label(item: Any) -> str | None:
    if isinstance(item, Mapping):
        for key in ("id", "name"):
            if isinstance(item.get(key), str):
                return item[key]
    return None


def path_str(tree: Any, path: TreePath) -> str:
    """Human-readable path: keys joined by '.', list elements by their id/name or ``[i]``."""
    out = ""
    node = tree
    for p in path:
        if isinstance(p, int):
            label = _item_label(node[p])
            out += f".{label}" if label is not None else f"[{p}]"
        else:
            out += f".{p}"
        node = _child(node, p)
    return out.lstrip(".") or "<root>"


_TOKEN_RE = re.compile(r"\[(\d+)\]|([^.\[\]]+)")


def _tokens(dotted: str) -> list[Any]:
    toks: list[Any] = []
    pos = 0
    for m in _TOKEN_RE.finditer(dotted):
        toks.append(int(m.group(1)) if m.group(1) is not None else m.group(2))
        pos = m.end()
    if pos != len(dotted) or not toks:
        raise ConfigError(f"malformed path {dotted!r}")
    return toks


def _match_item(items: Sequence[Any], toks: Sequence[Any]) -> tuple[int, int] | None:
    """(index, tokens consumed) of the element whose id/name is the longest dotted prefix of toks."""
    names = []
    for t in toks:
        if not isinstance(t, str):
            break
        names.append(t)
    for j in range(len(names), 0, -1):
        label = ".".join(names[:j])
        for i, item in enumerate(items):
            if _item_label(item) == label:
                return i, j
    return None


def _navigate(tree: Any, path: TreePath, toks: Sequence[Any]) -> TreePath | None:
    node = get_path(tree, path)
    path = tuple(path)
    i = 0
    while i < len(toks):
        t = toks[i]
        if isinstance(node, Mapping):
            if not isinstance(t, str) or t not in node:
                return None
            path, node, i = path + (t,), node[t], i + 1
        elif isinstance(node, list):
            if isinstance(t, int):
                if not 0 <= t < len(node):
                    return None
                path, node, i = path + (t,), node[t], i + 1
                continue
            hit = _match_item(node, toks[i:])
            if hit is None:
                return None
            path, node, i = path + (hit[0],), node[hit[0]], i + hit[1]
        elif (isinstance(node, Expr) and t == "where") or (isinstance(node, Fn) and t == "args"):
            path, node, i = path + (t,), _child(node, t), i + 1
        else:
            return None
    return path


def parse_path(tree: Any, dotted: str) -> TreePath:
    """Path of ``dotted`` (as printed by :func:`path_str`) from the root of ``tree``."""
    found = _navigate(tree, (), _tokens(dotted))
    if found is None:
        raise ConfigError(f"no config value at {dotted!r}")
    return found


def _lookup(tree: Any, scope: TreePath, parts: list[str], skip_where_at: TreePath | None = None) -> TreePath | None:
    """Scoped lookup of a dotted name (module docstring, step 3)."""
    for i in range(len(scope), -1, -1):
        q = tuple(scope[:i])
        node = get_path(tree, q)
        if isinstance(node, Mapping):
            where = node.get("where")
            if q != skip_where_at and isinstance(where, Mapping) and parts[0] in where:
                b = where[parts[0]]
                if isinstance(b, str):
                    base = _lookup(tree, q, _tokens(b), skip_where_at=q)
                    if base is None:
                        raise ConfigError(f"{path_str(tree, q + ('where', parts[0]))}: cannot resolve reference {b!r}")
                    return _navigate(tree, base, parts[1:])
                return _navigate(tree, q + ("where", parts[0]), parts[1:])
            if parts[0] in node and parts[0] != "where":
                return _navigate(tree, q + (parts[0],), parts[1:])
        elif isinstance(node, list):
            hit = _match_item(node, parts)
            if hit is not None:
                return _navigate(tree, q + (hit[0],), parts[hit[1]:])
    return None


def _walk(tree: Any, path: TreePath = ()):
    """Yield (path, node) depth-first; descends into an ``!expr``'s own ``where`` and a ``!fn``'s
    ``args`` (so sweep markers and expressions there are found), not into other markers."""
    yield path, tree
    if isinstance(tree, Mapping):
        for k, v in tree.items():
            yield from _walk(v, path + (k,))
    elif isinstance(tree, list):
        for i, v in enumerate(tree):
            yield from _walk(v, path + (i,))
    elif isinstance(tree, Expr):
        for k, v in tree.where.items():
            yield from _walk(v, path + ("where", k))
    elif isinstance(tree, Fn):
        for k, v in tree.args.items():
            yield from _walk(v, path + ("args", k))


def find_markers(tree: Any, types: tuple[type, ...] = SWEEP_MARKERS) -> list[tuple[TreePath, Any]]:
    """(path, marker) of every marker of ``types`` in tree order (default: sweep axes)."""
    return [(p, n) for p, n in _walk(tree) if isinstance(n, types)]


# ----------------------------------------------------------------------------------------- resolution


def _value_at(tree: Any, target: TreePath, where: str) -> Any:
    node = get_path(tree, target)
    if isinstance(node, SWEEP_MARKERS):
        return _DEFERRED
    if isinstance(node, PiFrac):
        return node.x * math.pi
    if isinstance(node, (Mapping, list)):
        raise ConfigError(f"{where}: reference to {path_str(tree, target)} is a mapping/list, not a value")
    if isinstance(node, Fn):
        raise ConfigError(f"{where}: reference to {path_str(tree, target)} is a !fn callable, not a value")
    return node


def _bind(tree: Any, path: TreePath, expr: Expr, where: str) -> dict[str, tuple[str, Any]]:
    container = path[:-1]
    out: dict[str, tuple[str, Any]] = {}
    for name in free_names(expr.source):
        parts = name.split(".")
        key = name if name in expr.where else (parts[0] if parts[0] in expr.where else None)
        if key is not None:
            b = expr.where[key]
            rest = parts[len(key.split(".")):]
            if isinstance(b, str):
                base = _lookup(tree, container, _tokens(b))
                target = None if base is None else _navigate(tree, base, rest)
                if target is None:
                    raise ConfigError(f"{where}: cannot resolve where-binding {key}: {b!r}")
                out[name] = ("path", target)
            else:
                if rest:
                    raise ConfigError(f"{where}: {name!r} navigates into a literal binding")
                # a literal (or marker / nested !expr) bound in the expression's own where: a path into
                # it, so a !sweep there is a sweep axis and a nested !expr is resolved first (WP-F2b)
                out[name] = ("path", path + ("where", key))
        elif name in expr.args:
            out[name] = ("arg", None)
        else:
            target = _lookup(tree, container, parts)
            if target is not None:
                out[name] = ("path", target)
            elif name not in CONSTANTS:
                raise ConfigError(
                    f"{where}: !expr {expr.source!r}: cannot resolve name {name!r}; bind it with where: or "
                    f"declare it as a callable argument (!expr {{source: ..., args: [{name}]}})"
                )
    return out


def resolve(tree: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve all callables / references in a merged config tree (sweep markers untouched)."""
    _check_size(tree)
    tree = copy.deepcopy(dict(tree))
    for p, node in _walk(tree):
        if isinstance(node, Mapping) and "where" in node and not isinstance(node["where"], Mapping):
            raise ConfigError(f"{path_str(tree, p + ('where',))}: where: must be a mapping of bindings")

    exprs: dict[TreePath, Expr] = {p: n for p, n in _walk(tree) if isinstance(n, Expr)}
    bindings: dict[TreePath, dict[str, tuple[str, Any]]] = {}
    graph: dict[TreePath, set[TreePath]] = {}
    for p, e in exprs.items():
        where = path_str(tree, p)
        bindings[p] = b = _bind(tree, p, e, where)
        deps: set[TreePath] = set()
        for kind, target in b.values():
            if kind != "path":
                continue
            if target in exprs:
                deps.add(target)
            elif any(x[: len(target)] == target for x in exprs):
                raise ConfigError(f"{where}: reference to {path_str(tree, target)} contains expressions, not a value")
        graph[p] = deps
    try:
        order = list(graphlib.TopologicalSorter(graph).static_order())
    except graphlib.CycleError as e:
        cycle = " -> ".join(path_str(tree, p) for p in e.args[1])
        raise ConfigError(f"circular !expr references: {cycle}") from None

    results: dict[TreePath, Any] = {}
    for p in order:
        e = exprs[p]
        where = path_str(tree, p)
        env: dict[str, Any] = {}
        deferred = False
        for name, (kind, target) in bindings[p].items():
            if kind == "arg":
                continue
            if target in results:
                v = results[target]
                if isinstance(v, Expr):
                    raise ConfigError(f"{where}: {name!r} refers to a callable !expr ({path_str(tree, target)})")
            else:
                v = _value_at(tree, target, where)
            if v is _DEFERRED:
                deferred = True
            env[name] = v
        if deferred:
            results[p] = _DEFERRED
            continue
        try:
            if e.args:
                for name, v in env.items():
                    if isinstance(v, str):
                        raise ConfigError(f"{where}: callable !expr cannot capture string value {name}={v!r}")
                compile_expr(e.source, e.args, env)
                results[p] = Expr(e.source, env, e.args)
            else:
                results[p] = safe_eval(e.source, env)
        except ExprError as err:
            raise ExprError(f"{where}: {err}") from None

    deferred_paths = [p for p, r in results.items() if r is _DEFERRED]
    for p, r in results.items():
        if r is not _DEFERRED:
            set_path(tree, p, r)

    tree = _convert_pi_frac(tree)
    _strip_where(tree, (), deferred_paths)
    if find_markers(tree):  # sweep axes left: keep dB names as written (config.db, "When")
        return tree
    return linearize(tree)


def _convert_pi_frac(node: Any) -> Any:
    if isinstance(node, PiFrac):
        return node.x * math.pi
    if isinstance(node, Mapping):
        return {k: _convert_pi_frac(v) for k, v in node.items()}
    if isinstance(node, list):
        return [_convert_pi_frac(v) for v in node]
    if isinstance(node, Fn):
        return Fn(node.name, _convert_pi_frac(node.args))
    return node


def _strip_where(node: Any, path: TreePath, deferred: list[TreePath]) -> None:
    if isinstance(node, dict):
        if "where" in node:
            n = len(path)
            if not any(d[:n] == path and len(d) > n and d[n] != "where" for d in deferred):
                del node["where"]
        for k, v in node.items():
            _strip_where(v, path + (k,), deferred)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            _strip_where(v, path + (i,), deferred)


# ---------------------------------------------------------------------------- callables and normalising


def as_callable(value: Any, arg_names: tuple[str, ...]) -> Callable[..., Any]:
    """Turn a callable-parameter value into a Python callable (G16 a: constant, callable or table).

    Accepts a Python callable (returned as is), ``Fn`` (registered function with keyword args), a
    resolved callable ``Expr`` (``args`` must equal ``arg_names``), or a number (constant function).
    """
    from cranberry.core.registry import get_fn

    if isinstance(value, Fn):
        return functools.partial(get_fn(value.name), **value.args)
    if isinstance(value, Expr):
        if value.args and tuple(value.args) != tuple(arg_names):
            raise ConfigError(f"!expr {value.source!r} has args {value.args}, expected {tuple(arg_names)}")
        return compile_expr(value.source, tuple(arg_names), value.where)
    if callable(value):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        const = value
        return lambda *args: const
    raise ConfigError(f"cannot use {value!r} as a callable parameter of {tuple(arg_names)}")


def normalize(obj: Any, path: str = "") -> Any:
    """Plain config data from Python values (``to_spec``): tuples → lists, numpy → Python, registered
    functions → ``!fn``, :class:`ExprFunction` → ``!expr``. Unregistered callables (lambdas) raise
    :class:`ConfigError` naming ``path`` (G12: they cannot round-trip to YAML)."""
    from cranberry.core.registry import fn_name

    if isinstance(obj, np.generic):  # before the float check: np.float64 subclasses float
        return obj.item()
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Mapping):
        out = {}
        for k, v in obj.items():
            if not isinstance(k, str):
                raise ConfigError(f"{path}: mapping key {k!r} must be a string")
            out[k] = normalize(v, f"{path}.{k}" if path else k)
        return out
    if isinstance(obj, (list, tuple)):
        return [normalize(v, f"{path}[{i}]") for i, v in enumerate(obj)]
    if isinstance(obj, Choice):
        return Choice(tuple(normalize(v, path) for v in obj.options))
    if isinstance(obj, Fn):
        return Fn(obj.name, normalize(obj.args, path))
    if isinstance(obj, Expr):
        return Expr(obj.source, normalize(obj.where, path), obj.args)
    if isinstance(obj, _MARKERS):
        return obj
    if isinstance(obj, ExprFunction):
        return Expr(obj.source, normalize(obj.env, path), obj.arg_names)
    if callable(obj):
        target, kwargs = obj, {}
        if isinstance(obj, functools.partial):
            if obj.args:
                raise ConfigError(f"{path}: functools.partial with positional args cannot be written to YAML")
            target, kwargs = obj.func, dict(obj.keywords)
        name = fn_name(target)
        if name is None:
            what = getattr(target, "__qualname__", repr(target))
            raise ConfigError(
                f"{path}: unnamed / unregistered callable {what!r} cannot be serialised to YAML (G12); "
                "register it with @register_fn(name) or use config.expr.compile_expr"
            )
        return Fn(name, normalize(kwargs, path))
    raise ConfigError(f"{path}: value of type {type(obj).__name__} cannot be represented in a config")
