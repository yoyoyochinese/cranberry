"""Restricted expression evaluator for ``!expr`` (G12). Owner: WP-F2.

YAML is untrusted input, so expressions are parsed with :mod:`ast` and interpreted node by node by
this module; ``eval``/``exec``/``compile`` of the source are never used and no Python attribute is
ever looked up on a value.

Grammar (anything else is rejected with :class:`ExprError` before evaluation)::

    expr     := literal | name | dotted | expr binop expr | unop expr | expr cmp expr [cmp expr ...]
              | expr 'if' expr 'else' expr | expr ('and'|'or') expr | func '(' expr, ... ')' | '(' expr ')'
    literal  := int | float | bool | str          (str only for comparisons / ternary results)
    binop    := '+' | '-' | '*' | '/' | '**'
    unop     := '-' | '+' | 'not'
    cmp      := '==' | '!=' | '<' | '<=' | '>' | '>='
    name     := identifier not starting with '_'
    dotted   := name '.' name ['.' name ...]      a *config reference* (e.g. ``arm1_915b.L_m``),
                                                  looked up as one key in ``env``; never getattr
    func     := a name of :data:`FUNCTIONS` (fixed arity, positional arguments only)

Constants ``pi``, ``e``, ``inf`` are available unless ``env`` binds the same name. Arithmetic needs
numeric operands (no string repetition/concatenation); integer powers with |exponent| > 1024 and
integers longer than 4096 bits are rejected (resource bound); source length ≤ 2000 characters,
≤ 500 AST nodes, nesting depth ≤ 64. In :func:`safe_eval`, floating-point errors (log(0), sqrt(-1),
overflow, division by zero) raise :class:`ExprError` instead of producing inf/nan.

``compile_expr`` returns an :class:`ExprFunction` (callable parameters such as Ref B's L(x), prompt
§5.5); it evaluates with numpy semantics so it also works element-wise on arrays.
"""

from __future__ import annotations

import ast
import math
import operator
from collections.abc import Mapping
from functools import lru_cache
from typing import Any

import numpy as np

from cranberry.config.spec import ConfigError

__all__ = [
    "ExprError",
    "FUNCTIONS",
    "CONSTANTS",
    "parse",
    "free_names",
    "safe_eval",
    "compile_expr",
    "ExprFunction",
]

MAX_SOURCE_LEN = 2000
MAX_NODES = 500
MAX_DEPTH = 64
MAX_INT_BITS = 4096
MAX_INT_EXPONENT = 1024


class ExprError(ConfigError):
    """An ``!expr`` is not in the allowed grammar, references an unknown name, or fails to evaluate."""


#: Whitelisted elementary functions: name -> (implementation, arity). numpy versions so that
#: :class:`ExprFunction` works on arrays; math-style aliases are provided.
FUNCTIONS: dict[str, tuple[Any, int]] = {
    "sqrt": (np.sqrt, 1),
    "exp": (np.exp, 1),
    "expm1": (np.expm1, 1),
    "log": (np.log, 1),
    "log10": (np.log10, 1),
    "log2": (np.log2, 1),
    "log1p": (np.log1p, 1),
    "sin": (np.sin, 1),
    "cos": (np.cos, 1),
    "tan": (np.tan, 1),
    "arcsin": (np.arcsin, 1),
    "arccos": (np.arccos, 1),
    "arctan": (np.arctan, 1),
    "arctan2": (np.arctan2, 2),
    "asin": (np.arcsin, 1),
    "acos": (np.arccos, 1),
    "atan": (np.arctan, 1),
    "atan2": (np.arctan2, 2),
    "sinh": (np.sinh, 1),
    "cosh": (np.cosh, 1),
    "tanh": (np.tanh, 1),
    "arcsinh": (np.arcsinh, 1),
    "arccosh": (np.arccosh, 1),
    "arctanh": (np.arctanh, 1),
    "hypot": (np.hypot, 2),
    "abs": (np.abs, 1),
    "sign": (np.sign, 1),
    "floor": (np.floor, 1),
    "ceil": (np.ceil, 1),
    "min": (np.minimum, 2),
    "max": (np.maximum, 2),
    "deg2rad": (np.deg2rad, 1),
    "rad2deg": (np.rad2deg, 1),
    "radians": (np.deg2rad, 1),
    "degrees": (np.rad2deg, 1),
}

#: Named constants (lowest lookup priority; ``env`` may shadow them).
CONSTANTS: dict[str, float] = {"pi": math.pi, "e": math.e, "inf": math.inf}

_BINOPS: dict[type, Any] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}
_UNARY: dict[type, Any] = {ast.UAdd: operator.pos, ast.USub: operator.neg, ast.Not: operator.not_}
_CMPOPS: dict[type, Any] = {
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
}
_OP_SYMBOL = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/", ast.Pow: "**"}


def _dotted(node: ast.AST) -> str | None:
    """'a.b.c' for a Name / chain of Attribute on a Name, else None."""
    parts: list[str] = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return ".".join(reversed(parts))


def _check(node: ast.AST, source: str, depth: int) -> None:
    if depth > MAX_DEPTH:
        raise ExprError(f"!expr {source!r}: nesting deeper than {MAX_DEPTH}")
    d = depth + 1
    if isinstance(node, ast.Constant):
        if type(node.value) not in (int, float, bool, str):
            raise ExprError(f"!expr {source!r}: literal {node.value!r} is not allowed (numbers, bools, strings only)")
    elif isinstance(node, (ast.Name, ast.Attribute)):
        name = _dotted(node)
        if name is None:
            raise ExprError(f"!expr {source!r}: attribute access is not allowed (only dotted config references)")
        if not isinstance(node.ctx, ast.Load):
            raise ExprError(f"!expr {source!r}: assignment is not allowed")
        for part in name.split("."):
            if part.startswith("_"):
                raise ExprError(f"!expr {source!r}: names starting with '_' (incl. dunder) are not allowed: {name!r}")
    elif isinstance(node, ast.BinOp):
        if type(node.op) not in _BINOPS:
            raise ExprError(f"!expr {source!r}: operator {type(node.op).__name__} is not allowed (+ - * / ** only)")
        _check(node.left, source, d)
        _check(node.right, source, d)
    elif isinstance(node, ast.UnaryOp):
        if type(node.op) not in _UNARY:
            raise ExprError(f"!expr {source!r}: unary operator {type(node.op).__name__} is not allowed")
        _check(node.operand, source, d)
    elif isinstance(node, ast.Compare):
        for op in node.ops:
            if type(op) not in _CMPOPS:
                raise ExprError(f"!expr {source!r}: comparison {type(op).__name__} is not allowed")
        _check(node.left, source, d)
        for c in node.comparators:
            _check(c, source, d)
    elif isinstance(node, ast.BoolOp):
        for v in node.values:
            _check(v, source, d)
    elif isinstance(node, ast.IfExp):
        _check(node.test, source, d)
        _check(node.body, source, d)
        _check(node.orelse, source, d)
    elif isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS:
            what = node.func.id if isinstance(node.func, ast.Name) else ast.unparse(node.func)
            raise ExprError(
                f"!expr {source!r}: call of {what!r} is not allowed; allowed functions: {sorted(FUNCTIONS)}"
            )
        if node.keywords or any(isinstance(a, ast.Starred) for a in node.args):
            raise ExprError(f"!expr {source!r}: only positional arguments are allowed in calls")
        arity = FUNCTIONS[node.func.id][1]
        if len(node.args) != arity:
            raise ExprError(f"!expr {source!r}: {node.func.id}() takes {arity} argument(s), got {len(node.args)}")
        for a in node.args:
            _check(a, source, d)
    else:
        raise ExprError(f"!expr {source!r}: {type(node).__name__} is not allowed in expressions")


@lru_cache(maxsize=4096)
def parse(source: str) -> ast.Expression:
    """Parse and validate ``source`` against the grammar; returns the AST (cached)."""
    if not isinstance(source, str):
        raise ExprError(f"!expr source must be a string, got {type(source).__name__}")
    if len(source) > MAX_SOURCE_LEN:
        raise ExprError(f"!expr longer than {MAX_SOURCE_LEN} characters")
    try:
        tree = ast.parse(source.strip(), mode="eval")
    except (SyntaxError, ValueError, RecursionError, MemoryError) as e:
        raise ExprError(f"!expr {source!r}: syntax error: {e}") from None
    if sum(1 for _ in ast.walk(tree)) > MAX_NODES:
        raise ExprError(f"!expr {source!r}: more than {MAX_NODES} syntax nodes")
    _check(tree.body, source, 0)
    return tree


def free_names(source: str) -> tuple[str, ...]:
    """Names / dotted references used by ``source`` (function names excluded), in first-use order."""
    out: dict[str, None] = {}

    def visit(node: ast.AST) -> None:
        name = _dotted(node) if isinstance(node, (ast.Name, ast.Attribute)) else None
        if name is not None:
            out[name] = None
            return
        if isinstance(node, ast.Call):
            for a in node.args:
                visit(a)
            return
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(parse(source).body)
    return tuple(out)


def _is_num(x: Any) -> bool:
    return isinstance(x, (int, float, np.number, np.ndarray)) and not isinstance(x, (str, np.str_))


def _check_int(x: Any, source: str) -> Any:
    if isinstance(x, int) and not isinstance(x, bool) and x.bit_length() > MAX_INT_BITS:
        raise ExprError(f"!expr {source!r}: integer result too large")
    if isinstance(x, (complex, np.complexfloating)):
        raise ExprError(f"!expr {source!r}: complex result (e.g. negative base with fractional exponent)")
    return x


class _Evaluator:
    def __init__(self, source: str, env: Mapping[str, Any]):
        self.source = source
        self.env = env

    def lookup(self, name: str) -> Any:
        if name in self.env:
            return self.env[name]
        if name in CONSTANTS:
            return CONSTANTS[name]
        raise ExprError(f"!expr {self.source!r}: unknown name {name!r}")

    def ev(self, node: ast.AST) -> Any:
        src = self.source
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, (ast.Name, ast.Attribute)):
            return self.lookup(_dotted(node))  # type: ignore[arg-type]
        if isinstance(node, ast.BinOp):
            left, right = self.ev(node.left), self.ev(node.right)
            if not (_is_num(left) and _is_num(right)):
                raise ExprError(
                    f"!expr {src!r}: operator {_OP_SYMBOL[type(node.op)]} needs numbers, got {left!r} and {right!r}"
                )
            if isinstance(node.op, ast.Pow) and isinstance(left, int) and isinstance(right, int):
                if abs(right) > MAX_INT_EXPONENT:
                    raise ExprError(f"!expr {src!r}: integer exponent {right} too large")
            return _check_int(_BINOPS[type(node.op)](left, right), src)
        if isinstance(node, ast.UnaryOp):
            v = self.ev(node.operand)
            if not isinstance(node.op, ast.Not) and not _is_num(v):
                raise ExprError(f"!expr {src!r}: unary minus/plus needs a number, got {v!r}")
            return _UNARY[type(node.op)](v)
        if isinstance(node, ast.Compare):
            left = self.ev(node.left)
            result: Any = True
            for op, comp in zip(node.ops, node.comparators):
                right = self.ev(comp)
                result = _CMPOPS[type(op)](left, right)
                if not result:
                    return result
                left = right
            return result
        if isinstance(node, ast.BoolOp):
            is_and = isinstance(node.op, ast.And)
            v: Any = None
            for sub in node.values:
                v = self.ev(sub)
                if (not v) if is_and else v:
                    return v
            return v
        if isinstance(node, ast.IfExp):
            return self.ev(node.body) if self.ev(node.test) else self.ev(node.orelse)
        if isinstance(node, ast.Call):
            fn = FUNCTIONS[node.func.id][0]  # type: ignore[attr-defined]
            args = [self.ev(a) for a in node.args]
            for a in args:
                if not _is_num(a):
                    raise ExprError(f"!expr {src!r}: {node.func.id}() needs numbers, got {a!r}")  # type: ignore[attr-defined]
            return fn(*args)
        raise ExprError(f"!expr {src!r}: {type(node).__name__} is not allowed")  # unreachable after parse


def _to_python(x: Any) -> Any:
    return x.item() if isinstance(x, np.generic) else x


def safe_eval(source: str, env: Mapping[str, Any]) -> Any:
    """Evaluate ``source`` with bindings ``env`` (name or dotted reference -> value) under the rules
    of the module docstring. Returns a Python scalar (int, float, bool or str)."""
    tree = parse(source)
    try:
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            return _to_python(_Evaluator(source, env).ev(tree.body))
    except ExprError:
        raise
    except (ArithmeticError, ValueError, TypeError) as e:
        raise ExprError(f"!expr {source!r}: evaluation failed: {type(e).__name__}: {e}") from None


class ExprFunction:
    """Callable compiled from an ``!expr`` with argument names (``compile_expr``).

    Equality compares (source, arg_names, env), so specs holding these compare by value. ``to_spec``
    turns it back into an ``Expr`` marker (``Expr(source, where=env, args=arg_names)``).
    """

    __slots__ = ("source", "arg_names", "env")

    def __init__(self, source: str, arg_names: tuple[str, ...], env: Mapping[str, Any]):
        self.source = source
        self.arg_names = tuple(arg_names)
        self.env = dict(env)

    def __call__(self, *args: Any) -> Any:
        if len(args) != len(self.arg_names):
            raise TypeError(f"!expr {self.source!r} takes {len(self.arg_names)} argument(s) {self.arg_names}, got {len(args)}")
        env = {**self.env, **dict(zip(self.arg_names, args))}
        try:
            return _to_python(_Evaluator(self.source, env).ev(parse(self.source).body))
        except ExprError:
            raise
        except (ArithmeticError, ValueError, TypeError) as e:
            raise ExprError(f"!expr {self.source!r}: evaluation failed: {type(e).__name__}: {e}") from None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ExprFunction):
            return NotImplemented
        return (self.source, self.arg_names, self.env) == (other.source, other.arg_names, other.env)

    def __hash__(self) -> int:
        return hash((self.source, self.arg_names))

    def __repr__(self) -> str:
        return f"ExprFunction({self.source!r}, args={self.arg_names}, env={self.env})"


def compile_expr(source: str, arg_names: tuple[str, ...], env: Mapping[str, Any] | None = None) -> ExprFunction:
    """Return a callable f(*args) for ``source`` (validated once; for callable params).

    Args:
        source: expression text (grammar above).
        arg_names: names bound to the positional call arguments.
        env: further constant bindings (e.g. resolved ``where:`` values). Every free name of
            ``source`` must be an argument, an ``env`` key or a constant.
    """
    arg_names = tuple(arg_names)
    for a in arg_names:
        if not isinstance(a, str) or not a.isidentifier() or a.startswith("_"):
            raise ExprError(f"!expr {source!r}: invalid argument name {a!r}")
    env = dict(env or {})
    unknown = [n for n in free_names(source) if n not in arg_names and n not in env and n not in CONSTANTS]
    if unknown:
        raise ExprError(f"!expr {source!r}: unknown name(s) {unknown}; arguments are {arg_names}")
    return ExprFunction(source, arg_names, env)
