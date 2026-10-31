"""Restricted expression evaluator for ``!expr`` (G12). Owner: WP-F2.

Allowed: numeric literals, names bound via ``where:`` / dotted parameter paths, + − * / ** unary −,
comparisons, the ternary ``a if c else b``, and calls of whitelisted elementary functions from
``math`` / ``numpy`` (sqrt, exp, log, log10, sin, cos, tan, arctan, ...). Forbidden: attribute access
on arbitrary objects, arbitrary calls, subscripts on non-arrays, comprehensions, lambdas, imports.
Implemented on ``ast``; never uses ``eval``/``exec`` (YAML is untrusted input).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

__all__ = ["safe_eval", "compile_expr"]


def safe_eval(source: str, env: Mapping[str, Any]) -> Any:
    """Evaluate ``source`` with bindings ``env`` under the rules above. WP-F2."""
    raise NotImplementedError("config.expr.safe_eval is implemented by WP-F2")


def compile_expr(source: str, arg_names: tuple[str, ...]):
    """Return a Python callable f(*args) for ``source`` (validated once; for callable params). WP-F2."""
    raise NotImplementedError("config.expr.compile_expr is implemented by WP-F2")
