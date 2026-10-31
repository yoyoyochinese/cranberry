"""Registries: block factories (``@register_block``) and named callables (``@register_fn``).

Blocks (PLAN §2.4): ``@register_block("pid")`` on ``make_pid(params, *, name=None) -> Block`` makes
YAML ``{block: pid}`` and ``make_pid(...)`` produce the same object. Factories register at import
time; :func:`ensure_loaded` imports every module that defines factories (framework modules and all
``cranberry.blocks.<family>`` subpackages) so the registry is complete before a spec is built.

Functions (G12): ``@register_fn("refB_loss_vs_x")`` registers a Python callable usable from YAML as
``!fn {name: refB_loss_vs_x, args: {...}}``; ``to_spec`` can only serialise registered functions.

The kinds ``seq``, ``parallel`` and ``stack`` are reserved for the pipeline combinators.
"""

from __future__ import annotations

import importlib
import pkgutil
from collections.abc import Callable, Mapping
from typing import TYPE_CHECKING, Any

from cranberry.core.block import KIND_RE

if TYPE_CHECKING:
    from cranberry.core.block import Block, BlockSpec

__all__ = [
    "RESERVED_KINDS",
    "BlockFactory",
    "register_block",
    "get_factory",
    "registered_kinds",
    "make_block",
    "register_fn",
    "get_fn",
    "fn_name",
    "registered_fns",
    "ensure_loaded",
]

#: Factory signature: ``make_<kind>(params, *, name=None) -> Block``.
BlockFactory = Callable[..., "Block"]

RESERVED_KINDS: frozenset[str] = frozenset({"seq", "parallel", "stack"})

#: Framework modules that register block factories (besides cranberry.blocks.*).
_FRAMEWORK_BLOCK_MODULES: tuple[str, ...] = (
    "cranberry.core.delay",
    "cranberry.path.segments",
    "cranberry.path.beams",
)

_BLOCKS: dict[str, BlockFactory] = {}
_FNS: dict[str, Callable[..., Any]] = {}


def _qualname(fn: Callable[..., Any]) -> str:
    return f"{getattr(fn, '__module__', '?')}.{getattr(fn, '__qualname__', repr(fn))}"


def register_block(kind: str) -> Callable[[BlockFactory], BlockFactory]:
    """Decorator registering a block factory under ``kind`` (unique; re-import of the same function is allowed)."""
    if not KIND_RE.match(kind):
        raise ValueError(f"invalid block kind {kind!r}")
    if kind in RESERVED_KINDS:
        raise ValueError(f"block kind {kind!r} is reserved for pipeline combinators")

    def deco(factory: BlockFactory) -> BlockFactory:
        old = _BLOCKS.get(kind)
        if old is not None and _qualname(old) != _qualname(factory):
            raise ValueError(f"block kind {kind!r} already registered by {_qualname(old)}")
        _BLOCKS[kind] = factory
        return factory

    return deco


def get_factory(kind: str) -> BlockFactory:
    try:
        return _BLOCKS[kind]
    except KeyError:
        raise KeyError(f"unknown block kind {kind!r}; registered: {sorted(_BLOCKS)}") from None


def registered_kinds() -> tuple[str, ...]:
    return tuple(sorted(_BLOCKS))


def make_block(spec: BlockSpec) -> Block:
    """Instantiate a leaf block from its spec: ``factory(spec.params, name=spec.name)`` + rate_div.

    The factory receives ``rate_div`` through the spec only; the returned block's ``spec`` must
    equal the input spec (round-trip contract, VT-F1).
    """
    factory = get_factory(spec.kind)
    block = factory(dict(spec.params), name=spec.name)
    if spec.rate_div != 1:
        from cranberry.core.block import BaseBlock

        if not isinstance(block, BaseBlock):
            raise TypeError(f"block {spec.name!r}: rate_div != 1 needs a BaseBlock-derived block")
        block.set_rate_div(spec.rate_div)
    return block


def register_fn(name: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator registering a callable parameter function for YAML ``!fn`` (G12)."""
    if not name or not isinstance(name, str):
        raise ValueError("function name must be a non-empty str")

    def deco(fn: Callable[..., Any]) -> Callable[..., Any]:
        old = _FNS.get(name)
        if old is not None and _qualname(old) != _qualname(fn):
            raise ValueError(f"function {name!r} already registered by {_qualname(old)}")
        _FNS[name] = fn
        return fn

    return deco


def get_fn(name: str) -> Callable[..., Any]:
    try:
        return _FNS[name]
    except KeyError:
        raise KeyError(f"unknown registered function {name!r}") from None


def fn_name(fn: Callable[..., Any]) -> str | None:
    """Registered name of ``fn`` or None (unregistered lambdas cannot round-trip to YAML, G12)."""
    for k, v in _FNS.items():
        if v is fn:
            return k
    return None


def registered_fns() -> Mapping[str, Callable[..., Any]]:
    return dict(_FNS)


def ensure_loaded() -> None:
    """Import all modules that register factories (framework modules + every cranberry.blocks subpackage)."""
    for mod in _FRAMEWORK_BLOCK_MODULES:
        importlib.import_module(mod)
    import cranberry.blocks as blocks_pkg

    for info in pkgutil.walk_packages(blocks_pkg.__path__, prefix="cranberry.blocks."):
        importlib.import_module(info.name)
