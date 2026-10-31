"""``loss_budget``: loss elements given directly in the config (PLAN §4.3 LOSS, G06). Owner: WP-I1 addition.

Offline (O) block: registers one ``LossElement`` per entry of ``elements`` with ``losses()``; the time
loop never steps it (no reads / writes, ``__call__`` is the identity). For losses that no physics
block owns, e.g. Ref A §4.2 "optical components (mirrors, lenses) 1 %" on the squeezed path. Blocks
that own a loss (OPA waveguide, EOM, tap, detector, mode overlap) register it themselves; listing it
here as well would count it twice (the LossRegistry only rejects repeated *overlap* kinds, G06).

Params (all required, D34)::

    elements: [{id, eta, kind, applies_to, source, segment, observer}, ...]
        eta        power transmission in [0, 1] (linear; a dB key would be renamed on load, G16-i)
        kind       one of LOSS_KINDS
        applies_to subset of LOSS_TARGETS (gbs / bhd / lo / probe)
        segment / observer   where it sits (exactly one non-null)
        source     provenance tag
    toggles: {enabled}   off ⇒ every element is recorded disabled (η treated as 1)

Element ids are prefixed with the block name (``<name>.<id>``).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["LossBudget", "make_loss_budget"]

KIND = "loss_budget"
_KEYS = ("id", "eta", "kind", "applies_to", "source", "segment", "observer")


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class LossBudget(BaseBlock):
    """The ``loss_budget`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "elements", "toggles", block=name)
        if not isinstance(params["elements"], (list, tuple)) or not params["elements"]:
            raise ValueError(f"block {name!r}: elements must be a non-empty list")
        ids = []
        for e in params["elements"]:
            require(e, *_KEYS, block=f"{name}.elements")
            if (e["segment"] is None) == (e["observer"] is None):
                raise ValueError(f"block {name!r}: element {e['id']!r} needs exactly one of segment / observer")
            ids.append(e["id"])
        if len(set(ids)) != len(ids):
            raise ValueError(f"block {name!r}: duplicate element ids {ids}")
        super().__init__(KIND, name, params)
        self._elements = [self._element(e) for e in params["elements"]]  # validates eta / kind / targets now

    def _element(self, e: Mapping[str, Any]) -> LossElement:
        return LossElement(id=f"{self.name}.{e['id']}", eta=float(e["eta"]), kind=str(e["kind"]),
                           applies_to=frozenset(e["applies_to"]), source=str(e["source"]),
                           enabled=self._toggles["enabled"],
                           segment=None if e["segment"] is None else str(e["segment"]),
                           observer=None if e["observer"] is None else str(e["observer"]))

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def losses(self) -> list[LossElement]:
        return list(self._elements)


@register_block(KIND)
def make_loss_budget(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``loss_budget`` block (module docstring)."""
    return LossBudget(params, name)
