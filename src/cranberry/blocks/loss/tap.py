"""``tap``: conventional tap after the squeezer (Ref A Fig. 1(a); T12 §6, T14). Owner: WP-B1.

O part: the squeezed (and probe) path loses the tapped fraction r_t, registered as one
``LossElement`` η = 1 − r_t on ``segment`` (G06; never multiplied in here).
T part: writes the tapped probe power P_tap = r_t · P_pr [W] to ``out`` every step, for the error
signal generators (T14). P_pr (``p_probe_w``, probe power leaving the OPA) is a runtime parameter;
r_t is structural (it is also the registered loss).

Toggle ``enabled`` (TAP_ENABLE, T12 §8.1): off ⇒ the element stays in the record disabled (η = 1)
and P_tap = 0.

Params (all required, D34): ``tap_ratio`` r_t ∈ [0, 1), ``p_probe_w`` [W], ``out`` (slot name),
``segment``, ``applies_to`` (subset of LOSS_TARGETS), ``source`` (provenance tag), ``toggles``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["Tap", "make_tap"]

KIND = "tap"


class Tap(BaseBlock):
    """The ``tap`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "tap_ratio", "p_probe_w", "out", "segment", "applies_to", "source", block=name)
        if not 0.0 <= float(params["tap_ratio"]) < 1.0:
            raise ValueError(f"block {name!r}: tap_ratio must be in [0, 1), got {params['tap_ratio']!r}")
        super().__init__(KIND, name, params, writes=(params["out"],), runtime_params=("p_probe_w",))
        self.slot_units = {params["out"]: "W"}

    @property
    def _r(self) -> float:
        return float(self.params["tap_ratio"]) if self._toggles["enabled"] else 0.0

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        layout.add(self.params["out"], (), np.float64, "W")

    def bind(self, layout) -> None:
        self._i_out = layout.slot(self.params["out"]).offset
        self._i_p = layout.param(self.param_name("p_probe_w")).offset
        i_out, i_p, r = self._i_out, self._i_p, self._r

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803, E741
            S[i_out] = r * P[i_p]

        self.kernel = kernel

    def __call__(self, state, params):
        state.S[self._i_out] = self._r * params.P[self._i_p]
        return state

    def losses(self) -> list[LossElement]:
        return [LossElement(id=f"{self.name}.tap", eta=1.0 - float(self.params["tap_ratio"]), kind="component",
                            applies_to=set(self.params["applies_to"]), source=str(self.params["source"]),
                            enabled=self._toggles["enabled"], segment=str(self.params["segment"]))]


@register_block(KIND)
def make_tap(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``tap`` block (module docstring)."""
    return Tap(params, name)
