"""``dispersion_comp``: per-port dispersion compensation and delay trim (D44, T03 §5.4, §9-1). Owner: WP-B11.

Offline (O) block. Port k's compensator (DCF, chirped FBG, delay stage) adds GDD, TOD and a delay to
the squeezed path, so the residuals that ``mode_overlap_table`` looks up are (D51-6: δL_k sets the
dispersion, δτ_k the timing)::

    φ2,res = β2 δL_k + gdd_k,   φ3,res = β3 δL_k + tod_k,   δτ_res = δτ_k + delay_k

and its insertion loss is one registered ``component`` element on the port's segment (G06).

Params (D34): ``ports`` {port: {gdd_s2, tod_s3, delay_s, segment, insertion_eta}}, ``applies_to``,
``source_tag``, ``toggles`` {enabled}: off ⇒ no compensation (residual = input) and the elements are
registered disabled.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["DispersionComp", "make_dispersion_comp"]

KIND = "dispersion_comp"


class DispersionComp(BaseBlock):
    """The ``dispersion_comp`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "ports", "applies_to", "source_tag", block=name)
        for k, v in params["ports"].items():
            require(v, "gdd_s2", "tod_s3", "delay_s", "segment", "insertion_eta", block=f"{name}.ports.{k}")
            if not 0.0 <= float(v["insertion_eta"]) <= 1.0:
                raise ValueError(f"block {name!r}: ports.{k}.insertion_eta must be in [0, 1]")
        super().__init__(KIND, name, params)

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def residual(self, port: str, phi2_s2: float, phi3_s3: float, dt_s: float) -> tuple[float, float, float]:
        """(φ2, φ3, δτ) after the compensator of ``port``."""
        if not self._toggles["enabled"]:
            return phi2_s2, phi3_s3, dt_s
        c = self.params["ports"][port]
        return phi2_s2 + float(c["gdd_s2"]), phi3_s3 + float(c["tod_s3"]), dt_s + float(c["delay_s"])

    def losses(self) -> list[LossElement]:
        return [LossElement(id=f"{self.name}.{k}", eta=float(v["insertion_eta"]), kind="component",
                            applies_to=set(self.params["applies_to"]), source=str(self.params["source_tag"]),
                            enabled=self._toggles["enabled"], segment=str(v["segment"]))
                for k, v in self.params["ports"].items()]


@register_block(KIND)
def make_dispersion_comp(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``dispersion_comp`` block (module docstring)."""
    return DispersionComp(params, name)
