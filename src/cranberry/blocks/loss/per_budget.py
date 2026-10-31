"""``per_budget``: PER accumulation over connectors, splices, couplers, polarizer, coil and fiber
(T09 §1.3, §3.5; D51-8). Owner: WP-B7. Offline (O) block: nothing is done in the time loop.

Every element is reduced to its fast-axis power fraction x_k:

    joint at misalignment θ (connector, splice, component port):  x = sin²θ,
        PER = −10 log10(tan²θ)       (3° → 25.61 dB, T09 test 3)
    element given by PER [dB] (coupler PER ≥ 25 dB, polarizer ER):  x = q/(1 + q), q = 10^(−PER/10)
    distributed fiber crosstalk (h-parameter):  x = tanh(h·L)   (h·L = 1e−3 → −30.00 dB, T09 test 4)
    element given directly as x (e.g. ``coil_birefringence.x_bend``)

The total uses the random-phase power sum x_tot = Σ x_k, PER_tot = −10 log10(x_tot/(1 − x_tot)). The
time loop does not use x_tot. It uses the coherent per-joint angles (``theta_rad(id)`` →
``segment_jones`` θ_in / θ_out), where the phases come from Δψ_b(t) (T09 §1.3). If a polarizing
element (polarizer, PBS, GBS input) follows, the fast power is lost: ``loss`` registers
η = 1 − x_tot (G06). ``loss: null`` registers nothing, e.g. Ref D, where the PD does not resolve
polarization.

Toggle ``pol_crosstalk``: off ⇒ x_tot = 0 and the loss element stays in the record disabled.

Params (all required, D34)::

    elements: [{id, type: connector|splice|coupler|polarizer|coil|fiber,
                one of theta_deg | per_db | x | {h_per_m, L_m}, source}]
    loss:     null | {segment, applies_to, source}
    toggles:  {pol_crosstalk}
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["PerBudget", "make_per_budget", "per_db_from_angle", "angle_from_per_db", "x_from_per_db",
           "per_db_from_x", "ELEMENT_TYPES"]

KIND = "per_budget"
ELEMENT_TYPES: tuple[str, ...] = ("connector", "splice", "coupler", "polarizer", "coil", "fiber")
_VALUE_KEYS = ("theta_deg", "per_db", "x", "h_per_m")


def per_db_from_angle(theta_rad: float) -> float:
    """PER of a joint with slow-axis misalignment θ: −10 log10(tan²θ) [dB]."""
    return -10.0 * math.log10(math.tan(theta_rad) ** 2)


def angle_from_per_db(per_db: float) -> float:
    """Misalignment angle θ = arctan(10^(−PER/20)) [rad] of a joint with the given PER."""
    return math.atan(10.0 ** (-per_db / 20.0))


def x_from_per_db(per_db: float) -> float:
    """Fast-axis power fraction x = q/(1 + q), q = 10^(−PER/10)."""
    q = 10.0 ** (-per_db / 10.0)
    return q / (1.0 + q)


def per_db_from_x(x: float) -> float:
    """PER = −10 log10(x/(1 − x)) [dB]."""
    return -10.0 * math.log10(x / (1.0 - x))


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class PerBudget(BaseBlock):
    """The ``per_budget`` block (module docstring)."""

    toggle_names = ("pol_crosstalk",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "elements", "loss", "toggles", block=name)
        super().__init__(KIND, name, params)
        self.x: dict[str, float] = {}
        for i, e in enumerate(params["elements"]):
            where = f"block {name!r}, elements[{i}]"
            require(e, "id", "type", "source", block=where)
            if e["type"] not in ELEMENT_TYPES:
                raise ValueError(f"{where}: type must be one of {ELEMENT_TYPES}")
            given = [k for k in _VALUE_KEYS if k in e]
            if len(given) != 1:
                raise ValueError(f"{where}: give exactly one of theta_deg, per_db, x, (h_per_m, L_m)")
            if e["id"] in self.x:
                raise ValueError(f"{where}: duplicate id {e['id']!r}")
            k = given[0]
            if k == "theta_deg":
                x = math.sin(math.radians(float(e["theta_deg"]))) ** 2
            elif k == "per_db":
                x = x_from_per_db(float(e["per_db"]))
            elif k == "x":
                x = float(e["x"])
            else:
                require(e, "L_m", block=where)
                x = math.tanh(float(e["h_per_m"]) * float(e["L_m"]))
            if not 0.0 <= x < 1.0:
                raise ValueError(f"{where}: fast-axis fraction must be in [0, 1), got {x}")
            self.x[e["id"]] = x
        if params["loss"] is not None:
            require(params["loss"], "segment", "applies_to", "source", block=f"{name}.loss")

    @property
    def x_total(self) -> float:
        """Σ x_k (random-phase power sum); 0 with pol_crosstalk off."""
        return math.fsum(self.x.values()) if self._toggles["pol_crosstalk"] else 0.0

    @property
    def per_total_db(self) -> float:
        """Total PER = −10 log10(x_tot/(1 − x_tot)) [dB] (inf when x_tot = 0)."""
        x = self.x_total
        return math.inf if x == 0.0 else per_db_from_x(x)

    def theta_rad(self, element_id: str) -> float:
        """Equivalent joint angle θ = arcsin √x of one element (for ``segment_jones``)."""
        return math.asin(math.sqrt(self.x[element_id])) if self._toggles["pol_crosstalk"] else 0.0

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def losses(self) -> list[LossElement]:
        loss = self.params["loss"]
        if loss is None:
            return []
        x = math.fsum(self.x.values())  # recorded value; disabled (η = 1) with the toggle off
        return [LossElement(id=f"{self.name}.per", eta=1.0 - x, kind="component",
                            applies_to=set(loss["applies_to"]), source=str(loss["source"]),
                            enabled=self._toggles["pol_crosstalk"], segment=str(loss["segment"]))]


@register_block(KIND)
def make_per_budget(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``per_budget`` block (module docstring)."""
    return PerBudget(params, name)
