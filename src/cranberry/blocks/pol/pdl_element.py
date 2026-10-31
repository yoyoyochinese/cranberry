"""``pdl_element``: polarization-dependent loss of a component (T09 §1.4, §6.2). Owner: WP-B7.

Transmissions T_max (along the axis at angle α from lab x) and T_min = q·T_max, q = 10^(−PDL/10),
with PDL_dB = 10 log10(T_max/T_min). The polarization-averaged transmission η̄ = (T_max + T_min)/2 =
10^(−IL/10) goes to the LossRegistry once through ``losses()`` (G06). The time loop gets only the
polarization dependence, normalised to η̄:

    M = R(α)·diag(√(2/(1 + q)), √(2q/(1 + q)))·R(−α),  R(α) = [[cos α, −sin α], [sin α, cos α]]

so ⟨|M e|²⟩ over the Poincaré sphere is 1 and |M e|² ranges over 2/(1+q) … 2q/(1+q)
(T_mean·(1 ± (T_max − T_min)/(T_max + T_min)), T09 §1.4). The block writes M to the J slot of its own
jones segment (topology ``pol_params: {J_slot: <out>}``, a zero-length component segment), so
``beam_assembler`` puts |M e|² into ``beam.<b>.power`` and arg into the phase. M is symmetric, so the
backward beam's Mᵀ = M (PDL is reciprocal). The absolute Jones form is ``jones()`` = √η̄·M.

Toggle ``pol_pdl``: off ⇒ M = 1 (the averaged IL stays registered: it is not the PDL).

Params (all required, D34): ``out`` (J slot), ``pdl_db``, ``axis_rad`` (runtime), ``il_db``
(polarization-averaged IL), ``segment`` (loss location), ``applies_to``, ``source``, ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.core.units import il_db_to_eta
from cranberry.offline.loss_registry import LossElement

__all__ = ["PdlElement", "make_pdl_element", "pdl_matrix"]

KIND = "pdl_element"


def _step(Z, P, o, p_pdl, p_ax, on):  # noqa: N803
    if not on:
        Z[o], Z[o + 1], Z[o + 2], Z[o + 3] = 1.0, 0.0, 0.0, 1.0
        return
    q = 10.0 ** (-P[p_pdl] / 10.0)
    tx, ty = math.sqrt(2.0 / (1.0 + q)), math.sqrt(2.0 * q / (1.0 + q))
    c, s = math.cos(P[p_ax]), math.sin(P[p_ax])
    Z[o] = tx * c * c + ty * s * s
    Z[o + 1] = (tx - ty) * c * s
    Z[o + 2] = (tx - ty) * c * s
    Z[o + 3] = tx * s * s + ty * c * c


_STEP_NB = numba.njit(_step)


def pdl_matrix(pdl_db: float, axis_rad: float) -> np.ndarray:
    """M normalised to the polarization-averaged transmission (module docstring)."""
    z = np.zeros(4, dtype=complex)
    _step(z, np.array([pdl_db, axis_rad]), 0, 0, 1, True)
    return z.reshape(2, 2)


class PdlElement(BaseBlock):
    """The ``pdl_element`` block (module docstring)."""

    toggle_names = ("pol_pdl",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "out", "pdl_db", "axis_rad", "il_db", "segment", "applies_to", "source", block=name)
        super().__init__(KIND, name, params, writes=(params["out"],), runtime_params=("pdl_db", "axis_rad"))
        q = 10.0 ** (-float(params["pdl_db"]) / 10.0)
        if float(params["pdl_db"]) < 0 or 2.0 * il_db_to_eta(float(params["il_db"])) / (1.0 + q) > 1.0:
            raise ValueError(f"block {name!r}: need pdl_db ≥ 0 and T_max = 2η̄/(1 + q) ≤ 1 (il_db too small for pdl_db)")

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        layout.add(self.params["out"], (2, 2), np.complex128, "1", init=(1.0, 0.0, 0.0, 1.0))

    def bind(self, layout) -> None:
        o = layout.slot(self.params["out"]).offset
        p_pdl, p_ax = layout.param(self.param_name("pdl_db")).offset, layout.param(self.param_name("axis_rad")).offset
        on = self.toggles["pol_pdl"]
        self._args = (o, p_pdl, p_ax, on)
        nb = _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(Z, P, o, p_pdl, p_ax, on)

        self.kernel = kernel

    def __call__(self, state, params):
        _step(state.Z, params.P, *self._args)
        return state

    def jones(self) -> np.ndarray:
        """Absolute Jones form √η̄·M at the spec values (M = 1 with pol_pdl off)."""
        p = self.params
        m = pdl_matrix(float(p["pdl_db"]), float(p["axis_rad"])) if self._toggles["pol_pdl"] else np.eye(2)
        return math.sqrt(il_db_to_eta(float(p["il_db"]))) * m

    def losses(self) -> list[LossElement]:
        p = self.params
        return [LossElement.from_il_db(f"{self.name}.il", float(p["il_db"]), kind="component",
                                       applies_to=set(p["applies_to"]), source=str(p["source"]), enabled=True,
                                       segment=str(p["segment"]))]


@register_block(KIND)
def make_pdl_element(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``pdl_element`` block (module docstring)."""
    return PdlElement(params, name)
