"""``faraday_nonreciprocal``: Faraday rotation, the non-reciprocal term that breaks J_back = Jᵀ
(T09 §1.6, §1.8; optional, default off). Owner: WP-B7.

A field B_∥ along the fiber rotates the linear polarization by ρ = V_F·B_∥·L (circular states get
phases ∓ρ). In the fixed lab basis (x, y) the rotation sense does not follow the propagation direction:

    N = R(ρ) = [[cos ρ, −sin ρ], [sin ρ, cos ρ]]   for both d = +1 and d = −1

whereas a reciprocal rotation would give Nᵀ = R(−ρ) backward. The block writes N to the segment's
``J_nr_slot`` (topology jones ``pol_params``, F6 addition): ``beam_assembler`` applies N·U forward and
Uᵀ·N backward, so Jᵀ is applied to the reciprocal part U only. A same-polarization counter-propagating
probe then differs from the squeezed beam by up to 2ρ: e.g. x → QWP → circular → N → QWP⁻¹ gives
arg t_s − arg t_p = −2ρ. T09 §1.8 [AE]: V_F ≈ 0.6 rad/(T·m) (silica, 1550 nm, UNVERIFIED), earth field
50 µT all parallel over 63 m ⇒ ρ ≈ 1.9×10⁻³ rad (worst case; a spool largely cancels, PM birefringence
suppresses it).

ponytail: lumped at the segment exit, with a static B_∥ (runtime param). A distributed Faraday term
inside a birefringent fiber does not commute with U; split the segment if that matters, and add a
B-field driver slot when a magnetic noise source exists.

Toggle ``pol_faraday``: off ⇒ N = 1.

Params (all required, D34): ``out`` (J_nr slot), ``verdet_rad_per_T_m``, ``B_par_T``, ``L_m``
(runtime), ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block

__all__ = ["FaradayNonreciprocal", "make_faraday_nonreciprocal"]

KIND = "faraday_nonreciprocal"
_RT = ("verdet_rad_per_T_m", "B_par_T", "L_m")


def _step(Z, P, o, pv, pb, pl, on):  # noqa: N803
    rho = P[pv] * P[pb] * P[pl] if on else 0.0
    c, s = math.cos(rho), math.sin(rho)
    Z[o], Z[o + 1], Z[o + 2], Z[o + 3] = c, -s, s, c


_STEP_NB = numba.njit(_step)


class FaradayNonreciprocal(BaseBlock):
    """The ``faraday_nonreciprocal`` block (module docstring)."""

    toggle_names = ("pol_faraday",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "out", block=name)
        super().__init__(KIND, name, params, writes=(params["out"],), runtime_params=_RT)

    @property
    def rho_rad(self) -> float:
        """Rotation ρ = V_F·B_∥·L at the spec values (0 with pol_faraday off)."""
        p = self.params
        return float(p["verdet_rad_per_T_m"]) * float(p["B_par_T"]) * float(p["L_m"]) if self._toggles["pol_faraday"] else 0.0

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        layout.add(self.params["out"], (2, 2), np.complex128, "1", init=(1.0, 0.0, 0.0, 1.0))

    def bind(self, layout) -> None:
        o = layout.slot(self.params["out"]).offset
        pv, pb, pl = (layout.param(self.param_name(k)).offset for k in _RT)
        on = self.toggles["pol_faraday"]
        self._args = (o, pv, pb, pl, on)
        nb = _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(Z, P, o, pv, pb, pl, on)

        self.kernel = kernel

    def __call__(self, state, params):
        _step(state.Z, params.P, *self._args)
        return state


@register_block(KIND)
def make_faraday_nonreciprocal(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``faraday_nonreciprocal`` block (module docstring)."""
    return FaradayNonreciprocal(params, name)
