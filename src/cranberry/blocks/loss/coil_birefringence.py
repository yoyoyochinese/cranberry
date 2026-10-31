"""``coil_birefringence``: fiber-stretcher / spool coil, here the Evanescent Optics 915B (T09 §3.4, §1.7;
D35, D53). Owner: WP-B7. Offline (O) block: nothing is done in the time loop.

    L_coil  = n_turns · turn_length          (915B: 38 × 0.22 m ≈ 8.4 m [manual × product page, UNVERIFIED])
    dΔψ_b/dT = K_T · L_coil                  (K_T < 0, D51-2: −2.05 → −17.1 rad/K; sweep −1.7…−2.9
                                              rad/(m·K) → −14…−24 rad/K)
    B        = λ / L_b                       (PM15-U25D: L_b 3.0–5.0 mm, nominal 4.2 mm)
    B_bend   = bend_coeff · (r/R)²           (bend_coeff = ¼ n³(p11 − p12)(1 + ν) ≈ 0.136, UNVERIFIED)
    δ0       = B_bend / (2B)                 (rotation of the principal axes in the bends)
    x_bend   = 2 · n_turns · δ0²             (two half-circles per turn, random bend orientation,
                                              incoherent [AE]; R = 15 mm → −31 dB, 20 mm → −36 dB)

``x_bend`` goes to ``per_budget`` as an element ``{type: coil, x: …}``. ``dpsi_dT_rad_per_K`` is the
coil's share of the birefringent drift. In the time loop the drift is produced by ``segment_jones``
(K_T · L · ΔT) on the coil segment. This block only gives the number for presets and checks. The
bend radius R is not in the documents (D36, sweep 15–20 mm), and neither is the bend loss (0.01–2 dB,
UNVERIFIED). The bend loss is registered with ``losses()`` on ``segment`` (G06).

ponytail: an aligned (non-random) bend orientation adds the amplitudes, which can be several times to
10× worse. Add an ``orientation`` choice when a measured coil needs it.

Toggles: ``pol_bend_coil`` (off ⇒ x_bend = 0), ``bend_loss`` (off ⇒ element kept, disabled).

Params (all required, D34): ``n_turns``, ``turn_length_m``, ``bend_radius_m``, ``cladding_radius_m``,
``bend_coeff``, ``beat_length_m``, ``wavelength_m``, ``K_T_rad_per_m_K``, ``bend_loss_db``,
``segment``, ``applies_to``, ``source``, ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from cranberry.blocks.loss.per_budget import _noop, per_db_from_x
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["CoilBirefringence", "make_coil_birefringence"]

KIND = "coil_birefringence"
_KEYS = ("n_turns", "turn_length_m", "bend_radius_m", "cladding_radius_m", "bend_coeff", "beat_length_m",
         "wavelength_m", "K_T_rad_per_m_K", "bend_loss_db", "segment", "applies_to", "source")


class CoilBirefringence(BaseBlock):
    """The ``coil_birefringence`` block (module docstring)."""

    toggle_names = ("pol_bend_coil", "bend_loss")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, *_KEYS, block=name)
        super().__init__(KIND, name, params)
        p = {k: float(params[k]) for k in _KEYS[:9]}
        if min(p["n_turns"], p["turn_length_m"], p["bend_radius_m"], p["beat_length_m"], p["wavelength_m"]) <= 0:
            raise ValueError(f"block {name!r}: lengths and n_turns must be > 0")
        self.length_m = p["n_turns"] * p["turn_length_m"]
        self.dpsi_dT_rad_per_K = p["K_T_rad_per_m_K"] * self.length_m
        self.B = p["wavelength_m"] / p["beat_length_m"]
        self.B_bend = p["bend_coeff"] * (p["cladding_radius_m"] / p["bend_radius_m"]) ** 2
        self.delta0 = self.B_bend / (2.0 * self.B)
        self._x = 2.0 * p["n_turns"] * self.delta0 ** 2

    @property
    def x_bend(self) -> float:
        """Fast-axis power fraction from the bends (0 with pol_bend_coil off)."""
        return self._x if self._toggles["pol_bend_coil"] else 0.0

    @property
    def per_bend_db(self) -> float:
        """Bend PER = −10 log10(x/(1 − x)) [dB] (inf when off)."""
        return per_db_from_x(self.x_bend) if self.x_bend > 0.0 else math.inf

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def losses(self) -> list[LossElement]:
        p = self.params
        return [LossElement.from_il_db(f"{self.name}.bend", float(p["bend_loss_db"]), kind="fiber",
                                       applies_to=set(p["applies_to"]), source=str(p["source"]),
                                       enabled=self._toggles["bend_loss"], segment=str(p["segment"]))]


@register_block(KIND)
def make_coil_birefringence(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``coil_birefringence`` block (module docstring)."""
    return CoilBirefringence(params, name)
