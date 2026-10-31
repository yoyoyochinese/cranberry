"""``coil_birefringence``: fiber-stretcher / spool coil, here the Evanescent Optics 915B (T09 §3.4, §1.7;
D35, D53). Owner: WP-B7. Offline (O) block: nothing is done in the time loop.

    L_coil  = n_turns · turn_length          (915B: 38 × 0.22 m ≈ 8.4 m [manual × product page, UNVERIFIED])
    dΔψ_b/dT = K_T · L_coil                  (K_T < 0, D51-2: −2.05 → −17.1 rad/K; sweep −1.7…−2.9
                                              rad/(m·K) → −14…−24 rad/K)
    B        = λ / L_b                       (PM15-U25D: L_b 3.0–5.0 mm, nominal 4.2 mm)
    B_bend   = bend_coeff · (r/R)²           (bend_coeff = ¼ n³(p11 − p12)(1 + ν) ≈ 0.136, UNVERIFIED)
    δ0       = B_bend / (2B)                 (rotation of the principal axes in the bends)
    x_bend   = 2 · n_turns · δ0²             (orientation: random — two half-circles per turn, random
                                              bend orientation, incoherent [AE]; R = 15 mm → −31 dB,
                                              20 mm → −36 dB)

orientation: aligned (fiber wound without twist, every bend at the same angle α to the slow axis): the
axis steps ε = δ0·sin 2α at the bend boundaries add coherently with the birefringent phase β z between
them (β = 2π/L_b + K_T·ΔT). Per turn: half-circle πR, straight s = L_turn/2 − πR, half-circle, straight,
boundaries at z = 0 (+ε), πR (−ε), πR + s (+ε), 2πR + s (−ε):

    x_bend(ΔT) = |ε·(1 − e^{−iβπR})·(1 + e^{−iβ(πR + s)})·Σ_{m<n} e^{−iβ m L_turn}|²

It depends on mm-level geometry through βπR ~ 70 rad and rotates with temperature (≈ K_T·L_turn
= 0.45 rad/K per turn, T09 §3.4). Mean over the phases 4·n·ε² (2× the random value at α = 45°); peak
(all phases matched) 16·n²·ε². ``x_bend`` is the value at ΔT = 0, ``x_bend_at(dT)`` the sweep.

``x_bend`` goes to ``per_budget`` as an element ``{type: coil, x: …}``. ``dpsi_dT_rad_per_K`` is the
coil's share of the birefringent drift. In the time loop the drift is produced by ``segment_jones``
(K_T · L · ΔT) on the coil segment. This block only gives the number for presets and checks. The
bend radius R is not in the documents (D36, sweep 15–20 mm), and neither is the bend loss (0.01–2 dB,
UNVERIFIED). The bend loss is registered with ``losses()`` on ``segment`` (G06).

Toggles: ``pol_bend_coil`` (off ⇒ x_bend = 0), ``bend_loss`` (off ⇒ element kept, disabled).

Params (all required, D34): ``n_turns``, ``turn_length_m``, ``bend_radius_m``, ``cladding_radius_m``,
``bend_coeff``, ``beat_length_m``, ``wavelength_m``, ``K_T_rad_per_m_K``, ``bend_loss_db``,
``orientation`` (random | aligned), ``bend_angle_rad`` (α, aligned only), ``segment``, ``applies_to``,
``source``, ``toggles``.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping
from typing import Any

from cranberry.blocks.loss.per_budget import _noop, per_db_from_x
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["CoilBirefringence", "make_coil_birefringence", "ORIENTATIONS"]

KIND = "coil_birefringence"
_KEYS = ("n_turns", "turn_length_m", "bend_radius_m", "cladding_radius_m", "bend_coeff", "beat_length_m",
         "wavelength_m", "K_T_rad_per_m_K", "bend_loss_db", "segment", "applies_to", "source", "orientation")
ORIENTATIONS: tuple[str, ...] = ("random", "aligned")


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
        self.orientation = params["orientation"]
        if self.orientation not in ORIENTATIONS:
            raise ValueError(f"block {name!r}: orientation must be one of {ORIENTATIONS}")
        self._p = p
        if self.orientation == "aligned":
            require(params, "bend_angle_rad", block=name)
            self.straight_m = p["turn_length_m"] / 2.0 - math.pi * p["bend_radius_m"]
            if self.straight_m < 0:
                raise ValueError(f"block {name!r}: turn_length_m < 2πR (no room for two half-circles)")
            self.eps = self.delta0 * math.sin(2.0 * float(params["bend_angle_rad"]))
            self._x = self._aligned(0.0)
        else:
            self._x = 2.0 * p["n_turns"] * self.delta0 ** 2

    def _aligned(self, dT: float) -> float:
        p = self._p
        beta = 2.0 * math.pi / p["beat_length_m"] + p["K_T_rad_per_m_K"] * dT
        a1 = (1 - cmath.exp(-1j * beta * math.pi * p["bend_radius_m"])) * (
            1 + cmath.exp(-1j * beta * (math.pi * p["bend_radius_m"] + self.straight_m)))
        ph = beta * p["turn_length_m"]
        turns = sum(cmath.exp(-1j * ph * m) for m in range(int(p["n_turns"])))
        return abs(self.eps * a1 * turns) ** 2

    def x_bend_at(self, dT: float) -> float:
        """Bend fast-axis fraction at coil temperature change dT [K] (random: constant; 0 when off)."""
        if not self._toggles["pol_bend_coil"]:
            return 0.0
        return self._aligned(dT) if self.orientation == "aligned" else self._x

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
