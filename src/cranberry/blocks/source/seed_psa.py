"""``seed_psa``: a weak coherent seed through the degenerate phase-sensitive amplifier of the squeezer (time loop, SRC;
Phase 6 fix of reviews/phase6_wp6c_p.md C-1; conv §3–§4). Owner: WP6-C fix.

A seed injected BEFORE the squeezer (lock/dual_band_slow, ``dbs_seed_inject: before_squeezer``) is not a θ-unit copy of
the squeezed light: the amplifier multiplies its quadrature along θ_amp = θ_sq − π/2 by e^{r} and the orthogonal one by
e^{−r} (S† â S = â cosh r − â† e^{iφ_ξ} sinh r, conv §3). With δ = φ_in − θ_amp the seed's input angle to the amplified
axis, the output field is α ∝ e^{r} cos δ + i e^{−r} sin δ::

    δ        = delta_rad + φ_seed − θ_sq                 (both beams at the same observer: the common downstream
                                                          path cancels, δ is the seed-vs-pump angle at the crystal)
    out      = θ_sq − π/2 + atan2(e^{−r} sin δ, e^{r} cos δ)
    out_i    = p_in_rel · (e^{2r} cos² δ + e^{−2r} sin² δ)  (output intensity, relative units)

``delta_rad`` is the preparation angle against θ_amp (0: on the amplified axis, π/2: on the squeezed axis);
``in_seed`` is the seed's own phase (laser through the seed delivery path, e.g. beam ``dsin<k>``), ``in_sqz`` the
squeezing angle θ_sq of the same port (beam ``sqz<k>``), so the seed's own input phase noise enters out with the gain
dψ/dδ = e^{−2r}/(cos²δ + e^{−4r} sin²δ): e^{2r} at δ = π/2 (36.6 at r = 1.8), e^{−2r} at δ = 0, 2e^{−2r}/(1 + e^{−4r})
at ±π/4. Two seeds prepared at (π/2, 0) keep the π/2 offset but their intensity ratio is e^{−4r}; at (π/4, −π/4) the
intensities are equal and the offset collapses to 2 atan(e^{−2r}).

Toggle ``amplification`` (off: r = 0, out = φ_seed + delta_rad − π/2, out_i = p_in_rel). Params (D34): ``in_seed``,
``in_sqz``, ``out``, ``out_i`` (slots); runtime ``r`` (≥ 0), ``delta_rad``, ``p_in_rel`` (≥ 0).

ponytail: single-mode amplifier gain e^{±r} on the seed's mode (the optimal mode of a pulsed OPO: r of that mode); the
pump-depletion, the OPO cavity's frequency response and the seed-pump relative-phase lock are not modelled.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block

__all__ = ["SeedPSA", "make_seed_psa", "psa_phase", "psa_gain"]

KIND = "seed_psa"


def psa_phase(delta: float, r: float) -> float:
    """Output angle to the amplified axis of a seed at input angle ``delta`` (module docstring)."""
    return math.atan2(math.exp(-r) * math.sin(delta), math.exp(r) * math.cos(delta))


def psa_gain(delta: float, r: float) -> float:
    """Intensity gain e^{2r} cos²δ + e^{−2r} sin²δ."""
    return math.exp(2.0 * r) * math.cos(delta) ** 2 + math.exp(-2.0 * r) * math.sin(delta) ** 2


def _step(S, I, P, i_seed, i_sqz, i_out, i_i, p_r, p_d, p_w, amp):  # noqa: N803, E741
    r = P[p_r] if amp else 0.0
    sq = S[i_sqz]
    d = P[p_d] + S[i_seed] - sq
    c, s = math.cos(d), math.sin(d)
    S[i_out] = sq - 0.5 * math.pi + math.atan2(math.exp(-r) * s, math.exp(r) * c)
    S[i_i] = P[p_w] * (math.exp(2.0 * r) * c * c + math.exp(-2.0 * r) * s * s)


class SeedPSA(KitBlock):
    """The ``seed_psa`` block (module docstring)."""

    KIND = KIND
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("amplification",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "in_seed", "in_sqz", "out", "out_i", "r", "delta_rad", "p_in_rel", block=name)
        if not float(params["r"]) >= 0.0 or not float(params["p_in_rel"]) >= 0.0:
            raise ValueError(f"block {name!r}: need r >= 0 and p_in_rel >= 0")
        super().__init__(params, name, reads=(params["in_seed"], params["in_sqz"]),
                         writes=(params["out"], params["out_i"]), runtime=("r", "delta_rad", "p_in_rel"))
        self.slot_units = {params["out"]: "rad", params["out_i"]: "1"}

    def slots(self):
        return [(self.params["out"], (), np.float64, "rad", 0, 0.0), (self.params["out_i"], (), np.float64, "1", 0, 0.0)]

    def make_args(self, layout):
        pr = self.params
        return (self.o(layout, pr["in_seed"]), self.o(layout, pr["in_sqz"]), self.o(layout, pr["out"]),
                self.o(layout, pr["out_i"]), self.p(layout, "r"), self.p(layout, "delta_rad"), self.p(layout, "p_in_rel"),
                self.toggles["amplification"])


@register_block(KIND)
def make_seed_psa(params: Mapping[str, Any], *, name: str | None = None) -> SeedPSA:
    """Factory of the ``seed_psa`` block (module docstring)."""
    return SeedPSA(params, name)
