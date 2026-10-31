"""``single_mode_reference``: prompt §7-1 single-mode model (Ref A; protected list). Owner: WP-B1.

R± = L + (1 − L)·exp(±2√(αP)) (shot-noise units, conv §2), phase jitter θ (rms, conv §6):

* ``average: exact``       R′ = (R− + R+)/2 ∓ (R+ − R−)/2 · e^{−2θ²}       (conv §6, T04 eq. (1))
* ``average: small_angle`` R′± ≈ R± cos²θ + R∓ sin²θ                          (prompt §7-1 form)

Offline (O) block: the time loop never steps it (no reads / writes, ``__call__`` is the identity).
Its product is :meth:`SingleModeReference.source_output` (``SourceOutput`` N = 1, F3 constructor,
D55 frame rule 4: cov at φ_p = 0, θ_off = 0). L is the *total* lumped loss of §7-1 (D48) and is
already inside R±, so ``losses()`` is empty (registering it would count it twice, G06).

Params (all required, D34): ``loss`` L, ``alpha_per_w`` α [1/W], ``pump_w`` P [W] (Ref A: OPA
*output* pump, D48), ``theta_rms_rad`` θ, ``phi_opa0_rad`` φ_OPA,0, ``average`` (exact | small_angle).
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.core.units import V_VAC, lin_to_db
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput

__all__ = ["AVERAGES", "SingleModeReference", "make_single_mode_reference", "optimal_pump_w", "r_prime"]

KIND = "single_mode_reference"
AVERAGES = ("exact", "small_angle")


def r_prime(r_minus, r_plus, theta_rms_rad, average: str):
    """(R′−, R′+) [linear, shot-noise units] from (R−, R+) and rms jitter θ (module docstring)."""
    if average == "exact":
        chi = gaussian.chi_gaussian(theta_rms_rad)
        return (gaussian.jitter_average(r_minus, r_plus, chi, delta0_rad=0.0),
                gaussian.jitter_average(r_minus, r_plus, chi, delta0_rad=np.pi / 2))
    if average == "small_angle":
        return (gaussian.jitter_average_small_angle(r_minus, r_plus, theta_rms_rad),
                gaussian.jitter_average_small_angle(r_plus, r_minus, theta_rms_rad))
    raise ValueError(f"average must be one of {AVERAGES}, got {average!r}")


def optimal_pump_w(alpha_per_w: float, theta_rms_rad: float) -> float:
    """P_opt = r_opt²/α, r_opt = ½ ln cot θ (small-angle model; T01 §2.1). R′−,min = L + (1 − L) sin 2θ."""
    return (0.5 * np.log(1.0 / np.tan(theta_rms_rad))) ** 2 / alpha_per_w


class SingleModeReference(BaseBlock):
    """The ``single_mode_reference`` block (module docstring)."""

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "loss", "alpha_per_w", "pump_w", "theta_rms_rad", "phi_opa0_rad", "average", block=name)
        if params["average"] not in AVERAGES:
            raise ValueError(f"block {name!r}: average must be one of {AVERAGES}, got {params['average']!r}")
        super().__init__(KIND, name, params, toggle_names=())

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    @cached_property
    def source_output(self) -> SourceOutput:
        p = self.params
        return SourceOutput.single_mode_reference(
            loss=float(p["loss"]), alpha_per_w=float(p["alpha_per_w"]), pump_w=float(p["pump_w"]),
            phi_opa0_rad=float(p["phi_opa0_rad"]), provenance={"block": self.name},
        )

    def r_pm(self) -> tuple[float, float]:
        """(R−, R+) without jitter [linear, shot-noise units]."""
        vm, vp, _ = self.source_output.project(None)
        return float(vm) / V_VAC, float(vp) / V_VAC

    def r_prime(self) -> tuple[float, float]:
        """(R′−, R′+) with the block's θ and averaging [linear]."""
        rm, rp = self.r_pm()
        a, b = r_prime(rm, rp, float(self.params["theta_rms_rad"]), self.params["average"])
        return float(a), float(b)

    def r_prime_db(self) -> tuple[float, float]:
        """(R′−, R′+) [dB re shot noise]."""
        a, b = self.r_prime()
        return float(lin_to_db(a)), float(lin_to_db(b))


@register_block(KIND)
def make_single_mode_reference(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``single_mode_reference`` block (module docstring)."""
    return SingleModeReference(params, name)
