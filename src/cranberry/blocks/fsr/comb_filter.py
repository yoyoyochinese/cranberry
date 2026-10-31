"""``fsr_comb_filter``: removal of the OPO's non-degenerate FSR comb pairs (T17 §2.7, §8.1; D50).
Owner: WP-B11.

Offline (O) block. The cavity comb makes the pairs ω0 ± m·FSR (m ≥ 1) two-mode squeezed: each
sideband alone is thermal with sinh² r_m photons (diagonal pairs only, pump bandwidth ≪ FSR,
[L-SI §1.2]). BHD (LO near ω0) does not see them; photon counting (GBS) sees what the filter leaks:

    n_FSR = η_ins Σ_{m=1}^{m_max} [T_f(+m·FSR) + T_f(−m·FSR)] sinh² r_m

T_f(ν) = Π_stages T_k(ν), ν relative to ω0/2π:

* ``flat_top`` (FBG used in reflection, band-pass): T = floor + (1 − floor)·exp(−|2(ν − c)/B|^{2·order})
  (Borealis: two ~5 GHz bands offset by ~1 GHz, ≥ 30 dB on ±1, ±2, … [M-SI p.19])
* ``mzi`` (unbalanced MZI, [L-Methods Eq. 2]): T = η[(1 − e)(½ + ½cos(2π(ν − ν′)/F)) + e], F = c/ΔL
  the MZI free spectral range, e the extinction floor (Liu: 3 stages, 99.8 %, > 40 dB)

The (0,0) passband transmission η_pass = η_ins·T_f(0) is registered as one ``component`` loss on
``segment`` (G06; never multiplied in here).

Params (D34): ``fsr_hz``, ``stages`` (list of {type: flat_top, center_hz, bandwidth_hz, order, floor}
or {type: mzi, period_hz, offset_hz, extinction, eta}), ``r_m`` (number or callable of ``m``),
``m_max``, ``insertion_eta``, ``segment``, ``applies_to``, ``source_tag``, ``toggles`` {fsr_modes}:
off ⇒ ideal comb removal, n_FSR = 0 (the passband loss stays registered).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.config.resolve import as_callable
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["FsrCombFilter", "make_fsr_comb_filter", "stage_transmission"]

KIND = "fsr_comb_filter"
_STAGE_KEYS = {"flat_top": ("center_hz", "bandwidth_hz", "order", "floor"),
               "mzi": ("period_hz", "offset_hz", "extinction", "eta")}


def stage_transmission(stage: Mapping[str, Any], nu_hz) -> np.ndarray:
    """Power transmission of one stage at ν [Hz] relative to the carrier (module docstring)."""
    nu = np.asarray(nu_hz, dtype=float)
    if stage["type"] == "flat_top":
        x = np.abs(2 * (nu - float(stage["center_hz"])) / float(stage["bandwidth_hz"]))
        fl = float(stage["floor"])
        return fl + (1 - fl) * np.exp(-(x ** (2 * float(stage["order"]))))
    e = float(stage["extinction"])
    return float(stage["eta"]) * ((1 - e) * (0.5 + 0.5 * np.cos(2 * np.pi * (nu - float(stage["offset_hz"])) / float(stage["period_hz"]))) + e)


class FsrCombFilter(BaseBlock):
    """The ``fsr_comb_filter`` block (module docstring)."""

    toggle_names = ("fsr_modes",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "fsr_hz", "stages", "r_m", "m_max", "insertion_eta", "segment", "applies_to", "source_tag", block=name)
        super().__init__(KIND, name, params)
        for i, st in enumerate(params["stages"]):
            if st.get("type") not in _STAGE_KEYS:
                raise ValueError(f"block {name!r}: stages[{i}].type must be one of {tuple(_STAGE_KEYS)}")
            require(st, *_STAGE_KEYS[st["type"]], block=f"{name}.stages[{i}]")
        if not 0.0 <= float(params["insertion_eta"]) <= 1.0 or int(params["m_max"]) < 0:
            raise ValueError(f"block {name!r}: need 0 <= insertion_eta <= 1 and m_max >= 0")

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def transmission(self, nu_hz) -> np.ndarray:
        """T_f(ν) of the stage cascade."""
        t = np.ones_like(np.asarray(nu_hz, dtype=float))
        for st in self.params["stages"]:
            t = t * stage_transmission(st, nu_hz)
        return t

    @property
    def eta_passband(self) -> float:
        return float(self.params["insertion_eta"]) * float(self.transmission(0.0))

    @property
    def n_fsr(self) -> float:
        """Thermal photons of the FSR pairs leaking through the filter (per pulse / mode set)."""
        if not self._toggles["fsr_modes"] or int(self.params["m_max"]) == 0:
            return 0.0
        m = np.arange(1, int(self.params["m_max"]) + 1)
        r = np.broadcast_to(np.asarray(as_callable(self.params["r_m"], ("m",))(m), dtype=float), m.shape)
        f = float(self.params["fsr_hz"])
        leak = self.transmission(m * f) + self.transmission(-m * f)
        return float(self.params["insertion_eta"]) * float(np.sum(leak * np.sinh(r) ** 2))

    def losses(self) -> list[LossElement]:
        return [LossElement(id=f"{self.name}.passband", eta=self.eta_passband, kind="component",
                            applies_to=set(self.params["applies_to"]), source=str(self.params["source_tag"]),
                            enabled=True, segment=str(self.params["segment"]))]


@register_block(KIND)
def make_fsr_comb_filter(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``fsr_comb_filter`` block (module docstring)."""
    return FsrCombFilter(params, name)
