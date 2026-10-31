"""``rep_rate_tail_check``: repetition-rate limit from the OPO mode tail (T17 §5.3, D50). Owner: WP-B11.

Offline CHK block (G16-j, D34: warns, never stops). The optimal mode of a pulse-pumped OPO decays as
e^{−κ_c t}; its energy after T from the pump turn-on is ε_tail(T) ≈ 8.8 e^{−2κ_c T} (κ_c T ≥ 2, both
presets, T17 §5.2). What leaks into the next slot is ε_tail(T_rep) — a coherent mixing with the next
pulse, not a loss (T11 §3.3 exception (a)). Highest rate with tail ≤ ε: f_rep,max = 2κ_c/ln(8.8/ε)
(ε = 10⁻³: Borealis-type ≲ 45 MHz, Liu-type ≲ 87 MHz). The check runs at construction (build time);
events go to ``chk_log`` (code ``rep_rate_tail``) and :attr:`result`.

``tail_model`` (D34, no default): ``analytic`` (the formula; ``fwhm_hz`` gives κ_c = πΔν_c) or
``source`` (``source``: nested ``pulsed_opo_source`` spec; its exact ε_tail(T_rep), f_rep,max still
from the formula with the source's κ_c).

Params: ``rep_rate_hz``, ``threshold`` ε, ``tail_model`` + its params, ``toggles`` {enabled}: off ⇒ no
check (identity).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.core.block import BaseBlock, Block, BlockSpec, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import make_block, register_block

__all__ = ["RepRateTailCheck", "make_rep_rate_tail_check", "tail_analytic", "max_rep_rate_hz"]

KIND = "rep_rate_tail_check"
TAIL_PREFACTOR = 8.8  # T17 §5.2 [sim], both presets, κ_c T ≥ 2 within 5 %


def tail_analytic(kappa_c: float, t_s: float) -> float:
    """ε_tail(T) ≈ 8.8 e^{−2κ_c T}, clipped to 1."""
    return float(min(1.0, TAIL_PREFACTOR * np.exp(-2.0 * kappa_c * t_s)))


def max_rep_rate_hz(kappa_c: float, threshold: float) -> float:
    """f_rep,max = 2κ_c / ln(8.8/ε)."""
    return float(2.0 * kappa_c / np.log(TAIL_PREFACTOR / threshold))


class RepRateTailCheck(BaseBlock):
    """The ``rep_rate_tail_check`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "rep_rate_hz", "threshold", "tail_model", block=name)
        if params["tail_model"] not in ("analytic", "source"):
            raise ValueError(f"block {name!r}: tail_model must be 'analytic' or 'source'")
        require(params, "fwhm_hz" if params["tail_model"] == "analytic" else "source", block=name)
        super().__init__(KIND, name, params)
        self.chk_log = ChkLog()
        self.result: dict[str, Any] | None = None
        if self._toggles["enabled"]:
            self.result = self._check()

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _check(self) -> dict[str, Any]:
        p = self.params
        trep, eps = 1.0 / float(p["rep_rate_hz"]), float(p["threshold"])
        if p["tail_model"] == "analytic":
            kc = np.pi * float(p["fwhm_hz"])
            tail = tail_analytic(kc, trep)
        else:
            src = make_block(BlockSpec.from_obj(p["source"]))
            kc, tail = src.kappa_c, src.tail(trep)
        res = {"tail": tail, "f_rep_max_hz": max_rep_rate_hz(kc, eps), "threshold": eps, "ok": tail <= eps}
        if not res["ok"]:
            self.chk_log.warn("rep_rate_tail", self.name,
                              f"mode tail into the next slot {tail:.3g} > {eps:g} at f_rep = {p['rep_rate_hz']:g} Hz "
                              f"(limit ≈ {res['f_rep_max_hz']:.3g} Hz)", **res)
        return res


@register_block(KIND)
def make_rep_rate_tail_check(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``rep_rate_tail_check`` block (module docstring)."""
    return RepRateTailCheck(params, name)
