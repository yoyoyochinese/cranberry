"""``polarizer_projection``: in-line polarizer, pass axis = lab x (T09 §1.4, §3.5; D51-8). Owner: WP-B7.

The block reads the Jones vector e arriving at the polarizer, for example ``beam.<b>.jones`` of a
monitor beam whose observer is the polarizer. It writes:

    e_out = diag(1, r)·e,  r = 10^(−ER/20)
    out_dp    = |e_out|²/|e|² − 1  relative power change after the polarizer, against an aligned
                                   input (= cos²(Θ/2) + 10^(−ER/10) sin²(Θ/2) − 1 ≈ −Θ²/4)
    out_jones = e_out/|e_out|      normalised output state: the input of the PM part downstream
                                   (the fast-axis leakage x = |e_out,y|²)

out_dp is the polarizer power drift that the DC fringe lock sees (T15 §2.4-4). It is ready for
``dc_fringe.in_dp``. Sum it with RIN in the preset if both are on. The insertion loss IL is
registered once through ``losses()`` and never applied here (G06). The axis misalignment of the
output pigtail is the next segment's ``theta_in_rad`` (``segment_jones``), so it is not repeated
here. Ref D (D51-8, Opneti): IL 0.3 dB typ, ER 30 dB typ / ≥ 28 dB, about −2 dB with connectors.

Toggles (all off ⇒ out_dp = 0, out_jones = (1, 0)):
    pol_drift_before_polarizer   off ⇒ the input is treated as aligned (e = (1, 0))
    pol_extinction               off ⇒ ideal polarizer (r = 0)

Params (all required, D34): ``in``, ``out_jones``, ``out_dp`` (slots), ``er_db`` (runtime),
``il_db``, ``segment`` (loss location), ``applies_to``, ``source``, ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["PolarizerProjection", "make_polarizer_projection"]

KIND = "polarizer_projection"


def _step(S, Z, P, i_in, i_out, i_dp, p_er, drift, ext):  # noqa: N803
    e0, e1 = (Z[i_in], Z[i_in + 1]) if drift else (1.0 + 0j, 0j)
    r = 10.0 ** (-P[p_er] / 20.0) if ext else 0.0
    n_in = abs(e0) ** 2 + abs(e1) ** 2
    a, b = e0, r * e1
    n_out = abs(a) ** 2 + abs(b) ** 2
    if n_out == 0.0:
        Z[i_out], Z[i_out + 1] = 0j, 0j
        S[i_dp] = -1.0
        return
    k = 1.0 / math.sqrt(n_out)
    Z[i_out], Z[i_out + 1] = a * k, b * k
    S[i_dp] = n_out / n_in - 1.0


_STEP_NB = numba.njit(_step)


class PolarizerProjection(BaseBlock):
    """The ``polarizer_projection`` block (module docstring)."""

    toggle_names = ("pol_drift_before_polarizer", "pol_extinction")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "in", "out_jones", "out_dp", "er_db", "il_db", "segment", "applies_to", "source", block=name)
        drift = bool(dict(params.get("toggles") or {}).get("pol_drift_before_polarizer", False))
        super().__init__(KIND, name, params, reads=(params["in"],) if drift else (),
                         writes=(params["out_jones"], params["out_dp"]), runtime_params=("er_db",))
        self.slot_units = {params["out_dp"]: "1"}

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        layout.add(self.params["out_jones"], (2,), np.complex128, "1", init=(1.0, 0.0))
        layout.add(self.params["out_dp"], (), np.float64, "1")

    def bind(self, layout) -> None:
        drift, ext = self.toggles["pol_drift_before_polarizer"], self.toggles["pol_extinction"]
        if drift:
            sl = layout.slot(self.params["in"])
            if sl.array != "Z" or sl.shape != (2,) or sl.ring_len:
                raise TypeError(f"block {self.name!r}: {sl.name!r} must be a complex (2,) Jones-vector slot")
        i_in = layout.slot(self.params["in"]).offset if drift else 0
        i_out, i_dp = layout.slot(self.params["out_jones"]).offset, layout.slot(self.params["out_dp"]).offset
        p_er = layout.param(self.param_name("er_db")).offset
        self._args = (i_in, i_out, i_dp, p_er, drift, ext)
        nb = _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, Z, P, i_in, i_out, i_dp, p_er, drift, ext)

        self.kernel = kernel

    def __call__(self, state, params):
        _step(state.S, state.Z, params.P, *self._args)
        return state

    def losses(self) -> list[LossElement]:
        p = self.params
        return [LossElement.from_il_db(f"{self.name}.il", float(p["il_db"]), kind="component",
                                       applies_to=set(p["applies_to"]), source=str(p["source"]), enabled=True,
                                       segment=str(p["segment"]))]


@register_block(KIND)
def make_polarizer_projection(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``polarizer_projection`` block (module docstring)."""
    return PolarizerProjection(params, name)
