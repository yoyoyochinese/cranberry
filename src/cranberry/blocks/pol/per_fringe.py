"""``per_fringe``: Jones → (x, Δχ) adapter for the Ref D DC fringe (T09 §3.5, T15 eq. (7)). Owner: WP-B7.

The PDB450C does not resolve polarization, so the fringe is the sum of the slow–slow and fast–fast
interference terms. ``dc_fringe`` (WP-B3) models that sum as μ = (1 − x) + x·e^{iΔχ}. This block
computes μ exactly from the Jones matrices of the two arms. It then writes (x, Δχ) so that
dc_fringe's μ equals it:

    a_arm = U_arm·e / |e|,  U_arm = J_n ⋯ J_1      (arm J slots in propagation order, d = +1)
    t_arm = g^†·U_arm·j_b                          (what beam_assembler puts into beam.<arm>.phi)
    μ = (a_pr^† a_lo)·(t_pr t_lo*)/|t_pr t_lo|      ⇒ Re(μ e^{iφ_D}) is the fringe term, with
                                                     φ_D = beam.lo.phi − beam.pr.phi (conv §5)
    x = |μ − 1|²/(2(1 − Re μ)),  e^{iΔχ} = 1 + (μ − 1)/x     (x = Δχ = 0 when μ = 1)

The inversion is exact: every μ in the unit disk lies on exactly one circle (1 − x) + x e^{iΔχ}.
Put dc_fringe's mu_abs = 1 and mu_phase_rad = 0 when it reads these slots. In the simplest case
(common input leakage x, no coupling inside the arms) the result is x itself and
Δχ = Δψ_b,pr − Δψ_b,lo = Δψ_b,arm1 − Δψ_b,arm2. T15 eq. (7) writes arm 2 − arm 1, which is a sign
slip there; T09 §3.5 (−arg(1 + ε e^{−iδ_b}) in φ_pr − φ_lo) agrees with this block. The lock-point
shift is −arg μ, with maximum arcsin(x/(1 − x)).

Inputs: ``in_jones`` is the state entering coupler 1 (complex (2,), e.g. ``polarizer_projection``
out_jones). ``j_pr`` / ``j_lo`` are lists of J slots (``segment_jones``). ``beam_jones_in`` /
``beam_jones_out`` are the arm beams' topology jones_in / jones_out (both arms must share them).

Toggle ``per_crosscoupling`` (same name as in dc_fringe): off ⇒ x = Δχ = 0.

Params (all required, D34): ``in_jones``, ``j_pr``, ``j_lo``, ``beam_jones_in``, ``beam_jones_out``,
``out_x``, ``out_dchi`` (slots, e.g. ``pol.<name>.x`` / ``pol.<name>.dchi``), ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block

__all__ = ["PerFringe", "make_per_fringe", "mu_to_x_dchi"]

KIND = "per_fringe"


def _build(jit: bool):
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def apply(Z, offs, v0, v1):
        for o in offs:
            v0, v1 = Z[o] * v0 + Z[o + 1] * v1, Z[o + 2] * v0 + Z[o + 3] * v1
        return v0, v1

    @dec
    def invert(mu):
        d = mu - 1.0
        den = 2.0 * (1.0 - mu.real)
        if den <= 0.0:
            return 0.0, 0.0
        x = (d.real * d.real + d.imag * d.imag) / den
        z = 1.0 + d / x
        return x, math.atan2(z.imag, z.real)

    @dec
    def step(S, Z, i_in, PR, LO, jb, gc, i_x, i_d, on):  # noqa: N803
        if not on:
            S[i_x], S[i_d] = 0.0, 0.0
            return
        e0, e1 = Z[i_in], Z[i_in + 1]
        p0, p1 = apply(Z, PR, e0, e1)
        l0, l1 = apply(Z, LO, e0, e1)
        m = p0.conjugate() * l0 + p1.conjugate() * l1
        m = m / math.sqrt((abs(p0) ** 2 + abs(p1) ** 2) * (abs(l0) ** 2 + abs(l1) ** 2))
        b0, b1 = apply(Z, PR, jb[0], jb[1])
        tp = gc[0] * b0 + gc[1] * b1
        b0, b1 = apply(Z, LO, jb[0], jb[1])
        tl = gc[0] * b0 + gc[1] * b1
        m = m * tp * tl.conjugate() / (abs(tp) * abs(tl))
        S[i_x], S[i_d] = invert(m)

    return step, invert


(_STEP_PY, _INVERT_PY), (_STEP_NB, _) = _build(False), _build(True)


def mu_to_x_dchi(mu: complex) -> tuple[float, float]:
    """(x, Δχ) with (1 − x) + x e^{iΔχ} = μ (|μ| ≤ 1)."""
    return _INVERT_PY(complex(mu))


class PerFringe(BaseBlock):
    """The ``per_fringe`` block (module docstring)."""

    toggle_names = ("per_crosscoupling",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "in_jones", "j_pr", "j_lo", "beam_jones_in", "beam_jones_out", "out_x", "out_dchi",
                "toggles", block=name)
        for k in ("j_pr", "j_lo"):
            if not isinstance(params[k], (list, tuple)) or not params[k]:
                raise ValueError(f"block {name!r}: {k} must be a non-empty list of J slot names")
        for k in ("beam_jones_in", "beam_jones_out"):
            if np.asarray(params[k]).shape != (2,):
                raise ValueError(f"block {name!r}: {k} must be a 2-vector")
        on = bool(dict(params.get("toggles") or {}).get("per_crosscoupling", False))
        reads = (params["in_jones"], *params["j_pr"], *params["j_lo"]) if on else ()
        super().__init__(KIND, name, params, reads=tuple(dict.fromkeys(reads)),
                         writes=(params["out_x"], params["out_dchi"]))
        self.slot_units = {params["out_x"]: "1", params["out_dchi"]: "rad"}

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        layout.add(self.params["out_x"], (), np.float64, "1")
        layout.add(self.params["out_dchi"], (), np.float64, "rad")

    def bind(self, layout) -> None:
        on = self.toggles["per_crosscoupling"]
        p = self.params

        def zoff(name, shape):
            sl = layout.slot(name)
            if sl.array != "Z" or sl.shape != shape or sl.ring_len:
                raise TypeError(f"block {self.name!r}: {name!r} must be a complex {shape} slot")
            return sl.offset

        i_in = zoff(p["in_jones"], (2,)) if on else 0
        PR = np.array([zoff(s, (2, 2)) for s in p["j_pr"]] if on else [0], dtype=np.int64)
        LO = np.array([zoff(s, (2, 2)) for s in p["j_lo"]] if on else [0], dtype=np.int64)
        jb = np.asarray(p["beam_jones_in"], dtype=np.complex128)
        gc = np.conj(np.asarray(p["beam_jones_out"], dtype=np.complex128))
        i_x, i_d = layout.slot(p["out_x"]).offset, layout.slot(p["out_dchi"]).offset
        self._args = (i_in, PR, LO, jb, gc, i_x, i_d, on)
        nb = _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, Z, i_in, PR, LO, jb, gc, i_x, i_d, on)

        self.kernel = kernel

    def __call__(self, state, params):
        _STEP_PY(state.S, state.Z, *self._args)
        return state


@register_block(KIND)
def make_per_fringe(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``per_fringe`` block (module docstring)."""
    return PerFringe(params, name)
