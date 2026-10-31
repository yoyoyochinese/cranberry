"""``opo_angle_map``: OPO detuning → segment × beam phase of the OPO segment (T17 §2.5, §8.3; T02).
Owner: WP-B11.

Time-loop (T) block. Per step, for every beam b crossing the OPO segment::

    psi.<segment>.<b> = c_s,b(x_ss)·Δ_s + c_p,b(x_ss)·Δ_p          [rad]

Δ_s [rad/s] from ``in_delta_s`` (e.g. ``opo_length_lock.out``); Δ_p from ``in_delta_p`` or, when that is
None, Δ_p = r·Δ_s with the runtime ratio r = ``delta_p_ratio`` (2 for a cavity-length change, T17 §1;
0 when only the signal resonance moves). x_ss from ``in_x_ss`` (e.g. ``pump_power_stabilization``) or
the runtime param ``x_ss``. (c_s, c_p) [rad per rad/s] are interpolated linearly in x_ss from
``coef`` on ``x_grid`` (x outside the grid is clamped to the end points: keep the grid wide enough);
``pulsed_opo_source.angle_map_params`` builds them from its detuning tables (roles ``squeezed``: the
squeezing-angle rotation as a phase of the squeezed beam; ``pump``: the transmitted pump phase
arg∫b dt, the beam a Liu-type pump-phase lock looks at). Beams not crossing the OPO have no entry.

D56-4: this block writes ``psi.<segment>.<beam>`` itself, so the OPO segment must be left out of
``segment_phase.segments``. The map is first order in Δ (linear up to Δ/κ_c ≈ 0.05, T17 V-P9); the
second-order loss of squeezing is not modelled (−7.431 → −7.405 dB at Δ/κ_c = 0.05).

Params (D34): ``segment``, ``beams`` (list), ``in_delta_s``, ``in_delta_p`` (slot or None),
``delta_p_ratio`` (required if ``in_delta_p`` is None, runtime), ``in_x_ss`` (slot or None), ``x_ss``
(required if ``in_x_ss`` is None, runtime), ``x_grid`` (ascending), ``coef`` {beam: [[c_s, c_p] per
x_grid point]}, ``toggles`` {map}: off ⇒ psi = 0 (identity of the OPO segment).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block

__all__ = ["OpoAngleMap", "make_opo_angle_map"]

KIND = "opo_angle_map"


def _make_step(xg, cs, cp, out, i_ds, i_dp, p_ratio, i_x, p_x, on):
    nx = xg.size

    def step(S, P):  # noqa: N803
        ds = S[i_ds]
        dp = S[i_dp] if i_dp >= 0 else P[p_ratio] * ds
        x = S[i_x] if i_x >= 0 else P[p_x]
        i, w = 0, 0.0
        if nx > 1:
            x = min(max(x, xg[0]), xg[nx - 1])
            i = 0
            while i < nx - 2 and x > xg[i + 1]:
                i += 1
            w = (x - xg[i]) / (xg[i + 1] - xg[i])
        for j in range(out.size):
            if on:
                a = cs[j, i] if nx == 1 else (1.0 - w) * cs[j, i] + w * cs[j, i + 1]
                b = cp[j, i] if nx == 1 else (1.0 - w) * cp[j, i] + w * cp[j, i + 1]
                S[out[j]] = a * ds + b * dp
            else:
                S[out[j]] = 0.0

    return step


class OpoAngleMap(BaseBlock):
    """The ``opo_angle_map`` block (module docstring)."""

    toggle_names = ("map",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "segment", "beams", "in_delta_s", "in_delta_p", "in_x_ss", "x_grid", "coef", block=name)
        runtime = []
        reads = [str(params["in_delta_s"])]
        if params["in_delta_p"] is None:
            runtime.append("delta_p_ratio")
        else:
            reads.append(str(params["in_delta_p"]))
        if params["in_x_ss"] is None:
            runtime.append("x_ss")
        else:
            reads.append(str(params["in_x_ss"]))
        seg = str(params["segment"])
        self.beams = [str(b) for b in params["beams"]]
        self.outs = [f"psi.{seg}.{b}" for b in self.beams]
        super().__init__(KIND, name, params, reads=reads, writes=self.outs, runtime_params=runtime)
        xg = np.asarray(params["x_grid"], dtype=float)
        if xg.ndim != 1 or xg.size < 1 or np.any(np.diff(xg) <= 0):
            raise ValueError(f"block {name!r}: x_grid must be a non-empty ascending list")
        missing = set(self.beams) - set(params["coef"])
        if missing:
            raise ValueError(f"block {name!r}: coef missing for beams {sorted(missing)}")
        coef = np.array([np.asarray(params["coef"][b], dtype=float) for b in self.beams])
        if coef.shape != (len(self.beams), xg.size, 2):
            raise ValueError(f"block {name!r}: coef[beam] must be (len(x_grid), 2) = ({xg.size}, 2) [c_s, c_p]")
        self._xg, self._cs, self._cp = xg, np.ascontiguousarray(coef[..., 0]), np.ascontiguousarray(coef[..., 1])
        self.slot_units = {str(params["in_delta_s"]): "rad/s", **{o: "rad" for o in self.outs}}
        if params["in_delta_p"] is not None:
            self.slot_units[str(params["in_delta_p"])] = "rad/s"

    def init(self, layout, rng):
        super().init(layout, rng)
        for o in self.outs:
            layout.add(o, (), np.float64, "rad")

    def bind(self, layout):
        p = self.params
        sl = lambda k: layout.slot(str(p[k])).offset if p[k] is not None else -1  # noqa: E731
        pr = lambda k: layout.param(self.param_name(k)).offset if k in self.runtime_params else -1  # noqa: E731
        out = np.array([layout.slot(o).offset for o in self.outs], dtype=np.int64)
        args = (self._xg, self._cs, self._cp, out, sl("in_delta_s"), sl("in_delta_p"), pr("delta_p_ratio"),
                sl("in_x_ss"), pr("x_ss"), self._toggles["map"])
        self._py = _make_step(*args)
        nb = numba.njit(_make_step(*args))

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803, E741
            nb(S, P)

        self.kernel = kernel

    def __call__(self, state, params):
        self._py(state.S, params.P)
        return state


@register_block(KIND)
def make_opo_angle_map(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opo_angle_map`` block (module docstring)."""
    return OpoAngleMap(params, name)
