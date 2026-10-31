"""Locks as linear combinations of beam phases (conv §9). Owner: WP-F6.

e(t) = Σ_b c_b φ_b(t) − setpoint. No wrapping here (controllers see the continuous phase error;
metrics wrap to (−π, π] themselves).

Block ``lock_error`` (WP-F6): writes the true phase error of each selected lock, e.g. for Ref D
lockD: e = φ_LO − φ_pr − π/2 (conv §5, §9). It is the *physical* error in the conv §9 sign; a sensor /
device block (``dc_fringe``, DSC1 with e = S − v) converts it into its own signal and sign inside
itself. The combo beams are all evaluated at the lock's observer (checked by ``Topology``).

    reads   beam.<b>.phi          rad   every beam of the selected locks' combos
    writes  lock.<id>.err         rad
    P       <block>.<id>.setpoint_rad   runtime parameter, initialised from the topology

Params (all required, D34): ``topology`` (mapping or Topology), ``locks`` ("all" | [lock id, ...]).
No toggles (a lock error has no on/off physics; switch the lock off in the controller).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.path.topology import Lock, as_topology

__all__ = ["lock_error", "LockError", "make_lock_error"]


def lock_error(lock: Lock, beam_phases_rad: Mapping[str, float]) -> float:
    """e = Σ_b combo[b]·φ_b − setpoint [rad] for one sample (conv §9 sign)."""
    missing = [b for b in lock.combo if b not in beam_phases_rad]
    if missing:
        raise KeyError(f"lock {lock.id!r}: no phase given for beam(s) {missing}")
    return float(sum(c * float(beam_phases_rad[b]) for b, c in lock.combo.items()) - lock.setpoint_rad)


class LockError(BaseBlock):
    """The ``lock_error`` block (module docstring)."""

    toggle_names = ()

    def __init__(self, params: Mapping[str, Any], name: str | None):
        from cranberry.path.segments import select_ids

        name = name or "lock_error"
        require(params, "topology", "locks", block=name)
        topo = as_topology(params["topology"])
        self.locks = [topo.locks[i] for i in select_ids(params["locks"], topo.locks, "locks", name)]
        reads = list(dict.fromkeys(f"beam.{b}.phi" for lk in self.locks for b in lk.combo))
        super().__init__("lock_error", name, params, reads=reads, writes=[f"lock.{lk.id}.err" for lk in self.locks])
        self.slot_units = {s: "rad" for s in (*self.reads, *self.writes)}

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        for lk in self.locks:
            layout.add(f"lock.{lk.id}.err", (), np.float64, "rad")
            layout.add_param(f"{self.name}.{lk.id}.setpoint_rad", lk.setpoint_rad, (), "rad")

    def bind(self, layout) -> None:
        out, sp, t0, n, src, c = [], [], [], [], [], []
        for lk in self.locks:
            out.append(layout.slot(f"lock.{lk.id}.err").offset)
            sp.append(layout.param(f"{self.name}.{lk.id}.setpoint_rad").offset)
            t0.append(len(src))
            n.append(len(lk.combo))
            for b, cb in lk.combo.items():
                src.append(layout.slot(f"beam.{b}.phi").offset)
                c.append(cb)
        self._sp = sp
        OUT, SP, T0, N = (np.array(x, dtype=np.int64) for x in (out, sp, t0, n))
        SRC = np.array(src if src else [0], dtype=np.int64)
        C = np.array(c if c else [0.0], dtype=np.float64)
        n_locks = len(self.locks)

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            for j in range(n_locks):
                e = 0.0
                for t in range(T0[j], T0[j] + N[j]):
                    e += C[t] * S[SRC[t]]
                S[OUT[j]] = e - P[SP[j]]

        self.kernel = kernel

    def __call__(self, state, params):
        for lk in self.locks:
            sp = float(params[f"{self.name}.{lk.id}.setpoint_rad"])
            phases = {b: float(state[f"beam.{b}.phi"]) for b in lk.combo}
            e = 0.0
            for b, c in lk.combo.items():
                e += c * phases[b]
            state[f"lock.{lk.id}.err"] = e - sp
        return state


@register_block("lock_error")
def make_lock_error(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``lock_error`` block (module docstring)."""
    return LockError(params, name)
