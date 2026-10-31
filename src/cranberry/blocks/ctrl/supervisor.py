"""``lock_supervisor``: lock-loss detection and re-acquisition (PLAN §4.3 CTRL; T16 §3.6, E83).

The supervisor sits in line between a controller and its actuator chain: out = u_in + offset, where
``offset`` is its own correction (so it never writes another block's state and needs no extra
feedback path). Per step (h = rate_div / f_sim; ``state`` 1 locked, 0 unlocked, 2 holding after a
reset / hop):

* lock loss (from state 1): the output u_in + offset lies outside [u_min_v, u_max_v] for ≥ t_sat_s,
  or |e_in| > e_th for ≥ t_e_s (``e_th = inf`` disables; t = 0 means "on the first sample").
  Counts ``n_loss``, stores the time of the first loss ``t_loss`` (−1 before).
* fringe slips (``true_in``, the true phase error [rad], may be null): k = round(δ / 2π); every change
  of k adds |Δk| to ``n_slip`` (T16 §3.6 item 2: counted on the unwrapped true phase).
* re-acquisition (structural ``relock``; toggle ``relock`` False = detection only, out = u_in):
  ``none``; ``reset_center`` (offset = centre − u_in: the output jumps to the range centre and the loop
  re-acquires by itself); ``fringe_hop`` (offset −= n·u_2pi_v with n the integer bringing the output
  closest to the centre: same fringe point modulo 2π, T16 §3.6 item 2); ``ramp_search`` (output sweeps
  [u_min, u_max] as a sawtooth of period ``ramp_period_s`` starting at the current output; it engages
  bumplessly (offset = u_ramp − u_in) when |e| < e_cap and sign(Δe·Δu) = ``slope_sign``, i.e. at a zero
  crossing of the correct slope, T16 §3.6 item 3). After reset / hop the state is 2 for ``t_hop_s``
  (data to discard), then 1. Every re-acquisition counts ``n_relock``.

Params: ``u_in``, ``e_in``, ``true_in`` (slot or null), ``out``, ``relock``; runtime ``u_min_v``,
``u_max_v``, ``t_sat_s``, ``e_th``, ``t_e_s``; ``reset_center`` / ``fringe_hop`` add ``t_hop_s``,
``fringe_hop`` adds ``u_2pi_v`` (915B: 2π/K_s = 0.968 V); ``ramp_search`` adds ``ramp_period_s``,
``e_cap``, ``slope_sign``. Slots: ``out`` (V); ``lock.<name>.state``, ``.n_loss``, ``.n_relock``,
``.n_slip`` (int); ``lock.<name>.t_loss`` (s); ``lock.<name>.x`` (internal).
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

__all__ = ["RELOCK_MODES", "make_lock_supervisor"]

RELOCK_MODES: dict[str, int] = {"none": 0, "reset_center": 1, "fringe_hop": 2, "ramp_search": 3}
_BASE = ("u_min_v", "u_max_v", "t_sat_s", "e_th", "t_e_s")
_EXTRA = {0: (), 1: ("t_hop_s",), 2: ("t_hop_s", "u_2pi_v"), 3: ("ramp_period_s", "e_cap", "slope_sign")}
_ALL = _BASE + ("t_hop_s", "u_2pi_v", "ramp_period_s", "e_cap", "slope_sign")
_TWO_PI = 2.0 * math.pi


def _step(S, I, P, i_t, i_u, i_e, i_true, i_out, x, st, nl, nr, ns, i_tl, po, h, mode):  # noqa: N803
    u_in = S[i_u]
    e = S[i_e]
    umin, umax = P[po[0]], P[po[1]]
    if i_true >= 0:
        k = math.floor(S[i_true] / _TWO_PI + 0.5)
        dk = k - S[x + 6]
        if dk != 0.0:
            I[ns] += int(abs(dk))
            S[x + 6] = k
    state = I[st]
    if state == 0 and mode == 3:
        T = P[po[7]]
        rt = S[x + 3] + h
        if rt >= T:
            rt -= T
        S[x + 3] = rt
        u_r = umin + (umax - umin) * rt / T
        de, du = e - S[x + 4], u_r - S[x + 5]
        S[x + 4], S[x + 5] = e, u_r
        if abs(e) < P[po[8]] and de * du * P[po[9]] > 0.0:
            S[x] = u_r - u_in
            I[st] = 1
            I[nr] += 1
            S[x + 1], S[x + 2] = 0.0, 0.0
        S[i_out] = u_r
        return
    u = u_in + S[x]
    S[x + 1] = S[x + 1] + 1.0 if (u > umax or u < umin) else 0.0
    S[x + 2] = S[x + 2] + 1.0 if abs(e) > P[po[3]] else 0.0
    if state == 2:
        S[x + 7] -= h
        if S[x + 7] <= 0.0:
            state = 1
            I[st] = 1
    if state == 1 and (S[x + 1] > 0.0 and S[x + 1] * h >= P[po[2]] or S[x + 2] > 0.0 and S[x + 2] * h >= P[po[4]]):
        I[nl] += 1
        if S[i_tl] < 0.0:
            S[i_tl] = S[i_t]
        S[x + 1], S[x + 2] = 0.0, 0.0
        c = 0.5 * (umin + umax)
        if mode == 1:
            S[x] = c - u_in
            I[st], I[nr] = 2, I[nr] + 1
            S[x + 7] = P[po[5]]
        elif mode == 2:
            S[x] -= math.floor((u - c) / P[po[6]] + 0.5) * P[po[6]]
            I[st], I[nr] = 2, I[nr] + 1
            S[x + 7] = P[po[5]]
        elif mode == 3:
            I[st] = 0
            r = (u - umin) / (umax - umin)
            r = 0.0 if r < 0.0 else (1.0 if r > 1.0 else r)
            S[x + 3] = r * P[po[7]]
            S[x + 4], S[x + 5] = e, u
        else:
            I[st] = 0
        u = u_in + S[x]
    S[i_out] = u


class LockSupervisor(KitBlock):
    KIND = "lock_supervisor"
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("relock",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "u_in", "e_in", "true_in", "out", "relock", block=name)
        if params["relock"] not in RELOCK_MODES:
            raise ValueError(f"block {name!r}: relock must be one of {sorted(RELOCK_MODES)}")
        mode = RELOCK_MODES[params["relock"]]
        reads = ("t", params["u_in"], params["e_in"]) + (() if params["true_in"] is None else (params["true_in"],))
        self.lk = {k: f"lock.{name}.{k}" for k in ("state", "n_loss", "n_relock", "n_slip", "t_loss", "x")}
        super().__init__(params, name, reads=reads, writes=(params["out"], *self.lk.values()),
                         runtime=_BASE + _EXTRA[mode])
        self.mode = mode if self.toggles["relock"] else 0
        self.slot_units = {params["out"]: "V", params["u_in"]: "V"}

    def slots(self):
        lk = self.lk
        return [(self.params["out"], (), np.float64, "V", 0, 0.0),
                (lk["state"], (), np.int64, "1", 0, 1), (lk["n_loss"], (), np.int64, "1", 0, 0),
                (lk["n_relock"], (), np.int64, "1", 0, 0), (lk["n_slip"], (), np.int64, "1", 0, 0),
                (lk["t_loss"], (), np.float64, "s", 0, -1.0), (lk["x"], (8,), np.float64, "1", 0, 0.0)]

    def make_args(self, layout):
        po = np.array([self.p(layout, k) if layout.has_param(self.param_name(k)) else -1 for k in _ALL], np.int64)
        lk = self.lk
        return (self.o(layout, "t"), self.o(layout, self.params["u_in"]), self.o(layout, self.params["e_in"]),
                self.o(layout, self.params["true_in"]), self.o(layout, self.params["out"]), self.o(layout, lk["x"]),
                self.o(layout, lk["state"]), self.o(layout, lk["n_loss"]), self.o(layout, lk["n_relock"]),
                self.o(layout, lk["n_slip"]), self.o(layout, lk["t_loss"]), po, self.h, self.mode)


@register_block("lock_supervisor")
def make_lock_supervisor(params: Mapping[str, Any], *, name: str | None = None) -> LockSupervisor:
    """Factory of the ``lock_supervisor`` block (module docstring)."""
    return LockSupervisor(params, name)
