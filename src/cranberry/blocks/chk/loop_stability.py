"""``loop_stability_check``: linear stability of a lock loop from its analytic open-loop model (Phase 5
M1D-02; refD's 915B-resonance grid at the D53 gains, D51-11, D65). Offline CHK block (G16-j, D34: warns,
never stops, never prunes a sweep point).

Continuous open loop, closed loop 1/(1 + L) (negative feedback):

    L(s) = k0 · (kp + ki/s + kd·s) · e^{−s·delay} · Π_p 1/(1 + s/(2π f_p)) · Π_m ω_m²/(s² + s ω_m/Q_m + ω_m²)

A pole at ``inf`` is absent (presets switch a pole off with an expression). Nyquist criterion for an
open loop without right-half-plane poles and n_int = 1 integrator (ki ≠ 0) or none: stable ⇔ the phase
of 1 + L(i2πf) changes by n_int·π/2 as f goes from the bottom to the top of ``f_grid_hz`` (the contour
closes at f → 0 through the right half plane, M1 review script ``m1_refD_stability_grid.py``). The
grid must resolve the narrowest mode (f_m/Q_m) and reach |L| ≪ 1 at the top; |L(f_hi)| ≥ 0.1 is
reported as inconclusive. Reported with ``max_l_at_180`` (largest |L| at a −180° crossing; > 1
unstable, 1/x = gain margin) and its frequency.

The result is in :attr:`result`; an unstable (or inconclusive) point goes to ``chk_log`` (code
``loop_unstable`` / ``loop_check_grid``) as a :class:`~cranberry.core.checks.ChkWarning`. Whether such
points may be dropped from a prediction band is not decided here (refD: open user decision under D04).

Params (all required, D34): ``k0`` (plant gain from controller output to controller input, DC),
``kp``, ``ki_per_s``, ``kd_s``, ``delay_s``, ``poles_hz`` (list, ``inf`` = absent), ``modes``
(list of {f_hz, q}), ``f_grid_hz`` [lo, hi, n] (log grid), ``toggles`` {enabled}: off ⇒ no check.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block

__all__ = ["LoopStabilityCheck", "make_loop_stability_check", "open_loop", "nyquist_stable"]

KIND = "loop_stability_check"


def open_loop(f: np.ndarray, p: Mapping[str, Any]) -> np.ndarray:
    """L(i2πf) of the module docstring for the params ``p``."""
    s = 2j * np.pi * np.asarray(f, float)
    lg = float(p["k0"]) * (float(p["kp"]) + float(p["ki_per_s"]) / s + float(p["kd_s"]) * s)
    lg = lg * np.exp(-s * float(p["delay_s"]))
    for fp in p["poles_hz"]:
        if math.isfinite(float(fp)):
            lg = lg / (1.0 + s / (2.0 * np.pi * float(fp)))
    for m in p["modes"]:
        w = 2.0 * np.pi * float(m["f_hz"])
        lg = lg * w * w / (s * s + s * w / float(m["q"]) + w * w)
    return lg


def nyquist_stable(f: np.ndarray, lg: np.ndarray, n_int: int) -> tuple[bool, float, float]:
    """(stable, max |L| at a −180° crossing, its frequency) on the grid ``f`` (module docstring)."""
    darg = float(np.unwrap(np.angle(1.0 + lg))[[0, -1]] @ [-1.0, 1.0])
    ph = np.unwrap(np.angle(lg))
    cross = np.nonzero(np.diff(np.floor((ph + np.pi) / (2.0 * np.pi))))[0]
    if cross.size:
        i = cross[np.argmax(np.abs(lg[cross]))]
        mx, f180 = float(abs(lg[i])), float(f[i])
    else:
        mx, f180 = 0.0, math.nan
    return abs(darg / (np.pi / 2) - n_int) < 0.2, mx, f180


class LoopStabilityCheck(BaseBlock):
    """The ``loop_stability_check`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "k0", "kp", "ki_per_s", "kd_s", "delay_s", "poles_hz", "modes", "f_grid_hz", block=name)
        super().__init__(KIND, name, params)
        lo, hi, n = params["f_grid_hz"]
        if not (0 < float(lo) < float(hi) and int(n) >= 2):
            raise ValueError(f"block {name!r}: f_grid_hz must be [lo, hi, n] with 0 < lo < hi, n >= 2")
        self.chk_log = ChkLog()
        self.result: dict[str, Any] | None = None
        if self._toggles["enabled"]:
            self.result = self._check()

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _check(self) -> dict[str, Any]:
        p = self.params
        lo, hi, n = p["f_grid_hz"]
        f = np.geomspace(float(lo), float(hi), int(n))
        lg = open_loop(f, p)
        n_int = 1 if float(p["ki_per_s"]) != 0.0 else 0
        stable, mx, f180 = nyquist_stable(f, lg, n_int)
        res = {"stable": bool(stable), "max_l_at_180": mx, "f_180_hz": f180, "l_at_f_hi": float(abs(lg[-1])),
               "modes": [dict(m) for m in p["modes"]]}
        if res["l_at_f_hi"] >= 0.1:
            self.chk_log.warn("loop_check_grid", self.name,
                              f"|L| = {res['l_at_f_hi']:.3g} at the top of the grid ({float(hi):g} Hz): "
                              "stability not decided", **res)
        elif not stable:
            self.chk_log.warn("loop_unstable", self.name,
                              f"linear loop unstable (max |L| at −180° = {mx:.3g} at {f180:.4g} Hz, modes "
                              f"{res['modes']}): the time loop will oscillate / slip instead of locking", **res)
        return res


@register_block(KIND)
def make_loop_stability_check(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``loop_stability_check`` block (module docstring)."""
    return LoopStabilityCheck(params, name)
