"""``loop_stability_check``: linear stability of a lock loop from its analytic open-loop model (Phase 5
M1D-02; refD's 915B-resonance grid at the D53 gains, D51-11, D65). Offline CHK block (G16-j, D34: warns,
never stops, never prunes a sweep point).

Continuous open loop, closed loop 1/(1 + L) (negative feedback):

    L(s) = k0 · C_1(s) · e^{−s·delay} · Π_p 1/(1 + s/(2π f_p)) · Π_m ω_m²/(s² + s ω_m/Q_m + ω_m²) · Σ(s)
    C_k(s) = kp_k + ki_k/s + kd_k·s                                  (discretisation continuous)
    C_k(s) = kp_k + ki_k·h/2 + ki_k/s + kd_k·(1 − e^{−s h})/h          (discretisation dsc1)

**Controller discretisation** (optional ``discretisation``: ``continuous`` (absent) | ``dsc1``, the latter
with ``h_s`` = h, the controller update period; per controller, so a stage may carry its own). ``dsc1`` is
the DSC1 law (DSC1 p.19, ``pid`` discretisation ``dsc1``): integral h·Σ K_i e including e[n] and the
backward difference (K_d/h)(e[n] − e[n−1]). Its integrator is h/(1 − z⁻¹) = h/2 + (h/2)coth(sh/2): the h/2
is an extra proportional gain ki·h/2, and (h/2)coth(sh/2) is replaced by 1/s (relative error (ωh)²/12 of
the I term; the exact form has poles at every multiple of 1/h that a continuous Nyquist test cannot
take). The backward difference is kept exact (entire in s; it is kd·s·e^{−sh/2}·sinc, i.e. half a period
of extra lag on the D path). ``delay_s`` stays the transport latency (ADC average, output hold): the
controller's own discretisation is not a delay. **Why it matters (Ref D, D53 gains):** the continuous
law has complex zeros at √(ki/kd)/2π = 28.3 kHz with ζ = kp/(2√(ki·kd)) = 0.03, a deep notch; the DSC1 has
kp + ki·h (ζ_eff ≈ 0.12, 3.9 × kp) there, so the continuous check under-states |L| by up to ×6 at 25–31 kHz
and mis-states the phase away from it (the I-dominated side gains lead, ≈ 9° at 14 kHz; the D-dominated
side loses it, ≈ 30° at 52–57 kHz). With the continuous law
the frozen Ref D prediction classified 6 of 128 LHS points opposite to its time loop (lhs_002, lhs_042
"stable" but oscillating; lhs_024, 034, 044 at 52–57 kHz "unstable" but locking); with ``dsc1`` all five
agree and the predicted closed-loop pole (growth rate, frequency) matches the noise-free time loop
(reviews/phase5_refD_stability_mismatch.md). A controller of another discrete law is still checked with
``continuous`` (Tustin / forward Euler are not modelled here; at MHz update rates the difference is small).

**What the verdict means (modelling limits).** The check is linear and small-signal:
* The plant gain ``k0`` is the fringe slope at the lock point (Ref D: γ_coh·|μ| ≈ 1, offsets neglected).
* The sampled-data loop is treated as continuous: aliasing of the plant response around multiples of
  1/h (and of the faster base rate) is neglected; for Ref D the exact multirate z-domain loop
  (tests/integration/refD_support.py) agrees with the ``dsc1`` form within ≈ 2° and ≤ 7 % within ±5 % of the mode
  (the ``dsc1`` form slightly over-states |L|, i.e. errs on the unstable side).
* ``unstable`` means the lock point is linearly unstable, not that the phase must run away. In a fringe
  lock the error sin(δ) has describing-function gain 2J₁(A)/A < 1 for an oscillation of amplitude A, so a
  point with max |L| at −180° only slightly above 1 (≲ 1.05) settles into a bounded limit cycle of a few
  tenths of a rad at the −180° frequency (Ref D: lhs_057, lhs_118, 10.7–11.45 kHz Q 10, std 0.25–0.29 rad,
  also with noise off) that a loose "locked" criterion (std < 0.3 rad) counts as a lock. Larger max |L|
  gives oscillations across the fringe (≈ 1.3 rad rms) or slips. Treat max |L| within ≈ ±10 % of 1 as
  marginal: the outcome there also depends on the details above.
* ``stable`` says nothing about large-signal effects (pull-in from a fringe top, saturation, the DAC
  limit cycle of D65-3).

Single loop (no ``stages``): Σ = 1. **Stacked (cascade / offload) loop** (Phase 5 M2C-02, Ref C's
EOM → piezo mirror → stretcher, T16 §4.4): ``stages`` lists the controllers 2..N; controller k reads
controller k − 1's output and drives actuator k, whose response relative to actuator 1 (whose gain is in
``k0``) is A_k = k_act,k · Π_p 1/(1 + s/(2π f_p)) · Π_m (modes as above). The phases of all actuators add:

    Σ(s) = 1 + Σ_{k=2..N} (Π_{j=2..k} C_j(s)) · A_k(s)

so the top-level ``poles_hz`` / ``modes`` / ``delay_s`` are the parts common to every path (sensor,
demodulator, latency) and of actuator 1. A pole at ``inf`` is absent (presets switch a pole off with an
expression). Nyquist criterion for an open loop without right-half-plane poles and n_int integrators
(the pole order of L at s = 0: the largest number of integrators ki ≠ 0 in C_1 … C_k over the paths k
with k_act ≠ 0): stable ⇔ the phase of 1 + L(i2πf) changes by n_int·π/2 as f goes from the bottom to the
top of ``f_grid_hz`` (the contour passes right of the origin; the indentation contributes −n_int·π, M1
review script ``m1_refD_stability_grid.py``). The grid must start where the n_int-fold integrator
dominates, resolve the narrowest mode (f_m/Q_m) and reach |L| ≪ 1 at the top; |L(f_hi)| ≥ 0.1 is
reported as inconclusive. Reported with ``max_l_at_180`` (largest |L| at a −180° crossing; > 1
unstable, 1/x = gain margin; for a stacked loop with several crossings it is not a margin) and its
frequency.

The result is in :attr:`result`; an unstable (or inconclusive) point goes to ``chk_log`` (code
``loop_unstable`` / ``loop_check_grid``) as a :class:`~cranberry.core.checks.ChkWarning`. Whether such
points may be dropped from a prediction band is not decided here (refD: open user decision under D04).

Params (all required, D34): ``k0`` (plant gain from controller output to controller input, DC; its
sign makes L positive for negative feedback), ``kp``, ``ki_per_s``, ``kd_s``, ``delay_s``, ``poles_hz``
(list, ``inf`` = absent), ``modes`` (list of {f_hz, q}), ``f_grid_hz`` [lo, hi, n] (log grid),
``toggles`` {enabled}: off ⇒ no check. Optional ``stages`` (absent = single loop, the M1 form): list of
{kp, ki_per_s, kd_s, k_act, poles_hz, modes} (all keys required per stage). Optional ``discretisation``
(+ ``h_s`` for ``dsc1``) at the top level and per stage (above).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block

__all__ = ["LoopStabilityCheck", "make_loop_stability_check", "open_loop", "n_integrators", "nyquist_stable"]

KIND = "loop_stability_check"


STAGE_KEYS = ("kp", "ki_per_s", "kd_s", "k_act", "poles_hz", "modes")
DISCRETISATIONS = ("continuous", "dsc1")


def _ctrl(s: np.ndarray, q: Mapping[str, Any]) -> np.ndarray:
    kp, ki, kd = float(q["kp"]), float(q["ki_per_s"]), float(q["kd_s"])
    if q.get("discretisation", "continuous") == "dsc1":
        h = float(q["h_s"])
        return kp + ki * h / 2.0 + ki / s + kd * (1.0 - np.exp(-s * h)) / h
    return kp + ki / s + kd * s


def _check_discretisation(q: Mapping[str, Any], where: str) -> None:
    d = q.get("discretisation", "continuous")
    if d not in DISCRETISATIONS:
        raise ValueError(f"block {where!r}: discretisation must be one of {DISCRETISATIONS}, got {d!r}")
    if d == "dsc1":
        require(q, "h_s", block=where)


def _dynamics(s: np.ndarray, q: Mapping[str, Any]) -> np.ndarray:
    out = np.ones_like(s)
    for fp in q["poles_hz"]:
        if math.isfinite(float(fp)):
            out = out / (1.0 + s / (2.0 * np.pi * float(fp)))
    for m in q["modes"]:
        w = 2.0 * np.pi * float(m["f_hz"])
        out = out * w * w / (s * s + s * w / float(m["q"]) + w * w)
    return out


def open_loop(f: np.ndarray, p: Mapping[str, Any]) -> np.ndarray:
    """L(i2πf) of the module docstring for the params ``p``."""
    s = 2j * np.pi * np.asarray(f, float)
    lg = float(p["k0"]) * _ctrl(s, p) * np.exp(-s * float(p["delay_s"])) * _dynamics(s, p)
    total, chain = np.ones_like(s), np.ones_like(s)
    for st in p.get("stages") or ():
        chain = chain * _ctrl(s, st)
        total = total + chain * float(st["k_act"]) * _dynamics(s, st)
    return lg * total


def n_integrators(p: Mapping[str, Any]) -> int:
    """Pole order of L at s = 0 (module docstring)."""
    n = cum = int(float(p["ki_per_s"]) != 0.0)
    for st in p.get("stages") or ():
        cum += int(float(st["ki_per_s"]) != 0.0)
        if float(st["k_act"]) != 0.0:
            n = max(n, cum)
    return n


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
        _check_discretisation(params, name)
        for i, st in enumerate(params.get("stages") or ()):
            require(st, *STAGE_KEYS, block=f"{name}.stages[{i}]")
            _check_discretisation(st, f"{name}.stages[{i}]")
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
        n_int = n_integrators(p)
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
