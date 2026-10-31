"""``pid``: canonical PID with selectable discretisation and anti-windup (PLAN §2.10, §4.3 CTRL; T16 §4.3).

Canonical form (``hardware.mapping``, conv §9; error e = v − S, input minus setpoint)::

    u = clip( K_p e + K_i/(s + ω_IS) e + K_i2/s² e + K_d s/(1 + s/ω_DS) e + O,  u_min, u_max )

ω_IS = 0 means a pure integrator, ω_DS = 0 an unfiltered derivative, K_i2 = 0 no double integrator.
A device with another error sign (DSC1: e_DSC1 = S − v, DSC1 p.18) is handled by the canonical gains
(``hardware.mapping.dsc1_to_canonical`` flips the signs): the block itself always uses e = v − S.

Discretisation (structural ``discretisation``; h = rate_div / f_sim). Each first-order term
(b1 s + b0)/(s + a0) becomes y[n] = −A1 y[n−1] + B0 x[n] + B1 x[n−1] with

* ``dsc1`` (backward Euler, s → (1 − z⁻¹)/h): reproduces DSC1 p.19 exactly for ω_IS = ω_DS = K_i2 = 0:
  u = K_p e + h Σ K_i e + (K_d/h)(e[n] − e[n−1]) (bumpless: K_i inside the sum, p.18);
* ``tustin`` (s → (2/h)(1 − z⁻¹)/(1 + z⁻¹), no prewarp; same as ``hardware.mapping.canonical_to_fpga``;
  an unfiltered derivative has its pole at z = −1);
* ``forward_euler`` (s → (z − 1)/h; an unfiltered derivative falls back to the backward difference).

K_i2/s² is two cascaded 1/s sections. ``discretisation`` comes from the spec or ``device.pid``.

Anti-windup (structural ``anti_windup``; DSC1's behaviour is UNVERIFIED, T16 §3.6 → a model option):
``none`` (integrators run on, output clipped), ``clamp`` (conditional integration: integrator updates
that push further into the active limit are discarded), ``back_calculation`` (I += h/T_t·(u − u_unclipped),
runtime ``aw_tracking_s`` = T_t).

Params: ``in`` (error-input slot), ``out`` (V), ``device``, ``discretisation``, ``anti_windup``;
runtime ``kp``, ``ki_per_s``, ``kd_s``, ``omega_is_rad_per_s``, ``omega_ds_rad_per_s``, ``ki2_per_s2``
(= ``hardware.mapping.CANONICAL_KEYS``, the SimController contract), ``setpoint_v``, ``offset_v``,
``u_min_v``, ``u_max_v`` (CV limits) and ``aw_tracking_s`` (back_calculation only).
Toggles: ``feedback`` (False: loop open, u = clip(O), states frozen; SimController engage/disengage),
``p_term``, ``i_term`` (also gates K_i2, as DSC1 p.18), ``d_term``.
Slots: ``out``; ``ctrl.<name>.x`` = [e[n−1], I, D, A (1st integral of e), J (K_i2 term)];
``ctrl.<name>.sat`` (int: +1 / −1 at the upper / lower CV limit, else 0).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock, device_section
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.hardware.mapping import CANONICAL_KEYS

__all__ = ["DISCRETISATIONS", "ANTI_WINDUP", "Pid", "make_pid", "first_order"]

DISCRETISATIONS: dict[str, int] = {"dsc1": 0, "tustin": 1, "forward_euler": 2}
ANTI_WINDUP: dict[str, int] = {"none": 0, "clamp": 1, "back_calculation": 2}
_RUNTIME = (*CANONICAL_KEYS, "setpoint_v", "offset_v", "u_min_v", "u_max_v")


def first_order(method, b1, b0, a0, h):
    """(B0, B1, A1) of (b1 s + b0)/(s + a0) for ``method`` 0 backward / 1 Tustin / 2 forward Euler."""
    if method == 0:
        d = 1.0 + a0 * h
        return (b1 + b0 * h) / d, -b1 / d, -1.0 / d
    if method == 1:
        k = 2.0 / h
        d = k + a0
        return (b1 * k + b0) / d, (b0 - b1 * k) / d, (a0 - k) / d
    return b1, -b1 + b0 * h, -(1.0 - a0 * h)


def _make_step(fo):
    def step(S, I, P, i_in, i_out, i_sat, x0, po, h, method, aw, feedback, p_on, i_on, d_on):  # noqa: N803
        e = S[i_in] - P[po[6]]
        off, umin, umax = P[po[7]], P[po[8]], P[po[9]]
        if not feedback:
            u = off
            u = umin if u < umin else (umax if u > umax else u)
            S[i_out] = u
            I[i_sat] = 1 if off > umax else (-1 if off < umin else 0)
            return
        e1, i_old, d_old, a_old, j_old = S[x0], S[x0 + 1], S[x0 + 2], S[x0 + 3], S[x0 + 4]
        kp, ki, kd, wis, wds, ki2 = P[po[0]], P[po[1]], P[po[2]], P[po[3]], P[po[4]], P[po[5]]
        i_new = 0.0
        a_new = 0.0
        j_new = 0.0
        if i_on:
            b0, b1, a1 = fo(method, 0.0, ki, wis, h)
            i_new = -a1 * i_old + b0 * e + b1 * e1
            if ki2 != 0.0:
                b0, b1, a1 = fo(method, 0.0, 1.0, 0.0, h)
                a_new = -a1 * a_old + b0 * e + b1 * e1
                j_new = -a1 * j_old + ki2 * (b0 * a_new + b1 * a_old)
        d_new = 0.0
        if d_on:
            if wds > 0.0:
                b0, b1, a1 = fo(method, kd * wds, 0.0, wds, h)
                d_new = -a1 * d_old + b0 * e + b1 * e1
            elif method == 1:
                d_new = -d_old + kd * (2.0 / h) * (e - e1)
            else:
                d_new = kd * (e - e1) / h
        pterm = kp * e if p_on else 0.0
        u_raw = pterm + i_new + j_new + d_new + off
        u = umin if u_raw < umin else (umax if u_raw > umax else u_raw)
        if aw == 1 and u != u_raw and (u_raw - u) * ((i_new - i_old) + (j_new - j_old)) > 0.0:
            i_new, a_new, j_new = i_old, a_old, j_old
            u_raw = pterm + i_new + j_new + d_new + off
            u = umin if u_raw < umin else (umax if u_raw > umax else u_raw)
        elif aw == 2:
            i_new += h / P[po[10]] * (u - u_raw)
        S[x0], S[x0 + 1], S[x0 + 2], S[x0 + 3], S[x0 + 4] = e, i_new, d_new, a_new, j_new
        S[i_out] = u
        I[i_sat] = 1 if u_raw > umax else (-1 if u_raw < umin else 0)

    return step


_first_order_nb = numba.njit(inline="always")(first_order)


class Pid(KitBlock):
    KIND = "pid"
    STEP = staticmethod(_make_step(first_order))
    STEP_NB = numba.njit(_make_step(_first_order_nb))
    toggle_names = ("feedback", "p_term", "i_term", "d_term")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        dev = device_section(params, "pid", name)
        dev = {"discretisation": dev["discretisation"]} if "discretisation" in dev else {}
        require(params, "in", "out", "anti_windup", block=name)
        runtime = _RUNTIME + (("aw_tracking_s",) if params["anti_windup"] == "back_calculation" else ())
        self.x_name, self.sat_name = f"ctrl.{name}.x", f"ctrl.{name}.sat"
        super().__init__(params, name, reads=(params["in"],), writes=(params["out"], self.x_name, self.sat_name),
                         runtime=runtime, device=dev)
        self.need("discretisation")
        if self.values["discretisation"] not in DISCRETISATIONS:
            raise ValueError(f"block {name!r}: discretisation must be one of {sorted(DISCRETISATIONS)}, got "
                             f"{self.values['discretisation']!r} (Moku's discrete law is undocumented: use tustin)")
        if params["anti_windup"] not in ANTI_WINDUP:
            raise ValueError(f"block {name!r}: anti_windup must be one of {sorted(ANTI_WINDUP)}")
        self.slot_units = {params["out"]: "V"}

    def slots(self):
        return [(self.params["out"], (), np.float64, "V", 0, 0.0),
                (self.x_name, (5,), np.float64, "1", 0, 0.0),
                (self.sat_name, (), np.int64, "1", 0, 0)]

    def make_args(self, layout):
        keys = (*_RUNTIME, "aw_tracking_s")
        po = np.array([self.p(layout, k) if layout.has_param(self.param_name(k)) else -1 for k in keys], np.int64)
        t = self.toggles
        return (self.o(layout, self.params["in"]), self.o(layout, self.params["out"]), self.o(layout, self.sat_name),
                self.o(layout, self.x_name), po, self.h, DISCRETISATIONS[self.values["discretisation"]],
                ANTI_WINDUP[self.params["anti_windup"]], t["feedback"], t["p_term"], t["i_term"], t["d_term"])

    def freq_response(self, f_hz: Any, params) -> np.ndarray:
        """Discrete C(e^{iωh}) of the unclipped linear law (feedback on, current toggles and P values),
        Laplace convention z = e^{+i2πfh}; for analytic loop checks (VT-13 z-domain)."""
        g = {k: float(params[self.param_name(k)]) for k in CANONICAL_KEYS}
        m = DISCRETISATIONS[self.values["discretisation"]]
        h, t = self.h, self.toggles
        zi = np.exp(-2j * np.pi * np.atleast_1d(np.asarray(f_hz, dtype=float)) * h)

        def sec(b1, b0, a0):
            B0, B1, A1 = first_order(m, b1, b0, a0, h)
            return (B0 + B1 * zi) / (1.0 + A1 * zi)

        c = g["kp"] * np.ones_like(zi) if t["p_term"] else np.zeros_like(zi)
        if t["i_term"]:
            c = c + sec(0.0, g["ki_per_s"], g["omega_is_rad_per_s"])
            if g["ki2_per_s2"]:
                c = c + g["ki2_per_s2"] * sec(0.0, 1.0, 0.0) ** 2
        if t["d_term"]:
            w = g["omega_ds_rad_per_s"]
            if w > 0:
                c = c + sec(g["kd_s"] * w, 0.0, w)
            elif m == 1:
                c = c + g["kd_s"] * (2.0 / h) * (1 - zi) / (1 + zi)
            else:
                c = c + g["kd_s"] * (1 - zi) / h
        return c


@register_block("pid")
def make_pid(params: Mapping[str, Any], *, name: str | None = None) -> Pid:
    """Factory of the ``pid`` block (module docstring)."""
    return Pid(params, name)
