"""``dc_fringe``: DC fringe detector #5 (T15 eq. (1)–(7) superset, T14 §7, D06, G16-d). Owner: WP-B3.

Two arms meet at coupler 2 (power ratio κ_2) and a balanced detector (PDB450C for Ref D). With
φ_D = φ_LO − φ_pr (conv §5 Ref D mapping: arm without 915B = "LO" = ``in_lo``, arm with 915B =
"pr" = ``in_pr``; T15 §2.2: φ_D = φ_2 − φ_1)::

    P_1 = κ_1 t_1 P_in,  P_2 = (1 − κ_1) t_2 P_in,  P_det = P_1 + P_2
    v = V_0 [ ε + s·𝒱 |μ| cos(φ_D + φ_μ) ] + V_off + n(t)                                  (T15 (3))
    V_0 = G ℜ P_det (1 + δp)
    ε   = s (2κ_2 − 1)(P_1 − P_2)/P_det + ε_R                       (toggle optical_imbalance)
    𝒱   = 4√(κ_2(1 − κ_2) P_1 P_2)/P_det
    μ   = |μ_0| e^{iφ_μ0} · γ_coh · [(1 − x) + x e^{iΔχ}]           (T15 (2), (7); γ_coh, PER toggles)

s = ``coupler_sign`` ∈ {+1, −1}: which coupler output goes to INPUT+ (conv §5); swapping the outputs
flips the coupler terms but not the PD responsivity imbalance ε_R. γ_coh = exp(−σ²_HF/2) is the
visibility loss from phase noise above the modelled band, taken from the ``BandLedger``
(:func:`gamma_coh_from_ledger`); the modelled band must not be put into γ_coh again (T15 §2.6).
Δχ: differential birefringent phase (runtime ``per_dchi_rad``). Optional slot inputs ``in_per_x`` /
``in_per_dchi`` (WP-I2 addition, B7 request): when both are given (not null) and ``per_crosscoupling``
is on, x and Δχ are read from these slots (e.g. ``per_fringe`` out_x / out_dchi, computed from the
arm Jones matrices) instead of the runtime params, and they are declared in ``reads`` only then. Put
``mu_abs = 1``, ``mu_phase_rad = 0`` in that case (per_fringe docstring).

Noise: shot noise of both PDs S_I = 2eℜP_det(1 + δp), electronic i_n²; RIN enters through
``in_dp`` (relative power δp, toggle ``rin``), offset drift through ``in_voff`` (toggle
``offset_drift``). Saturation clips ±``v_sat_v`` after the noise (toggle ``saturation``).

Outputs: ``out`` = v (the PV that e.g. the DSC1 reads; its e = S − v convention lives in the DSC1
block, conv §9) and, if ``out_err`` is a slot, e_v = loop_sign·(S − v) with S = ``setpoint_v``: for
loop_sign = +1, s = +1 the slope de_v/dφ_D is positive at the +π/2 lock point (conv §9 sign), so a
controller with positive gain gives negative feedback there; flipping loop_sign makes +π/2 unstable
and −π/2 stable (T15 V11).
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba

from cranberry.blocks.err._common import E_CHARGE_C, ErrBlock, helpers, register
from cranberry.core.block import require

__all__ = ["DcFringe", "make_dc_fringe", "fringe_lock_point", "gamma_coh_from_ledger"]

_P = ("p_in_w", "kappa1", "kappa2", "t_arm1", "t_arm2", "responsivity_a_per_w", "eps_r", "gain_v_per_a",
      "v_off_v", "mu_abs", "mu_phase_rad", "gamma_coh", "per_x", "per_dchi_rad", "coupler_sign",
      "i_n_a_rthz", "v_sat_v", "setpoint_v", "loop_sign")
_T = ("shot_noise", "electronic_noise", "rin", "offset_drift", "optical_imbalance", "per_crosscoupling",
      "visibility_hf_averaging", "saturation")
C_DCF = namedtuple("C_DCF", ("s_lo", "s_pr", "s_dp", "s_voff", "s_perx", "s_perd", "s_out", "s_err", "t", "buf",
                             "idx", "fb", "mode",
                             "law", *(f"f_{k}" for k in _T), *(f"p_{k}" for k in _P)))


def _build(jit: bool):
    draw = helpers(jit).draw
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def step(S, Z, I, P, c):  # noqa: N803
        p_in, k1, k2 = P[c.p_p_in_w], P[c.p_kappa1], P[c.p_kappa2]
        p1, p2 = k1 * P[c.p_t_arm1] * p_in, (1.0 - k1) * P[c.p_t_arm2] * p_in
        pdet = p1 + p2
        sgn = P[c.p_coupler_sign]
        dp = S[c.s_dp] if c.f_rin else 0.0
        voff = P[c.p_v_off_v] + (S[c.s_voff] if c.f_offset_drift else 0.0)
        eps = sgn * (2.0 * k2 - 1.0) * (p1 - p2) / pdet + P[c.p_eps_r] if c.f_optical_imbalance else 0.0
        vis = 4.0 * math.sqrt(k2 * (1.0 - k2) * p1 * p2) / pdet
        mu_re = P[c.p_mu_abs] * math.cos(P[c.p_mu_phase_rad])
        mu_im = P[c.p_mu_abs] * math.sin(P[c.p_mu_phase_rad])
        if c.f_visibility_hf_averaging:
            mu_re *= P[c.p_gamma_coh]
            mu_im *= P[c.p_gamma_coh]
        if c.f_per_crosscoupling:
            x = S[c.s_perx] if c.s_perx >= 0 else P[c.p_per_x]
            dchi = S[c.s_perd] if c.s_perd >= 0 else P[c.p_per_dchi_rad]
            pr, pi = 1.0 - x + x * math.cos(dchi), x * math.sin(dchi)
            mu_re, mu_im = mu_re * pr - mu_im * pi, mu_re * pi + mu_im * pr
        phid = S[c.s_lo] - S[c.s_pr]
        v0 = P[c.p_gain_v_per_a] * P[c.p_responsivity_a_per_w] * pdet * (1.0 + dp)
        # |μ| cos(φ_D + φ_μ) = Re(μ e^{iφ_D})
        v = v0 * (eps + sgn * vis * (mu_re * math.cos(phid) - mu_im * math.sin(phid))) + voff
        if c.buf >= 0:
            s_i = 0.0
            if c.f_shot_noise:
                s_i += 2.0 * E_CHARGE_C * P[c.p_responsivity_a_per_w] * pdet * (1.0 + dp)
            if c.f_electronic_noise:
                s_i += P[c.p_i_n_a_rthz] ** 2
            v += P[c.p_gain_v_per_a] * math.sqrt(s_i * c.fb * 0.5) * draw(S, I, c.buf, c.idx)
        if c.f_saturation:
            v = min(max(v, -P[c.p_v_sat_v]), P[c.p_v_sat_v])
        S[c.s_out] = v
        if c.s_err >= 0:
            S[c.s_err] = P[c.p_loop_sign] * (P[c.p_setpoint_v] - v)

    return step


class DcFringe(ErrBlock):
    KIND, CONST, STEP = "dc_fringe", C_DCF, (_build(False), _build(True))

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in_lo", "in_pr", "in_dp", "in_voff", "out", "out_err", "toggles", block=name)
        tg = params["toggles"]
        for tog, key in (("rin", "in_dp"), ("offset_drift", "in_voff")):
            if tg.get(tog) and params[key] is None:
                raise ValueError(f"block {name!r}: toggle {tog!r} is on but {key} is null")
        per_on = bool(tg.get("per_crosscoupling"))
        per_slots = per_on and params.get("in_per_x") is not None
        if per_slots != (per_on and params.get("in_per_dchi") is not None):
            raise ValueError(f"block {name!r}: give both in_per_x and in_per_dchi (or neither)")
        self.setup(
            self.KIND, name, params,
            inputs={"lo": (params["in_lo"], "rad"), "pr": (params["in_pr"], "rad"),
                    "dp": (params["in_dp"] if tg.get("rin") else None, "1"),
                    "voff": (params["in_voff"] if tg.get("offset_drift") else None, "V"),
                    "perx": (params["in_per_x"] if per_slots else None, "1"),
                    "perd": (params["in_per_dchi"] if per_slots else None, "rad")},
            outputs={"out": (params["out"], (), "V"), "err": (params["out_err"], (), "V")},
            runtime=_P, toggles=_T,
            n_draws=int(bool(tg.get("shot_noise")) or bool(tg.get("electronic_noise"))),
        )


make_dc_fringe = register("dc_fringe", DcFringe)


def fringe_lock_point(v0: float, eps: float, vis: float, setpoint_v: float, v_off_v: float,
                      phi_mu_rad: float = 0.0) -> tuple[float, float]:
    """Static lock point and slope on the +π/2 branch, T15 eq. (4)–(5) (s = +1, |μ| folded into ``vis``).

    Returns (φ_lock [rad], K_d = ∂v/∂φ_D [V/rad]). Raises if the setpoint is outside the fringe.
    """
    c_lock = (setpoint_v - v_off_v) / (v0 * vis) - eps / vis
    if abs(c_lock) >= 1.0:
        raise ValueError(f"no lock point: |c_lock| = {abs(c_lock)} >= 1")
    return math.acos(c_lock) - phi_mu_rad, -v0 * vis * math.sqrt(1.0 - c_lock * c_lock)


def gamma_coh_from_ledger(ledger: Any, metric: str, beams: Any) -> float:
    """γ_coh = exp(−σ²_HF/2) with σ²_HF the out-of-band phase variance of the beam combination (T15 §2.6, G05)."""
    return math.exp(-0.5 * ledger.hf_variance(metric, beams))
