"""Coherent-probe error-signal generators #1, #2, #3, #6 (T14 §3–§5, §8). Owner: WP-B3.

All phases are beam phases at the named interaction point (conv §9), read from slots. Signal
representation, sensor noise and dither harmonics: ``blocks.err._common``. A± = amplitude gains of
the amplified / de-amplified quadrature (``amplitude_gains_opa`` / ``_opo`` / ``_from_table``);
g_c = (A+ + A−)/2, g_s = (A+ − A−)/2, Ḡ = (A+² + A−²)/2, ΔG = (A+² − A−²)/2 (T14 §1, §4).

``parametric_beat_detector`` (#1a phase-detection OPA, #1b conventional tap, #2 Ref B lock 1, Liu
scheme B detection OPA)
    Φ1 = 2φ_pr − φ_p − φ_OPA (pump–probe phase seen by the parametric process, conv §4, T14 §3.2).
    PD power P = η P_in [Ḡ + ΔG cos(·)]; v = G ℜ P.

    * ``mode: heterodyne`` (probe shifted by f_shift_hz, conv §4: A ∝ e^{i(φ_pr − 2π f_sh t)}):
      v = GℜηP_in [Ḡ + ΔG cos(2π·2f_sh·t − Φ1)], carrier 2·f_sh (T14 §3.2: demodulate at 2 f_sh).
      Baseband (a, b) = GℜηP_in ΔG (cos Φ1, sin Φ1), so with demod phase φ_d1:
      I = K cos(Φ1 + φ_d1), K = GℜηP_in sinh 2r for the OPA. Toggle ``f_sh_leak`` adds
      f_sh_leak_v·cos(2π f_sh t) in RF-rate mode only (non-ideal leak; outside the 2 f_sh band).
    * ``mode: dither`` (Ref B lock 1, DC extremum): Φ1 + β sin(2π f_m t), β = ``dither_m_rad`` is the
      dither amplitude **in Φ1** (2m for a probe dithered with amplitude m). Baseband b = −2GℜηP_in
      ΔG J1(β) sin Φ1 (T14 §4); zeros at Φ1 = 0 (max gain), π (max de-gain).
    * Conventional tap (#1b) = same block with P_in = T_tap·P_pr and the squeezer's A±; the tap loss
      on the squeezed path is the ``tap`` loss element of WP-B1 (not registered here, no double count).

    Shot noise S_I = 2eℜP_DC (instantaneous / slow DC power), electronic i_n²; S_v = G² S_I.

``bhd_probe_beat`` (#1c Ref A lock 2, #3 Ref B lock 2)
    u = φ_pr − (φ_p + φ_OPA)/2, v = φ_LO − (φ_p + φ_OPA)/2 = Δ + π/2 (T14 §3.4).
    BHD signal ∝ K2 [g_c cos(Ω t − (u − v)) + g_s cos(Ω t − (u + v))], K2 = G·2ℜ√(P_LO P_pr).

    * ``mode: heterodyne`` (carrier f_shift_hz): I = K2 [g_c cos(u′ − v) + g_s cos(u′ + v)],
      u′ = u + φ_d2; zero at Δ = 0 for u′ = 0; a lock-1 error δu′ moves it by
      tan ε = (A−/A+) tan δu′.
    * ``mode: dither`` (probe dithered, amplitude ``dither_m_rad`` in u, T14 §5).
    * ``mode: dc``: the DC interference K2 [g_c cos(u − v) + g_s cos(u + v)] (scalar).

    Noise: S_I = 2eℜP_LO·R(f_c) (``noise_rel_shot`` = quantum noise at the carrier relative to shot
    noise, includes squeezing) + i_n².

``liu_probe_interferometer`` (#6, T14 §8)
    X = φ_pr − φ_ref (counter-propagating probe of port k at the detector vs the forward reference).
    Amplitude K = s·G·2ℜ√(P_pk P_ref)·V_pr,ref with probe duty D: gated detection (toggle ``gating``)
    samples during the pulse, P_pk = P̄/D, noise PSD ×1/D; ungated averages, K ×D → 2ℜ√(D P̄ P_ref).
    Port sampling (f_rep/32, T14 §8.4) is the block's ``rate_div``.

    * ``mode: homodyne_dc``: v = K cos X + V_off (DC fringe form of #5 with a strong reference).
    * ``mode: heterodyne``: carrier f_shift_hz, (a, b) = K (cos X, sin X).
    * ``mode: dither``: amplitude ``dither_m_rad`` in X.

    Noise: S_I = 2eℜ(P_ref + P_pr,DC) (+ i_n²), ×1/D when gated. Same / different laser for probe
    and reference (D40) is a property of the beams (their phases carry the laser phase at their own
    emission times, F6); the block only takes the difference. Schemes A/B of the Liu ⊕ Ref-A
    integration (T14 §8.3, D42) are wirings of this block and ``parametric_beat_detector``.
    Leak information for T10 (``leak.*``) is written by ``probe_pulser`` (WP-B8), not here.
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba

from cranberry.core.block import require
from cranberry.blocks.err._common import E_CHARGE_C, ErrBlock, choice, helpers, register

__all__ = [
    "ParametricBeatDetector",
    "BhdProbeBeat",
    "LiuProbeInterferometer",
    "make_parametric_beat_detector",
    "make_bhd_probe_beat",
    "make_liu_probe_interferometer",
]

TWO_PI = 2.0 * math.pi
DITHER_KEYS = ("dither_m_rad", "f_dither_hz", "ram_a", "ram_phase_rad")
DET_KEYS = ("responsivity_a_per_w", "gain_v_per_a", "i_n_a_rthz")

_PBD_P = ("a_plus", "a_minus", "phi_opa_rad", "p_in_w", "eta", *DET_KEYS, "f_shift_hz", "f_sh_leak_v", *DITHER_KEYS)
_BHD_P = ("a_plus", "a_minus", "phi_opa_rad", "p_lo_w", "p_pr_w", "noise_rel_shot", *DET_KEYS, "f_shift_hz", *DITHER_KEYS)
_LIU_P = ("p_pr_avg_w", "p_ref_w", "duty", "visibility", "coupler_sign", "v_off_v", *DET_KEYS, "f_shift_hz", *DITHER_KEYS)
_COMMON_C = ("t", "buf", "idx", "fb", "mode", "law", "f_rf", "f_shot_noise", "f_electronic_noise", "f_ram_offset")
C_PBD = namedtuple("C_PBD", ("s_pr", "s_p", "s_out", *_COMMON_C, "f_f_sh_leak", *(f"p_{k}" for k in _PBD_P)))
C_BHD = namedtuple("C_BHD", ("s_pr", "s_p", "s_lo", "s_out", *_COMMON_C, *(f"p_{k}" for k in _BHD_P)))
C_LIU = namedtuple("C_LIU", ("s_pr", "s_ref", "s_out", *_COMMON_C, "f_gating", *(f"p_{k}" for k in _LIU_P)))

HET, DITHER, DC = 0, 1, 2


def _build(jit: bool):
    h = helpers(jit)
    harmonic, emit = h.harmonic, h.emit
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def ram(P, c, t):
        if c.f_ram_offset:
            return 1.0 + P[c.p_ram_a] * math.sin(TWO_PI * P[c.p_f_dither_hz] * t + P[c.p_ram_phase_rad])
        return 1.0

    @dec
    def ram_a(P, c):
        return P[c.p_ram_a] if c.f_ram_offset else 0.0

    @dec
    def det_coefs(P, c, popt):
        """(n0, n1) with S_v = n0 + n1·p for PD optical power popt·p (electronic, shot)."""
        g2 = P[c.p_gain_v_per_a] ** 2
        n0 = g2 * P[c.p_i_n_a_rthz] ** 2 if c.f_electronic_noise else 0.0
        n1 = g2 * 2.0 * E_CHARGE_C * P[c.p_responsivity_a_per_w] * popt if c.f_shot_noise else 0.0
        return n0, n1

    @dec
    def pbd_step(S, Z, I, P, c):  # noqa: N803
        ap, am = P[c.p_a_plus], P[c.p_a_minus]
        gbar, dg = 0.5 * (ap * ap + am * am), 0.5 * (ap * ap - am * am)
        phi1 = 2.0 * S[c.s_pr] - S[c.s_p] - P[c.p_phi_opa_rad]
        popt = P[c.p_eta] * P[c.p_p_in_w]
        kv = P[c.p_gain_v_per_a] * P[c.p_responsivity_a_per_w] * popt
        n0, n1 = det_coefs(P, c, popt)
        t = S[c.t]
        b = 0.0
        if c.mode == HET:
            sa = sb = n0 + n1 * gbar
            if c.f_rf:
                a = kv * (gbar + dg * math.cos(TWO_PI * 2.0 * P[c.p_f_shift_hz] * t - phi1))
                if c.f_f_sh_leak:
                    a += P[c.p_f_sh_leak_v] * math.cos(TWO_PI * P[c.p_f_shift_hz] * t)
            else:
                a, b = kv * dg * math.cos(phi1), kv * dg * math.sin(phi1)
        elif c.f_rf:
            p = (gbar + dg * math.cos(phi1)) * ram(P, c, t)
            a = kv * p
            sa = sb = n0 + n1 * p
        else:
            _, a, b, sa, sb = harmonic(0, kv, gbar, dg, phi1, 0.0, 0.0, P[c.p_dither_m_rad], ram_a(P, c),
                                       P[c.p_ram_phase_rad], n0, n1, 0.0, 0.0, 0.0, 0.0)
        emit(S, I, c, c.s_out, a, b, sa, sb)

    @dec
    def bhd_step(S, Z, I, P, c):  # noqa: N803
        ap, am = P[c.p_a_plus], P[c.p_a_minus]
        gc, gs = 0.5 * (ap + am), 0.5 * (ap - am)
        ref = 0.5 * (S[c.s_p] + P[c.p_phi_opa_rad])
        u, v = S[c.s_pr] - ref, S[c.s_lo] - ref
        k2 = P[c.p_gain_v_per_a] * 2.0 * P[c.p_responsivity_a_per_w] * math.sqrt(P[c.p_p_lo_w] * P[c.p_p_pr_w])
        t = S[c.t]
        b = 0.0
        if c.mode == HET:
            if c.f_rf:
                w = TWO_PI * P[c.p_f_shift_hz] * t
                a = k2 * (gc * math.cos(w - (u - v)) + gs * math.cos(w - (u + v)))
            else:
                a = k2 * (gc * math.cos(u - v) + gs * math.cos(u + v))
                b = k2 * (gc * math.sin(u - v) + gs * math.sin(u + v))
        elif c.mode == DITHER and not c.f_rf:
            _, a, b, _, _ = harmonic(0, k2, 0.0, gc, u - v, gs, u + v, P[c.p_dither_m_rad], ram_a(P, c),
                                     P[c.p_ram_phase_rad], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        else:  # dc, or dither at RF rate
            a = k2 * (gc * math.cos(u - v) + gs * math.cos(u + v))
            if c.mode == DITHER:
                a *= ram(P, c, t)
        s_v = 0.0
        if c.f_shot_noise:
            s_v += 2.0 * E_CHARGE_C * P[c.p_responsivity_a_per_w] * P[c.p_p_lo_w] * P[c.p_noise_rel_shot]
        if c.f_electronic_noise:
            s_v += P[c.p_i_n_a_rthz] ** 2
        s_v *= P[c.p_gain_v_per_a] ** 2  # stationary (LO shot noise): no cyclostationary weighting
        emit(S, I, c, c.s_out, a, b, s_v, s_v)

    @dec
    def liu_step(S, Z, I, P, c):  # noqa: N803
        x = S[c.s_pr] - S[c.s_ref]
        d, pbar, pref = P[c.p_duty], P[c.p_p_pr_avg_w], P[c.p_p_ref_w]
        if c.f_gating:
            p_eff, p_dc, nscale = pbar / d, pref + pbar / d, 1.0 / d
        else:
            p_eff, p_dc, nscale = d * pbar, pref + pbar, 1.0
        k = (P[c.p_coupler_sign] * P[c.p_gain_v_per_a] * 2.0 * P[c.p_responsivity_a_per_w]
             * math.sqrt(p_eff * pref) * P[c.p_visibility])
        t = S[c.t]
        b = 0.0
        if c.mode == HET:
            if c.f_rf:
                a = k * math.cos(TWO_PI * P[c.p_f_shift_hz] * t - x)
            else:
                a, b = k * math.cos(x), k * math.sin(x)
        elif c.mode == DITHER:
            if c.f_rf:
                a = k * math.cos(x) * ram(P, c, t)
            else:
                _, a, b, _, _ = harmonic(0, k, 0.0, 1.0, x, 0.0, 0.0, P[c.p_dither_m_rad], ram_a(P, c),
                                         P[c.p_ram_phase_rad], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        else:  # homodyne_dc
            a = k * math.cos(x) + P[c.p_v_off_v]
        n0, n1 = det_coefs(P, c, 1.0)
        s_v = (n0 + n1 * p_dc) * nscale  # stationary (reference shot noise dominates)
        emit(S, I, c, c.s_out, a, b, s_v, s_v)

    return pbd_step, bhd_step, liu_step


_PY = _build(False)
_NB = _build(True)


class _Coherent(ErrBlock):
    """Common constructor: mode, slots, per-mode runtime keys and toggles."""

    KIND = ""
    MODES: tuple[str, ...] = ()
    INPUTS: tuple[str, ...] = ()
    BASE_P: tuple[str, ...] = ()
    MODE_P: Mapping[str, tuple[str, ...]] = {}
    MODE_T: Mapping[str, tuple[str, ...]] = {}
    CARRIER_MODES: tuple[str, ...] = ("heterodyne", "dither")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        mode = choice(params, "mode", self.MODES, name)
        mname = self.MODES[mode]
        require(params, *self.INPUTS, "out", "toggles", block=name)
        carrier = mname in self.CARRIER_MODES
        toggles = ("shot_noise", "electronic_noise", *(("rf_rate_mode",) if carrier else ()), *self.MODE_T.get(mname, ()))
        tg = params["toggles"]
        rf = bool(tg.get("rf_rate_mode", True)) if carrier else True  # DC outputs are scalar
        noisy = bool(tg.get("shot_noise")) or bool(tg.get("electronic_noise"))
        self._rf = rf
        self.setup(
            self.KIND, name, params,
            inputs={r.removeprefix("in_"): (params[r], "rad") for r in self.INPUTS},
            outputs={"out": (params["out"], () if rf else (2,), "V")},
            runtime=(*self.BASE_P, *self.MODE_P.get(mname, ())),
            toggles=toggles,
            n_draws=(0 if not noisy else (1 if rf else 2)),
            mode={"heterodyne": HET, "dither": DITHER}.get(mname, DC),
        )

    @property
    def flags(self) -> dict[str, bool]:
        f = {k: v for k, v in self.toggles.items() if k != "rf_rate_mode"}
        f["rf"] = self._rf
        return f


class ParametricBeatDetector(_Coherent):
    KIND, MODES, CONST, STEP = "parametric_beat_detector", ("heterodyne", "dither"), C_PBD, (_PY[0], _NB[0])
    INPUTS = ("in_pr", "in_p")
    BASE_P = ("a_plus", "a_minus", "phi_opa_rad", "p_in_w", "eta", *DET_KEYS)
    MODE_P = {"heterodyne": ("f_shift_hz", "f_sh_leak_v"), "dither": DITHER_KEYS}
    MODE_T = {"heterodyne": ("f_sh_leak",), "dither": ("ram_offset",)}


class BhdProbeBeat(_Coherent):
    KIND, MODES, CONST, STEP = "bhd_probe_beat", ("heterodyne", "dither", "dc"), C_BHD, (_PY[1], _NB[1])
    INPUTS = ("in_pr", "in_p", "in_lo")
    BASE_P = ("a_plus", "a_minus", "phi_opa_rad", "p_lo_w", "p_pr_w", "noise_rel_shot", *DET_KEYS)
    MODE_P = {"heterodyne": ("f_shift_hz",), "dither": DITHER_KEYS}
    MODE_T = {"dither": ("ram_offset",)}


class LiuProbeInterferometer(_Coherent):
    KIND, MODES, CONST, STEP = "liu_probe_interferometer", ("heterodyne", "dither", "homodyne_dc"), C_LIU, (_PY[2], _NB[2])
    INPUTS = ("in_pr", "in_ref")
    BASE_P = ("p_pr_avg_w", "p_ref_w", "duty", "visibility", "coupler_sign", *DET_KEYS)
    MODE_P = {"heterodyne": ("f_shift_hz",), "dither": DITHER_KEYS, "homodyne_dc": ("v_off_v",)}
    MODE_T = {"heterodyne": ("gating",), "dither": ("gating", "ram_offset"), "homodyne_dc": ("gating",)}


make_parametric_beat_detector = register("parametric_beat_detector", ParametricBeatDetector)
make_bhd_probe_beat = register("bhd_probe_beat", BhdProbeBeat)
make_liu_probe_interferometer = register("liu_probe_interferometer", LiuProbeInterferometer)
