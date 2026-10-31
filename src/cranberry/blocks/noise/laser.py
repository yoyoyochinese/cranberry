"""``laser_source``: laser frequency noise, phase and RIN (T05 §2, §3.5, §4; T15; G16-d). Owner: WP-B2.

One block per laser; every beam of that source reads the same ``laser.<s>.phi`` through
``beam_assembler`` (emission time, fractional delay; F6 keeps the history), so inter-beam correlation
follows from the path delays (T05 §3.5 item 4).

Model (each term toggled independently; off ⇒ zero contribution; T05 §2.1)::

    δν(t) = δν_white + δν_flicker + Σ peaks_freq + ν̇·t + Δν_d sin(2π f_d t) + Δν_SBS sin(2π f_SBS t)
    φ_L(t) = −2π ∫ δν dt                                   (conv §4: E ∝ A e^{−iω0 t})
    δp(t)  = δp_white + Σ peaks_rin,   power = P0 (1 + δp)

* white FM h0 [Hz²/Hz] (one-sided, Δν = π h0, conv §7): per-step sample of std √(h0 f_b/2) (f_b = block
  rate); φ is its exact random walk (Var[φ(t) − φ(t − τ)] = 2π² h0 τ for τ = k/f_b).
* flicker FM h₋₁/f [Hz²]: OU sum (``_lgs.flicker``), corners from ``flicker_f_min_hz`` up to the block
  Nyquist frequency, ``flicker_per_decade`` OUs per decade (G11).
* technical peaks: band-pass sections (centre f, quality q, area [Hz²] for ``freq`` / [1] for ``rin``).
* deterministic drift / dither / SBS: closed-form phase (no accumulation error).
* RIN: in-loop white at the level of the first segment of ``rin_psd`` (piecewise-constant one-sided
  S_RIN [1/Hz], rows [f_upper_hz, s_per_hz] ascending); the first segment must cover the loop band.
  DX1 SC (D35, T05 §2.2): [[2e7, 1e−14], [1e9, 1e−15]] with D49's flat −140 dB/Hz below 100 kHz.

Band ledger (G05, F3 contract): the loop represents [0, f_N], f_N = f_b/2. For every entry of
``oob_pairs`` ({target, tau_s}: a beam pair whose laser phases differ by the delay τ) the block
declares a ``freq`` band on ``target`` with psd_outside = 4 sin²(πfτ)·S_ν(f) (the decorrelation folded
in; the ledger converts with S_φ = S_ν/f²). RIN above f_N is declared as ``rin`` on ``rin_targets``.
:func:`laser_oob_variance` integrates the same quantity accurately (offline table, T05 §3.5 item 3).

Filter-cavity pair {target, dnu_fc_hz, tau_fc_s, tau_ref_s} (T05 §6; WP-I1, corrected per D65 / M1-A-01):
one beam passes ``filter_cavity`` (``laser_lowpass``, FWHM Δν, fed with φ_L(t − tau_fc_s), tau_fc_s =
that beam's τ_b = the cavity's ``in_delay_s``), the other beam reads φ_L(t − tau_ref_s); both delays
default to 0. The white-FM random walk is an exact sample of the continuous Brownian phase, so the
loop's discrete chain (linear interpolation at both delays, the cavity's low pass) already contains
the aliased part above f_N; it differs from the continuous truth by O(h) terms of either sign (cavity
discretisation +, interpolation loss −, D57-7 / D65-4). The pair therefore declares, for white FM, a
``phase`` row with variance_outside = (continuous total − expected value of the discrete model)
(:func:`fc_pair_white_variances`), not the PSD above f_N. The other S_ν terms (flicker, peaks: band-
limited generators) keep the ``freq`` row |H_FC e^{−i2πfδ} − 1|²·S_ν above f_N, δ = tau_fc_s − tau_ref_s.

Slots: writes ``laser.<s>.phi`` [rad], ``.dnu`` [Hz], ``.dp`` [1], ``.power`` [W], ``.acc`` [rad]
(random-phase accumulator), plus the generator states / draw buffers under ``rng.<name>.*``.

Params (D34: required when the toggle using them is on)::

    source, power_w, oob_pairs, rin_targets, toggles{fm_white, fm_flicker, peaks, drift, dither, sbs, rin}
    fm_white: h0_hz2_per_hz        fm_flicker: h_m1_hz2, flicker_f_min_hz, flicker_per_decade
    peaks: peaks [{f_hz, q, area, quantity: freq|rin}]      drift: nu_dot_hz_per_s
    dither: dither_amp_hz, dither_f_hz                     sbs: sbs_amp_hz, sbs_f_hz (sinusoid, [AE])
    rin: rin_psd
"""

from __future__ import annotations

import math
import warnings
from collections import namedtuple
from collections.abc import Callable, Mapping
from typing import Any

import numba
import numpy as np
import scipy.signal

from cranberry.blocks.noise._lgs import EMPTY_BANK, Bank, bank_step, bank_step_nb, flicker, resonance
from cranberry.core.block import BaseBlock, Block, MissingParameterError, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer, substream
from cranberry.offline.band_ledger import BandDecl, integrate_psd

__all__ = ["LaserSource", "make_laser_source", "s_nu_model", "s_rin_model", "laser_oob_variance", "rin_cmrr_excess",
           "resonance_psd", "fc_pair_white_variances", "TOGGLES"]

TOGGLES = ("fm_white", "fm_flicker", "peaks", "drift", "dither", "sbs", "rin")
_NEED = {"fm_white": ("h0_hz2_per_hz",), "fm_flicker": ("h_m1_hz2", "flicker_f_min_hz", "flicker_per_decade"),
         "peaks": ("peaks",), "drift": ("nu_dot_hz_per_s",), "dither": ("dither_amp_hz", "dither_f_hz"),
         "sbs": ("sbs_amp_hz", "sbs_f_hz"), "rin": ("rin_psd",)}
_E = 1.602176634e-19
_H = 6.62607015e-34
_C = 299_792_458.0


def resonance_psd(f, f0_hz: float, q: float, area: float) -> np.ndarray:
    """One-sided PSD of ``_lgs.resonance`` (∫ df = area)."""
    f = np.asarray(f, dtype=float)
    k2 = 2.0 * area * q / (math.pi * f0_hz)
    g = f * f0_hz / q
    return k2 * g * g / ((f0_hz**2 - f * f) ** 2 + g * g)


def _peaks(params, quantity):
    return [p for p in params["peaks"] if p["quantity"] == quantity]


def s_nu_model(params: Mapping[str, Any]) -> Callable[[np.ndarray], np.ndarray] | None:
    """S_ν(f) [Hz²/Hz] of the stochastic FM terms switched on (None if none)."""
    tg = params["toggles"]
    h0 = float(params["h0_hz2_per_hz"]) if tg["fm_white"] else 0.0
    hm1 = float(params["h_m1_hz2"]) if tg["fm_flicker"] else 0.0
    pk = _peaks(params, "freq") if tg["peaks"] else []
    if h0 == 0.0 and hm1 == 0.0 and not pk:
        return None

    def s_nu(f):
        f = np.asarray(f, dtype=float)
        s = h0 + hm1 / f
        for p in pk:
            s = s + resonance_psd(f, float(p["f_hz"]), float(p["q"]), float(p["area"]))
        return s

    return s_nu


def s_rin_model(params: Mapping[str, Any]) -> Callable[[np.ndarray], np.ndarray] | None:
    """S_RIN(f) [1/Hz]: piecewise-constant ``rin_psd`` (+ rin peaks), None if RIN and rin peaks are off."""
    tg = params["toggles"]
    rows = [(float(a), float(b)) for a, b in params["rin_psd"]] if tg["rin"] else []
    pk = _peaks(params, "rin") if tg["peaks"] else []
    if not rows and not pk:
        return None
    edges = np.array([r[0] for r in rows])
    levels = np.array([r[1] for r in rows] + [0.0])

    def s_rin(f):
        f = np.asarray(f, dtype=float)
        s = levels[np.searchsorted(edges, f, side="left")] if rows else np.zeros_like(f)
        for p in pk:
            s = s + resonance_psd(f, float(p["f_hz"]), float(p["q"]), float(p["area"]))
        return s

    return s_rin


def _pair_psd(pair: Mapping[str, Any], s_nu: Callable[[np.ndarray], np.ndarray], block: str):
    """Decorrelated S_ν of one ``oob_pairs`` entry: ``tau_s`` → 4 sin²(πfτ)·S_ν (delay, T05 §3.1);
    ``dnu_fc_hz`` → |H_FC(f) e^{−i2πfδ} − 1|²·S_ν, H_FC = 1/(1 + 2if/Δν), δ = tau_fc_s − tau_ref_s (one
    beam through a filter cavity of FWHM Δν, the other not, T05 §6; module docstring)."""
    given = [k for k in ("tau_s", "dnu_fc_hz") if k in pair]
    if len(given) != 1:
        raise MissingParameterError(f"block {block!r}: each oob_pairs entry needs exactly one of tau_s / dnu_fc_hz, "
                                    f"got {given}")
    if given[0] == "tau_s":
        tau = float(pair["tau_s"])
        return lambda f: 4.0 * np.sin(np.pi * np.asarray(f) * tau) ** 2 * s_nu(f)
    dnu, tau_fc, tau_ref = _fc_pair(pair, block)
    delta = tau_fc - tau_ref

    def psd(f):
        f = np.asarray(f, dtype=float)
        h = np.exp(-2j * np.pi * f * delta) / (1.0 + 2j * f / dnu)
        return np.abs(h - 1.0) ** 2 * s_nu(f)

    return psd


def _fc_pair(pair: Mapping[str, Any], block: str) -> tuple[float, float, float]:
    dnu = float(pair["dnu_fc_hz"])
    tau_fc, tau_ref = float(pair.get("tau_fc_s", 0.0)), float(pair.get("tau_ref_s", 0.0))
    if not dnu > 0.0:
        raise ValueError(f"block {block!r}: oob_pairs dnu_fc_hz must be > 0")
    if tau_fc < 0.0 or tau_ref < 0.0:
        raise ValueError(f"block {block!r}: oob_pairs tau_fc_s / tau_ref_s must be >= 0")
    return dnu, tau_fc, tau_ref


def fc_pair_white_variances(h0: float, dnu_fc_hz: float, tau_fc_s: float, tau_ref_s: float,
                            h_s: float) -> tuple[float, float]:
    """(continuous, discrete) variance [rad²] of Δ = (H_FC φ_L)(t − τ_fc) − φ_L(t − τ_ref) for white FM h0.

    Continuous (all f): Var = 2π²h0 ∫K(t)² dt with K the response to a unit phase step, p = πΔν:
    δ = τ_fc − τ_ref ≥ 0 → δ + 1/(2p); δ = −u < 0 → u − 2(1 − e^{−pu})/p + 1/(2p) (δ = 0: πh0/Δν).
    Discrete (expected value of the loop's chain at step h, no lock): φ_n is an exact random walk with
    increments of variance 2π²h0h; LO = low pass (y_n = a y_{n−1} + (1 − a)(x_n + x_{n−1})/2, a = e^{−ph},
    ``filter_cavity``) of the linear interpolation at τ_fc, reference = linear interpolation at τ_ref
    (``core.delay`` / ``beam_assembler``). Var = 2π²h0h Σ K_n², K = step response; once both inputs are
    past, K_n decays geometrically with ratio a, so the tail is summed in closed form.
    """
    from cranberry.path.segments import frac_delay

    p = math.pi * dnu_fc_hz
    d = tau_fc_s - tau_ref_s
    if d >= 0.0:
        k_int = d + 0.5 / p
    else:
        u = -d
        k_int = u - 2.0 * (-math.expm1(-p * u)) / p + 0.5 / p
    v_cont = 2.0 * math.pi**2 * h0 * k_int

    a = math.exp(-p * h_s)
    kf, mf = frac_delay(tau_fc_s, h_s)
    kr, mr = frac_delay(tau_ref_s, h_s)
    n0 = max(kf, kr) + 3
    x = np.zeros(n0)
    r = np.zeros(n0)
    x[kf] += 1.0 - mf
    x[kf + 1] += mf
    r[kr] += 1.0 - mr
    r[kr + 1] += mr
    b = 0.5 * (1.0 - a)
    y = scipy.signal.lfilter([b, b], [1.0, -a], x)
    k = np.cumsum(y - r)
    v_disc = 2.0 * math.pi**2 * h0 * h_s * (float(np.sum(k * k)) + k[-1] ** 2 * a * a / (1.0 - a * a))
    return v_cont, v_disc


def laser_oob_variance(s_nu: Callable[[np.ndarray], np.ndarray], tau_s: float, f_lo_hz: float, f_hi_hz: float,
                       *, pts_per_period: int = 64, max_periods: float = 2e4) -> float:
    """σ² = ∫_{f_lo}^{f_hi} 4 sin²(πfτ) S_ν(f)/f² df [rad²] (T05 §3.5 item 3), resolving the sin² oscillation.

    Log grid below 0.1/τ (smooth), linear grid with ``pts_per_period`` points per period 1/τ up to
    ``max_periods`` periods, beyond that the period average 2·S_ν/f² (error ~ the oscillating residual).
    """
    if f_hi_hz <= f_lo_hz or tau_s == 0.0:
        return 0.0
    tau = abs(tau_s)
    w = lambda f: 4.0 * np.sin(np.pi * f * tau) ** 2  # noqa: E731
    f1 = min(f_hi_hz, max(f_lo_hz, 0.1 / tau))
    total = integrate_psd(s_nu, f_lo_hz, f1, weight=w, divide_f2=True, n_per_decade=2000) if f1 > f_lo_hz else 0.0
    f2 = min(f_hi_hz, f1 + max_periods / tau)
    if f2 > f1:
        n = int(math.ceil((f2 - f1) * tau * pts_per_period)) + 1
        f = np.linspace(f1, f2, n)
        total += float(np.trapezoid(w(f) * np.asarray(s_nu(f), dtype=float) / (f * f), f))
    if f_hi_hz > f2:
        total += integrate_psd(lambda f: 2.0 * np.asarray(s_nu(f), dtype=float), f2, f_hi_hz, divide_f2=True)
    return total


def rin_cmrr_excess(p_lo_w: float, eta: float, s_rin_per_hz: float, cm_leak: float, wavelength_m: float) -> float:
    """Excess noise / shot noise of LO RIN leaking through a balanced detector: ε² ℜ P_LO S_RIN / (2e),
    ℜ = η e λ/(h c) (T05 §4.1; ε = common-mode leakage, CMRR = −20 log10 ε)."""
    resp = eta * _E * wavelength_m / (_H * _C)
    return cm_leak**2 * resp * p_lo_w * s_rin_per_hz / (2.0 * _E)


LC = namedtuple("LC", "t o_phi o_dnu o_dp o_pow o_acc f_white f_flicker f_peaks f_drift f_dither f_sbs f_rin "
                      "sig_w sig_rin nu_dot da df sa sf p0 h wbuf widx rbuf ridx bfl tfl bpk tpk")


def _build(jit: bool):
    bs = bank_step_nb if jit else bank_step

    def step(S, I, c):  # noqa: N803
        t = S[c.t]
        S[c.o_dnu] = 0.0
        S[c.o_dp] = 0.0
        if c.f_white:
            S[c.o_dnu] += c.sig_w * S[c.wbuf + I[c.widx]]
            I[c.widx] += 1
        if c.f_flicker:
            bs(S, I, c.bfl, c.tfl)
        if c.f_peaks:
            bs(S, I, c.bpk, c.tpk)
        if c.f_rin:
            S[c.o_dp] += c.sig_rin * S[c.rbuf + I[c.ridx]]
            I[c.ridx] += 1
        rnd = S[c.o_dnu]
        phi_det = 0.0
        if c.f_drift:
            S[c.o_dnu] += c.nu_dot * t
            phi_det -= math.pi * c.nu_dot * t * t
        if c.f_dither:
            S[c.o_dnu] += c.da * math.sin(2.0 * math.pi * c.df * t)
            phi_det -= c.da * (1.0 - math.cos(2.0 * math.pi * c.df * t)) / c.df
        if c.f_sbs:
            S[c.o_dnu] += c.sa * math.sin(2.0 * math.pi * c.sf * t)
            phi_det -= c.sa * (1.0 - math.cos(2.0 * math.pi * c.sf * t)) / c.sf
        S[c.o_phi] = S[c.o_acc] + phi_det
        S[c.o_acc] -= 2.0 * math.pi * c.h * rnd
        S[c.o_pow] = c.p0 * (1.0 + S[c.o_dp])

    return numba.njit(step) if jit else step


_STEP_PY, _STEP_NB = _build(False), _build(True)


class LaserSource(BaseBlock):
    """The ``laser_source`` block (module docstring)."""

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or "laser_source"
        require(params, "source", "power_w", "oob_pairs", "rin_targets", "toggles", block=name)
        for t, keys in _NEED.items():
            if params["toggles"].get(t):
                require(params, *keys, block=name)
        s = params["source"]
        self.o = {k: f"laser.{s}.{k}" for k in ("phi", "dnu", "dp", "power", "acc")}
        super().__init__("laser_source", name, params, reads=("t",), writes=tuple(self.o.values()))
        units = {"phi": "rad", "dnu": "Hz", "dp": "1", "power": "W", "acc": "rad"}
        self.slot_units = {self.o[k]: u for k, u in units.items()}
        self._fb: float | None = None
        self.banks: dict[str, Bank] = {}
        self.white: dict[str, DrawBuffer] = {}

    def init(self, layout, rng):
        super().init(layout, rng)
        tg, p = self.toggles, self.params
        h = block_dt_s(layout.f_sim_hz, self.rate_div)
        self._fb, f_n = 1.0 / h, 0.5 / h
        for slot, unit in self.slot_units.items():
            layout.add(slot, (), np.float64, unit, init=float(p["power_w"]) if slot == self.o["power"] else 0.0)
        for key in ("dither", "sbs"):
            if tg[key] and not 0.0 < float(p[f"{key}_f_hz"]) < f_n:
                raise ValueError(f"block {self.name!r}: {key}_f_hz must be in (0, f_N = {f_n} Hz)")
        if tg["rin"] and float(p["rin_psd"][0][0]) < f_n:
            raise ValueError(f"block {self.name!r}: first rin_psd segment must cover the loop band [0, {f_n}] Hz")
        if tg["fm_flicker"]:
            b = self.banks["flicker"] = Bank(self.name, "flicker")
            for sec in flicker(float(p["h_m1_hz2"]), float(p["flicker_f_min_hz"]), f_n, int(p["flicker_per_decade"])):
                b.add(sec, [(self.o["dnu"], 0, 1.0)])
        if tg["peaks"] and p["peaks"]:
            b = self.banks["peaks"] = Bank(self.name, "peaks")
            for pk in p["peaks"]:
                if not 0.0 < float(pk["f_hz"]) < f_n or pk["quantity"] not in ("freq", "rin"):
                    raise ValueError(f"block {self.name!r}: peak {pk} needs 0 < f_hz < f_N and quantity freq|rin")
                b.add(resonance(float(pk["f_hz"]), float(pk["q"]), float(pk["area"])),
                      [(self.o["dnu" if pk["quantity"] == "freq" else "dp"], 0, 1.0)])
        extra: list[str] = []
        for b in self.banks.values():
            extra += b.init(layout, rng, self.rate_div, h)
        for key, on in (("white", tg["fm_white"]), ("rin", tg["rin"])):
            if on:
                d = self.white[key] = DrawBuffer(self.name, 1, "standard_normal", label=key)
                d.init(layout, substream(rng, key), self.rate_div)
                extra += d.slot_names
        self.writes = (*self.writes, *extra)

    def bind(self, layout):
        tg, p = self.toggles, self.params
        for d in self.white.values():
            d.bind(layout)
        off = lambda s: layout.slot(s).offset  # noqa: E731
        bank = {k: (b.bind(layout), b.scratch_off(layout)) for k, b in self.banks.items()}
        fl, pk = bank.get("flicker", (EMPTY_BANK, -1)), bank.get("peaks", (EMPTY_BANK, -1))
        w, r = self.white.get("white"), self.white.get("rin")
        g = lambda k: float(p[k]) if k in p else 0.0  # noqa: E731
        fb = self._fb
        self._c = LC(
            off("t"), *(off(self.o[k]) for k in ("phi", "dnu", "dp", "power", "acc")),
            tg["fm_white"], tg["fm_flicker"], bool(tg["peaks"] and p["peaks"]), tg["drift"], tg["dither"],
            tg["sbs"], tg["rin"],
            math.sqrt(g("h0_hz2_per_hz") * fb / 2.0) if tg["fm_white"] else 0.0,
            math.sqrt(float(p["rin_psd"][0][1]) * fb / 2.0) if tg["rin"] else 0.0,
            g("nu_dot_hz_per_s"), g("dither_amp_hz"), g("dither_f_hz") or 1.0, g("sbs_amp_hz"), g("sbs_f_hz") or 1.0,
            float(p["power_w"]), 1.0 / fb,
            w.buf_off if w else -1, w.idx_off if w else -1, r.buf_off if r else -1, r.idx_off if r else -1,
            fl[0], fl[1], pk[0], pk[1],
        )
        c, nb = self._c, _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        for b in self.banks.values():
            b.refill(state, n0, n_samples)
        for d in self.white.values():
            d.refill(state, n0, n_samples)

    def __call__(self, state, params):
        _STEP_PY(state.S, state.I, self._c)
        return state

    def bands(self) -> list[BandDecl]:
        if self._fb is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        f_n = 0.5 * self._fb
        out = []
        p, tg = self.params, self.toggles
        s_nu = s_nu_model(p)
        # S_ν without white FM: the band-limited generators (flicker, peaks) of the FC pairs' freq row
        s_nw = s_nu_model({**p, "toggles": {**tg, "fm_white": False}})
        if s_nu is not None:
            for pair in p["oob_pairs"]:
                tgt = (str(pair["target"]),)
                if "dnu_fc_hz" not in pair:
                    out.append(BandDecl(source=self.name, quantity="freq", targets=tgt, f_lo_hz=0.0, f_hi_hz=f_n,
                                        unit="Hz", psd_outside=_pair_psd(pair, s_nu, self.name)))
                    continue
                if s_nw is not None:
                    out.append(BandDecl(source=self.name, quantity="freq", targets=tgt, f_lo_hz=0.0, f_hi_hz=f_n,
                                        unit="Hz", psd_outside=_pair_psd(pair, s_nw, self.name)))
                if tg["fm_white"]:
                    dnu, tau_fc, tau_ref = _fc_pair(pair, self.name)
                    v_cont, v_disc = fc_pair_white_variances(float(p["h0_hz2_per_hz"]), dnu, tau_fc, tau_ref,
                                                             1.0 / self._fb)
                    corr = v_cont - v_disc
                    if corr < 0.0:
                        # ponytail: the ledger has no negative rows (phase_to_quadrature needs σ²_HF ≥ 0), so a
                        # discrete chain that over-represents the continuous variance (cavity discretisation
                        # ≈ ε²/12 of πh0/Δν, ε = πΔν h, when the delays do not dominate) is left as is; upgrade:
                        # signed correction rows in BandLedger + phase_to_quadrature.
                        warnings.warn(f"block {self.name!r}: FC pair on {tgt[0]!r}: the discrete loop exceeds the "
                                      f"continuous white-FM variance by {-corr:.3e} rad² ({-corr / v_cont:.2%}); "
                                      "ledger row set to 0", stacklevel=2)
                        corr = 0.0
                    out.append(BandDecl(source=self.name, quantity="phase", targets=tgt, f_lo_hz=0.0, f_hi_hz=f_n,
                                        unit="rad", variance_outside=corr))
        s_rin = s_rin_model(p)
        if s_rin is not None and p["rin_targets"]:
            out.append(BandDecl(source=self.name, quantity="rin", targets=tuple(p["rin_targets"]),
                                f_lo_hz=0.0, f_hi_hz=f_n, unit="1", psd_outside=s_rin))
        return out


@register_block("laser_source")
def make_laser_source(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``laser_source`` block (module docstring)."""
    return LaserSource(params, name)
