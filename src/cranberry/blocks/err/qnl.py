"""``quantum_noise_power``: quantum noise lock #4 (Ref C, T14 §6). Owner: WP-B3.

Δ = φ_LO − θ_sq, θ_sq = (φ_p + φ_OPA)/2 + π/2 (conv §4–§5). Band-integrated BHD noise power over
[f1, f2] (offline, :func:`qnl_band_integrals`)::

    P(Δ) = p_mean − p_delta cos 2Δ                       (P_∓ = P(0), P(π/2); p_delta = ΔP)
    Q(Δ) = ∫ S_x(f; Δ)² df = q0 + q1 cos 2Δ + q2 cos² 2Δ   (flat band: Q = P²/B)

Detector (``detector_law``, structural): square v = k P, linear envelope v = k √P, log
v = k·10 log10 P. LO dither Δ0 + m sin(2π f_m t) (``dither_m_rad`` = m, the amplitude in Δ); the
baseband first harmonic is computed exactly (``_common`` quadrature); square law: b = 2kΔP J1(2m)
sin 2Δ0 (T14 §6.2). Optional LO RAM (toggle ``ram_offset``): the LO-proportional part only,
P → (P − p_dark)·(1 + a sin(2π f_m t + φ_a)) + p_dark, p_dark = ∫|H|²S_el df (the BHD dark / circuit
noise is not modulated by the LO's RAM, T14 §6.4; D65-1, M1-A-03: Ref C shift 2.54 mrad at a = 10⁻²).

Radiometer sensor noise (toggle ``radiometer_noise``, T14 §6.3), equivalent Gaussian model: the
detector output has low-frequency one-sided PSD S_v = κ_det · v′(P)² · 2Q(Δ) (square law:
S_P̂ = 2∫S_x² df; flat band: Var(P̂)/P² = 1/(BT)), κ_det = detector-law factor (T14 MC: square 0.97,
linear 1.02, log 1.32; validated by :func:`radiometer_mc`). Valid for f_m ≪ B, loop bandwidth ≪ f_m,
video bandwidth ≫ f_m, BT ≫ 1. The noise is cyclostationary over the dither cycle; the baseband mode
uses the demodulation-weighted cycle average (``_common``), equal to the RF-rate mode. The weight is
law-dependent (D65-1, M1-A-04): Ref C at the squeezing point 1.2257 for the square law (3.58 mrad,
continuous loop b_loop 10 kHz), exactly 1 for the log law, whose slope is 0.905 × the stationary one
instead (4.11 mrad with κ_det = 1.32). Ponytail: the noise is white up to the block Nyquist frequency
(no video-filter shape).

RF-rate Monte Carlo (cross-validation, not a time-loop mode): :func:`radiometer_mc` draws
band-limited Gaussian BHD noise at f_s ≥ 1 GHz, applies the detector law and a boxcar of length T,
and returns the measured κ_det.
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.err._common import DETECTOR_LAWS, ErrBlock, choice, helpers, register
from cranberry.core.block import require
from cranberry.core.fft import fft_agrawal, ifft_agrawal, omega_grid

__all__ = ["QuantumNoisePower", "make_quantum_noise_power", "qnl_band_integrals", "radiometer_mc"]

TWO_PI = 2.0 * math.pi
_P = ("phi_opa_rad", "p_mean", "p_delta", "p_dark", "q0", "q1", "q2", "k_det", "kappa_det", "dither_m_rad",
      "f_dither_hz", "ram_a", "ram_phase_rad")
C_QNL = namedtuple("C_QNL", ("s_lo", "s_p", "s_out", "t", "buf", "idx", "fb", "mode", "law", "f_rf",
                             "f_radiometer_noise", "f_ram_offset", *(f"p_{k}" for k in _P)))


def _build(jit: bool):
    h = helpers(jit)
    law, noise, harmonic, emit = h.law, h.noise, h.harmonic, h.emit
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def step(S, Z, I, P, c):  # noqa: N803
        delta = S[c.s_lo] - 0.5 * (S[c.s_p] + P[c.p_phi_opa_rad]) - 0.5 * math.pi
        c2 = math.cos(2.0 * delta)
        pw = P[c.p_p_mean] - P[c.p_p_delta] * c2
        k = P[c.p_k_det]
        ra = P[c.p_ram_a] if c.f_ram_offset else 0.0
        kq = 2.0 * P[c.p_kappa_det] if c.f_radiometer_noise else 0.0  # S_v = κ v′(P)² · 2Q
        nq0, nq1, nq2 = kq * P[c.p_q0], kq * P[c.p_q1], kq * P[c.p_q2]
        b = 0.0
        dark = P[c.p_p_dark]  # not modulated by the LO RAM (M1-A-03)
        if c.f_rf:
            p = (pw - dark) * (1.0 + ra * math.sin(TWO_PI * P[c.p_f_dither_hz] * S[c.t] + P[c.p_ram_phase_rad])) + dark
            a = law(c.law, k, p)
            sa = sb = noise(c.law, k, p, c2, 0.0, 0.0, nq0, nq1, nq2)
        else:
            _, a, b, sa, sb = harmonic(c.law, k, P[c.p_p_mean] - dark, -P[c.p_p_delta], 2.0 * delta, 0.0, 0.0,
                                       2.0 * P[c.p_dither_m_rad], ra, P[c.p_ram_phase_rad], 0.0, 0.0, nq0, nq1, nq2,
                                       dark)
        emit(S, I, c, c.s_out, a, b, sa, sb)

    return step


class QuantumNoisePower(ErrBlock):
    KIND, CONST, STEP = "quantum_noise_power", C_QNL, (_build(False), _build(True))

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        law_code = choice(params, "detector_law", DETECTOR_LAWS, name)
        require(params, "in_lo", "in_p", "out", "toggles", block=name)
        tg = params["toggles"]
        self._rf = bool(tg.get("rf_rate_mode"))
        self.setup(
            self.KIND, name, params,
            inputs={"lo": (params["in_lo"], "rad"), "p": (params["in_p"], "rad")},
            outputs={"out": (params["out"], () if self._rf else (2,), "V")},
            runtime=_P, toggles=("radiometer_noise", "ram_offset", "rf_rate_mode"),
            n_draws=(0 if not tg.get("radiometer_noise") else (1 if self._rf else 2)),
            law=law_code,
        )

    @property
    def flags(self) -> dict[str, bool]:
        f = {k: v for k, v in self.toggles.items() if k != "rf_rate_mode"}
        f["rf"] = self._rf
        return f


make_quantum_noise_power = register("quantum_noise_power", QuantumNoisePower)


def qnl_band_integrals(f_hz: np.ndarray, s_sn: np.ndarray, r_minus: np.ndarray, r_plus: np.ndarray,
                       s_el: np.ndarray, h2: np.ndarray) -> dict[str, float]:
    """Runtime parameters (p_mean, p_delta, p_dark, q0, q1, q2) of ``quantum_noise_power`` on a band grid.

    S_x(f; Δ) = |H|² [S_SN (R_− cos²Δ + R_+ sin²Δ) + S_el] = A(f) + C(f) cos 2Δ with
    A = |H|²[S_SN (R_+ + R_−)/2 + S_el], C = −|H|² S_SN (R_+ − R_−)/2 (T14 §6.1); trapezoid rule
    over the given grid (it spans [f1, f2]). p_dark = ∫|H|² S_el df is the part of p_mean that the LO RAM
    does not modulate.
    """
    f = np.asarray(f_hz, dtype=float)
    a = h2 * (s_sn * (r_plus + r_minus) / 2.0 + s_el)
    c = -h2 * s_sn * (r_plus - r_minus) / 2.0
    integ = lambda y: float(np.trapezoid(y, f))  # noqa: E731
    return {"p_mean": integ(a), "p_delta": -integ(c), "p_dark": integ(h2 * s_el), "q0": integ(a * a),
            "q1": integ(2.0 * a * c), "q2": integ(c * c)}


def radiometer_mc(rng: np.random.Generator, *, f1_hz: float, f2_hz: float, fs_hz: float, t_avg_s: float,
                  n_blocks: int, law: str, batch_blocks: int = 512) -> float:
    """RF-rate Monte Carlo of the radiometer noise; returns κ_det = (measured relative variance)·B·T.

    Band-limited Gaussian noise (ideal flat band [f1, f2], conv §8 transforms) sampled at ``fs_hz``;
    per block of length T: square law P̂ = ⟨x²⟩ (statistic Var(P̂)/⟨P̂⟩²), linear envelope
    E = ⟨|z|⟩ (4 Var(E)/⟨E⟩²), log L = ⟨ln |z|²⟩ (Var(L)); z = analytic signal (T14 §6.3).
    """
    if law not in DETECTOR_LAWS:
        raise ValueError(f"law must be one of {DETECTOR_LAWS}")
    m = int(round(t_avg_s * fs_hz))
    dt = 1.0 / fs_hz
    stats: list[np.ndarray] = []
    done = 0
    while done < n_blocks:
        nb = min(batch_blocks, n_blocks - done)
        n = nb * m
        om = omega_grid(n, dt)
        spec = fft_agrawal(rng.standard_normal(n).astype(complex), dt)
        band = (np.abs(om) >= TWO_PI * f1_hz) & (np.abs(om) <= TWO_PI * f2_hz)
        if law == "square":
            y = ifft_agrawal(np.where(band, spec, 0.0), dt).real ** 2
        else:  # analytic signal: one side of the band, doubled (which side is irrelevant for |z|)
            z = ifft_agrawal(np.where(band & (om > 0), 2.0 * spec, 0.0), dt)
            y = np.abs(z) if law == "linear" else np.log(np.abs(z) ** 2)
        stats.append(y.reshape(nb, m).mean(axis=1))
        done += nb
    s = np.concatenate(stats)
    bt = (f2_hz - f1_hz) * m * dt
    if law == "square":
        rel = s.var() / s.mean() ** 2
    elif law == "linear":
        rel = 4.0 * s.var() / s.mean() ** 2
    else:
        rel = s.var()
    return float(rel * bt)
