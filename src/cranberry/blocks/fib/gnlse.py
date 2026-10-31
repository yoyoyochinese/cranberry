"""``fiber_gnlse``: classical GNLSE by RK4IP with per-term toggles and the CW shortcut (WP-B9, T06).

Equation (T06 eq. (1), Agrawal convention, conv §4/§8; FFTs only through ``core.fft``)::

    ∂Ã/∂z = D̂(Ω) Ã + N̂(Ã)
    D̂(Ω)  = −α/2 + i Σ_{k≥2} β_k Ω^k / k!                         (= −α/2 + i[β(ω0+Ω) − β0 − β1Ω])
    N̂(Ã)  = iγ (1 + Ω/ω0) F{ A · [c_K |A|² + c_R (h_R ⊛ |A|²)] }  + i b γ P_x(z) Ã

* ``h_R`` is given analytically in frequency, h̃_R(Ω) = ∫ h_R e^{+iΩt} dt (conv §8, T06 §2.3/§3.4):
  Blow–Wood (τ1 = 12.2 fs, τ2 = 32 fs) or Lin–Agrawal (f_a = 0.75, f_b = 0.21, f_c = 0.04,
  τ_b = 96 fs). The model constants belong to the named model; ``f_R`` is a spec parameter (D34).
* Toggle weights (D49 "total γ kept when Raman off"): kerr on/raman on → c_K = 1 − f_R, c_R = f_R;
  kerr on/raman off → c_K = 1, c_R = 0; kerr off/raman on → c_K = 0, c_R = f_R (Raman part only,
  for decompositions); both off → no nonlinearity. ``self_steepening`` off → factor 1.
* ``xpm_external``: counter-propagating field as an averaged intensity term (T06 §2.5-3, eq. (3) with
  the probe power averaged over the 2n_gL/c window): P_x(z) = P̄ e^{−α(L−z)} (P̄ enters at z = L).
  Co-propagating same-polarisation fields are already inside A (no extra term, T06 §2.5-1).
* RK4IP (Hult 2007, T06 §3.1) on the spectrum (FFT order): the linear step e^{hD̂/2} is exact,
  the global error is O(h⁴) from N̂ only. ``control = halving``: step-doubling local error
  (Sinkin 2003) with the fine solution kept (no Richardson extrapolation, so VT-09 stays 4th order).
* CW (``laser_mode = cw`` with ``cw_shortcut`` on): exact solution T06 eq. (5),
  A(L) = A(0) e^{−αL/2} e^{iγ_eff P L_eff} (+ XPM bγP̄ L_eff); quantum sidebands rotate by
  θ_disp(Ω) = L Σ_{even k} β_k Ω^k/k! (T06 §8.7). Same result as RK4IP on a grid (VT-03).

Tangent-linear layer for WP-B10 (T07 §2.2 (a)): :meth:`GnlseModel.step` and :func:`propagate`
also carry tangent fields δÃ through the *same* RK4IP stages (the variational equation
δÃ' = D̂ δÃ + dN̂(Ã)[δÃ], not complex-linear: δ|A|² = 2 Re(A* δA)). The tangent of the discrete
map is exact (tested against finite differences).

Grid: N samples, Δt = T_win/N, time grid centred (t0 = −(N//2)Δt). Spectra returned to the user are
``fft_agrawal(A, dt, t0=t0)`` (absolute spectral phase); the solver works with t0 = 0 internally (all
operators are diagonal in Ω, so the t0 phase cancels).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from math import factorial
from typing import Any

import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.fft import fft_agrawal, fftshift, ifft_agrawal, omega_grid
from cranberry.core.registry import register_block

__all__ = [
    "C_LIGHT",
    "RAMAN_MODELS",
    "TOGGLES",
    "EDGE_FRACTION",
    "EDGE_TOL",
    "raman_response",
    "raman_first_moment_s",
    "db_per_km_to_alpha_per_m",
    "l_eff",
    "GnlseModel",
    "Propagation",
    "propagate",
    "grid_edge_fractions",
    "grid_edge_monitor",
    "cw_shortcut",
    "cw_sideband_rotation",
    "CwResult",
    "GnlseResult",
    "FiberGnlse",
    "make_fiber_gnlse",
    "model_from_params",
    "pulse_field",
    "solve_pulsed",
    "time_axis",
]

KIND = "fiber_gnlse"
C_LIGHT = 299_792_458.0
RAMAN_MODELS = ("blow_wood", "lin_agrawal")
TOGGLES = ("loss", "beta2", "beta3", "beta_higher", "kerr", "raman", "self_steepening", "xpm_external",
           "cw_shortcut")
#: grid-edge monitor (T06 §3.3, V11): energy in the outer 5 % of points on each side ≤ 1e−6
EDGE_FRACTION = 0.05
EDGE_TOL = 1e-6

# Raman model constants (T06 §2.3; values cross-checked in gnlse-python 2.0.0 and PyNLO 0.1.2,
# UNVERIFIED against the primary papers)
_BW_TAU1, _BW_TAU2 = 12.2e-15, 32e-15
_LA_FA, _LA_FB, _LA_FC, _LA_TAUB = 0.75, 0.21, 0.04, 96e-15


def raman_response(model: str, omega: np.ndarray) -> np.ndarray:
    """h̃_R(Ω) = ∫ h_R(t) e^{+iΩt} dt (h̃_R(0) = 1), analytic (T06 §2.3, conv §8)."""
    w = np.asarray(omega, dtype=float)
    t1, t2 = _BW_TAU1, _BW_TAU2
    ha = ((t1**2 + t2**2) / (t1**2 * t2**2)) / ((1.0 / t2 - 1j * w) ** 2 + 1.0 / t1**2)
    if model == "blow_wood":
        return ha
    if model == "lin_agrawal":
        x = w * _LA_TAUB
        hb = (1.0 - 2j * x) / (1.0 - 1j * x) ** 2
        return (_LA_FA + _LA_FC) * ha + _LA_FB * hb
    raise ValueError(f"raman_model must be one of {RAMAN_MODELS}, got {model!r}")


def raman_first_moment_s(model: str) -> float:
    """∫ t h_R(t) dt [s] (T_R = f_R × this, Gordon SSFS). BW: 2τ1²τ2/(τ1²+τ2²); h_b has zero moment."""
    t1, t2 = _BW_TAU1, _BW_TAU2
    m_a = 2 * t1**2 * t2 / (t1**2 + t2**2)
    if model == "blow_wood":
        return m_a
    if model == "lin_agrawal":
        return (_LA_FA + _LA_FC) * m_a
    raise ValueError(f"raman_model must be one of {RAMAN_MODELS}, got {model!r}")


def db_per_km_to_alpha_per_m(x_db_per_km: float) -> float:
    """Power attenuation α [1/m] from α_dB [dB/km]: α = α_dB·ln(10)/10 / 1000 (T06 §1)."""
    return float(x_db_per_km) * np.log(10.0) / 10.0 / 1000.0


def l_eff(alpha_per_m: float, length_m: float) -> float:
    """L_eff = (1 − e^{−αL})/α (→ L for α → 0)."""
    a = float(alpha_per_m) * float(length_m)
    return float(length_m) if a == 0.0 else float(-np.expm1(-a) / alpha_per_m)


# --------------------------------------------------------------------------------------------------
# model + integrator
# --------------------------------------------------------------------------------------------------


class GnlseModel:
    """Operators of one fiber on one grid (toggles already folded in; see module docstring).

    Args:
        n_t, dt_s: grid size and step.
        lambda0_m: carrier wavelength.
        alpha_per_m: power attenuation (0 = loss off).
        betas: (β2, β3, β4, ...) [s^k/m] (zeros for terms switched off).
        gamma_per_w_m: total nonlinear coefficient γ (Kerr + Raman, T06 §2.3).
        c_kerr, c_raman: weights of |A|² and h_R ⊛ |A|² (module docstring).
        raman_model: name or None (only used if ``c_raman`` != 0).
        self_steepening: (1 + Ω/ω0) factor on/off.
        xpm_rate_per_m: b·γ·P̄ of the external counter-propagating field (0 = off).
        length_m: fiber length (z reference for the XPM loss profile).
    """

    def __init__(self, n_t: int, dt_s: float, *, lambda0_m: float, alpha_per_m: float, betas, gamma_per_w_m: float,
                 c_kerr: float, c_raman: float, raman_model: str | None, self_steepening: bool,
                 xpm_rate_per_m: float = 0.0, length_m: float = 0.0):
        self.n_t, self.dt = int(n_t), float(dt_s)
        self.omega = omega_grid(self.n_t, self.dt)
        self.omega0 = 2 * np.pi * C_LIGHT / float(lambda0_m)
        self.alpha = float(alpha_per_m)
        disp = sum(float(b) * self.omega ** (k + 2) / factorial(k + 2) for k, b in enumerate(betas))
        self.D = -self.alpha / 2 + 1j * np.asarray(disp + 0.0 * self.omega)
        self.gamma = float(gamma_per_w_m)
        self.c_kerr, self.c_raman = float(c_kerr), float(c_raman)
        self.hR = raman_response(raman_model, self.omega) if self.c_raman != 0.0 else None
        self.shock = (1.0 + self.omega / self.omega0) if self_steepening else np.ones_like(self.omega)
        self.xpm_rate = float(xpm_rate_per_m)
        self.length = float(length_m)
        self.has_nl = self.gamma != 0.0 and (self.c_kerr != 0.0 or self.c_raman != 0.0)

    # --- transforms (internal t0 = 0) ---
    def to_t(self, at: np.ndarray) -> np.ndarray:
        return ifft_agrawal(at, self.dt)

    def to_w(self, a: np.ndarray) -> np.ndarray:
        return fft_agrawal(a, self.dt)

    def _resp(self, intensity: np.ndarray) -> np.ndarray:
        """c_K I + c_R h_R ⊛ I (real)."""
        out = self.c_kerr * intensity
        if self.hR is not None:
            out = out + self.c_raman * self.to_t(self.hR * self.to_w(intensity)).real
        return out

    def _xpm(self, z: float) -> float:
        return self.xpm_rate * np.exp(-self.alpha * (self.length - z)) if self.xpm_rate else 0.0

    def nonlinear(self, at: np.ndarray, z: float) -> np.ndarray:
        """N̂(Ã) (spectrum, FFT order)."""
        out = 1j * self._xpm(z) * at
        if self.has_nl:
            a = self.to_t(at)
            out = out + 1j * self.gamma * self.shock * self.to_w(a * self._resp(np.abs(a) ** 2))
        return out

    def nonlinear_tangent(self, at: np.ndarray, dat: np.ndarray, z: float) -> np.ndarray:
        """dN̂(Ã)[δÃ]; ``dat`` may carry leading batch axes (several tangent directions)."""
        out = 1j * self._xpm(z) * dat
        if self.has_nl:
            a, da = self.to_t(at), self.to_t(dat)
            dI = 2.0 * (a.conj() * da).real
            out = out + 1j * self.gamma * self.shock * self.to_w(da * self._resp(np.abs(a) ** 2) + a * self._resp(dI))
        return out

    def step(self, at: np.ndarray, h: float, z: float, dat: np.ndarray | None = None):
        """One RK4IP step z → z + h (T06 §3.1). Returns (Ã, δÃ or None)."""
        e = np.exp(0.5 * h * self.D)
        tan = dat is not None

        def f(zz, a, da):
            return self.nonlinear(a, zz), (self.nonlinear_tangent(a, da, zz) if tan else None)

        ai = e * at
        di = e * dat if tan else None
        n1, m1 = f(z, at, dat)
        k1 = e * h * n1
        l1 = e * h * m1 if tan else None
        n2, m2 = f(z + h / 2, ai + k1 / 2, di + l1 / 2 if tan else None)
        k2 = h * n2
        n3, m3 = f(z + h / 2, ai + k2 / 2, di + h * m2 / 2 if tan else None)
        k3 = h * n3
        n4, m4 = f(z + h, e * (ai + k3), e * (di + h * m3) if tan else None)
        out = e * (ai + k1 / 6 + k2 / 3 + k3 / 3) + h * n4 / 6
        if not tan:
            return out, None
        dout = e * (di + l1 / 6 + h * m2 / 3 + h * m3 / 3) + h * m4 / 6
        return out, dout

    def linear(self, at: np.ndarray, z: float) -> np.ndarray:
        """Linear-only propagation e^{zD̂} Ã (the 'linear mode' reference)."""
        return np.exp(z * self.D) * at

    def energy(self, at: np.ndarray) -> float:
        """E = ∫|A|² dT = (1/2π)∫|Ã|² dΩ [J] (discrete Parseval)."""
        return float(np.sum(np.abs(at) ** 2) / (self.n_t * self.dt))

    def photon_number(self, at: np.ndarray) -> float:
        """∝ ∫|Ã|²/(ω0 + Ω) dΩ (arbitrary scale; T06 §2.4)."""
        return float(np.sum(np.abs(at) ** 2 / (self.omega0 + self.omega)))

    def centroid(self, at: np.ndarray) -> float:
        """Spectral centroid ⟨Ω⟩ [rad/s] of |Ã|²."""
        p = np.abs(at) ** 2
        return float(np.sum(self.omega * p) / np.sum(p))

    def rms_width(self, at: np.ndarray) -> float:
        """rms spectral width [rad/s]."""
        p = np.abs(at) ** 2
        c = np.sum(self.omega * p) / np.sum(p)
        return float(np.sqrt(np.sum((self.omega - c) ** 2 * p) / np.sum(p)))


@dataclass(frozen=True)
class Propagation:
    """Raw integrator output. ``at_z``/``dat_z``: spectra (FFT order, t0 = 0) at ``z`` (n_save + 1 points)."""

    z: np.ndarray
    at_z: np.ndarray
    dat_z: np.ndarray | None
    phi_nl_peak_rad: np.ndarray  # γ(c_K + c_R) ∫ max_T |A|² dz up to each saved z (trapezoid over steps)
    n_steps: int                 # accepted steps
    n_rejected: int


def propagate(model: GnlseModel, at0: np.ndarray, length_m: float, *, n_steps: int, n_save: int,
              control: str = "fixed", rtol: float | None = None, dat0: np.ndarray | None = None) -> Propagation:
    """Integrate z = 0 → L. ``fixed``: ``n_steps`` equal steps (multiple of ``n_save``).
    ``halving``: initial h = L/n_steps, step-doubling local error ≤ ``rtol`` (Sinkin 2003)."""
    if control not in ("fixed", "halving"):
        raise ValueError(f"control must be 'fixed' or 'halving', got {control!r}")
    n_steps, n_save = int(n_steps), int(n_save)
    if n_steps < 1 or n_save < 1:
        raise ValueError("n_steps and n_save must be >= 1")
    if control == "fixed" and n_steps % n_save:
        raise ValueError(f"n_steps ({n_steps}) must be a multiple of n_save ({n_save})")
    if control == "halving" and not (rtol and rtol > 0):
        raise ValueError("control 'halving' needs rtol > 0")
    z_save = np.linspace(0.0, float(length_m), n_save + 1)
    at = np.asarray(at0, dtype=complex)
    dat = None if dat0 is None else np.asarray(dat0, dtype=complex)
    ats, dats, phis = [at], [dat], [0.0]
    g_nl = model.gamma * (model.c_kerr + model.c_raman) if model.has_nl else 0.0  # CW-equivalent rate
    peak = lambda a: g_nl * float(np.max(np.abs(model.to_t(a)) ** 2))  # noqa: E731
    phi, p_prev = 0.0, peak(at)
    z, h, n_acc, n_rej = 0.0, float(length_m) / n_steps, 0, 0
    for z_next in z_save[1:]:
        if control == "fixed":
            for _ in range(n_steps // n_save):
                at, dat = model.step(at, h, z, dat)
                z += h
                p = peak(at)
                phi, p_prev, n_acc = phi + 0.5 * h * (p + p_prev), p, n_acc + 1
            z = z_next
        else:
            while z_next - z > 1e-12 * length_m:
                ht = min(h, z_next - z)
                coarse, _ = model.step(at, ht, z, None)
                mid, dmid = model.step(at, ht / 2, z, dat)
                fine, dfine = model.step(mid, ht / 2, z + ht / 2, dmid)
                err = np.linalg.norm(fine - coarse) / max(np.linalg.norm(fine), 1e-300)
                if err > 2 * rtol:
                    h, n_rej = ht / 2, n_rej + 1
                    continue
                at, dat, z = fine, dfine, z + ht
                p = peak(at)
                phi, p_prev, n_acc = phi + 0.5 * ht * (p + p_prev), p, n_acc + 1
                if err > rtol:
                    h = ht / 2**0.2
                elif err < rtol / 2:
                    h = max(h, ht * 2**0.2)
            z = z_next
        ats.append(at)
        dats.append(dat)
        phis.append(phi)
    return Propagation(z_save, np.array(ats), None if dat0 is None else np.array(dats), np.array(phis), n_acc, n_rej)


def grid_edge_fractions(model: GnlseModel, at: np.ndarray) -> tuple[float, float]:
    """(time, frequency) energy fraction in the outer ``EDGE_FRACTION`` of points on each side (T06 §3.3)."""
    m = max(1, int(round(EDGE_FRACTION * model.n_t)))
    # internal order has index 0 at T = 0 / Ω = 0, so fftshift puts the window edges at both ends
    it = np.abs(fftshift(model.to_t(at))) ** 2
    iw = np.abs(fftshift(at)) ** 2
    frac = lambda p: float((p[:m].sum() + p[-m:].sum()) / p.sum())  # noqa: E731
    return frac(it), frac(iw)


def grid_edge_monitor(model: GnlseModel, at: np.ndarray, source: str, log: ChkLog | None = None) -> tuple[float, float]:
    """CHK ``grid_edge`` (PLAN §4.3 CHK ``grid_edge_monitor``, G16-j): warn (never stop) if an edge fraction > EDGE_TOL."""
    ft, fw = grid_edge_fractions(model, at)
    if max(ft, fw) > EDGE_TOL:
        (log or ChkLog()).warn("grid_edge", source,
                               f"GNLSE energy near window edges: time {ft:.2e}, frequency {fw:.2e} (> {EDGE_TOL:g}); "
                               "widen T_win or refine dt", time_fraction=ft, freq_fraction=fw, tol=EDGE_TOL)
    return ft, fw


# --------------------------------------------------------------------------------------------------
# CW shortcut
# --------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class CwResult:
    """CW shortcut output (T06 eq. (5)). ``a_out = a_in·√η·e^{iφ_NL}``."""

    a_in: complex
    a_out: complex
    eta: float
    phi_nl_rad: float
    length_m: float
    betas: tuple[float, ...]

    def sideband_rotation_rad(self, omega_rad_s) -> np.ndarray:
        """Quadrature-angle rotation θ_disp(Ω) = L Σ_{even k} β_k Ω^k/k! (T06 §8.7)."""
        return cw_sideband_rotation(omega_rad_s, self.length_m, self.betas)


def cw_shortcut(power_w: float, *, length_m: float, alpha_per_m: float, gamma_eff_per_w_m: float,
                xpm_rate_per_m: float = 0.0, betas=()) -> CwResult:
    """Exact CW solution: A(L) = √P e^{−αL/2} e^{i(γ_eff P + bγP̄) L_eff}, A(0) = √P (real)."""
    le = l_eff(alpha_per_m, length_m)
    eta = float(np.exp(-alpha_per_m * length_m))
    phi = (gamma_eff_per_w_m * power_w + xpm_rate_per_m) * le
    a0 = complex(np.sqrt(power_w))
    return CwResult(a0, a0 * np.sqrt(eta) * np.exp(1j * phi), eta, float(phi), float(length_m), tuple(map(float, betas)))


def cw_sideband_rotation(omega_rad_s, length_m: float, betas) -> np.ndarray:
    """θ_disp(Ω) = L·[β(ω0+Ω) + β(ω0−Ω) − 2β0]/2 = L Σ_{even k} β_k Ω^k/k! (odd orders are a delay)."""
    w = np.asarray(omega_rad_s, dtype=float)
    return length_m * sum(float(b) * w ** (k + 2) / factorial(k + 2) for k, b in enumerate(betas) if k % 2 == 0) + 0.0 * w


# --------------------------------------------------------------------------------------------------
# block
# --------------------------------------------------------------------------------------------------


def pulse_field(shape: str, peak_power_w: float, t0_s: float, t: np.ndarray) -> np.ndarray:
    """A(0, T): ``sech`` √P0 sech(T/T0) or ``gaussian`` √P0 e^{−T²/(2T0²)} (T06 §1)."""
    x = np.asarray(t) / t0_s
    if shape == "sech":
        return np.sqrt(peak_power_w) / np.cosh(x)
    if shape == "gaussian":
        return np.sqrt(peak_power_w) * np.exp(-0.5 * x**2)
    raise ValueError(f"pulse shape must be 'sech' or 'gaussian', got {shape!r}")


_FIBER_KEYS = ("L_m", "alpha_per_m", "beta2_s2_per_m", "beta3_s3_per_m", "beta_higher", "gamma_per_w_m",
               "lambda0_m", "raman_model", "f_R", "xpm_b", "xpm_power_w")


def _check_fiber(params: Mapping[str, Any], name: str) -> None:
    require(params, *_FIBER_KEYS, block=name)
    if params["raman_model"] not in RAMAN_MODELS:
        raise ValueError(f"block {name!r}: raman_model must be one of {RAMAN_MODELS}, got {params['raman_model']!r}")
    if not 0.0 <= float(params["f_R"]) < 1.0:
        raise ValueError(f"block {name!r}: f_R must be in [0, 1), got {params['f_R']!r}")


def _betas(params: Mapping[str, Any], tg: Mapping[str, bool]) -> tuple[float, ...]:
    hi = [float(b) for b in params["beta_higher"]] if tg["beta_higher"] else [0.0] * len(params["beta_higher"])
    return (float(params["beta2_s2_per_m"]) if tg["beta2"] else 0.0,
            float(params["beta3_s3_per_m"]) if tg["beta3"] else 0.0, *hi)


def _weights(params: Mapping[str, Any], tg: Mapping[str, bool]) -> tuple[float, float]:
    """(c_K, c_R), D49: total γ kept when Raman is off."""
    f_r = float(params["f_R"])
    return ((1.0 - f_r if tg["raman"] else 1.0) if tg["kerr"] else 0.0, f_r if tg["raman"] else 0.0)


def _alpha(params, tg) -> float:
    return float(params["alpha_per_m"]) if tg["loss"] else 0.0


def _xpm_rate(params, tg) -> float:
    return float(params["xpm_b"]) * float(params["gamma_per_w_m"]) * float(params["xpm_power_w"]) if tg["xpm_external"] else 0.0


def model_from_params(params: Mapping[str, Any], toggles: Mapping[str, bool], n_t: int, dt_s: float) -> GnlseModel:
    """:class:`GnlseModel` from a fiber param mapping (``_FIBER_KEYS``) and the toggles (shared by the FIB tables)."""
    c_k, c_r = _weights(params, toggles)
    return GnlseModel(n_t, dt_s, lambda0_m=float(params["lambda0_m"]), alpha_per_m=_alpha(params, toggles),
                      betas=_betas(params, toggles), gamma_per_w_m=float(params["gamma_per_w_m"]), c_kerr=c_k,
                      c_raman=c_r, raman_model=params["raman_model"], self_steepening=toggles["self_steepening"],
                      xpm_rate_per_m=_xpm_rate(params, toggles), length_m=float(params["L_m"]))


def time_axis(n_t: int, t_window_s: float) -> tuple[np.ndarray, float]:
    """Centred time grid T_n = (n − N//2)·Δt, Δt = T_win/N."""
    dt = float(t_window_s) / int(n_t)
    return (np.arange(int(n_t)) - int(n_t) // 2) * dt, dt


@dataclass(frozen=True)
class GnlseResult:
    """RK4IP output. Time arrays use the centred grid ``t_s``; spectra are Agrawal (t0 = t_s[0]) in FFT order."""

    t_s: np.ndarray
    omega_rad_s: np.ndarray
    z_m: np.ndarray
    a_tz: np.ndarray             # A(z_k, T), shape (n_save + 1, N)
    phi_nl_peak_rad: np.ndarray  # per z_k
    energy_j: tuple[float, float]
    photon_number: tuple[float, float]  # arbitrary scale (in, out)
    edge_fractions: tuple[float, float]  # (time, frequency) at z = L
    n_steps: int
    n_rejected: int
    model: GnlseModel

    @property
    def a_in(self) -> np.ndarray:
        return self.a_tz[0]

    @property
    def a_out(self) -> np.ndarray:
        return self.a_tz[-1]

    def spectrum(self, k: int = -1) -> np.ndarray:
        """Ã(z_k, Ω) with absolute spectral phase (t0 = t_s[0])."""
        return fft_agrawal(self.a_tz[k], self.model.dt, t0=float(self.t_s[0]))


def _to_centered(model: GnlseModel, at: np.ndarray, n_shift: int) -> np.ndarray:
    return np.roll(model.to_t(at), n_shift, axis=-1)


def solve_pulsed(model: GnlseModel, a0_t: np.ndarray, t_s: np.ndarray, length_m: float, *, n_steps: int, n_save: int,
                 control: str, rtol: float | None, source: str = KIND, log: ChkLog | None = None,
                 monitor_edges: bool = True) -> GnlseResult:
    """RK4IP from a field on the centred grid ``t_s`` (t_s[N//2] = 0); grid-edge CHK at z = 0 and L
    (``monitor_edges = False`` for a CW field, which fills the window by definition)."""
    n_shift = model.n_t // 2  # internal sample order: index 0 ↔ T = 0
    at0 = model.to_w(np.roll(np.asarray(a0_t, dtype=complex), -n_shift))
    if monitor_edges:
        grid_edge_monitor(model, at0, source, log)
    pr = propagate(model, at0, length_m, n_steps=n_steps, n_save=n_save, control=control, rtol=rtol)
    edge = grid_edge_monitor(model, pr.at_z[-1], source, log) if monitor_edges else grid_edge_fractions(model, pr.at_z[-1])
    return GnlseResult(
        t_s=np.asarray(t_s), omega_rad_s=model.omega, z_m=pr.z, a_tz=_to_centered(model, pr.at_z, n_shift),
        phi_nl_peak_rad=pr.phi_nl_peak_rad,
        energy_j=(model.energy(pr.at_z[0]), model.energy(pr.at_z[-1])),
        photon_number=(model.photon_number(pr.at_z[0]), model.photon_number(pr.at_z[-1])),
        edge_fractions=edge, n_steps=pr.n_steps, n_rejected=pr.n_rejected, model=model)


class FiberGnlse(BaseBlock):
    """The ``fiber_gnlse`` offline block (module docstring). Product: :attr:`result`.

    Params (all required, D34; parameters of every term are required even when its toggle is off,
    so ablation can switch any toggle):

    * ``toggles``: ``loss``, ``beta2``, ``beta3``, ``beta_higher``, ``kerr``, ``raman``,
      ``self_steepening``, ``xpm_external``, ``cw_shortcut``
    * fiber: ``L_m``, ``alpha_per_m`` (power, 1/m), ``beta2_s2_per_m``, ``beta3_s3_per_m``,
      ``beta_higher`` (list β4, β5, ... [s^k/m], may be empty), ``gamma_per_w_m`` (total γ),
      ``lambda0_m``, ``raman_model`` (blow_wood | lin_agrawal), ``f_R``, ``xpm_b``, ``xpm_power_w``
    * ``laser_mode``: ``cw`` → ``power_w``; ``pulsed`` → ``pulse: {shape: sech|gaussian, peak_power_w, t0_s}``
    * ``grid: {n_t, t_window_s}`` and ``steps: {n_steps, n_save, control: fixed|halving[, rtol]}``
      (pulsed, or CW with ``cw_shortcut`` off)
    """

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        _check_fiber(params, name)
        require(params, "laser_mode", block=name)
        mode = params["laser_mode"]
        if mode not in ("cw", "pulsed"):
            raise ValueError(f"block {name!r}: laser_mode must be 'cw' or 'pulsed', got {mode!r}")
        super().__init__(KIND, name, params)
        if mode == "cw":
            require(params, "power_w", block=name)
        else:
            require(params, "pulse", block=name)
            require(params["pulse"], "shape", "peak_power_w", "t0_s", block=f"{name}.pulse")
        if mode == "pulsed" or not self._toggles["cw_shortcut"]:
            require(params, "grid", "steps", block=name)
            require(params["grid"], "n_t", "t_window_s", block=f"{name}.grid")
            require(params["steps"], "n_steps", "n_save", "control", block=f"{name}.steps")
        self.chk_log = ChkLog()

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    @cached_property
    def result(self) -> CwResult | GnlseResult:
        p, tg = self.params, self._toggles
        if p["laser_mode"] == "cw" and tg["cw_shortcut"]:
            c_k, c_r = _weights(p, tg)
            return cw_shortcut(float(p["power_w"]), length_m=float(p["L_m"]), alpha_per_m=_alpha(p, tg),
                               gamma_eff_per_w_m=float(p["gamma_per_w_m"]) * (c_k + c_r),
                               xpm_rate_per_m=_xpm_rate(p, tg), betas=_betas(p, tg))
        g, s = p["grid"], p["steps"]
        t, dt = time_axis(g["n_t"], g["t_window_s"])
        model = model_from_params(p, tg, int(g["n_t"]), dt)
        if p["laser_mode"] == "cw":
            a0 = np.full(t.shape, np.sqrt(float(p["power_w"])), dtype=complex)
        else:
            q = p["pulse"]
            a0 = pulse_field(q["shape"], float(q["peak_power_w"]), float(q["t0_s"]), t)
        return solve_pulsed(model, a0, t, float(p["L_m"]), n_steps=int(s["n_steps"]), n_save=int(s["n_save"]),
                            control=s["control"], rtol=s.get("rtol"), source=self.name, log=self.chk_log,
                            monitor_edges=p["laser_mode"] == "pulsed")


@register_block(KIND)
def make_fiber_gnlse(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``fiber_gnlse`` block (:class:`FiberGnlse`)."""
    return FiberGnlse(params, name)

