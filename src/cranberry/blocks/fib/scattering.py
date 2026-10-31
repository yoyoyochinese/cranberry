"""Scattering-noise tables and the SBS check (WP-B10; T07 §3–§6, T10 §3.4, T17 §6; PLAN §4.3 FIB / CHK).

Blocks (offline, never stepped): ``raman_noise_table``, ``gawbs_table``, ``brillouin_table``,
``stray_light_table`` (cached through ``TableCache``) and the CHK ``sbs_check``. The functions are
shared with the linearized quantum layer (``blocks.fib.quantum``).

Formulas (conv §2: added phase-insensitive noise N photons/mode ⇒ V + N; N = ∫ n(f)|ũ(f)|² df;
ũ(f) = ∫u e^{+i2πft} dt, ∫|ũ|² df = 1, |ũ(0)|² = τ_eff; T17 §6)

* Raman (T07 §2.2, §3.1; Drummond & Corney 2001 Eq. (77) in conv units). g_R(Ω) = 2γ f_R |Im h̃_R(Ω)|
  with h̃_R from B9 (``raman_response``). Ω = (bright) − (scattered) angular frequency, Stokes Ω > 0.
  Scattered photon-flux density per unit bright power and length [photons s⁻¹ Hz⁻¹ / (W m)]:
  n_R/(PL) = g_R(|Ω|)·(n_th(|Ω|) + 1) (Stokes), g_R(|Ω|)·n_th(|Ω|) (anti-Stokes), n_th = 1/(e^{ħΩ/k_BT} − 1).
  Ω → 0: g_R n_th → g_R′ k_BT/ħ with g_R′ = 2γ f_R ∫t h_R dt (finite, T07 §3.1).
  Phase-noise field (two-sided, per m) S_R(Ω) = hν0 g_R(|Ω|)(n_th + 1/2) [rad² s/m].
  D48 §7-6: only the thermal part (∝ n_th) vanishes as T → 0; spontaneous Stokes g_R(Ω)PL remains.
* GAWBS (T07 §4.1). R0,m peaks: (1 − α_v²)J0(y) − α_v² J2(y) = 0, α_v = V_S/V_L, f_m = y_m V_L/(π d)
  (m = 1 ↔ 30.5 MHz at 125 µm, London et al. numbering: 321.75 MHz measured = R0,7). Two-sided
  phase-noise field S_G(Ω) = hν0 Σ_m 2γ_OM,m L_m(|Ω|)(n_th + 1/2), L_m = Lorentzian of peak 1 and FWHM
  Γ_m = 2π Δf_m (the factor 2 in "2γ_OM" is UNVERIFIED, T07). Accumulated over L: one-sided
  S_φ(f) = 2L·S_G(2πf); per-mode variance σ_m² = 2π L hν0 γ_OM,m (n_th(Ω_m) + 1/2) Δf_m (Lorentz area,
  n_th at the peak); correlation C(τ) = Σ_m hν0 2γ_OM,m (n_th,m + 1/2)(Γ_m/2) cos(Ω_m τ) e^{−Γ_m|τ|/2}
  [rad²/m]. ``mode_sum_factor`` scales every term for the modes not modelled (TR2,m, PM asymmetric
  peaks; T07: total 1–5 × R0,7 [AE]). PM vs SMF enters through the cladding diameter (80 / 125 µm)
  and the per-mode (γ_OM, Δf) lists. Scattered photons from a co-propagating bright field P:
  N_G = Σ_m (P σ_m²/(2hν0))(|ũ(f_m)|² + |ũ(−f_m)|²) (T17 §6).
* Brillouin (T07 §5, T10 §3.4.3). Smith threshold P_th = 21 A_eff/(g_B L_eff)·(1 + Δν_L/Δν_B);
  gain exponent G_B = g_B P L_eff/A_eff (L_eff = (1 − e^{−αL})/α). Backscatter coefficient (both
  sidebands) b_B = (g_B/A_eff) hν0 (2n_th(ν_B) + 1)(π Δν_B/2) [1/m] (4.3e−9 /m for the T07 values);
  Stokes share (n_th + 1)/(2n_th + 1) at ν_pr − ν_B, anti-Stokes n_th/(2n_th + 1) at ν_pr + ν_B.
  Toggle ``sbs_amplification`` multiplies by (e^{G_B} − 1)/G_B (T07 §5.2 [AE]).
* Rayleigh (T10 §3.4.3, 00_summary C13): b_R from the datasheet (6.2e−8 /m adopted), same frequency
  as the probe, line width = probe line width.
* Counter-propagating probe → squeezed direction (T10 eq. (9)): flux Φ = ε_pol b P̄ L_bs/hν0 with
  P̄ = P_peak·duty, L_bs = ∫_0^L e^{−2αz} dz. Distributed scattering spreads over 2n_gL/c, so only
  P̄ matters (probe pulse width and timing do not, T10 §3.4.4; gating per port is WP-B8). ε_pol:
  ``same`` 1, ``pm_orthogonal`` ε_xt + T_⊥, ``smf_orthogonal`` 1/3 + (2/3)T_⊥ (T10 §3.4.3).
  Backward Raman uses the forward n_R (capture assumed equal, T07 §3.1 [AE, UNVERIFIED]).
* Stray light (T10 eq. (9), T17 §6): clicks per gate = T_gate·[Σ_lines Φ_k F(f_k) I_k + ∫ n(f) F(f) df]
  (F: GBS-side filter, I_WDM on Rayleigh and reflections), photons per mode = Σ Φ_k I_k |ũ(f_k)|² +
  ∫ n(f)|ũ(f)|² df (lines evaluated at their centre). Offsets f are relative to the squeezed carrier.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv

from cranberry.blocks.fib.gnlse import C_LIGHT, RAMAN_MODELS, l_eff, raman_first_moment_s, raman_response
from cranberry.blocks.fib.tables import _ascending, _cached
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.offline.tables import Table, TableCache

__all__ = [
    "H_PLANCK", "HBAR", "K_B", "POLARIZATIONS", "MODE_SHAPES", "STRAY_COLUMNS", "RAMAN_COLUMNS",
    "n_thermal", "raman_gain", "raman_scatter_density", "raman_phase_psd",
    "gawbs_peak_frequencies", "gawbs_field_psd", "gawbs_phase_psd", "gawbs_variances", "gawbs_correlation", "gawbs_mode_photons",
    "sbs_threshold_w", "brillouin_gain_exponent", "brillouin_backscatter_coeff", "backscatter_length",
    "eps_pol", "mode_weight", "mode_tau_eff", "mode_average", "backscatter_sources", "stray_photons",
    "make_raman_noise_table", "make_gawbs_table", "make_brillouin_table", "make_stray_light_table", "make_sbs_check",
]

H_PLANCK = 6.62607015e-34
HBAR = H_PLANCK / (2 * np.pi)
K_B = 1.380649e-23
SMITH_FACTOR = 21.0  # G_B at the Smith threshold (T07 §5.1, UNVERIFIED)
POLARIZATIONS = ("same", "pm_orthogonal", "smf_orthogonal")
MODE_SHAPES = ("gaussian", "exponential")
RAMAN_COLUMNS = ("g_r_per_w_m", "n_scatter_per_w_m", "n_thermal_per_w_m", "s_phase_rad2_s_per_m")
STRAY_COLUMNS = ("clicks_per_gate", "photons_per_mode", "rayleigh_gate", "brillouin_gate", "raman_gate",
                 "reflections_gate", "rayleigh_mode", "brillouin_mode", "raman_mode", "reflections_mode")
_REV = 1


def _hnu(lambda0_m: float) -> float:
    return H_PLANCK * C_LIGHT / float(lambda0_m)


def n_thermal(omega_rad_s, temperature_k: float) -> np.ndarray:
    """Bose occupation n_th(|Ω|, T) (T = 0 → 0; Ω = 0 → inf)."""
    w = np.abs(np.asarray(omega_rad_s, dtype=float))
    if temperature_k <= 0.0:
        return np.zeros_like(w)
    with np.errstate(divide="ignore", over="ignore"):
        return 1.0 / np.expm1(HBAR * w / (K_B * float(temperature_k)))


# --- Raman --------------------------------------------------------------------------------------------

def raman_gain(omega_rad_s, *, gamma_per_w_m: float, f_R: float, raman_model: str) -> np.ndarray:
    """g_R(|Ω|) = 2γ f_R |Im h̃_R(Ω)| [1/(W m)]."""
    return 2.0 * gamma_per_w_m * f_R * np.abs(raman_response(raman_model, np.abs(np.asarray(omega_rad_s, float))).imag)


def _g_nth(omega, temperature_k, gamma_per_w_m, f_R, raman_model):
    """g_R(|Ω|)·n_th(|Ω|) with the finite Ω → 0 limit g_R′ k_BT/ħ."""
    w = np.abs(np.asarray(omega, dtype=float))
    g = raman_gain(w, gamma_per_w_m=gamma_per_w_m, f_R=f_R, raman_model=raman_model)
    if temperature_k <= 0.0:
        return np.zeros_like(w)
    zero = w == 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        out = g * n_thermal(w, temperature_k)
    g1 = 2.0 * gamma_per_w_m * f_R * raman_first_moment_s(raman_model)
    return np.where(zero, g1 * K_B * temperature_k / HBAR, out)


def raman_scatter_density(omega_rad_s, *, gamma_per_w_m, f_R, raman_model, temperature_k, thermal_only=False):
    """n_R/(PL) [photons s⁻¹ Hz⁻¹ per W per m]; Ω = bright − scattered (Stokes Ω > 0). ``thermal_only``:
    the g_R n_th part only (D48 §7-6)."""
    w = np.asarray(omega_rad_s, dtype=float)
    th = _g_nth(w, temperature_k, gamma_per_w_m, f_R, raman_model)
    if thermal_only:
        return th
    return th + np.where(w > 0, raman_gain(w, gamma_per_w_m=gamma_per_w_m, f_R=f_R, raman_model=raman_model), 0.0)


def raman_phase_psd(omega_rad_s, *, gamma_per_w_m, f_R, raman_model, temperature_k, lambda0_m):
    """Two-sided Raman phase-noise field spectrum S_R(Ω) = hν0 g_R(|Ω|)(n_th + 1/2) [rad² s/m]."""
    w = np.asarray(omega_rad_s, dtype=float)
    g = raman_gain(w, gamma_per_w_m=gamma_per_w_m, f_R=f_R, raman_model=raman_model)
    return _hnu(lambda0_m) * (_g_nth(w, temperature_k, gamma_per_w_m, f_R, raman_model) + 0.5 * g)


# --- GAWBS --------------------------------------------------------------------------------------------

def gawbs_peak_frequencies(cladding_diameter_m: float, v_l_m_s: float, v_s_m_s: float, n_modes: int) -> np.ndarray:
    """R0,m peak frequencies f_m [Hz], m = 1..n_modes (module docstring)."""
    a2 = (float(v_s_m_s) / float(v_l_m_s)) ** 2
    f = lambda y: (1 - a2) * jv(0, y) - a2 * jv(2, y)  # noqa: E731
    roots, y, dy = [], 1e-3, 0.05  # roots are ~π apart: a 0.05 scan never skips one
    while len(roots) < int(n_modes):
        if f(y) * f(y + dy) < 0:
            roots.append(brentq(f, y, y + dy, xtol=1e-14))
        y += dy
    return np.array(roots) * float(v_l_m_s) / (np.pi * float(cladding_diameter_m))


def _gawbs_modes(g: Mapping[str, Any]):
    """(f_m, γ_OM,m·mode_sum_factor, Δf_m) of a ``gawbs`` param mapping."""
    n = int(g["n_modes"])
    fm = gawbs_peak_frequencies(g["cladding_diameter_m"], g["v_l_m_s"], g["v_s_m_s"], n)
    gam = np.broadcast_to(np.asarray(g["gamma_om_per_w_m"], float), (n,)) * float(g["mode_sum_factor"])
    df = np.broadcast_to(np.asarray(g["linewidth_hz"], float), (n,))
    return fm, gam, df


def _occ(omega, temperature_k, thermal_only):
    return n_thermal(omega, temperature_k) + (0.0 if thermal_only else 0.5)


def gawbs_field_psd(f_hz, g: Mapping[str, Any], *, temperature_k, lambda0_m, thermal_only=False):
    """Two-sided GAWBS phase-noise field spectrum S_G(2π|f|) per unit length [rad² s/m]."""
    fm, gam, df = _gawbs_modes(g)
    f = np.abs(np.asarray(f_hz, dtype=float))[..., None]
    lor = (df / 2) ** 2 / ((f - fm) ** 2 + (df / 2) ** 2)
    return _hnu(lambda0_m) * np.sum(2 * gam * lor, axis=-1) * _occ(2 * np.pi * f[..., 0], temperature_k, thermal_only)


def gawbs_phase_psd(f_hz, g: Mapping[str, Any], *, temperature_k, lambda0_m, length_m, thermal_only=False):
    """One-sided accumulated GAWBS phase PSD S_φ(f) = 2L·S_G(2πf) [rad²/Hz] (f > 0)."""
    return 2.0 * float(length_m) * gawbs_field_psd(f_hz, g, temperature_k=temperature_k, lambda0_m=lambda0_m,
                                                   thermal_only=thermal_only)


def gawbs_variances(g: Mapping[str, Any], *, temperature_k, lambda0_m, length_m, thermal_only=False) -> np.ndarray:
    """Per-mode accumulated phase variance σ_m² [rad²] (Lorentz area, n_th at the peak)."""
    fm, gam, df = _gawbs_modes(g)
    occ = _occ(2 * np.pi * fm, temperature_k, thermal_only)
    return 2 * np.pi * float(length_m) * _hnu(lambda0_m) * gam * occ * df


def gawbs_correlation(tau_s, g: Mapping[str, Any], *, temperature_k, lambda0_m) -> np.ndarray:
    """C(τ) of the GAWBS phase-noise field per unit length [rad²/m] (module docstring)."""
    fm, gam, df = _gawbs_modes(g)
    w, gw = 2 * np.pi * fm, 2 * np.pi * df
    t = np.abs(np.asarray(tau_s, dtype=float))[..., None]
    amp = _hnu(lambda0_m) * 2 * gam * _occ(w, temperature_k, False) * gw / 2
    return np.sum(amp * np.cos(w * t) * np.exp(-gw * t / 2), axis=-1)


def gawbs_mode_photons(power_w, g, mode, *, temperature_k, lambda0_m, length_m) -> float:
    """Photons per mode scattered from a co-propagating bright field: Σ_m P σ_m² (|ũ(f_m)|² + |ũ(−f_m)|²)/(2hν0)."""
    fm = _gawbs_modes(g)[0]
    var = gawbs_variances(g, temperature_k=temperature_k, lambda0_m=lambda0_m, length_m=length_m)
    return float(np.sum(power_w * var * 2 * mode_weight(fm, mode) / (2 * _hnu(lambda0_m))))


# --- Brillouin, Rayleigh, polarization ---------------------------------------------------------------

def sbs_threshold_w(*, g_b_m_per_w, a_eff_m2, alpha_per_m, length_m, laser_linewidth_hz, dnu_b_hz) -> float:
    """Smith threshold P_th = 21 A_eff/(g_B L_eff)·(1 + Δν_L/Δν_B) [W]."""
    return SMITH_FACTOR * a_eff_m2 / (g_b_m_per_w * l_eff(alpha_per_m, length_m)) * (1 + laser_linewidth_hz / dnu_b_hz)


def brillouin_gain_exponent(power_w, *, g_b_m_per_w, a_eff_m2, alpha_per_m, length_m, laser_linewidth_hz, dnu_b_hz):
    """G_B = g_B P L_eff/A_eff (g_B reduced by Δν_B/(Δν_B + Δν_L))."""
    g = g_b_m_per_w * dnu_b_hz / (dnu_b_hz + laser_linewidth_hz)
    return g * np.asarray(power_w, dtype=float) * l_eff(alpha_per_m, length_m) / a_eff_m2


def brillouin_backscatter_coeff(*, g_b_m_per_w, a_eff_m2, nu_b_hz, dnu_b_hz, temperature_k, lambda0_m):
    """(b_B [1/m] both sidebands, Stokes share, anti-Stokes share)."""
    n = float(n_thermal(2 * np.pi * nu_b_hz, temperature_k))
    b = g_b_m_per_w / a_eff_m2 * _hnu(lambda0_m) * (2 * n + 1) * np.pi * dnu_b_hz / 2
    return b, (n + 1) / (2 * n + 1), n / (2 * n + 1)


def backscatter_length(alpha_per_m: float, length_m: float) -> float:
    """L_bs = ∫_0^L e^{−2αz} dz (probe in, scattered light out, T10 eq. (8))."""
    return l_eff(2.0 * alpha_per_m, length_m)


def eps_pol(polarization: str, eps_xt: float, t_perp: float) -> float:
    """Fraction of the backscatter in the squeezed polarization (T10 §3.4.3)."""
    if polarization == "same":
        return 1.0
    if polarization == "pm_orthogonal":
        return float(eps_xt) + float(t_perp)
    if polarization == "smf_orthogonal":
        return 1.0 / 3.0 + 2.0 / 3.0 * float(t_perp)
    raise ValueError(f"polarization must be one of {POLARIZATIONS}, got {polarization!r}")


# --- temporal mode -----------------------------------------------------------------------------------

def _check_mode(mode: Mapping[str, Any], name: str) -> None:
    require(mode, "shape", block=f"{name}.mode")
    if mode["shape"] == "gaussian":
        require(mode, "t_fwhm_s", block=f"{name}.mode")
    elif mode["shape"] == "exponential":
        require(mode, "tau_eff_s", block=f"{name}.mode")
    else:
        raise ValueError(f"{name}: mode.shape must be one of {MODE_SHAPES}, got {mode['shape']!r}")


def mode_tau_eff(mode: Mapping[str, Any]) -> float:
    """τ_eff = |ũ(0)|² = (∫u dt)²: Gaussian 2.129·τ_FWHM (intensity), exponential = given."""
    if mode["shape"] == "gaussian":
        return float(mode["t_fwhm_s"]) * np.sqrt(np.pi / np.log(2.0))
    return float(mode["tau_eff_s"])


def mode_weight(f_hz, mode: Mapping[str, Any]) -> np.ndarray:
    """|ũ(f)|² [s]: ``gaussian`` (intensity FWHM τ) τ_eff e^{−4 ln2 (f τ/0.441)²}-type, ``exponential``
    (u = √(2κ) e^{−κt}, τ_eff = 2/κ) τ_eff/(1 + (π τ_eff f)²). Both integrate to 1."""
    f = np.asarray(f_hz, dtype=float)
    te = mode_tau_eff(mode)
    if mode["shape"] == "gaussian":
        return te * np.exp(-np.pi * (te * f) ** 2)
    return te / (1.0 + (np.pi * te * f) ** 2)


def mode_average(n_of_f: Callable[[np.ndarray], np.ndarray], mode: Mapping[str, Any], n_pts: int = 4001) -> float:
    """∫ n(f)|ũ(f)|² df for a smooth density n (Gaussian: ±6σ trapezoid; exponential: f = tan θ/(πτ))."""
    te = mode_tau_eff(mode)
    if mode["shape"] == "gaussian":
        f = np.linspace(-6, 6, n_pts) / (np.sqrt(2 * np.pi) * te)
        return float(np.trapezoid(n_of_f(f) * mode_weight(f, mode), f))
    th = (np.arange(n_pts) + 0.5) * np.pi / n_pts - np.pi / 2  # midpoint rule in θ
    return float(np.mean(n_of_f(np.tan(th) / (np.pi * te))))


# --- counter-propagating probe → squeezed direction ---------------------------------------------------

def backscatter_sources(p: Mapping[str, Any], toggles: Mapping[str, bool], p_avg_w: float):
    """Lines [(name, f0 [Hz], flux [1/s], FWHM [Hz])] and the broadband density n(f) [1/(s Hz)] that a
    counter-propagating probe of average power ``p_avg_w`` sends into the squeezed direction. Offsets f
    are scattered − squeezed-carrier frequency. ``p``: the ``probe``/``brillouin``/``rayleigh`` groups and
    the fiber keys (module docstring); ``toggles``: ``rayleigh``, ``brillouin``, ``raman``."""
    pr, br = p["probe"], p["brillouin"]
    hnu = _hnu(p["lambda0_m"])
    eps = eps_pol(pr["polarization"], pr["eps_xt"], pr["t_perp"])
    flux = eps * p_avg_w * backscatter_length(p["alpha_per_m"], p["L_m"]) / hnu  # × b → photons/s
    d0 = float(pr["freq_offset_hz"])
    lines = []
    if toggles["rayleigh"]:
        lines.append(("rayleigh", d0, flux * float(p["rayleigh"]["b_r_per_m"]), float(pr["linewidth_hz"])))
    if toggles["brillouin"]:
        b, s, a = brillouin_backscatter_coeff(g_b_m_per_w=br["g_b_m_per_w"], a_eff_m2=br["a_eff_m2"],
                                              nu_b_hz=br["nu_b_hz"], dnu_b_hz=br["dnu_b_hz"],
                                              temperature_k=p["temperature_k"], lambda0_m=p["lambda0_m"])
        if toggles.get("sbs_amplification", False):
            gb = float(brillouin_gain_exponent(p_avg_w, g_b_m_per_w=br["g_b_m_per_w"], a_eff_m2=br["a_eff_m2"],
                                               alpha_per_m=p["alpha_per_m"], length_m=p["L_m"],
                                               laser_linewidth_hz=pr["linewidth_hz"], dnu_b_hz=br["dnu_b_hz"]))
            b *= np.expm1(gb) / gb if gb > 0 else 1.0
        lines += [("brillouin", d0 - br["nu_b_hz"], flux * b * s, float(br["dnu_b_hz"])),
                  ("brillouin", d0 + br["nu_b_hz"], flux * b * a, float(br["dnu_b_hz"]))]

    def n_raman(f):
        if not toggles["raman"]:
            return np.zeros_like(np.asarray(f, dtype=float))
        w = 2 * np.pi * (d0 - np.asarray(f, dtype=float))  # Ω = probe − scattered
        return flux * hnu * raman_scatter_density(w, gamma_per_w_m=p["gamma_per_w_m"], f_R=p["f_R"],
                                                  raman_model=p["raman_model"], temperature_k=p["temperature_k"])

    return lines, n_raman


def _filter(f, filt: Mapping[str, Any]) -> np.ndarray:
    f = np.abs(np.asarray(f, dtype=float))
    if filt["shape"] == "flat":
        return (f <= filt["fwhm_hz"] / 2).astype(float)
    if filt["shape"] == "gaussian":
        return np.exp(-4 * np.log(2) * (f / filt["fwhm_hz"]) ** 2)
    raise ValueError(f"gbs_filter.shape must be 'flat' or 'gaussian', got {filt['shape']!r}")


def stray_photons(p: Mapping[str, Any], toggles: Mapping[str, bool], p_avg_w: float) -> dict[str, float]:
    """{<mechanism>_gate, <mechanism>_mode, clicks_per_gate, photons_per_mode} for one probe power
    (module docstring; ``toggles``: stray_rayleigh, stray_brillouin, stray_raman, stray_reflections)."""
    tg = {"rayleigh": toggles["stray_rayleigh"], "brillouin": toggles["stray_brillouin"],
          "raman": toggles["stray_raman"], "sbs_amplification": False}
    lines, n_raman = backscatter_sources(p, tg, p_avg_w)
    d0, iso, gate, filt, mode = (float(p["probe"]["freq_offset_hz"]), float(p["wdm_isolation_lin"]),
                                 float(p["t_gate_s"]), p["gbs_filter"], p["mode"])
    if toggles["stray_reflections"]:
        refl = sum(float(r) * float(e) for r, e in p["reflections"])
        lines.append(("reflections", d0, p_avg_w * refl / _hnu(p["lambda0_m"]), float(p["probe"]["linewidth_hz"])))
    out = dict.fromkeys(STRAY_COLUMNS[2:], 0.0)
    for name, f0, phi, _ in lines:
        i = iso if name in ("rayleigh", "reflections") else 1.0
        out[f"{name}_gate"] += gate * phi * i * float(_filter(f0, filt))
        out[f"{name}_mode"] += phi * i * float(mode_weight(f0, mode))
    if toggles["stray_raman"]:
        fb = np.linspace(-3, 3, 6001) * filt["fwhm_hz"]
        out["raman_gate"] = gate * float(np.trapezoid(n_raman(fb) * _filter(fb, filt), fb))
        out["raman_mode"] = mode_average(n_raman, mode)
    out["clicks_per_gate"] = sum(v for k, v in out.items() if k.endswith("_gate"))
    out["photons_per_mode"] = sum(v for k, v in out.items() if k.endswith("_mode"))
    return out


# --- blocks -------------------------------------------------------------------------------------------

_RAMAN_KEYS = ("gamma_per_w_m", "f_R", "raman_model", "lambda0_m", "temperature_k")
_BRILLOUIN_KEYS = ("g_b_m_per_w", "a_eff_m2", "nu_b_hz", "dnu_b_hz")
_PROBE_KEYS = ("freq_offset_hz", "linewidth_hz", "polarization", "eps_xt", "t_perp")
_GAWBS_KEYS = ("cladding_diameter_m", "v_l_m_s", "v_s_m_s", "n_modes", "gamma_om_per_w_m", "linewidth_hz",
               "mode_sum_factor")


def _check_raman(params, name):
    require(params, *_RAMAN_KEYS, block=name)
    if params["raman_model"] not in RAMAN_MODELS:
        raise ValueError(f"block {name!r}: raman_model must be one of {RAMAN_MODELS}, got {params['raman_model']!r}")


def check_probe_groups(params: Mapping[str, Any], name: str) -> None:
    """Validate the ``probe`` / ``brillouin`` / ``rayleigh`` groups (shared with ``fiber_quantum_linearized``)."""
    require(params, "probe", "brillouin", "rayleigh", block=name)
    require(params["probe"], *_PROBE_KEYS, block=f"{name}.probe")
    require(params["brillouin"], *_BRILLOUIN_KEYS, block=f"{name}.brillouin")
    require(params["rayleigh"], "b_r_per_m", block=f"{name}.rayleigh")
    eps_pol(params["probe"]["polarization"], params["probe"]["eps_xt"], params["probe"]["t_perp"])


def check_gawbs_group(params: Mapping[str, Any], name: str) -> None:
    require(params, "gawbs", block=name)
    require(params["gawbs"], *_GAWBS_KEYS, block=f"{name}.gawbs")


class _TableBlock(BaseBlock):
    kind_name = ""

    def __call__(self, state, params):
        return state  # offline block, not a time-loop stage

    def tables(self, cache: TableCache | None = None) -> dict[str, Table]:
        return _cached(self.kind_name, self.name, self.params, cache, self.table_names, self._build)


class RamanNoiseTable(_TableBlock):
    """``raman_noise_table``: Table ``raman_noise`` (axis ``f_hz`` = (bright − scattered)/2π, ascending,
    two-sided allowed) with columns :data:`RAMAN_COLUMNS`. Params (required, D34): ``gamma_per_w_m``,
    ``f_R``, ``raman_model``, ``lambda0_m``, ``temperature_k``, ``f_hz``."""

    kind_name, table_names = "raman_noise_table", ("raman_noise",)

    def __init__(self, params, name):
        name = name or self.kind_name
        _check_raman(params, name)
        require(params, "f_hz", block=name)
        _ascending(name, "f_hz", params["f_hz"])
        super().__init__(self.kind_name, name, params)

    def _build(self):
        p = self.params
        w = 2 * np.pi * np.asarray(p["f_hz"], dtype=float)
        kw = {k: p[k] for k in ("gamma_per_w_m", "f_R", "raman_model")}
        cols = [raman_gain(w, **kw),
                raman_scatter_density(w, temperature_k=p["temperature_k"], **kw),
                raman_scatter_density(w, temperature_k=p["temperature_k"], thermal_only=True, **kw),
                raman_phase_psd(w, temperature_k=p["temperature_k"], lambda0_m=p["lambda0_m"], **kw)]
        return {"raman_noise": Table("raman_noise", {"f_hz": w / (2 * np.pi)}, np.stack(cols, -1), "mixed",
                                     {"columns": list(RAMAN_COLUMNS), "sources": ["T07 §2.2", "T07 §3.1", "D48"]})}


class GawbsTable(_TableBlock):
    """``gawbs_table``: Tables ``gawbs_peaks`` (axis ``m`` = 1..n_modes; columns f_m [Hz], σ_m² [rad²]
    accumulated over ``L_m``) and ``gawbs_phase_psd`` (axis ``f_hz``; one-sided S_φ [rad²/Hz]).
    Params (required, D34): ``gawbs`` {cladding_diameter_m, v_l_m_s, v_s_m_s, n_modes,
    gamma_om_per_w_m (scalar or per mode), linewidth_hz (scalar or per mode), mode_sum_factor},
    ``lambda0_m``, ``temperature_k``, ``L_m``, ``f_hz`` (ascending, > 0)."""

    kind_name, table_names = "gawbs_table", ("gawbs_peaks", "gawbs_phase_psd")

    def __init__(self, params, name):
        name = name or self.kind_name
        check_gawbs_group(params, name)
        require(params, "lambda0_m", "temperature_k", "L_m", "f_hz", block=name)
        _ascending(name, "f_hz", params["f_hz"])
        super().__init__(self.kind_name, name, params)

    def _build(self):
        p, g = self.params, self.params["gawbs"]
        kw = dict(temperature_k=p["temperature_k"], lambda0_m=p["lambda0_m"], length_m=p["L_m"])
        fm = _gawbs_modes(g)[0]
        m = np.arange(1, fm.size + 1, dtype=float)
        f = np.asarray(p["f_hz"], dtype=float)
        src = {"sources": ["T07 §4.1", "London et al. 2018 (R0,7)"]}
        return {"gawbs_peaks": Table("gawbs_peaks", {"m": m}, np.stack([fm, gawbs_variances(g, **kw)], -1), "mixed",
                                     {**src, "columns": ["f_hz", "sigma2_rad2"]}),
                "gawbs_phase_psd": Table("gawbs_phase_psd", {"f_hz": f}, gawbs_phase_psd(f, g, **kw), "rad^2/Hz", src)}


class BrillouinTable(_TableBlock):
    """``brillouin_table``: Table ``brillouin`` (axis ``L_m``; columns L_eff [m], P_th [W], G_B at
    ``power_w``, b_B [1/m], Φ_B per W of average probe power [1/(s W)] with ε_pol = 1 and loss).
    Params (required, D34): ``brillouin`` group, ``lengths_m``, ``alpha_per_m``, ``power_w``,
    ``laser_linewidth_hz``, ``temperature_k``, ``lambda0_m``."""

    kind_name, table_names = "brillouin_table", ("brillouin",)

    def __init__(self, params, name):
        name = name or self.kind_name
        require(params, "brillouin", "lengths_m", "alpha_per_m", "power_w", "laser_linewidth_hz", "temperature_k",
                "lambda0_m", block=name)
        require(params["brillouin"], *_BRILLOUIN_KEYS, block=f"{name}.brillouin")
        _ascending(name, "lengths_m", params["lengths_m"])
        super().__init__(self.kind_name, name, params)

    def _build(self):
        p, br = self.params, self.params["brillouin"]
        b = brillouin_backscatter_coeff(**br, temperature_k=p["temperature_k"], lambda0_m=p["lambda0_m"])[0]
        rows = []
        for L in np.asarray(p["lengths_m"], dtype=float):
            kw = dict(g_b_m_per_w=br["g_b_m_per_w"], a_eff_m2=br["a_eff_m2"], alpha_per_m=p["alpha_per_m"],
                      length_m=L, laser_linewidth_hz=p["laser_linewidth_hz"], dnu_b_hz=br["dnu_b_hz"])
            rows.append([l_eff(p["alpha_per_m"], L), sbs_threshold_w(**kw),
                         float(brillouin_gain_exponent(p["power_w"], **kw)), b,
                         b * backscatter_length(p["alpha_per_m"], L) / _hnu(p["lambda0_m"])])
        return {"brillouin": Table("brillouin", {"L_m": p["lengths_m"]}, np.array(rows), "mixed",
                                   {"columns": ["l_eff_m", "p_th_w", "g_b", "b_b_per_m", "flux_per_w"],
                                    "sources": ["T07 §5", "T10 §3.4.3"]})}


class StrayLightTable(_TableBlock):
    """``stray_light_table``: Table ``stray_light`` (axis ``p_peak_w``; columns :data:`STRAY_COLUMNS`),
    stray photons of a counter-propagating probe per detection gate and per squeezed mode (module
    docstring). ``clicks_per_gate`` / ``photons_per_mode`` feed ``measure.stray_light`` (D46).

    Params (required, D34): ``toggles`` {stray_rayleigh, stray_brillouin, stray_raman, stray_reflections};
    ``p_peak_w`` (ascending), ``duty``, ``probe`` {freq_offset_hz, linewidth_hz, polarization, eps_xt,
    t_perp}, ``brillouin`` {g_b_m_per_w, a_eff_m2, nu_b_hz, dnu_b_hz}, ``rayleigh`` {b_r_per_m};
    ``L_m``, ``alpha_per_m``, ``gamma_per_w_m``, ``f_R``, ``raman_model``, ``lambda0_m``, ``temperature_k``;
    ``t_gate_s``, ``reflections`` ([[R_i, ε_i], ...], may be empty), ``gbs_filter`` {shape: flat|gaussian,
    fwhm_hz}, ``wdm_isolation_lin`` (1 = none), ``mode`` {shape: gaussian, t_fwhm_s | exponential, tau_eff_s}.
    """

    kind_name, table_names = "stray_light_table", ("stray_light",)
    toggle_names = ("stray_rayleigh", "stray_brillouin", "stray_raman", "stray_reflections")

    def __init__(self, params, name):
        name = name or self.kind_name
        _check_raman(params, name)
        check_probe_groups(params, name)
        require(params, "p_peak_w", "duty", "L_m", "alpha_per_m", "t_gate_s", "reflections", "gbs_filter",
                "wdm_isolation_lin", "mode", block=name)
        require(params["gbs_filter"], "shape", "fwhm_hz", block=f"{name}.gbs_filter")
        _check_mode(params["mode"], name)
        _ascending(name, "p_peak_w", params["p_peak_w"])
        super().__init__(self.kind_name, name, params)

    def _build(self):
        p = self.params
        rows = [[stray_photons(p, self._toggles, pk * float(p["duty"]))[c] for c in STRAY_COLUMNS]
                for pk in np.asarray(p["p_peak_w"], dtype=float)]
        return {"stray_light": Table("stray_light", {"p_peak_w": p["p_peak_w"]}, np.array(rows), "photons",
                                     {"columns": list(STRAY_COLUMNS), "sources": ["T10 eq. (9)", "T07 §5.2", "T17 §6"]})}


class SbsCheck(BaseBlock):
    """CHK ``sbs_check`` (T07 §5.1, V8): warns ``sbs_gain`` if G_B > ``gain_warn`` and ``sbs_threshold``
    if P ≥ P_th (never stops a run). Params (required, D34): ``brillouin`` group, ``alpha_per_m``,
    ``L_m``, ``power_w`` (CW or the power per Brillouin linewidth of a pulse train), ``laser_linewidth_hz``,
    ``gain_warn`` (5 in T07). Product: :attr:`result` = {p_th_w, g_b} after :meth:`check`."""

    def __init__(self, params, name):
        name = name or "sbs_check"
        require(params, "brillouin", "alpha_per_m", "L_m", "power_w", "laser_linewidth_hz", "gain_warn", block=name)
        require(params["brillouin"], *_BRILLOUIN_KEYS, block=f"{name}.brillouin")
        super().__init__("sbs_check", name, params)
        self.chk_log = ChkLog()

    def __call__(self, state, params):
        return state

    def check(self, log: ChkLog | None = None) -> dict[str, float]:
        p, br = self.params, self.params["brillouin"]
        kw = dict(g_b_m_per_w=br["g_b_m_per_w"], a_eff_m2=br["a_eff_m2"], alpha_per_m=p["alpha_per_m"],
                  length_m=p["L_m"], laser_linewidth_hz=p["laser_linewidth_hz"], dnu_b_hz=br["dnu_b_hz"])
        pth, gb = sbs_threshold_w(**kw), float(brillouin_gain_exponent(p["power_w"], **kw))
        log = log or self.chk_log
        if gb > float(p["gain_warn"]):
            log.warn("sbs_gain", self.name, f"Brillouin gain exponent G_B = {gb:.3g} > {p['gain_warn']} "
                     f"(P = {p['power_w']:.3g} W, L = {p['L_m']} m): amplified spontaneous Stokes / RIN",
                     g_b=gb, gain_warn=float(p["gain_warn"]))
        if float(p["power_w"]) >= pth:
            log.warn("sbs_threshold", self.name, f"P = {p['power_w']:.3g} W ≥ SBS threshold {pth:.3g} W", p_th_w=pth)
        self.result = {"p_th_w": pth, "g_b": gb}
        return self.result


@register_block("raman_noise_table")
def make_raman_noise_table(params, *, name=None) -> Block:
    """Factory of ``raman_noise_table`` (:class:`RamanNoiseTable`)."""
    return RamanNoiseTable(params, name)


@register_block("gawbs_table")
def make_gawbs_table(params, *, name=None) -> Block:
    """Factory of ``gawbs_table`` (:class:`GawbsTable`)."""
    return GawbsTable(params, name)


@register_block("brillouin_table")
def make_brillouin_table(params, *, name=None) -> Block:
    """Factory of ``brillouin_table`` (:class:`BrillouinTable`)."""
    return BrillouinTable(params, name)


@register_block("stray_light_table")
def make_stray_light_table(params, *, name=None) -> Block:
    """Factory of ``stray_light_table`` (:class:`StrayLightTable`)."""
    return StrayLightTable(params, name)


@register_block("sbs_check")
def make_sbs_check(params, *, name=None) -> Block:
    """Factory of the CHK ``sbs_check`` (:class:`SbsCheck`)."""
    return SbsCheck(params, name)
