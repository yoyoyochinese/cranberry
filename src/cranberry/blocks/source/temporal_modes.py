"""Temporal-mode helpers shared by the pulsed blocks (WP-B11; T03 §2.4–§2.6, T17 §2.2, §5.2, §6.3).

Not a block. Everything acts on a ``temporal`` :class:`SourceOutput` (G07) or on a mode function u(T)
sampled on a uniform grid (∫|u|² dT = 1, Riemann weights Δt as in ``SourceOutput.lo_mode_cov``).

LO-mode convention (conv §2, §5; ``SourceOutput.lo_mode_cov``): â_LO = ∫u* â dT = Σ c_n â_n with
c_n = ∫u* f_n dT. A real quadrature combination w·X = Σ (w_q,n q_n + w_p,n p_n) is the quadrature
θ = π/2 of the mode with coefficients u = Σ_n (w_p,n − i w_q,n) f_n (then p_LO = w·X). Multiplying u by
e^{−iα} moves every measured angle by +α, so :func:`optimal_lo_mode` returns u·e^{−iφ_OPA,0/2}: its
squeezed quadrature sits at θ_sq,0 = φ_OPA,0/2 + π/2 and ``project`` gives θ_off = 0 for it (D55-4).

Delay / dispersion (conv §8, T03 §2.5): ũ_s(Ω) = ũ(Ω) e^{i(ΩΔt + φ2Ω²/2 + φ3Ω³/6)}, i.e. Δt > 0 is a
later pulse, applied with ``core.fft`` on a zero-padded grid (no wrap-around for |Δt| < pad window).
"""

from __future__ import annotations

import numpy as np

from cranberry.core import fft

__all__ = ["optimal_lo_mode", "tau_eff_s", "tail_energy", "mode_spectrum", "transform_mode", "overlap"]


def _dt(grid: np.ndarray) -> float:
    d = np.diff(grid)
    if not np.allclose(d, d[0], rtol=1e-9, atol=0.0):
        raise ValueError("temporal-mode helpers need a uniform grid")
    return float(d[0])


def optimal_lo_mode(so) -> tuple[np.ndarray, float]:
    """(u, V_min): LO mode on ``so.grid`` that sees the smallest quadrature variance (vacuum 1/2) and
    that variance (min eigenvalue of ``so.cov``; T17 §2.2 "best squeezing R_− = min eigenvalue")."""
    if so.basis != "temporal":
        raise ValueError("optimal_lo_mode needs a temporal SourceOutput")
    lam, vec = np.linalg.eigh(np.asarray(so.cov))
    n = so.n_modes
    w = vec[:, 0]
    d = (w[n:] - 1j * w[:n]) * np.exp(-0.5j * so.phi_opa0_rad)
    u = d @ so.modes
    w_t = np.gradient(so.grid)
    u = u / np.sqrt(np.sum(np.abs(u) ** 2 * w_t))
    if np.sum(u.real * w_t) < 0:  # sign is free (a mode and −mode are the same LO up to θ → θ + π)
        u = -u
    return u, float(lam[0])


def tau_eff_s(u: np.ndarray, grid: np.ndarray) -> float:
    """Effective mode time τ_eff = |∫u dT|² (conv §2) [s]."""
    return float(abs(np.sum(u) * _dt(grid)) ** 2)


def tail_energy(u: np.ndarray, grid: np.ndarray, t_s: float) -> float:
    """ε_tail(T) = ∫_T^∞ |u|² dT for a piecewise-constant u on bins centred on ``grid`` (exact)."""
    dt = _dt(grid)
    lo = grid - 0.5 * dt
    frac = np.clip((lo + dt - t_s) / dt, 0.0, 1.0)
    return float(np.sum(np.abs(u) ** 2 * frac) * dt)


def mode_spectrum(u: np.ndarray, grid: np.ndarray, pad: int) -> tuple[np.ndarray, np.ndarray]:
    """(f [Hz] ascending, |ũ(f)|² [s]) on a grid zero-padded ``pad``× (conv §8 sign; discrete Parseval
    Σ|ũ|²Δf = ∫|u|²dT holds exactly, and |ũ(0)|² = τ_eff)."""
    dt = _dt(grid)
    n = u.size * int(pad)
    a = np.zeros(n, dtype=complex)
    a[: u.size] = u
    spec = fft.fftshift(np.abs(fft.fft_agrawal(a, dt, t0=float(grid[0]))) ** 2)
    return fft.fftshift(fft.freq_grid(n, dt)), spec


def transform_mode(u: np.ndarray, grid: np.ndarray, *, dt_s: float = 0.0, phi2_s2: float = 0.0,
                   phi3_s3: float = 0.0, pad: int = 2) -> np.ndarray:
    """u_s(T) on ``grid``: u delayed by Δt and given the GDD/TOD φ2, φ3 (module docstring). The grid is
    zero-padded ``pad``× on both sides so that the shifted mode does not wrap around."""
    dt = _dt(grid)
    n, m = u.size, u.size * int(pad)
    a = np.zeros(n + 2 * m, dtype=complex)
    a[m : m + n] = u
    t0 = float(grid[0]) - m * dt
    om = fft.omega_grid(a.size, dt)
    spec = fft.fft_agrawal(a, dt, t0=t0) * np.exp(1j * (om * dt_s + 0.5 * phi2_s2 * om**2 + phi3_s3 * om**3 / 6.0))
    return fft.ifft_agrawal(spec, dt, t0=t0)[m : m + n]


def overlap(u_lo: np.ndarray, u_s: np.ndarray, grid: np.ndarray) -> complex:
    """κ = ∫ u_LO* u_s dT (T03 §1); η_t = |κ|²."""
    return complex(np.sum(np.conj(u_lo) * u_s) * _dt(grid))
