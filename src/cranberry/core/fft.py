"""Fourier transforms in the project convention (notes/conventions.md §8).

This is the ONLY module allowed to call ``numpy.fft`` / ``scipy.fft`` / ``scipy.signal`` spectral
routines directly (PLAN §7.3). Every other module goes through these wrappers so that the sign and
normalization are fixed in one place.

Convention (Agrawal):

    continuous:  Ã(Ω) = ∫ A(T) e^{+iΩT} dT,      A(T) = (1/2π) ∫ Ã(Ω) e^{−iΩT} dΩ
                 ∂/∂T  →  −iΩ
    numerical:   grid T_n = t0 + n Δt, n = 0..N−1   (conv §8 uses t0 = 0)
                 Ã(Ω_k) = N Δt · numpy.fft.ifft(A)[k] · e^{+iΩ_k t0}
                 Ω_k    = 2π · numpy.fft.fftfreq(N, Δt)[k]          [rad/s], FFT (unshifted) order
                 A      = numpy.fft.fft(Ã · e^{−iΩ t0}) / (N Δt)

With t0 = 0 these are exactly the conv §8 formulas. ``t0`` (time of sample 0) is an extension so
that a grid centred on T = 0 (t0 = −(N//2)·Δt) gives the absolute spectral phase of the continuous
transform instead of a linear-phase-shifted one.

Ω is the optical sideband angular frequency (offset from the carrier ω0, conv §1). A field shifted to
ω0 + Ω_c has envelope ∝ e^{−iΩ_c T} (conv §4) and appears at Ω = +Ω_c in Ã.

Discrete Parseval (exact for the DFT):  Σ_n |A_n|² Δt = (1/2π) Σ_k |Ã_k|² ΔΩ,  ΔΩ = 2π/(N Δt).

PSDs are one-sided (conv §7): ∫_0^∞ S_x(f) df = Var(x); see :func:`psd_onesided`.
"""

from __future__ import annotations

import numpy as np
import scipy.signal

__all__ = [
    "omega_grid",
    "freq_grid",
    "domega",
    "time_grid",
    "fft_agrawal",
    "ifft_agrawal",
    "fftshift",
    "ifftshift",
    "psd_onesided",
]


def _check_grid(n: int, dt: float) -> None:
    if int(n) != n or n < 1:
        raise ValueError(f"number of samples must be a positive integer, got {n!r}")
    if not (np.isfinite(dt) and dt > 0):
        raise ValueError(f"time step dt must be finite and > 0 [s], got {dt!r}")


def omega_grid(n: int, dt: float) -> np.ndarray:
    """Sideband angular-frequency grid Ω_k = 2π·fftfreq(n, dt) [rad/s], in FFT (unshifted) order."""
    _check_grid(n, dt)
    return 2.0 * np.pi * np.fft.fftfreq(int(n), dt)


def freq_grid(n: int, dt: float) -> np.ndarray:
    """Sideband frequency grid f_k = Ω_k/2π = fftfreq(n, dt) [Hz], in FFT (unshifted) order."""
    _check_grid(n, dt)
    return np.fft.fftfreq(int(n), dt)


def domega(n: int, dt: float) -> float:
    """Angular-frequency spacing ΔΩ = 2π/(n·dt) [rad/s]."""
    _check_grid(n, dt)
    return 2.0 * np.pi / (int(n) * dt)


def time_grid(n: int, dt: float, t0: float = 0.0) -> np.ndarray:
    """Time grid T_k = t0 + k·dt [s], k = 0..n−1 (conv §8 uses t0 = 0)."""
    _check_grid(n, dt)
    return t0 + dt * np.arange(int(n))


def _phase_along(axis_len: int, dt: float, t0: float, ndim: int, axis: int, sign: float) -> np.ndarray:
    shape = [1] * ndim
    shape[axis] = axis_len
    return np.exp(sign * 1j * omega_grid(axis_len, dt) * t0).reshape(shape)


def fft_agrawal(a: np.ndarray, dt: float, *, t0: float = 0.0, axis: int = -1) -> np.ndarray:
    """Forward transform Ã(Ω_k) = N·Δt·ifft(A)[k]·e^{iΩ_k t0} (conv §8), along ``axis``.

    Args:
        a: samples A(T_n) on the grid T_n = t0 + n·dt (any unit U).
        dt: time step [s].
        t0: time of sample 0 [s]. Default 0 reproduces conv §8 exactly.
        axis: transform axis.

    Returns:
        Complex array Ã(Ω_k) [U·s] in FFT order; pair with :func:`omega_grid`.
    """
    a = np.asarray(a)
    n = a.shape[axis]
    _check_grid(n, dt)
    out = (n * dt) * np.fft.ifft(a, axis=axis)
    if t0 != 0.0:
        out = out * _phase_along(n, dt, t0, out.ndim, axis, +1.0)
    return out


def ifft_agrawal(a_omega: np.ndarray, dt: float, *, t0: float = 0.0, axis: int = -1) -> np.ndarray:
    """Inverse transform A(T_n) = fft(Ã·e^{−iΩ t0})[n]/(N·Δt) (conv §8); exact inverse of fft_agrawal."""
    a_omega = np.asarray(a_omega)
    n = a_omega.shape[axis]
    _check_grid(n, dt)
    if t0 != 0.0:
        a_omega = a_omega * _phase_along(n, dt, t0, a_omega.ndim, axis, -1.0)
    return np.fft.fft(a_omega, axis=axis) / (n * dt)


def fftshift(x: np.ndarray, axes: int | tuple[int, ...] | None = None) -> np.ndarray:
    """``numpy.fft.fftshift`` (FFT order → ascending Ω), re-exported so callers avoid numpy.fft."""
    return np.fft.fftshift(x, axes=axes)


def ifftshift(x: np.ndarray, axes: int | tuple[int, ...] | None = None) -> np.ndarray:
    """``numpy.fft.ifftshift`` (ascending Ω → FFT order)."""
    return np.fft.ifftshift(x, axes=axes)


def psd_onesided(
    x: np.ndarray,
    fs_hz: float,
    *,
    nperseg: int | None = None,
    window: str = "hann",
    detrend: str | bool = "constant",
) -> tuple[np.ndarray, np.ndarray]:
    """One-sided power spectral density by Welch's method (conv §7: ∫_0^{fs/2} S df ≈ Var(x)).

    Thin wrapper around ``scipy.signal.welch(scaling="density", return_onesided=True)``; the
    window power normalization makes Σ S·Δf equal the (detrended) mean power up to window leakage.

    Args:
        x: real time series (unit U), sampled at ``fs_hz``. Complex input is rejected because a
            one-sided PSD is not defined for it.
        fs_hz: sample rate [Hz].
        nperseg: Welch segment length (default: min(len(x), 256), scipy's default).
        window: scipy window name.
        detrend: scipy detrend option ("constant" removes the mean of each segment).

    Returns:
        (f [Hz], S [U²/Hz]).
    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError("psd_onesided requires a real signal (one-sided PSD, conv §7)")
    if not (np.isfinite(fs_hz) and fs_hz > 0):
        raise ValueError(f"fs_hz must be finite and > 0, got {fs_hz!r}")
    f, s = scipy.signal.welch(
        x,
        fs=fs_hz,
        window=window,
        nperseg=nperseg,
        detrend=detrend,
        return_onesided=True,
        scaling="density",
    )
    return f, s
