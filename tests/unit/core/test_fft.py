"""core.fft against conv §8: analytic Gaussian transform (incl. sign), round trip, Parseval."""

import numpy as np
import pytest

from cranberry.core import fft

N = 4096
DT = 0.01  # arbitrary time unit; transforms are unit-agnostic
T0 = 1.0  # Gaussian width


def gauss_omega(omega, center=0.0, delay=0.0):
    """∫ exp(−(T−τ)²/2T0²) e^{−iΩc T} e^{+iΩT} dT = √(2π) T0 exp(−(Ω−Ωc)²T0²/2) e^{i(Ω−Ωc)τ}."""
    w = omega - center
    return np.sqrt(2 * np.pi) * T0 * np.exp(-(w**2) * T0**2 / 2) * np.exp(1j * w * delay)


def test_grids_match_numpy_definitions():
    np.testing.assert_array_equal(fft.omega_grid(N, DT), 2 * np.pi * np.fft.fftfreq(N, DT))
    np.testing.assert_array_equal(fft.freq_grid(N, DT), np.fft.fftfreq(N, DT))
    assert fft.domega(N, DT) == pytest.approx(2 * np.pi / (N * DT))
    np.testing.assert_allclose(fft.time_grid(N, DT), DT * np.arange(N))


def test_gaussian_conv8_grid_starting_at_zero():
    """conv §8 exactly (t0 = 0): pulse centred at τ > 0 → e^{+iΩτ} (sign of the +iΩT kernel)."""
    tau = N * DT / 2
    t = fft.time_grid(N, DT)
    a = np.exp(-((t - tau) ** 2) / (2 * T0**2))
    omega = fft.omega_grid(N, DT)
    np.testing.assert_allclose(fft.fft_agrawal(a, DT), gauss_omega(omega, delay=tau), atol=1e-12, rtol=0)


def test_gaussian_centered_grid_and_carrier_sign():
    """Envelope e^{−iΩc T} (field shifted to ω0 + Ωc, conv §4) appears at Ω = +Ωc."""
    t0 = -(N // 2) * DT
    t = fft.time_grid(N, DT, t0)
    omega_c = 5.0
    a = np.exp(-(t**2) / (2 * T0**2)) * np.exp(-1j * omega_c * t)
    omega = fft.omega_grid(N, DT)
    a_w = fft.fft_agrawal(a, DT, t0=t0)
    np.testing.assert_allclose(a_w, gauss_omega(omega, center=omega_c), atol=1e-12, rtol=0)
    assert omega[np.argmax(np.abs(a_w))] == pytest.approx(omega_c, abs=fft.domega(N, DT))


def test_derivative_is_minus_i_omega():
    """∂/∂T ↔ −iΩ (conv §8)."""
    t0 = -(N // 2) * DT
    t = fft.time_grid(N, DT, t0)
    a = np.exp(-(t**2) / (2 * T0**2))
    omega = fft.omega_grid(N, DT)
    da = fft.ifft_agrawal(-1j * omega * fft.fft_agrawal(a, DT, t0=t0), DT, t0=t0)
    np.testing.assert_allclose(da, -(t / T0**2) * a, atol=1e-10)


@pytest.mark.parametrize("t0", [0.0, -(N // 2) * DT, 3.7])
def test_round_trip(t0):
    rng = np.random.default_rng(1)
    a = rng.normal(size=N) + 1j * rng.normal(size=N)
    np.testing.assert_allclose(fft.ifft_agrawal(fft.fft_agrawal(a, DT, t0=t0), DT, t0=t0), a, atol=1e-12, rtol=0)
    a_w = rng.normal(size=N) + 1j * rng.normal(size=N)
    np.testing.assert_allclose(fft.fft_agrawal(fft.ifft_agrawal(a_w, DT, t0=t0), DT, t0=t0), a_w, rtol=1e-12, atol=1e-12)


def test_parseval():
    rng = np.random.default_rng(2)
    a = rng.normal(size=N) + 1j * rng.normal(size=N)
    a_w = fft.fft_agrawal(a, DT)
    lhs = np.sum(np.abs(a) ** 2) * DT
    rhs = np.sum(np.abs(a_w) ** 2) * fft.domega(N, DT) / (2 * np.pi)
    assert rhs == pytest.approx(lhs, rel=1e-12)


def test_axis_argument():
    rng = np.random.default_rng(3)
    a = rng.normal(size=(64, 3)) + 1j * rng.normal(size=(64, 3))
    out = fft.fft_agrawal(a, DT, t0=0.3, axis=0)
    for j in range(3):
        np.testing.assert_allclose(out[:, j], fft.fft_agrawal(a[:, j], DT, t0=0.3), rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(fft.ifft_agrawal(out, DT, t0=0.3, axis=0), a, atol=1e-13)


def test_shift_helpers_and_validation():
    x = np.arange(8)
    np.testing.assert_array_equal(fft.ifftshift(fft.fftshift(x)), x)
    with pytest.raises(ValueError):
        fft.omega_grid(0, DT)
    with pytest.raises(ValueError):
        fft.fft_agrawal(np.ones(4), -1.0)


def test_psd_onesided_white_noise_level_and_variance():
    """One-sided (conv §7): white noise of variance σ² at fs → S = 2σ²/fs, ∫S df = σ²."""
    rng = np.random.default_rng(4)
    fs, sigma = 1.0e3, 0.3
    x = rng.normal(scale=sigma, size=2**18)
    f, s = fft.psd_onesided(x, fs, nperseg=1024)
    assert np.mean(s[1:-1]) == pytest.approx(2 * sigma**2 / fs, rel=0.02)
    assert np.sum(s) * (f[1] - f[0]) == pytest.approx(sigma**2, rel=0.02)


def test_psd_onesided_sine_power():
    fs, amp, f0 = 1.0e4, 0.5, 1234.5
    t = np.arange(2**16) / fs
    f, s = fft.psd_onesided(amp * np.sin(2 * np.pi * f0 * t), fs, nperseg=4096)
    assert np.sum(s) * (f[1] - f[0]) == pytest.approx(amp**2 / 2, rel=0.01)
    with pytest.raises(ValueError):
        fft.psd_onesided(np.ones(16) * 1j, fs)
