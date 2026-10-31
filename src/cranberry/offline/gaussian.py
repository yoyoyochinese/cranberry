"""Gaussian-state covariance operations (PLAN §4.2, conv §2–§3, §5–§6). Owner: WP-F3.

Units and ordering (conv §2, §9, D54)
    ħ = 1, â = (q̂ + i p̂)/√2, [q̂, p̂] = i, vacuum covariance (1/2)·𝟙. Covariance matrices
    σ_ij = ⟨{Δx̂_i, Δx̂_j}⟩/2 of N modes use the ordering ``xxpp``: (q_1, …, q_N, p_1, …, p_N).
    Every function accepts leading batch dimensions (``cov`` of shape (..., 2N, 2N)), e.g. the
    per-sideband (n_f, 2, 2) covariances of a ``cw_sideband`` ``SourceOutput``.

Symplectic matrices act on the quadrature vector in the Heisenberg picture, X → M X, so a state
transforms as σ → M σ Mᵀ (:func:`apply_symplectic`). M is symplectic if M Ω Mᵀ = Ω with
Ω = [[0, 𝟙], [−𝟙, 0]] (xxpp, :func:`symplectic_form`).

Channels and gates (conv references in each docstring)
    * pure loss η:                σ → η σ + (1 − η)/2 · 𝟙                         (conv §2)
    * phase-insensitive noise N:  σ → σ + N · 𝟙   (N photons/mode, ΔV/V_SN = 2N)  (conv §2)
    * phase shift ψ (A → A e^{iψ}, conv §4): the squeezing angle rotates by +ψ   (:func:`rotation`)
    * squeezer S(ξ), ξ = r e^{iφ_ξ}: S†âS = â cosh r − â† e^{iφ_ξ} sinh r, θ_sq = φ_ξ/2 (conv §3)
    * beam splitter, real form of conv §5 (b̂_i = √η â_i + √(1−η) â_j, b̂_j = √(1−η) â_i − √η â_j)

Homodyne variance at angle θ: x̂_θ = q̂ cos θ + p̂ sin θ (conv §2), θ = φ_LO in the signal frame
(conv §5). Quasi-static phase jitter (conv §6, T04 eq. (1)):
V′ = (V_− + V_+)/2 − (V_+ − V_−)/2 · cos(2Δ0) · χ, χ = ⟨cos 2δ⟩ (Gaussian e^{−2σ²}, sinusoidal
dither J0(2m), independent components multiply).
"""

from __future__ import annotations

import numpy as np
from scipy.special import j0

__all__ = [
    "vacuum",
    "symplectic_form",
    "n_modes_of",
    "apply_symplectic",
    "rotation",
    "rotate",
    "squeezer",
    "squeeze",
    "squeezed_vacuum",
    "beam_splitter_matrix",
    "beam_splitter",
    "pure_loss",
    "add_noise",
    "homodyne_variance",
    "quadrature_variance",
    "is_symplectic",
    "is_physical",
    "symplectic_eigenvalues",
    "cov_from_noise_levels",
    "noise_levels_from_cov",
    "chi_gaussian",
    "chi_sinusoid",
    "chi_timeseries",
    "jitter_average",
    "jitter_average_small_angle",
]

_V_VAC = 0.5  # conv §2 (same value as core.units.V_VAC; kept local to avoid a hard import cycle)


# --------------------------------------------------------------------------------------------------
# basics
# --------------------------------------------------------------------------------------------------


def vacuum(n_modes: int = 1) -> np.ndarray:
    """Vacuum covariance (1/2)·𝟙 of ``n_modes`` modes (2N×2N, xxpp)."""
    _check_n(n_modes)
    return _V_VAC * np.eye(2 * n_modes)


def symplectic_form(n_modes: int) -> np.ndarray:
    """Ω = [[0, 𝟙_N], [−𝟙_N, 0]] for the xxpp ordering ([q_i, p_j] = i δ_ij)."""
    _check_n(n_modes)
    n = n_modes
    om = np.zeros((2 * n, 2 * n))
    om[:n, n:] = np.eye(n)
    om[n:, :n] = -np.eye(n)
    return om


def n_modes_of(cov: np.ndarray) -> int:
    """Number of modes N of a (..., 2N, 2N) covariance (validates squareness and even size)."""
    cov = np.asarray(cov)
    if cov.ndim < 2 or cov.shape[-1] != cov.shape[-2] or cov.shape[-1] % 2:
        raise ValueError(f"covariance must have shape (..., 2N, 2N), got {cov.shape}")
    return cov.shape[-1] // 2


def apply_symplectic(cov: np.ndarray, s: np.ndarray) -> np.ndarray:
    """σ → S σ Sᵀ (Heisenberg-picture matrix S on the quadrature vector, batch dims allowed)."""
    cov = np.asarray(cov, dtype=float)
    s = np.asarray(s, dtype=float)
    n_modes_of(cov)
    if s.shape[-2:] != cov.shape[-2:]:
        raise ValueError(f"symplectic matrix shape {s.shape} does not match covariance {cov.shape}")
    return s @ cov @ np.swapaxes(s, -1, -2)


def _check_n(n_modes: int) -> None:
    if int(n_modes) != n_modes or n_modes < 1:
        raise ValueError(f"n_modes must be a positive int, got {n_modes!r}")


def _check_mode(mode: int, n_modes: int) -> None:
    if not (0 <= mode < n_modes):
        raise ValueError(f"mode index {mode} out of range for {n_modes} modes")


def _embed(m2: np.ndarray, n_modes: int, mode: int) -> np.ndarray:
    """Embed a single-mode 2×2 (q, p) matrix acting on ``mode`` into a 2N×2N xxpp identity."""
    _check_mode(mode, n_modes)
    out = np.eye(2 * n_modes)
    k, kp = mode, mode + n_modes
    out[k, k], out[k, kp] = m2[0, 0], m2[0, 1]
    out[kp, k], out[kp, kp] = m2[1, 0], m2[1, 1]
    return out


# --------------------------------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------------------------------


def rotation(theta_rad: float | np.ndarray, n_modes: int = 1, mode: int | None = None) -> np.ndarray:
    """Phase-shift matrix for â → â e^{iθ} (conv §4 propagation A → A e^{iψ}).

    q' = q cos θ − p sin θ, p' = q sin θ + p cos θ; a squeezing angle θ_sq becomes θ_sq + θ.
    ``mode=None``: every mode is rotated (``theta_rad`` scalar or length-N array); otherwise only
    ``mode`` is rotated by the scalar ``theta_rad``.
    """
    _check_n(n_modes)
    if mode is not None:
        c, s = np.cos(float(theta_rad)), np.sin(float(theta_rad))
        return _embed(np.array([[c, -s], [s, c]]), n_modes, mode)
    th = np.broadcast_to(np.asarray(theta_rad, dtype=float), (n_modes,))
    c, s = np.diag(np.cos(th)), np.diag(np.sin(th))
    return np.block([[c, -s], [s, c]])


def rotate(cov: np.ndarray, theta_rad: float | np.ndarray, mode: int | None = None) -> np.ndarray:
    """Apply :func:`rotation` to ``cov``."""
    return apply_symplectic(cov, rotation(theta_rad, n_modes_of(cov), mode))


def squeezer(r: float, phi_xi_rad: float, n_modes: int = 1, mode: int = 0) -> np.ndarray:
    """Heisenberg matrix of S(ξ), ξ = r e^{iφ_ξ} (conv §3): S†âS = â cosh r − â† e^{iφ_ξ} sinh r.

    q' = q (cosh r − sinh r cos φ) − p sinh r sin φ, p' = −q sinh r sin φ + p (cosh r + sinh r cos φ).
    On vacuum: Var(x_θ) = ½[e^{−2r} cos²(θ − φ_ξ/2) + e^{2r} sin²(θ − φ_ξ/2)], θ_sq = φ_ξ/2.
    """
    if r < 0:
        raise ValueError(f"squeezing parameter r must be >= 0 (conv §3), got {r}")
    ch, sh = np.cosh(r), np.sinh(r)
    c, s = np.cos(phi_xi_rad), np.sin(phi_xi_rad)
    m2 = np.array([[ch - sh * c, -sh * s], [-sh * s, ch + sh * c]])
    _check_n(n_modes)
    return _embed(m2, n_modes, mode)


def squeeze(cov: np.ndarray, r: float, phi_xi_rad: float, mode: int = 0) -> np.ndarray:
    """Apply :func:`squeezer` to ``mode`` of ``cov``."""
    return apply_symplectic(cov, squeezer(r, phi_xi_rad, n_modes_of(cov), mode))


def squeezed_vacuum(r: float, theta_sq_rad: float) -> np.ndarray:
    """Single-mode squeezed vacuum S(r e^{2iθ_sq})|0⟩: V_− = e^{−2r}/2 at θ_sq (conv §3)."""
    return squeeze(vacuum(1), r, 2.0 * theta_sq_rad)


def beam_splitter_matrix(eta: float, n_modes: int, i: int, j: int) -> np.ndarray:
    """Beam splitter of power transmission η between modes i and j, real form of conv §5.

    b̂_i = √η â_i + √(1−η) â_j, b̂_j = √(1−η) â_i − √η â_j (η = 1/2 gives conv §5:
    b̂_1 = (â_s + â_LO)/√2, b̂_2 = (â_s − â_LO)/√2). Acts identically on q and p (passive).
    """
    if not (0.0 <= eta <= 1.0):
        raise ValueError(f"beam splitter transmission must be in [0, 1], got {eta}")
    _check_n(n_modes)
    _check_mode(i, n_modes)
    _check_mode(j, n_modes)
    if i == j:
        raise ValueError("beam splitter needs two distinct modes")
    t, rr = np.sqrt(eta), np.sqrt(1.0 - eta)
    o = np.eye(n_modes)
    o[i, i], o[i, j], o[j, i], o[j, j] = t, rr, rr, -t
    z = np.zeros_like(o)
    return np.block([[o, z], [z, o]])


def beam_splitter(cov: np.ndarray, eta: float, i: int, j: int) -> np.ndarray:
    """Beam splitter of power transmission ``eta`` between modes i and j (conv §5 real form)."""
    return apply_symplectic(cov, beam_splitter_matrix(eta, n_modes_of(cov), i, j))


# --------------------------------------------------------------------------------------------------
# channels
# --------------------------------------------------------------------------------------------------


def _per_mode(value: float | np.ndarray, n_modes: int, what: str) -> tuple[np.ndarray, bool]:
    v = np.asarray(value, dtype=float)
    if v.ndim == 0:
        return v, True
    if v.shape != (n_modes,):
        raise ValueError(f"{what} must be a scalar or a length-{n_modes} array (per mode), got {v.shape}")
    return v, False


def pure_loss(cov: np.ndarray, eta: float | np.ndarray) -> np.ndarray:
    """σ → η σ + (1 − η)/2 𝟙 (conv §2); per mode if ``eta`` is an array of length N.

    Per mode: σ_ij → √(η_i η_j) σ_ij + δ_ij (1 − η_i)/2 (beam splitter with vacuum, T12 §1.2).
    Applying the registry's total η (``LossRegistry.eta``) once is equivalent to applying each
    series element in turn (transmissions multiply, T12 V2).
    """
    cov = np.asarray(cov, dtype=float)
    n = n_modes_of(cov)
    e, scalar = _per_mode(eta, n, "eta")
    if np.any((e < 0.0) | (e > 1.0)):
        raise ValueError(f"eta must be in [0, 1], got {eta!r}")
    if scalar:
        return e * cov + (1.0 - e) * _V_VAC * np.eye(2 * n)
    g = np.sqrt(np.concatenate([e, e]))
    return cov * np.outer(g, g) + np.diag(np.concatenate([1.0 - e, 1.0 - e]) * _V_VAC)


def add_noise(cov: np.ndarray, n_photons: float | np.ndarray) -> np.ndarray:
    """Phase-insensitive added noise: σ → σ + N 𝟙 (N photons per mode, ΔV/V_SN = 2N, conv §2).

    ``n_photons`` is a scalar or a length-N array (per mode); N = ∫ n(f)|ũ(f)|² df (conv §2).
    """
    cov = np.asarray(cov, dtype=float)
    n = n_modes_of(cov)
    v, scalar = _per_mode(n_photons, n, "n_photons")
    if np.any(v < 0.0):
        raise ValueError(f"added noise must be >= 0, got {n_photons!r}")
    if scalar:
        return cov + v * np.eye(2 * n)
    return cov + np.diag(np.concatenate([v, v]))


# --------------------------------------------------------------------------------------------------
# measurement
# --------------------------------------------------------------------------------------------------


def homodyne_variance(cov: np.ndarray, theta_rad: float | np.ndarray, mode: int = 0) -> np.ndarray:
    """Var(x̂_θ) of ``mode``, x̂_θ = q̂ cos θ + p̂ sin θ (conv §2); broadcasts over batch dims and θ."""
    cov = np.asarray(cov, dtype=float)
    n = n_modes_of(cov)
    _check_mode(mode, n)
    k, kp = mode, mode + n
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return c * c * cov[..., k, k] + s * s * cov[..., kp, kp] + 2.0 * c * s * cov[..., k, kp]


def quadrature_variance(cov: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Var(wᵀX) = wᵀ σ w for a real weight vector w (length 2N, xxpp; Σw² = 1 for a unit LO mode)."""
    cov = np.asarray(cov, dtype=float)
    w = np.asarray(w, dtype=float)
    if w.shape != cov.shape[-1:]:
        raise ValueError(f"weight vector shape {w.shape} does not match covariance {cov.shape}")
    return np.einsum("i,...ij,j->...", w, cov, w)


# --------------------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------------------


def is_symplectic(s: np.ndarray, tol: float = 1e-9) -> bool:
    """True if S Ω Sᵀ = Ω within ``tol`` (max abs deviation, xxpp ordering)."""
    s = np.asarray(s, dtype=float)
    if s.ndim != 2 or s.shape[0] != s.shape[1] or s.shape[0] % 2:
        raise ValueError(f"symplectic matrix must be 2N×2N, got {s.shape}")
    om = symplectic_form(s.shape[0] // 2)
    return bool(np.max(np.abs(s @ om @ s.T - om)) <= tol)


def is_physical(cov: np.ndarray, tol: float = 1e-12) -> bool:
    """Robertson–Schrödinger condition σ + iΩ/2 ≥ 0 (min eigenvalue ≥ −tol), all batch entries.

    Also requires σ symmetric (within ``tol``·max|σ|).
    """
    cov = np.asarray(cov, dtype=float)
    n = n_modes_of(cov)
    scale = max(1.0, float(np.max(np.abs(cov))))
    if np.max(np.abs(cov - np.swapaxes(cov, -1, -2))) > tol * scale:
        return False
    m = cov + 0.5j * symplectic_form(n)
    return bool(np.min(np.linalg.eigvalsh(m)) >= -tol * scale)


def symplectic_eigenvalues(cov: np.ndarray) -> np.ndarray:
    """Symplectic eigenvalues ν_1 ≤ … ≤ ν_N of σ (moduli of the eigenvalues of iΩσ); vacuum → 1/2.

    Invariant under symplectic transformations; ν_k ≥ 1/2 for physical states (purity: all = 1/2).
    """
    cov = np.asarray(cov, dtype=float)
    n = n_modes_of(cov)
    ev = np.abs(np.linalg.eigvals(1j * symplectic_form(n) @ cov))
    ev = np.sort(ev, axis=-1)
    return ev[..., ::2]  # eigenvalues come in ±ν pairs; sorted moduli are (ν1, ν1, ν2, ν2, ...)


# --------------------------------------------------------------------------------------------------
# noise levels (V_−, V_+, θ_sq) ↔ 2×2 covariance
# --------------------------------------------------------------------------------------------------


def cov_from_noise_levels(v_minus, v_plus, theta_sq_rad) -> np.ndarray:
    """2×2 covariance with minimum variance V_− at θ_sq and V_+ at θ_sq + π/2 (conv §2–§3).

    Variances in conv §2 units (vacuum 1/2). Broadcasts: output shape = broadcast shape + (2, 2).
    """
    vm, vp, th = np.broadcast_arrays(
        np.asarray(v_minus, dtype=float), np.asarray(v_plus, dtype=float), np.asarray(theta_sq_rad, dtype=float)
    )
    if np.any(vm > vp):
        raise ValueError("cov_from_noise_levels: need V_- <= V_+")
    c, s = np.cos(th), np.sin(th)
    out = np.empty(vm.shape + (2, 2))
    out[..., 0, 0] = vm * c * c + vp * s * s
    out[..., 1, 1] = vm * s * s + vp * c * c
    out[..., 0, 1] = out[..., 1, 0] = (vm - vp) * c * s
    return out


def noise_levels_from_cov(cov: np.ndarray) -> tuple[float, float, float]:
    """Inverse of :func:`cov_from_noise_levels`: (V_−, V_+, θ_sq mod π) of a (..., 2, 2) covariance.

    Var(θ) = a + b cos(2θ − c) with a = (σ_qq + σ_pp)/2; V_∓ = a ∓ |b|; θ_sq ∈ [0, π). For an
    isotropic state (V_− = V_+) θ_sq is undefined and returned as 0. Batch input gives arrays.
    """
    cov = np.asarray(cov, dtype=float)
    if cov.shape[-2:] != (2, 2):
        raise ValueError(f"noise_levels_from_cov needs (..., 2, 2), got {cov.shape}")
    qq, pp, qp = cov[..., 0, 0], cov[..., 1, 1], 0.5 * (cov[..., 0, 1] + cov[..., 1, 0])
    a = 0.5 * (qq + pp)
    d = 0.5 * (qq - pp)
    b = np.hypot(d, qp)
    theta = np.mod(0.5 * (np.arctan2(qp, d) + np.pi), np.pi)
    theta = np.where(b <= 1e-15 * np.abs(a), 0.0, theta)
    if cov.ndim == 2:
        return float(a - b), float(a + b), float(theta)
    return a - b, a + b, theta


# --------------------------------------------------------------------------------------------------
# phase jitter → measured variance (conv §6, T04 §3.1)
# --------------------------------------------------------------------------------------------------


def chi_gaussian(sigma_rad):
    """χ = ⟨cos 2δ⟩ = e^{−2σ²} for zero-mean Gaussian δ of rms σ (conv §6)."""
    s = np.asarray(sigma_rad, dtype=float)
    return np.exp(-2.0 * s * s)


def chi_sinusoid(m_rad):
    """χ = J0(2m) for a sinusoidal dither δ = m sin(2πf_m t) of peak amplitude m (conv §6)."""
    return j0(2.0 * np.asarray(m_rad, dtype=float))


def chi_timeseries(delta_rad: np.ndarray, axis: int = -1):
    """χ = time average ⟨cos 2δ(t)⟩ of a residual-error series (any distribution, T04 §13)."""
    return np.mean(np.cos(2.0 * np.asarray(delta_rad, dtype=float)), axis=axis)


def jitter_average(v_minus, v_plus, chi, *, delta0_rad):
    """V′ = (V_− + V_+)/2 − (V_+ − V_−)/2 · cos(2Δ0) · χ (conv §6, exact for symmetric δ).

    ``delta0_rad`` is the set-point offset Δ0 = Δ_set + δ0 from the squeezed quadrature
    (0 → squeezing, π/2 → anti-squeezing; conv §5). χ of independent components multiply.
    Works for absolute variances (vacuum 1/2) and for shot-noise-normalised R± alike (linear).
    """
    vm = np.asarray(v_minus, dtype=float)
    vp = np.asarray(v_plus, dtype=float)
    return 0.5 * (vm + vp) - 0.5 * (vp - vm) * np.cos(2.0 * np.asarray(delta0_rad, dtype=float)) * np.asarray(chi)


def jitter_average_small_angle(v_main, v_other, sigma_rad):
    """Prompt §7-1 small-angle form V′ ≈ V_main cos²σ + V_other sin²σ (Ref A/B, conv §6).

    Squeezing: (V_−, V_+); anti-squeezing: (V_+, V_−). Differs from :func:`jitter_average` with
    χ = e^{−2σ²} by O(σ⁴) (T04 §3.1); kept for the §7-1 cross-check (``average: small_angle``).
    """
    s = np.asarray(sigma_rad, dtype=float)
    return np.asarray(v_main, dtype=float) * np.cos(s) ** 2 + np.asarray(v_other, dtype=float) * np.sin(s) ** 2
