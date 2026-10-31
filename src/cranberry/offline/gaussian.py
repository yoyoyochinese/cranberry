"""Gaussian-state covariance operations (PLAN §4.2, conv §2–§3). Owner: WP-F3.

Units: ħ = 1, vacuum covariance (1/2)·𝟙, ordering xxpp (see ``offline.source_output``).
Pure loss: σ → η σ + (1 − η)/2 · 𝟙 (conv §2). Squeezed vacuum with r, angle θ_sq:
Var(x_θ) = ½[e^{−2r} cos²(θ − θ_sq) + e^{2r} sin²(θ − θ_sq)] (conv §3).
"""

from __future__ import annotations

import numpy as np

__all__ = ["pure_loss", "beam_splitter", "is_symplectic", "cov_from_noise_levels", "noise_levels_from_cov"]


def pure_loss(cov: np.ndarray, eta: float | np.ndarray) -> np.ndarray:
    """σ → η σ + (1 − η)/2 𝟙 (per mode if ``eta`` is an array of length N). WP-F3."""
    raise NotImplementedError("offline.gaussian.pure_loss is implemented by WP-F3")


def beam_splitter(cov: np.ndarray, eta: float, i: int, j: int) -> np.ndarray:
    """Beam splitter of power transmission ``eta`` between modes i and j (conv §5 real form). WP-F3."""
    raise NotImplementedError("offline.gaussian.beam_splitter is implemented by WP-F3")


def is_symplectic(s: np.ndarray, tol: float = 1e-9) -> bool:
    """True if S Ω Sᵀ = Ω within ``tol`` (xxpp ordering). WP-F3."""
    raise NotImplementedError("offline.gaussian.is_symplectic is implemented by WP-F3")


def cov_from_noise_levels(v_minus: float, v_plus: float, theta_sq_rad: float) -> np.ndarray:
    """2×2 covariance with minimum variance V_− at θ_sq and V_+ at θ_sq + π/2. WP-F3."""
    raise NotImplementedError("offline.gaussian.cov_from_noise_levels is implemented by WP-F3")


def noise_levels_from_cov(cov: np.ndarray) -> tuple[float, float, float]:
    """Inverse of :func:`cov_from_noise_levels`: (V_−, V_+, θ_sq mod π). WP-F3."""
    raise NotImplementedError("offline.gaussian.noise_levels_from_cov is implemented by WP-F3")
