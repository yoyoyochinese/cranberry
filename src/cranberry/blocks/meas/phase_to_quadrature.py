"""Phase error → measured quadrature variance (conv §6, T04 §13, G05). Owner: WP-F4.

V′ = (V_− + V_+)/2 − (V_+ − V_−)/2 · cos(2Δ0) · χ, χ = ⟨cos 2δ⟩ (Gaussian: e^{−2σ²}; sinusoidal
dither of amplitude m: J0(2m); independent parts multiply). The out-of-band phase variance from the
``BandLedger`` is added here and only here (G05). ``dither_mode`` (timeseries | analytic) must not be
counted twice. Integration limits per conv §6 (lower: 1/measurement time; upper: flat-spectrum range
of the source, ~κ_c for OPO modes, ~1/τ_p for pulses).
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["AVERAGE_MODES", "phase_to_quadrature"]

AVERAGE_MODES: tuple[str, ...] = ("exact", "small_angle")


def phase_to_quadrature(
    delta_rad: np.ndarray,
    v_minus: np.ndarray,
    v_plus: np.ndarray,
    *,
    hf_variance_rad2: float,
    average: str,
    **options: Any,
) -> dict[str, np.ndarray]:
    """Map δ(t) (and V± from the loss-aware tables) to V′(t), R′± and per-band σ. WP-F4."""
    raise NotImplementedError("blocks.meas.phase_to_quadrature is implemented by WP-F4")
