"""Test fixtures of the CW OPO (WP-B6): Ref B and Ref C parameter sets for ``opo_cw_spectrum`` and the
measured-noise helper of the T02 validation tests. Not presets (those are WP-I3 / WP-I4 with full
provenance); values carry their source tags here because the VTs assert against them.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from cranberry.blocks.source.opo_cw import make_opo_cw_spectrum
from cranberry.offline import gaussian

__all__ = ["TOGGLES_ON", "REF_B", "REF_B_ETA_D", "REF_C", "REF_C_ETA_D", "opo", "measured_db", "db"]

TOGGLES_ON = {"frequency_dependence": True, "pump_dependent_loss": True, "detuning": False}


def ref_b_loss_ic(x):
    """Ref B Eq. (3): L_ic(x) = 0.00249 + 0.00222x, x = √(P/P_th) (T02 §3.4, read with √ as in the paper)."""
    return 0.00249 + 0.00222 * x


#: Ref B (Takeno 2007) at 100 mW, f = 1 MHz (T02 §3, D48): T = 0.123 [Ref B §2], l ≈ 0.5 m [Ref B §2],
#: P_th = 180 mW [Ref B §4]; G = 18.7 deliberately not used (D51-10).
REF_B: dict[str, Any] = {
    "cavity_model": "mirror", "T_oc": 0.123, "round_trip_m": 0.5, "loss_ic": ref_b_loss_ic,
    "pump_mode": "power", "pump_w": 0.100, "p_th_w": 0.180, "detuning_rad_s": 0.0, "f_s_hz": [1.0e6],
    "phi_opa0_rad": 0.0, "segment": "opo", "gain_table_grid": None, "coef_table_x": None, "toggles": TOGGLES_ON,
}
#: ζ ξ² η = 0.99 · 0.988² · 0.998 (propagation, visibility², QE) [Ref B §3]; applied downstream (G06).
REF_B_ETA_D = 0.99 * 0.988**2 * 0.998

#: Ref C (Verclas 2026) monolithic ppKTP: fitted linewidth 1.75 GHz FWHM (T02 §0-6, Tohermes 2024),
#: η_esc ≈ 0.99 [Tohermes 2024]; x is solved from the 30 MHz anti-squeezing in V-C1.
REF_C: dict[str, Any] = {
    "cavity_model": "linewidth", "fwhm_hz": 1.75e9, "escape_efficiency": 0.99,
    "pump_mode": "x", "x": 0.56, "detuning_rad_s": 0.0, "f_s_hz": [30e6, 0.25e9, 0.69e9, 1.21e9],
    "phi_opa0_rad": 0.0, "segment": "opo", "gain_table_grid": None, "coef_table_x": None, "toggles": TOGGLES_ON,
}
#: total η_tot = 0.79 (fitted loss 21 %, building B after the fibre, T02 §6.2) = ρ · REF_C_ETA_D
REF_C_ETA_D = 0.79 / 0.99


def db(r):
    return 10.0 * np.log10(r)


def opo(base: dict[str, Any], **kw: Any):
    """``opo_cw_spectrum`` block from ``base`` with overrides (``toggles`` merged)."""
    tg = dict(base["toggles"], **kw.pop("toggles", {}))
    return make_opo_cw_spectrum({**base, **kw, "toggles": tg})


def measured_db(block, eta_d: float, theta_rad: float, *, average: str = "small_angle") -> tuple[np.ndarray, np.ndarray]:
    """(R′−, R′+) [dB] per f_s: project, pure loss η_d (as the LossRegistry would), phase jitter θ̃
    (``small_angle`` = Ref B Eq. (1); ``gaussian`` = conv §6 exact)."""
    vm, vp, _ = block.source_output.project(None)
    rm, rp = (2.0 * (eta_d * v + 0.5 * (1.0 - eta_d)) for v in (vm, vp))
    if average == "small_angle":
        return (db(gaussian.jitter_average_small_angle(rm, rp, theta_rad)),
                db(gaussian.jitter_average_small_angle(rp, rm, theta_rad)))
    chi = gaussian.chi_gaussian(theta_rad)
    return (db(gaussian.jitter_average(rm, rp, chi, delta0_rad=0.0)),
            db(gaussian.jitter_average(rm, rp, chi, delta0_rad=np.pi / 2)))
