"""Test fixtures of the pulsed OPO (WP-B11): Borealis-type and Liu-type parameter sets for
``pulsed_opo_source`` (T17 §5.1). Not presets (WP-I5 owns ``source.borealis_opo`` /
``source.liu2026_opo`` with full provenance); source tags are kept here because the VTs assert on them.
"""

from __future__ import annotations

from typing import Any

from cranberry.blocks.source.pulsed_opo import make_pulsed_opo_source

__all__ = ["TOGGLES_ON", "BOREALIS", "LIU", "pulsed_opo"]

TOGGLES_ON = {"pulsed": True, "pump_resonant": True, "detuning": False, "include_initial_vacuum": True}

#: Madsen 2022 SI: signal FWHM 66 MHz, finesse 57; pump FWHM 127 MHz [M-SI Fig. S7]; 3 ns rectangular
#: pump at 6 MHz [M-SI p.15]; ρ = 0.95 [T17 §5.1, AE]; x_ss = 1.297 ↔ r_0 = 0.533 [M-SI Table S2, T17 §2.6].
BOREALIS: dict[str, Any] = {
    "cavity_model": "linewidth", "fwhm_hz": 66e6, "escape_efficiency": 0.95, "pump_fwhm_hz": 127e6,
    "pump_shape": "rect", "pump_duration_s": 3e-9, "x_ss": 1.297, "phi_opa0_rad": 0.0, "delta_s_rad_s": 0.0,
    "delta_p_rad_s": 0.0, "dt_kappa": 1 / 40, "window_kappa": 12.0, "segment": "opo", "rep_rate_hz": 6e6,
    "angle_table_x_ss": None, "angle_probe_d": 0.01, "toggles": TOGGLES_ON,
}

#: Liu 2026: FSR 7.4 GHz, finesse 59 (signal) / 29 (pump) [L-SI §1.1]; 1.6 ns rectangular pump at 20 MHz
#: [L-SI §1.1]; ρ = 0.99 [T17 §5.1, AE]; x_ss = 4.129 ↔ r_0 = 1.8 [L-本文, T17 §2.6].
LIU: dict[str, Any] = {
    "cavity_model": "fsr_finesse", "fsr_hz": 7.4e9, "finesse": 59.0, "escape_efficiency": 0.99,
    "pump_fwhm_hz": 7.4e9 / 29.0, "pump_shape": "rect", "pump_duration_s": 1.6e-9, "x_ss": 4.129, "phi_opa0_rad": 0.0,
    "delta_s_rad_s": 0.0, "delta_p_rad_s": 0.0, "dt_kappa": 1 / 40, "window_kappa": 12.0, "segment": "opo",
    "rep_rate_hz": 20e6, "angle_table_x_ss": None, "angle_probe_d": 0.01, "toggles": TOGGLES_ON,
}


def pulsed_opo(base: dict[str, Any], **kw: Any):
    """``pulsed_opo_source`` block from ``base`` with overrides (``toggles`` merged)."""
    tg = dict(base["toggles"], **kw.pop("toggles", {}))
    return make_pulsed_opo_source({**base, **kw, "toggles": tg})
