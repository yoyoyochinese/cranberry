"""Test support for WP-B2 blocks: pipeline-element builders with every parameter given (test values only;
presets give their own values with provenance, D34). Registers nothing.

Values: DX1 SC (T05 §2.2, D35, D49, D51-9): h0 = 8e3 Hz²/Hz, h₋₁ = 3e8 Hz², RIN −140 / −150 dB/Hz;
T08 §2.5–§2.6 defaults (h_ref = 0.04 rad²·Hz at 1 km, f_c = 1 kHz, f_lo = 0.1 Hz; Dong 2016 Table I).
"""

from __future__ import annotations

from typing import Any

import cranberry.path.testing as pt

__all__ = ["ALL_LASER", "laser", "thermal", "crystal", "pd_offset", "mech", "thermo", "gawbs", "DONG", "K_EPS",
           "seg_topo"]

ALL_LASER = ("fm_white", "fm_flicker", "peaks", "drift", "dither", "sbs", "rin")
#: k_ε = ξ_eo·β = 0.78 · 2π·1.468/1548.8 nm [rad/m per unit strain] (conv §4, T16)
K_EPS = 0.78 * 2 * 3.141592653589793 * 1.468 / 1548.8e-9
#: Dong 2016 Table I (T08 §2.6)
DONG = {"T_k": 300.0, "wavelength_m": 1550e-9, "dn_dT_per_k": 9.2e-6, "n": 1.468, "alpha_l_per_k": 5.5e-7,
        "kappa_th_w_per_m_k": 1.37, "d_th_m2_per_s": 8.2e-7, "k_max_per_m": 3.846e5, "k_min_per_m": 3.848e4,
        "e0_pa": 1.9e10, "area_m2": 4.91e-8, "phi0": 0.01}


def _blk(kind: str, name: str, params: dict[str, Any], rate_div: int = 1) -> dict[str, Any]:
    obj: dict[str, Any] = {"block": kind, "name": name, "params": params}
    if rate_div != 1:
        obj["rate_div"] = rate_div
    return obj


def laser(name: str = "las", source: str = "dx1", on=ALL_LASER, rate_div: int = 1, **over: Any) -> dict[str, Any]:
    p = {"source": source, "power_w": 3.33e-3, "oob_pairs": [], "rin_targets": [],
         "h0_hz2_per_hz": 8e3, "h_m1_hz2": 3e8, "flicker_f_min_hz": 0.1, "flicker_per_decade": 2,
         "peaks": [{"f_hz": 89e3, "q": 10.0, "area": 1e6, "quantity": "freq"},
                   {"f_hz": 177e3, "q": 10.0, "area": 1e-9, "quantity": "rin"}],
         "nu_dot_hz_per_s": 1e5, "dither_amp_hz": 1e6, "dither_f_hz": 1e3, "sbs_amp_hz": 1e8, "sbs_f_hz": 3e3,
         "rin_psd": [[2e7, 1e-14], [1e9, 1e-15]], "toggles": {t: t in on for t in ALL_LASER}}
    p.update(over)
    return _blk("laser_source", name, p, rate_div)


def seg_topo(segs: dict[str, tuple[float, int, str]], *, coef: dict[str, float] | None = None, beams=None):
    """Topology with segments {id: (L_m, n_sub, placement)} carrying the drivers of ``coef`` (same placement),
    one beam per segment ``b_<id>`` (unless ``beams``), coefficients {driver: value per metre} × L."""
    coef = coef or {"dT": 1.0, "strain": 1.0}
    seg_list, beam_list, obs, cf = [], [], [], []
    for s, (L, n, place) in segs.items():
        spatial = "lumped" if place == "uniform" else "distributed"
        seg_list.append(pt.seg_spec(s, L, n_sub=n, spatial=spatial, drivers=dict.fromkeys(coef, place)))
    for s in segs:
        if beams is None:
            beam_list.append(pt.beam_spec(f"b_{s}", [pt.step_spec(s, 1)], f"o_{s}"))
            obs.append(f"o_{s}")
    if beams is not None:
        beam_list, obs = beams
    for s, (L, _, _) in segs.items():
        for drv, k in coef.items():
            cf += pt.coef_spec(s, drv, "*", k * L)
    return pt.topo_spec(seg_list, beam_list, obs, cf)


def thermal(topo, segments, tau: dict[str, float], name: str = "th", on=("drift", "random", "periodic"),
            rate_div: int = 1, **over: Any) -> dict[str, Any]:
    p = {"topology": topo, "segments": segments, "tau_th_s": tau, "r_T_k_per_s": 1e-4, "b_T_k2_hz": 1e-8,
         "periodic_amp_k": 0.05, "periodic_period_s": 900.0,
         "toggles": {t: t in on for t in ("drift", "random", "periodic")}}
    p.update(over)
    return _blk("segment_thermal", name, p, rate_div)


def crystal(name: str = "xt", segment: str = "xtal", on=("drift", "random", "periodic"), **over: Any) -> dict[str, Any]:
    p = {"segment": segment, "tau_th_s": 10.0, "crystal_length_m": 0.045, "dT_fwhm_k": 2.7, "detune_k": 0.01,
         "r_T_k_per_s": 1e-5, "b_T_k2_hz": 1e-10, "periodic_amp_k": 1e-3, "periodic_period_s": 600.0,
         "toggles": {t: t in on for t in ("drift", "random", "periodic")}}
    p.update(over)
    return _blk("crystal_thermal", name, p)


def pd_offset(name: str = "pdo", out: str = "v.voff", on=("offset", "drift", "random"), **over: Any) -> dict[str, Any]:
    p = {"out": out, "offset_v": 5e-3, "r_v_per_s": 1e-6, "b_v2_hz": 1e-12,
         "toggles": {t: t in on for t in ("offset", "drift", "random")}}
    p.update(over)
    return _blk("pd_offset_drift", name, p)


def mech(topo, segments, name: str = "mn", on=("mech", "peaks", "correlation"), rate_div: int = 1,
         **over: Any) -> dict[str, Any]:
    p = {"topology": topo, "segments": segments, "h_ref_rad2_hz": 0.04, "L_ref_m": 1000.0, "gamma": 1.0,
         "f_c_hz": 1e3, "f_lo_hz": 0.1, "k_eps_rad_per_m": K_EPS,
         "peaks": [{"f_hz": 40.0, "q": 5.0, "area_rad2": 1e-6}], "rho_0": 0.9, "f_rho_hz": 10.0,
         "spatial_coherence": {}, "toggles": {t: t in on for t in ("mech", "peaks", "correlation")}}
    p.update(over)
    return _blk("segment_mech_noise", name, p, rate_div)


def thermo(topo, segments, name: str = "td", driver: str = "strain", accumulate: bool = True,
           on=("conduction", "thermomechanical"), **over: Any) -> dict[str, Any]:
    p = {"topology": topo, "segments": segments, "driver": driver, "accumulate": accumulate,
         "k_rad_per_m_unit": K_EPS, "material": dict(DONG), "tm_f_min_hz": 0.1, "tm_per_decade": 2,
         "toggles": {t: t in on for t in ("conduction", "thermomechanical")}}
    p.update(over)
    return _blk("segment_thermodynamic", name, p)


def gawbs(name: str = "gw", targets=("sqz",), sigma2: float = 1e-4, on: bool = True) -> dict[str, Any]:
    return _blk("gawbs_phase_jitter", name, {"targets": list(targets), "sigma2_rad2": sigma2, "toggles": {"gawbs": on}})
