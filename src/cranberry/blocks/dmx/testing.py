"""Test specs for the WP-B8 blocks (T11 §7, §9 numbers; T10 §3.4.5, §7). Not presets: presets give every
value with its source tag (D34, PLAN §7.3)."""

from __future__ import annotations

import math
from typing import Any

from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

__all__ = ["N_G", "LAMBDA0", "design", "switch", "delay", "pulser", "stray_params", "make", "compiled_matches_reference"]

N_G = 1.468  # T11 §7 / T10 (survey value, UNVERIFIED, T06)
LAMBDA0 = 1545.32e-9  # D07 placeholder (Ref A)


def design(K=32, layout="tree", switching="active", toggles=None, **over) -> dict[str, Any]:
    """``demux_design``: T11 test 1 tree (t = 0.99, ER = 23.01 dB → ε = 0.005), 100 MHz, N_b = 1."""
    p = {"K": K, "layout": layout, "switching": switching, "f_rep_hz": 100e6, "N_b": 1, "n_guard": 0, "L0_m": 0.0,
         "n_g": N_G, "alpha_db_per_m": 0.2e-3, "t_stage": 0.99, "stage_gdd_s2": 0.0,
         "segments": [f"D{k}" for k in range(K)], "loss_applies_to": ["gbs", "bhd"], "loss_source": "[T11 §9 test 1]",
         "er_db": -10.0 * math.log10(0.005), "t_rise_s": 10e-9, "sigma_j_s": 0.0, "t_ring_s": 0.0,
         "tau_mode_s": 1e-12, "f_tog_max_hz": 1e6,
         "toggles": {"enabled": True, "insertion_loss": True, "crosstalk": True, "rise_time_clipping": False}}
    p.update(over)
    p["toggles"] = {**p["toggles"], **(toggles or {})}
    return {"block": "demux_design", "name": "dmx", "params": p}


def switch(name="sw", n_stages=5, mode="timeseries", toggles=None, rate_div=1, **over) -> dict[str, Any]:
    """``demux_switch`` on segment ``dmx`` (driver ``eom``: a Pockels / LN switch is an EO phase): T11 test 5 numbers (κ_pm = 1, δV/V_π = 1.4235e−3)."""
    p = {"segment": "dmx", "driver": "eom", "n_sub": 1, "sub": [0, 1], "ring_len": 0, "n_stages": n_stages,
         "drive_noise_mode": mode, "kappa_pm": 1.0, "v_pi_v": 1000.0, "drive_noise_v_rms": 1.4235,
         "bias_offset_rad": 0.0, "bias_rw_rad_per_rts": 0.0,
         "toggles": {"enabled": True, "drive_phase_noise": True, "bias_drift": False}}
    p.update(over)
    p["toggles"] = {**p["toggles"], **(toggles or {})}
    out = {"block": "demux_switch", "name": name, "params": p}
    if rate_div != 1:
        out["rate_div"] = rate_div
    return out


def delay(name="dl", port=31, toggles=None, **over) -> dict[str, Any]:
    """``delay_line_k`` (T11 test 7: 100 MHz, n_g = 1.468, N_b = 1, L_0 = 0)."""
    p = {"port": port, "f_rep_hz": 100e6, "N_b": 1, "n_guard": 0, "L0_m": 0.0, "n_g": N_G,
         "kappa_t_s_per_m_k": 39e-15, "dL_m": 0.0, "beta2_s2_per_m": -21.7e-27, "alpha_db_per_m": 0.2e-3,
         "in_dT": None, "segment": f"F{port}", "loss_applies_to": ["gbs", "bhd"], "loss_source": "[T11 §7 placeholder]",
         "dtau_s": 0.0, "rep_rel_drift": 0.0,
         "toggles": {"static_mismatch": True, "delay_thermal": True, "rep_rate_drift": True, "fiber_loss": True}}
    p.update(over)
    p["toggles"] = {**p["toggles"], **(toggles or {})}
    return {"block": "delay_line_k", "name": name, "params": p}


def pulser(name="pp", toggles=None, **over) -> dict[str, Any]:
    """``probe_pulser`` (T10 test 7: T_f = 320 ns, τ_pr = 10 ns, ε_AOM = 0)."""
    p = {"period_s": 320e-9, "p_avg_w": 1e-3 * 10 / 320, "duty": 10 / 320, "t_offset_s": 5e-9, "eps_aom": 0.0,
         "toggles": {"gating": True, "aom_leak": True}}
    p.update(over)
    p["toggles"] = {**p["toggles"], **(toggles or {})}
    return {"block": "probe_pulser", "name": name, "params": p}


def stray_params(L_m: float, **over) -> dict[str, Any]:
    """B10 ``stray_photons`` parameter set of one port (T10 §3.4.5: same polarization and frequency,
    b_R = 6.18e−8 /m, α = 0, T_gate 1 ns, flat GBS filter)."""
    p = {"probe": {"freq_offset_hz": 0.0, "linewidth_hz": 1e3, "polarization": "same", "eps_xt": 0.0, "t_perp": 0.0},
         "brillouin": {"g_b_m_per_w": 5e-11, "a_eff_m2": 80e-12, "nu_b_hz": 10.8e9, "dnu_b_hz": 30e6},
         "rayleigh": {"b_r_per_m": 6.18e-8}, "L_m": L_m, "alpha_per_m": 0.0, "gamma_per_w_m": 1.3e-3, "f_R": 0.18,
         "raman_model": "blow_wood", "lambda0_m": LAMBDA0, "temperature_k": 295.0, "t_gate_s": 1e-9,
         "reflections": [], "gbs_filter": {"shape": "flat", "fwhm_hz": 1e9}, "wdm_isolation_lin": 1.0,
         "mode": {"shape": "gaussian", "t_fwhm_s": 1e-12}}
    p.update(over)
    return p


def make(obj, fs: float, seed: int = 3):
    """Build ``obj``; returns (block, State, Params)."""
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def compiled_matches_reference(obj, fs: float, n: int, seed: int = 3) -> float:
    """max |S_compiled − S_reference| after ``n`` samples from the same seed (VT-44 rule, 1e−12)."""
    b1, s1, p1 = make(obj, fs, seed)
    b2, s2, p2 = make(obj, fs, seed)
    run_reference(b1, s1, p1, n)
    run_compiled(b2, s2, p2, n)
    return float(abs(s1.S - s2.S).max()) if s1.S.size else 0.0
