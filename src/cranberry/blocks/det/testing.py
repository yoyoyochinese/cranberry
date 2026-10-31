"""Spec builders for the WP-B5 tests (DET blocks and the B5 MEAS blocks). Not physics; registers no kinds.

Every value is still given explicitly (D34); these helpers only shorten tests. Keyword overrides
replace keys; ``toggles`` overrides are merged.
"""

from __future__ import annotations

from typing import Any

__all__ = ["det_eff", "frontend", "pulsed_kernel", "zero_span", "swept", "pulse_sample", "overlap", "ALL_ON"]

ALL_ON = True


def _spec(kind: str, name: str, params: dict[str, Any], toggles: dict[str, bool], over: dict[str, Any]) -> dict:
    tg = {**toggles, **over.pop("toggles", {})}
    return {"block": kind, "name": name, "params": {**params, **over, "toggles": tg}}


def det_eff(name: str = "det", on: bool = ALL_ON, **over: Any) -> dict:
    """``detector_efficiency`` with the generic high-QE values written out (T13 §4.2, Ref A)."""
    p = dict(device=None, observer="bhd", applies_to=["bhd"], source="[T13 §4.2]", f_s_hz=3.0e6, p_lo_w=0.0154,
             qe=0.99, clearance_lin=10 ** 2.8, p_lo_ref_w=0.0154, imbalance_abs2=2.5e-4, n_ex=60.0,
             tia={"model": "first_order", "f3db_hz": 35.0e6}, dip_r=0.9, dip_tau_s=1.0 / (2 * 0.96e9))
    tg = {k: on for k in ("qe", "circuit_noise", "imbalance_leak", "freq_response", "dip")} | {"subtract_dark": False}
    return _spec("detector_efficiency", name, p, tg, over)


def frontend(name: str = "fe", on: bool = ALL_ON, **over: Any) -> dict:
    """``bhd_frontend`` with PDB450C-like values (T13 §4.1) and a LO of 15.4 mW."""
    p = dict(device=None, p_lo_w=0.0154, responsivity_a_per_w=1.0, gain_v_per_a=1.0e3, load="high_z",
             load_divisor_50ohm=2.0, coupler_sign=1, v_off_v=0.0, v_sat_v=4.6, clearance_lin=10 ** 2.8,
             p_lo_ref_w=0.0154, imbalance_abs2=2.5e-4, n_ex=60.0, tia={"model": "first_order", "f3db_hz": 150.0e6},
             dip_r=0.9, dip_tau_s=1.0 / (2 * 0.96e9), monitor_v_per_w=1.0e4, monitor_v_max_v=10.0)
    tg = {k: on for k in ("circuit_noise", "imbalance_leak", "freq_response", "dip", "saturation")}
    tg["subtract_dark"] = False
    return _spec("bhd_frontend", name, p, tg, over)


def pulsed_kernel(name: str = "pk", on: bool = ALL_ON, **over: Any) -> dict:
    p = dict(f_rep_hz=100.0e6, response={"model": "first_order", "bandwidth_hz": 50.0e6}, window_s=1.0e-9,
             n_taps=8, sigma_el_e=730.0, n_lo_photons=1.3e7, qe=1.0)
    return _spec("pulsed_bhd_kernel", name, p, {"pulsed_crosstalk": on, "circuit_noise": on}, over)


def zero_span(name: str = "esa", noise: bool = True, **over: Any) -> dict:
    """Ref A zero span: f_s 3 MHz, RBW 1 MHz, VBW 100 Hz, 0.15 s (T04 §5.1)."""
    p = dict(f_s_hz=3.0e6, rbw_hz=1.0e6, vbw_hz=100.0, enbw_factor=1.0, rbw_shape="rect", average="power",
             detector="sample", oversample=2, sweep_time_s=0.15, n_points=1000, n_traces=1)
    return _spec("esa_zero_span", name, p, {"radiometer_noise": noise}, over)


def swept(name: str = "sa", noise: bool = False, subtract: bool = True, **over: Any) -> dict:
    p = dict(f_start_hz=1.0e6, f_stop_hz=1.2e9, n_points=400, rbw_hz=1.0e6, vbw_hz=100.0, enbw_factor=1.0,
             rbw_shape="rect", average="power", n_sub=61)
    return _spec("esa_swept_spectrum", name, p, {"radiometer_noise": noise, "subtract_circuit_noise": subtract}, over)


def pulse_sample(name: str = "ps", mode: str = "pulse_resolved", on: bool = ALL_ON, **over: Any) -> dict:
    p = dict(mode=mode, k_slots=32, target_slots=[0], n_frames=4000)
    return _spec("pulse_sample", name, p, {"crosstalk": on, "electronic_noise": on, "sampling_noise": on}, over)


def overlap(name: str = "mm", on: bool = True, **over: Any) -> dict:
    """Spatial Gaussian overlap w₂/w₁ = 1.1 at the BHD observer (T12 V7)."""
    p = dict(overlap="spatial", model="gaussian_beam", w1_m=5.25e-6, w2_m=1.1 * 5.25e-6, d_m=0.0, theta_rad=0.0,
             wavelength_m=1.55e-6, n=1.0, segment=None, observer="bhd", applies_to=["bhd"], source="[T12 V7]")
    return _spec("mode_overlap", name, p, {"enabled": on}, over)
