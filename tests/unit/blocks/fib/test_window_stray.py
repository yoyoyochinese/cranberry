"""Phase 6 (D85-4) offline stray tables ``window_stray_light`` and ``forward_ref_stray_table``: analytic limits
(zero-length fibre ⇒ no tail, tail ends after 2τ_k, window mean W/(2 T_meas), infinite extinction ⇒ no leak),
placement bookkeeping (guard slot L_k doubling, D78), CHKs, D34 contract."""

import math
import warnings

import numpy as np
import pytest

from cranberry.blocks.fib import window_stray as ws
from cranberry.blocks.fib.quantum_testing import BRILLOUIN, GAMMA, GAWBS_SMF, LAMBDA0, PROBE
from cranberry.blocks.fib.scattering import stray_photons
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.registry import make_block
from cranberry.offline.tables import TableCache

C = 299792458.0
STRAY_TG = dict(stray_rayleigh=True, stray_brillouin=True, stray_raman=True, stray_reflections=True)
REF_TG = dict(mod_leak=True, edges=True, demux_crosstalk=True, forward_scatter=True, afterpulse=True)


def win_params(**over):
    p = dict(toggles=dict(STRAY_TG), probe={k: v for k, v in PROBE.items() if k not in ("p_peak_w", "duty")},
             brillouin=dict(BRILLOUIN), rayleigh=dict(b_r_per_m=6.2e-8), L_m=1000.0, alpha_per_m=0.0,
             gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood", lambda0_m=LAMBDA0, temperature_k=295.0,
             t_gate_s=10e-9, reflections=[[1e-6, 1.0]], reflection_z_m=[10.0],
             gbs_filter=dict(shape="flat", fwhm_hz=100e9), wdm_isolation_lin=1.0,
             mode=dict(shape="gaussian", t_fwhm_s=1e-9), n_g=1.468, p_on_avg_w=1e-6, eps_off=0.0,
             period_s=100e-6, duty_lock=0.5, phase_s=0.0, edge_s=0.0, slot_period_s=10e-9)
    p.update(over)
    return p


def win(**over):
    return make_block(BlockSpec("window_stray_light", "ws", win_params(**over)))


def test_window_factor_zero_length_fibre_has_no_tail(tmp_path):
    b = win(L_m=0.0, reflections=[], reflection_z_m=[])
    v = b.tables(TableCache(tmp_path))["window_stray"].values
    assert np.all(v == 0.0)                            # factor and every count: the probe is off, nothing returns
    t = np.linspace(0, 100e-6, 1001)
    a = ws.window_factor(t, period_s=100e-6, duty_lock=0.5, phase_s=0.0, edge_s=0.0, roundtrip_s=0.0)
    np.testing.assert_array_equal(a, (t % 100e-6 < 50e-6).astype(float))


def test_window_tail_ends_after_two_tau_and_mean_is_w_over_2tmeas(tmp_path):
    L = 1000.0
    w = 2 * 1.468 * L / C                              # 9.79 µs
    b = win(L_m=L, reflections=[], reflection_z_m=[])
    tab = b.tables(TableCache(tmp_path))["window_stray"]
    t, f = tab.grid["t_s"] - 50e-6, tab.values[:, 0]   # time since the probe switched off
    np.testing.assert_allclose(f, np.clip(1 - t / w, 0, None), atol=1e-9)
    assert np.all(f[t > w] == 0.0)
    r = b.check()
    assert r["mean_factor"] == pytest.approx(w / (2 * 50e-6), rel=1e-4)
    assert r["tail_s"] == pytest.approx(w) and r["tail_time_fraction"] == pytest.approx(w / 50e-6)
    cw = stray_photons(b.params, {**STRAY_TG, "stray_reflections": False}, 1e-6)
    assert r["photons_per_mode_mean"] == pytest.approx(r["mean_factor"] * cw["photons_per_mode"], rel=1e-9)
    assert r["photons_per_mode_max"] == pytest.approx(f[0] * cw["photons_per_mode"], rel=1e-9)


def test_window_extinction_floor_edges_and_reflection(tmp_path):
    b = win(eps_off=1e-3, L_m=10.0, edge_s=1e-6)
    r = b.check()
    # long after the tail the floor is ε_off times the CW value (reflection included: its delay is ≪ the window)
    v = b.tables(TableCache(tmp_path))["window_stray"].values
    cw = stray_photons(b.params, STRAY_TG, 1e-6)
    assert v[-1, ws.WINDOW_COLUMNS.index("photons_per_mode")] == pytest.approx(1e-3 * cw["photons_per_mode"], rel=1e-9)
    assert r["t_meas_s"] == pytest.approx(49e-6)       # the ramp halves are not measurement time
    # periodic: a lock window longer than W and a CW limit (duty → 1 is not allowed; factor ≤ 1 everywhere)
    t = np.linspace(0, 300e-6, 3001)
    h = ws.window_factor(t, period_s=100e-6, duty_lock=0.9, phase_s=3e-6, edge_s=2e-6, roundtrip_s=150e-6)
    assert np.all((h >= 0) & (h <= 1 + 1e-12))
    assert h.mean() == pytest.approx(0.9, rel=1e-3)    # the box average keeps the duty


def test_window_no_slot_warns_and_limit():
    b = win(slot_period_s=60e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = b.check()
    assert r["n_slots"] == 0 and [e.code for e in b.chk_log.events] == ["window_no_slot"]
    lim = dict(photons_per_mode_max=1e-30, source="test", statistic="max")
    b = win(limit=lim)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b.check()
    assert [e.code for e in b.chk_log.events] == ["stray_light_limit"]
    b = win(limit=dict(lim, photons_per_mode_max=1e3))
    b.check()
    assert b.chk_log.events == []


def test_window_contract():
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("window_stray_light", "ws", {k: v for k, v in win_params().items() if k != "eps_off"}))
    with pytest.raises(ValueError, match="duty_lock"):
        win(duty_lock=1.0)
    with pytest.raises(ValueError, match="reflection_z_m"):
        win(reflection_z_m=[])


# --- forward_ref_stray_table ---------------------------------------------------------------------------------------

def ref_params(**over):
    p = dict(toggles=dict(REF_TG), placement="cal_frames", p_peak_w=1e-3, t_ref_s=1e-9, freq_offset_hz=0.0,
             eps_mod=1e-5, edge_tau_s=1e-10, er_dmx_lin=10 ** 2.3, qe_det=0.93, p_afterpulse=1e-3, tau_afterpulse_s=1e-7,
             det_recovery_s=43e-9, t_rep_s=50e-9, N_b=4, n_guard=0, K=32, t_cal_s=30.0, t_gate_s=10e-9,
             gbs_filter=dict(shape="flat", fwhm_hz=1e8), mode=dict(shape="gaussian", t_fwhm_s=1e-8), L_m=1270.0,
             alpha_per_m=1e-4, gawbs=dict(GAWBS_SMF), gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood",
             lambda0_m=LAMBDA0, temperature_k=295.0, sigma_target_rad=0.01, lock_ratio_R=1.0)
    p.update(over)
    return p


def ref(**over):
    return make_block(BlockSpec("forward_ref_stray_table", "fr", ref_params(**over)))


def _quiet_check(b):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return b.check()


def test_forward_ref_infinite_extinction_has_no_leak():
    r = _quiet_check(ref(eps_mod=0.0, er_dmx_lin=math.inf))
    for c in ("mod_leak_mode", "mod_leak_gate", "forward_mode", "forward_gate", "crosstalk_mode", "crosstalk_gate"):
        assert r[f"{c}_mean"] == 0.0
    # the leak is linear in ε_mod·P_pk and equals the CW line of the GBS mode: Φ_leak·τ_eff
    r1 = _quiet_check(ref(edge_tau_s=0.0, p_afterpulse=0.0))
    hnu = 6.62607015e-34 * C / LAMBDA0
    tau_eff = 1e-8 * math.sqrt(math.pi / math.log(2))
    assert r1["mod_leak_mode_mean"] == pytest.approx(1e-5 * 1e-3 / hnu * tau_eff, rel=1e-9)
    assert r1["mod_leak_gate_mean"] == pytest.approx(1e-5 * 1e-3 / hnu * 10e-9, rel=1e-12)
    r2 = _quiet_check(ref(edge_tau_s=0.0, p_afterpulse=0.0, eps_mod=2e-5))
    assert r2["forward_mode_mean"] == pytest.approx(2 * r1["forward_mode_mean"], rel=1e-9)


def test_forward_ref_edges_afterpulse_and_toggles_off():
    off = _quiet_check(ref(toggles=dict.fromkeys(REF_TG, False)))
    assert off["photons_per_mode_mean"] == 0.0 and off["clicks_per_gate_mean"] == 0.0
    # guard slot, N_b = 1: own pulse one slot (50 ns) before the data slot, exponential tail of 5 ns
    p = dict(placement="guard_slot", n_guard=1, N_b=1, eps_mod=0.0, er_dmx_lin=math.inf, edge_tau_s=5e-9)
    r = _quiet_check(ref(**p))
    hnu = 6.62607015e-34 * C / LAMBDA0
    x = 50e-9 - 0.5e-9
    phi = 1e-3 / hnu
    assert r["edges_gate_mean"] == pytest.approx(phi * 5e-9 * (math.exp(-(x - 5e-9) / 5e-9) - math.exp(-(x + 5e-9) / 5e-9)),
                                                 rel=1e-9)
    p_click = -math.expm1(-0.93 * phi * 1e-9)
    assert r["afterpulse_gate_mean"] == pytest.approx(p_click * 1e-3 * (math.exp(-45e-9 / 1e-7) - math.exp(-55e-9 / 1e-7)),
                                                      rel=1e-9)
    assert r["blind_prob_mean"] == 0.0                 # 50 ns > 43 ns recovery
    assert _quiet_check(ref(**p, det_recovery_s=1e-5))["blind_prob_mean"] == pytest.approx(p_click)


def test_forward_ref_placements_and_d78():
    b = ref(placement="guard_slot", n_guard=1, N_b=1, t_rep_s=1 / 6e6, edge_tau_s=20e-9)
    r = _quiet_check(b)
    assert r["lk_factor"] == pytest.approx(2.0)        # N_b = 1: the guard slot doubles L_k (D62-3)
    assert r["dnu_max_hz"] == pytest.approx(r["dnu_max_without_ref_hz"] / 2)
    assert r["dtau_s"] == pytest.approx(31 * 2 / 6e6)
    # D78 table: F_wg6 / F_bor6 (Δτ_31 = 5.17 µs) at 10 mrad, R = 1: 3.08 Hz
    assert r["dnu_max_without_ref_hz"] == pytest.approx(3.08, rel=2e-3)
    assert [e.code for e in b.chk_log.events] == ["forward_ref_guard_slot"]
    c = ref(t_rep_s=1 / 6e6, N_b=1, edge_tau_s=20e-9)
    rc = _quiet_check(c)
    assert rc["lk_factor"] == 1.0 and c.chk_log.events == []
    # calibration frames: only the frames next to a calibration frame, each with weight T_f/t_cal
    assert 0.0 < r["edges_mode_mean"] and rc["edges_mode_mean"] < r["edges_mode_mean"] * 1e-6
    with pytest.raises(ValueError, match="n_guard"):
        ref(placement="guard_slot", n_guard=0)
    with pytest.raises(ValueError, match="placement"):
        ref(placement="intra_slot")


def test_forward_ref_table_and_limit(tmp_path):
    b = ref(limit=dict(photons_per_mode_max=1e-30, source="test"))
    tab = b.tables(TableCache(tmp_path))["forward_ref_stray"]
    assert tab.values.shape == (4, len(ws.FORWARD_REF_COLUMNS))
    _quiet_check(b)
    assert [e.code for e in b.chk_log.events] == ["stray_light_limit"]
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("forward_ref_stray_table", "fr", {k: v for k, v in ref_params().items() if k != "eps_mod"}))
