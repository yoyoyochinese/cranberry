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
             eps_mod=1e-5, eps_gate=1.0, edge_tau_s=1e-10, er_dmx_lin=10 ** 2.3, qe_det=0.93, p_afterpulse=1e-3,
             tau_afterpulse_s=1e-7, det_recovery_s=43e-9, t_rep_s=50e-9, N_b=4, n_guard=0, K=32, t_cal_s=30.0,
             t_gate_s=10e-9, gbs_filter=dict(shape="flat", fwhm_hz=1e8), mode=dict(shape="gaussian", t_fwhm_s=1e-8),
             L_m=1270.0, gawbs=dict(GAWBS_SMF), gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood",
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


# --- Phase 5 review WS-02 … WS-09 ----------------------------------------------------------------------------------

def _tail(b, tmp_path, col="factor"):
    tab = b.tables(TableCache(tmp_path))["window_stray"]
    return tab.grid["t_s"] - 50e-6, tab.values[:, ws.WINDOW_COLUMNS.index(col)]


def test_window_tail_attenuation_weight(tmp_path):
    """WS-03: the tail is ∫_τ^W e^{−κs} ds / ∫_0^W e^{−κs} ds (κ = αc/n_g), τ = time since switch-off."""
    L, alpha = 3000.0, 0.5e-3 * math.log(10) / 10    # 0.5 dB/km, 2αL = 0.69
    w, k = 2 * 1.468 * L / C, alpha * C / 1.468
    t, f = _tail(win(L_m=L, alpha_per_m=alpha, reflections=[], reflection_z_m=[]), tmp_path)
    np.testing.assert_allclose(f, np.where(t < w, (np.exp(-k * t) - np.exp(-k * w)) / (1 - np.exp(-k * w)), 0.0),
                               atol=1e-12)
    # with ramps: brute-force quadrature of the envelope (the W = 0 path) with the weight
    timing = dict(period_s=20e-6, duty_lock=0.4, phase_s=1e-6, edge_s=2e-6)
    tq = np.array([0.5e-6, 2.3e-6, 8.9e-6, 9.6e-6, 10.4e-6, 14e-6, 19e-6, 27e-6])
    s = np.linspace(0.0, w, 200001)
    wt = np.exp(-k * s)
    brute = [np.trapezoid(ws.window_factor(x - s, roundtrip_s=0.0, **timing) * wt, s) / np.trapezoid(wt, s) for x in tq]
    np.testing.assert_allclose(ws.window_factor(tq, roundtrip_s=w, kappa_per_s=k, **timing), brute, atol=2e-6)
    # κW → 0 (series branch) reproduces the plain box
    np.testing.assert_allclose(ws.window_factor(tq, roundtrip_s=w, kappa_per_s=1e-3, **timing),
                               ws.window_factor(tq, roundtrip_s=w, **timing), atol=1e-8)


def test_reflection_at_far_end_is_delayed_by_w(tmp_path):
    """WS-07 (b) / WS-08: a reflection at z = L returns the full on level until W after switch-off, then nothing."""
    L = 1000.0
    w = 2 * 1.468 * L / C
    b = win(L_m=L, reflections=[[1e-4, 1.0]], reflection_z_m=[L])
    t, v = _tail(b, tmp_path, "reflections_mode")
    only = dict.fromkeys(("stray_rayleigh", "stray_brillouin", "stray_raman"), False) | {"stray_reflections": True}
    one = stray_photons(b.params, only, 1e-6)["reflections_mode"]
    assert one > 0
    np.testing.assert_allclose(v, np.where(t < w, one, 0.0), rtol=1e-12)


def test_window_tail_longer_than_period_gives_the_duty():
    """WS-07 (g): W = 3T (F_wg100-like, T = 10 µs): the box spans whole periods, the factor is the duty everywhere."""
    T = 10e-6
    b = win(L_m=3 * T * C / (2 * 1.468), period_s=T, duty_lock=0.65, edge_s=1e-7, slot_period_s=0.5e-6,
            reflections=[], reflection_z_m=[])
    r = b.check()
    assert r["n_slots"] > 1 and r["mean_factor"] == pytest.approx(0.65, rel=1e-9)
    assert r["photons_per_mode_max"] == pytest.approx(r["photons_per_mode_mean"], rel=1e-9)


HNU = 6.62607015e-34 * C / LAMBDA0


def test_forward_ref_eps_gate_cal_frames_only():
    """WS-02: a slow gate closes the reference path in data frames (calibration frames only)."""
    p0 = dict(edge_tau_s=0.0, p_afterpulse=0.0)
    r1, r2 = _quiet_check(ref(**p0)), _quiet_check(ref(**p0, eps_gate=1e-4))
    for c in ("mod_leak_mode", "mod_leak_gate", "forward_mode", "forward_gate"):
        assert r2[f"{c}_mean"] == pytest.approx(1e-4 * r1[f"{c}_mean"], rel=1e-9)
    with pytest.raises(ValueError, match="eps_gate"):
        ref(placement="guard_slot", n_guard=1, eps_gate=1e-4)
    ref(placement="guard_slot", n_guard=1, eps_gate=1.0)


def test_forward_ref_cal_frame_positions_and_weights():
    """WS-07 (c): own pulse of the calibration frame at Δ = T_f − (N_b − 1 − j)·T_rep, weight T_f/t_cal."""
    tr, nb, K, tcal, tau = 50e-9, 4, 32, 30.0, 3e-6
    tf = K * nb * tr
    tg = dict(REF_TG, edges=False, demux_crosstalk=False)
    b = ref(toggles=tg, eps_mod=0.0, tau_afterpulse_s=tau, det_recovery_s=1e-5, p_peak_w=1e-6)
    rows = ws.forward_ref_stray(b.params, b._toggles)
    p_click = -math.expm1(-0.93 * 1e-6 / HNU * 1e-9)
    for j in range(nb):
        d = tf - (nb - 1 - j) * tr
        ap = tf / tcal * p_click * 1e-3 * (math.exp(-(d - 5e-9) / tau) - math.exp(-(d + 5e-9) / tau))
        assert rows[j, ws.FORWARD_REF_COLUMNS.index("afterpulse_gate")] == pytest.approx(ap, rel=1e-9)
        assert rows[j, ws.FORWARD_REF_COLUMNS.index("blind_prob")] == pytest.approx(tf / tcal * p_click, rel=1e-9)
    assert _quiet_check(b)["data_slot_fraction"] == pytest.approx(1 - tf / tcal)


def _edge_gate(x, tau_e):
    return tau_e * (math.exp(-max(x - 5e-9, 0) / tau_e) - math.exp(-max(x + 5e-9, 0) / tau_e))


def test_forward_ref_crosstalk_positions():
    """WS-06 / WS-07 (d): leaked pulses of the other ports (edges only, ε_dmx = 1/ER)."""
    tr, tau_e, eps = 50e-9, 5e-9, 1e-2
    phi = 1e-3 / HNU
    tg = dict.fromkeys(REF_TG, False) | {"edges": True, "demux_crosstalk": True}
    col = ws.FORWARD_REF_COLUMNS.index("crosstalk_gate")
    # guard slot with n_guard = 2 (B = 6): the next port's guard pulse is B − 1 − j slots after slot j
    g = ref(toggles=tg, placement="guard_slot", n_guard=2, er_dmx_lin=1 / eps, edge_tau_s=tau_e)
    for j, v in enumerate(ws.forward_ref_stray(g.params, g._toggles)[:, col]):
        assert v == pytest.approx(eps * phi * _edge_gate((5 - j) * tr - 0.5e-9, tau_e), rel=1e-6)
    # calibration frames (worst port): Δ = (1 + n_guard + j)·T_rep and −(n_guard + N_b − j)·T_rep, weight T_f/t_cal
    c = ref(toggles=tg, er_dmx_lin=1 / eps, edge_tau_s=tau_e)
    wf = 32 * 4 * tr / 30.0
    for j, v in enumerate(ws.forward_ref_stray(c.params, c._toggles)[:, col]):
        x1, x2 = (1 + j) * tr - 0.5e-9, (4 - j) * tr - 0.5e-9
        assert v == pytest.approx(wf * eps * phi * (_edge_gate(x1, tau_e) + _edge_gate(x2, tau_e)), rel=1e-9)


def test_forward_ref_guard_slot_counts_every_leaked_pulse():
    """WS-05: in the guard slot every other port's leaked guard pulse clicks and afterpulses (F_bor6-like timing)."""
    tr, K, B, tau, er = 1 / 6e6, 32, 2, 1e-6, 10 ** 2.3
    tg = dict.fromkeys(REF_TG, False) | {"afterpulse": True, "demux_crosstalk": True}
    common = dict(toggles=tg, placement="guard_slot", N_b=1, n_guard=1, t_rep_s=tr, er_dmx_lin=er, p_peak_w=1e-6,
                  eps_mod=0.0)
    r = _quiet_check(ref(**common, tau_afterpulse_s=tau))
    n1 = 0.93 * 1e-6 / HNU * 1e-9

    def ap(d, n):
        return -math.expm1(-n) * 1e-3 * (math.exp(-(d - 5e-9) / tau) - math.exp(-(d + 5e-9) / tau))

    own = ap(tr, n1)
    expect = own + sum(ap((m * B + 1) * tr, n1 / er) for m in range(1, K)) + ap((K * B + 1) * tr, n1)
    assert r["afterpulse_gate_mean"] == pytest.approx(expect, rel=1e-9)
    assert r["afterpulse_gate_mean"] > 3 * own          # Phase 5 review: ×3.5 at τ_ap = 1 µs
    rb = _quiet_check(ref(**common, det_recovery_s=2e-5))
    assert 0.999 < rb["blind_prob_mean"] <= 1.0          # a probability, not a sum of click probabilities


def test_forward_scatter_magnitude():
    """WS-07 (f) / WS-09: forward GAWBS = gawbs_mode_photons, forward Raman = P_out·L·n_R (not L_eff)."""
    from cranberry.blocks.fib.scattering import gawbs_mode_photons, mode_average, raman_scatter_density

    p0 = dict(edge_tau_s=0.0, p_afterpulse=0.0)
    full = _quiet_check(ref(**p0))
    nog = _quiet_check(ref(**p0, gawbs=dict(GAWBS_SMF, gamma_om_per_w_m=[0.0] * 12)))
    P, mode, L = 1e-5 * 1e-3, dict(shape="gaussian", t_fwhm_s=1e-8), 1270.0
    g = gawbs_mode_photons(P, GAWBS_SMF, mode, temperature_k=295.0, lambda0_m=LAMBDA0, length_m=L)
    assert full["forward_mode_mean"] - nog["forward_mode_mean"] == pytest.approx(g, rel=1e-9)
    raman = mode_average(lambda f: P * L * raman_scatter_density(-2 * np.pi * np.asarray(f), gamma_per_w_m=GAMMA,
                                                                 f_R=0.18, raman_model="blow_wood", temperature_k=295.0),
                         mode)
    assert nog["forward_mode_mean"] == pytest.approx(raman, rel=1e-9) and raman > 0
