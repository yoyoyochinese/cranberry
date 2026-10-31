"""WP-I5: the deliberately bad configurations of 00_summary §2.5 (D34) degrade in the expected direction.

Static (offline-stage) checks run by default: stray light of a mW same-polarization probe (T10, T17), switch-tree and
all-fibre losses (T11, T12), integrating-BHD dilution −10 → −0.12 dB (T11 §3.7), fibre EOM on the squeezed path
(T12), LO shape overlaps (T17 §6.3), FSR comb without filter (T17 V-P15), long pump K ≈ 2 (T17 §5.2), seed photons
in the GBS (L-A-seed, T10), forward GAWBS of a co-propagating probe (T07). Time-loop checks (``slow``): orthogonal PM
probe without TEC (T09, T10 test 4), 1 THz probe drift (T10), Moku ±0.5 V range (T16), SMF orthogonal probe (T09
test 2), scheme B pump following (T14 §8.3).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]


def bad(name: str) -> str:
    return f"presets/bad/{name}.yaml"


@pytest.fixture(scope="module")
def built():
    cache = {}

    def get(preset, overlays=()):
        key = (preset, tuple(overlays))
        if key not in cache:
            cache[key] = fs.build_check(preset, tuple(overlays))
        return cache[key]

    return get


def test_every_bad_preset_states_what_it_demonstrates():
    for b in fs.manifest()["bad"]:
        head = (fs.ROOT / bad(b)).read_text().splitlines()[0]
        assert head.startswith("# BAD (D34"), b


def test_mw_same_polarization_probe_floods_the_squeezed_mode(built):
    base, mw = built(fs.nominal("F_liu20")), built(bad("bad_mw_probe_samepol"))
    s0, s1 = fs.stray_light(base), fs.stray_light(mw)
    # 1 µW orthogonal (ε_pol = ε_xt + T_⊥ = 2e−3) → 1 mW same polarization: ×1e3 power × 500 polarization
    assert s1["photons_per_mode"] / s0["photons_per_mode"] == pytest.approx(1e3 / 2e-3, rel=0.05)
    assert s1["photons_per_mode"] > 1e3              # T17 §6: ≈ 435 per OPO mode at 63 m; ×20 length here
    assert "sbs_gain" in fs.chk_codes(mw)            # 10 mW peak in 1.27 km (T07 §5.1)


def test_switch_tree_and_all_fibre_losses(built):
    wg6 = built(fs.nominal("F_wg6")).metrics.efficiency[1].eta_gbs
    ln = built(bad("bad_ln_switch_tree")).metrics.efficiency[1].eta_gbs
    nf = built(bad("bad_all_fiber_NF")).metrics.efficiency[1].eta_gbs
    p = fs.params(built(fs.nominal("F_wg6")))
    t_pc = p["dmx_t_stage"] * (1 - 1 / p["dmx_er_lin"])            # dB names load linearised (D56-1)
    t_ln = 0.537 * (1 - 1e-3)
    assert ln / wg6 == pytest.approx((t_ln / t_pc) ** 5, rel=1e-9)   # five tree stages (T11 §3.2)
    assert nf < ln < 0.05 < wg6                                      # T11 §7.2 "squeezing はほぼ残らない"
    assert 10 * math.log10(1 - nf) > -0.1                            # bound on squeezing ≈ −10 log(L) → ~0 dB


def test_integrating_bhd_dilution(built):
    base, integ = built(fs.nominal("F_wg100")), built(bad("bad_integrating_bhd"))
    k = fs.params(integ)["K_ports"]
    assert integ.metrics.efficiency[1].eta_bhd == pytest.approx(base.metrics.efficiency[1].eta_bhd / k, rel=1e-12)
    assert integ.metrics.efficiency[1].eta_gbs == base.metrics.efficiency[1].eta_gbs   # GBS unaffected
    # T11 §7.7: K = 32, η = 1, −10 dB → (0.1 + 31)/32 → −0.124 dB
    assert 10 * math.log10((10 ** -1.0 + (k - 1)) / k) == pytest.approx(-0.124, abs=0.001)
    assert integ.metrics.noise_levels["i.port0"]["eta"] == pytest.approx(
        base.metrics.noise_levels["i.port0"]["eta"] / k, rel=1e-12)


def test_eom_on_squeezed_path(built):
    base, eom = built(fs.nominal("F_wg100")), built(bad("bad_eom_on_squeezed"))
    ratio = eom.metrics.efficiency[1].eta_gbs / base.metrics.efficiency[1].eta_gbs
    assert ratio == pytest.approx(10 ** -0.32, rel=1e-12)            # 2.7 dB + 2 × 0.25 dB (T12 §3.2)
    # T12: a single 2.7 dB EOM with the Ref A upstream caps squeezing at ≈ 3.1 dB; the bound of this port:
    bound = -10 * math.log10(1 - eom.metrics.efficiency[1].eta_gbs)
    assert bound < 1.5 < -10 * math.log10(1 - base.metrics.efficiency[1].eta_gbs)


def test_lo_shape_overlaps(built):
    gauss = fs.block(built(fs.nominal("F_bor6")), "lo_shape").eta_t
    rect = fs.block(built(bad("bad_gauss_rect_lo")), "lo_shape").eta_t
    assert gauss == pytest.approx(0.955, abs=0.005)                  # T17 §6.3 as corrected by D60-5 / D70-2
    assert rect == pytest.approx(0.866, abs=0.005)                   # D70-2 (T17: 0.864)
    e_g = built(fs.nominal("F_bor6")).metrics.efficiency[1].eta_bhd
    e_r = built(bad("bad_gauss_rect_lo")).metrics.efficiency[1].eta_bhd
    assert e_r / e_g == pytest.approx(rect / gauss, rel=1e-9)


def test_no_fsr_filter(built):
    n_base = fs.block(built(fs.nominal("F_bor6")), "fsr").n_fsr
    blk = fs.block(built(bad("bad_no_fsr_filter")), "fsr")
    r = fs.params(built(bad("bad_no_fsr_filter")))["r_fsr_pairs"]
    assert blk.n_fsr == pytest.approx(2 * 40 * math.sinh(r) ** 2, rel=1e-9)   # T17 V-P15: T_f ≡ 1
    assert blk.n_fsr > 1e5 * n_base


def test_long_pump_is_two_mode(built):
    base = fs.block(built(fs.nominal("F_bor6")), "src")
    long = fs.block(built(bad("bad_long_pump_opo")), "src")
    assert base.schmidt_number == pytest.approx(1.10, abs=0.01)      # T17 V-P5
    assert long.schmidt_number == pytest.approx(2.1, abs=0.15)       # T17 §5.2 "K ≈ 2"
    assert float(long.supermode_r[0]) == pytest.approx(float(base.supermode_r[0]), abs=0.08)   # same r_0


def test_rep_rate_tail_mixing(built):
    tail = built(bad("bad_bor_100MHz")).offline.source_outputs["src"].tail_into_next_slot
    assert tail > 0.1 and built(fs.nominal("F_bor6")).offline.source_outputs["src"].tail_into_next_slot < 1e-6


def test_seed_in_gbs_photons(built):
    res = built(bad("bad_refA_seed_in_gbs"))
    p = fs.params(res)
    n = p["seed_photons_per_gate"]
    assert n == pytest.approx(p["p_seed_w"] * p["t_gate_s"] * p["lambda0_m"] / (fs.H_PLANCK * fs.C_LIGHT), rel=1e-12)
    # T10: 1 µW → ≈ 8×10⁴ photons per 10 ns gate
    assert n * 1e-8 / p["t_gate_s"] == pytest.approx(7.8e4, rel=0.05)
    assert n > 1e5 * fs.stray_light(built(fs.nominal("F_wg100")))["clicks_per_gate"]


def test_copropagating_probe_gawbs(built):
    res = built(bad("bad_copropagating_probe"))
    n_g = fs.gawbs_photons_per_mode(res, 1e-3)
    base = fs.stray_light(built(fs.nominal("F_wg100")))["photons_per_mode"]
    assert n_g > 100 * base                         # forward GAWBS lands in the squeezed mode (T07 §4)


# ------------------------------------------------------------------------------------------------ time loop (slow)

DRIFT_RUN = {"f_sim_hz": 1.0e5}
DRIFT_P = {"latency_s": 1.0e-5, "tec_rate_div": 1000, "tau_th_spool_s": 1.0}
QUIET = {"laser.fm_white": False, "env.mech.ports.mech": False, "env.mech.lo.mech": False, "env.mech.dist.mech": False,
         "env.mech.fs_lock.mech": False, "env.mech.fs_sq.mech": False, "dmx0.drive_phase_noise": False,
         "dmx1.drive_phase_noise": False, "dmx0.bias_drift": False, "dmx1.bias_drift": False,
         "opo.noise.length_noise": False}


def drift_run(preset, overlays=(), duration=3.0, params=None):
    return run(fs.load(preset, params={**DRIFT_P, **(params or {})}, run={**DRIFT_RUN, "duration_s": duration},
                       overlays=overlays, toggles=QUIET))


def theta_slope(res, t_from):
    return fs.slope(res.t_s, res.series["metric.mii.port1"], t_from)


@pytest.mark.slow
def test_orthogonal_pm_probe_without_tec_drifts():
    fast_tec = {"tec_tau_th_s": 1.0, "f_12_hz": 0.1, "r_T_spool_k_per_s": 3.0e-5, "dr_T_k_per_s": 0.0}   # Θ stays < π/2
    good = drift_run(fs.nominal("F_liu20"), duration=10.0, params=fast_tec)
    bad_ = drift_run(bad("bad_orth_pm_no_tec"), duration=10.0, params=fast_tec)
    p = fs.params(bad_)
    expect = p["K_T_pm_rad_per_m_K"] * p["L_k_m"] * (p["r_T_spool_k_per_s"] + 0.5 * p["dr_T_k_per_s"])
    s_bad, s_good = theta_slope(bad_, 6.0), theta_slope(good, 6.0)
    assert s_bad == pytest.approx(expect, rel=0.1) and s_bad < 0     # K_T·L·dT/dt, K_T < 0 (conv §9)
    assert abs(s_good) < 0.2 * abs(s_bad)                            # the TEC removes the spool's ΔT (G08)


@pytest.mark.slow
def test_wdm_1thz_probe_differential_drift():
    same = drift_run(fs.nominal("F_liu20"), overlays=("probe/pol_same",))
    wdm = drift_run(bad("bad_wdm_1THz_probe"))
    p = fs.params(wdm)
    expect = -2 * math.pi * 1e12 * p["kappa_tau_s_per_m_K"] * p["L_k_m"] * (p["r_T_spool_k_per_s"] + 0.5 * p["dr_T_k_per_s"])
    s_wdm = theta_slope(wdm, 1.0)
    assert abs(theta_slope(same, 1.0)) < 0.02 * abs(s_wdm)           # same pol / freq: 0 (T10 test 4)
    assert s_wdm == pytest.approx(expect * (1 - math.exp(-2.0)), rel=0.25)   # τ_th = 1 s lag over [1, 3] s


@pytest.mark.slow
def test_moku_05v_range_loses_lock():
    good = drift_run(fs.nominal("F_liu20"))
    moku = drift_run(bad("bad_moku_pid_05V"))
    assert int(good.final_state["lock.sup1.n_loss"]) == 0
    assert int(moku.final_state["lock.sup1.n_loss"]) >= 1            # 10 % of the 915B range (T16 §4.5)
    assert np.max(np.abs(moku.series["act.pzt1.dpsi"])) == pytest.approx(3.25, rel=0.01)


@pytest.mark.slow
def test_smf_orthogonal_probe_unpredictable():
    good = run(fs.load(fs.nominal("F_liu20"), run={"duration_s": 0.2}))
    smf = run(fs.load(bad("bad_smf_orth"), run={"duration_s": 0.2}))
    assert smf.metrics.metric_ii[1].stats.std > 3 * good.metrics.metric_ii[1].stats.std


@pytest.mark.slow
def test_scheme_b_actuators_follow_the_pump():
    tg = {**QUIET, "env.mech.fs_sq.mech": True, "env.thermal.F0.drift": False, "env.thermal.F1.drift": False}
    loud = {"h_ref_fs_rad2_hz": 1.0e-1}                                # free-space noise well above the DAC limit cycle
    a = run(fs.load(fs.nominal("F_liu20"), params=loud, run={"duration_s": 0.1}, toggles=tg))
    b = run(fs.load(bad("bad_liuPD_B"), params=loud, run={"duration_s": 0.1}, toggles=tg))
    half = len(a.series["act.pzt1.dpsi"]) // 2
    s_a, s_b = (np.std(r.series["act.pzt1.dpsi"][half:]) for r in (a, b))
    assert s_b > 10 * s_a                                            # pump-path noise now drives every port actuator
    assert b.metrics.metric_ii[1].stats.var == pytest.approx(a.metrics.metric_ii[1].stats.var, rel=0.3)  # ψ_k still cancels
