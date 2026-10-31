"""Unit tests of the WP-B10 scattering tables and ``sbs_check``: contract (D34), toggles off ⇒ zero,
linearity, mode weights, polarization factors, cache."""

import warnings

import numpy as np
import pytest

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.quantum_testing import BRILLOUIN, GAMMA, GAWBS_SMF, LAMBDA0, PROBE
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import make_block
from cranberry.offline.tables import TableCache

STRAY_TG = dict(stray_rayleigh=True, stray_brillouin=True, stray_raman=True, stray_reflections=True)


def stray_params(**over):
    p = dict(toggles=dict(STRAY_TG), p_peak_w=[1e-6, 1e-3, 2e-3], duty=1.0, probe={k: v for k, v in PROBE.items()
                                                                                   if k not in ("p_peak_w", "duty")},
             brillouin=dict(BRILLOUIN), rayleigh=dict(b_r_per_m=6.2e-8), L_m=63.0, alpha_per_m=0.0,
             gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood", lambda0_m=LAMBDA0, temperature_k=295.0,
             t_gate_s=10e-9, reflections=[[1e-6, 1.0]], gbs_filter=dict(shape="flat", fwhm_hz=100e9),
             wdm_isolation_lin=1.0, mode=dict(shape="gaussian", t_fwhm_s=1e-12))
    p.update(over)
    return p


def test_mode_weight_normalised_and_tau_eff():
    for mode in (dict(shape="gaussian", t_fwhm_s=1e-12), dict(shape="exponential", tau_eff_s=14.5e-9)):
        assert sc.mode_average(lambda f: np.ones_like(f), mode) == pytest.approx(1.0, rel=1e-4)
        assert sc.mode_weight(0.0, mode) == pytest.approx(sc.mode_tau_eff(mode), rel=1e-15)
    assert sc.mode_tau_eff(dict(shape="gaussian", t_fwhm_s=1e-12)) == pytest.approx(2.129e-12, rel=1e-3)
    with pytest.raises(ValueError, match="mode.shape"):
        make_block(BlockSpec("stray_light_table", "st", stray_params(mode=dict(shape="sech"))))


def test_eps_pol():
    assert sc.eps_pol("same", 0.1, 0.2) == 1.0
    assert sc.eps_pol("pm_orthogonal", 1e-3, 1e-3) == pytest.approx(2e-3)
    assert sc.eps_pol("smf_orthogonal", 0.0, 0.0) == pytest.approx(1 / 3)
    with pytest.raises(ValueError):
        sc.eps_pol("x", 0, 0)


def test_stray_table_linear_toggles_and_cache(tmp_path):
    cache = TableCache(tmp_path)
    t = make_block(BlockSpec("stray_light_table", "st", stray_params())).tables(cache)["stray_light"]
    col = {c: i for i, c in enumerate(sc.STRAY_COLUMNS)}
    v = t.values
    np.testing.assert_allclose(v[2], 2 * v[1], rtol=1e-12)  # ∝ P̄
    assert v[1, col["clicks_per_gate"]] == pytest.approx(sum(v[1, col[f"{m}_gate"]] for m in
                                                            ("rayleigh", "brillouin", "raman", "reflections")))
    t2 = make_block(BlockSpec("stray_light_table", "st", stray_params())).tables(cache)["stray_light"]
    assert t2.provenance["spec_hash"] == t.provenance["spec_hash"]
    for tg in STRAY_TG:
        off = make_block(BlockSpec("stray_light_table", "st", stray_params(toggles=dict(STRAY_TG, **{tg: False}))))
        vo = off.tables(TableCache(tmp_path / tg))["stray_light"].values
        m = tg.removeprefix("stray_")
        assert vo[1, col[f"{m}_gate"]] == 0.0 and vo[1, col[f"{m}_mode"]] == 0.0
        assert v[1, col[f"{m}_gate"]] > 0.0
    none = make_block(BlockSpec("stray_light_table", "st", stray_params(toggles=dict.fromkeys(STRAY_TG, False))))
    assert not np.any(none.tables(TableCache(tmp_path / "none"))["stray_light"].values)
    # duty: P̄ = P_peak·D; orthogonal PM polarization scales by ε_xt + T_⊥
    half = sc.stray_photons(stray_params(), STRAY_TG, 0.5e-3)
    assert half["photons_per_mode"] == pytest.approx(0.5 * v[1, col["photons_per_mode"]], rel=1e-12)
    pm = stray_params(probe=dict(stray_params()["probe"], polarization="pm_orthogonal", eps_xt=1e-3, t_perp=1e-3),
                      toggles=dict(STRAY_TG, stray_reflections=False))
    same = sc.stray_photons(stray_params(), dict(STRAY_TG, stray_reflections=False), 1e-3)
    assert sc.stray_photons(pm, pm["toggles"], 1e-3)["clicks_per_gate"] == pytest.approx(
        2e-3 * same["clicks_per_gate"], rel=1e-12)


def test_frequency_offset_and_wdm():
    p = stray_params(probe=dict(stray_params()["probe"], freq_offset_hz=7e9), wdm_isolation_lin=1e-3,
                     mode=dict(shape="exponential", tau_eff_s=14.5e-9))
    out = sc.stray_photons(p, STRAY_TG, 1e-3)
    base = sc.stray_photons(stray_params(mode=dict(shape="exponential", tau_eff_s=14.5e-9)), STRAY_TG, 1e-3)
    assert out["rayleigh_gate"] == pytest.approx(1e-3 * base["rayleigh_gate"], rel=1e-12)  # inside the 100 GHz filter
    w = sc.mode_weight(7e9, p["mode"]) / sc.mode_weight(0.0, p["mode"])
    assert out["rayleigh_mode"] == pytest.approx(1e-3 * w * base["rayleigh_mode"], rel=1e-12)
    assert w < 2.2e-5  # T17 §6: |ũ(7 GHz)|²/τ_eff ≤ 2.2e-5 (Borealis, one-sided exponential bound)


def test_raman_gawbs_brillouin_tables(tmp_path):
    rt = make_block(BlockSpec("raman_noise_table", "r", dict(gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood",
                                                            lambda0_m=LAMBDA0, temperature_k=295.0,
                                                            f_hz=[-13.2e12, 0.0, 1e11, 13.2e12])))
    v = rt.tables(TableCache(tmp_path))["raman_noise"].values
    assert v[1, 0] == 0.0 and v[1, 1] == v[1, 2] > 0  # Ω = 0: finite thermal limit, no spontaneous part
    assert v[3, 1] - v[3, 2] == pytest.approx(v[3, 0])  # Stokes: n_th + 1
    assert v[0, 1] == v[0, 2]  # anti-Stokes: n_th only
    gt = make_block(BlockSpec("gawbs_table", "g", dict(gawbs=dict(GAWBS_SMF), lambda0_m=LAMBDA0, temperature_k=295.0,
                                                      L_m=63.0, f_hz=[1e6, 323.3e6])))
    tabs = gt.tables(TableCache(tmp_path))
    assert tabs["gawbs_peaks"].values.shape == (12, 2) and tabs["gawbs_phase_psd"].values[1] > 0
    bt = make_block(BlockSpec("brillouin_table", "b", dict(brillouin=dict(BRILLOUIN), lengths_m=[63.0, 630.0],
                                                          alpha_per_m=0.0, power_w=1e-3, laser_linewidth_hz=0.0,
                                                          temperature_k=295.0, lambda0_m=LAMBDA0)))
    bv = bt.tables(TableCache(tmp_path))["brillouin"].values
    assert bv[1, 1] == pytest.approx(bv[0, 1] / 10) and bv[1, 2] == pytest.approx(10 * bv[0, 2])
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("gawbs_table", "g", dict(gawbs={k: v for k, v in GAWBS_SMF.items() if k != "v_s_m_s"},
                                                      lambda0_m=LAMBDA0, temperature_k=295.0, L_m=63.0, f_hz=[1e6])))


def test_sbs_check_warns_and_never_raises():
    p = dict(brillouin=dict(BRILLOUIN), alpha_per_m=0.2 * np.log(10) / 1e4, L_m=630.0, power_w=15.4e-3,
             laser_linewidth_hz=0.0, gain_warn=5.0)
    b = make_block(BlockSpec("sbs_check", "sbs", p))
    with pytest.warns(ChkWarning, match="sbs_gain"):
        r = b.check()
    assert r["g_b"] == pytest.approx(5.98, rel=0.01) and [e.code for e in b.chk_log.events] == ["sbs_gain"]
    quiet = make_block(BlockSpec("sbs_check", "sbs", dict(p, L_m=63.0)))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        quiet.check()
    over = make_block(BlockSpec("sbs_check", "sbs", dict(p, power_w=0.1)))
    with pytest.warns(ChkWarning):
        over.check()
    assert {e.code for e in over.chk_log.events} == {"sbs_gain", "sbs_threshold"}
