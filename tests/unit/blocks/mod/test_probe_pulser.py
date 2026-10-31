"""``probe_pulser`` (WP-B8): leak slots, toggles, compiled == reference, per-port gating of B10's stray table."""

import numpy as np
import pytest

import cranberry.blocks.dmx.testing as dt
from cranberry.blocks.fib.scattering import resolve_mode, stray_photons
from cranberry.blocks.mod.probe_pulser import gate_factors, gated_stray_photons
from cranberry.blocks.source.pulsed_testing import BOREALIS, TOGGLES_ON
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_reference
from cranberry.path.topology import C_M_PER_S

TG = {"stray_rayleigh": True, "stray_brillouin": True, "stray_raman": False, "stray_reflections": True}


def leak(obj):
    b, st, pr = dt.make(obj, 1e6)
    run_reference(b, st, pr, 1)
    return {k: float(st[f"leak.pp.{k}"]) for k in ("p_avg_w", "p_peak_w", "duty", "t_offset_s", "eps_aom")}


def test_leak_slots_and_toggles():
    d = 10 / 320
    out = leak(dt.pulser(p_avg_w=1e-6, eps_aom=1e-3))
    assert out["p_peak_w"] == pytest.approx(1e-6 / d) and out["duty"] == d and out["eps_aom"] == 1e-3
    assert out["p_avg_w"] == pytest.approx(1e-6 + 1e-3 * (1 - d) * 1e-6 / d, rel=1e-15)
    assert leak(dt.pulser(p_avg_w=1e-6, eps_aom=1e-3, toggles={"aom_leak": False}))["p_avg_w"] == pytest.approx(1e-6)
    cw = leak(dt.pulser(p_avg_w=1e-6, eps_aom=1e-3, toggles={"gating": False}))
    assert cw == {"p_avg_w": 1e-6, "p_peak_w": 1e-6, "duty": 1.0, "t_offset_s": 5e-9, "eps_aom": 0.0}
    assert dt.compiled_matches_reference(dt.pulser(eps_aom=1e-3), 1e6, 5) <= 1e-12
    obj = dt.pulser()
    del obj["params"]["t_offset_s"]
    with pytest.raises(MissingParameterError):
        dt.make(obj, 1e6)


def test_gate_factor_limits():
    # CW-like train (duty → 1 of the period): factor 1
    f, _ = gate_factors(period_s=100e-9, tau_pr_s=100e-9, t_offset_s=3e-9, t_gate_s=1e-9, roundtrip_s=37e-9)
    assert f == pytest.approx(1.0, rel=1e-12)
    # the gate average over all offsets is the duty (energy conservation of the backscatter)
    offs = np.linspace(0, 320e-9, 641)[:-1]
    fs = [gate_factors(period_s=320e-9, tau_pr_s=10e-9, t_offset_s=o, t_gate_s=1e-9, roundtrip_s=617e-9)[0] for o in offs]
    assert np.mean(fs) == pytest.approx(10 / 320, rel=1e-3)
    # a reflection: rectangle overlap with the gate
    _, (r,) = gate_factors(period_s=320e-9, tau_pr_s=10e-9, t_offset_s=-5e-9, t_gate_s=1e-9, roundtrip_s=1e-9,
                           delays_s=[0.0])
    assert r == pytest.approx(1.0)


def test_gated_stray_cw_and_reflections():
    p = dt.stray_params(30.0, reflections=[[1e-6, 1.0]])
    kw = dict(p_avg_w=1e-6, duty=10 / 320, period_s=320e-9, eps_aom=0.0, n_g=dt.N_G, reflection_z_m=[2.0])
    cw = gated_stray_photons(p, TG, gating=False, t_offset_s=0.0, **kw)
    assert cw == stray_photons(p, TG, 1e-6)
    # reflection at z = 2 m returns 19.6 ns later: offset −24.6 ns puts it on the gate, 5 ns misses it
    d = 2 * dt.N_G * 2.0 / C_M_PER_S
    hit = gated_stray_photons(p, TG, gating=True, t_offset_s=-d - 5e-9, **kw)
    miss = gated_stray_photons(p, TG, gating=True, t_offset_s=5e-9, **kw)
    assert hit["reflections_gate"] == pytest.approx(cw["reflections_gate"] * 320 / 10, rel=1e-12)
    assert miss["reflections_gate"] == 0.0
    assert miss["clicks_per_gate"] == pytest.approx(sum(miss[k] for k in miss if k.endswith("_gate") and k != "clicks_per_gate"))
    # AOM leak keeps a CW floor ε·P_pk
    leaky = gated_stray_photons(p, TG, gating=True, t_offset_s=5e-9, **{**kw, "eps_aom": 1e-3})
    assert leaky["reflections_gate"] == pytest.approx(1e-3 * cw["reflections_gate"] * 320 / 10, rel=1e-12)


def test_photons_per_mode_weighted_by_the_mode_not_the_gate():
    """M5-02 (conv §2): Borealis-type ns mode (x_ss = 2.4, τ_eff 14.5 ns), 25 ns gate on the mode centroid,
    10 ns reflection pulse (no delay). The gate factor is 0.40 for every start; the mode factor
    |∫u·1_pulse|²/τ_eff is 0.035–0.47 (reviews/scripts/m5_stray_gate_vs_mode.py)."""
    src = {"block": "pulsed_opo_source", "name": "opo", "params": dict(BOREALIS, x_ss=2.4, toggles=TOGGLES_ON)}
    tg = dict.fromkeys(TG, False) | {"stray_reflections": True}
    p = dt.stray_params(30.0, reflections=[[1e-6, 1.0]], t_gate_s=25e-9, mode={"source": src})
    p["mode"] = resolve_mode(p["mode"])  # build the nested source once
    cw = stray_photons(p, tg, 1e-6)
    tau_eff = cw["reflections_mode"] / (cw["reflections_gate"] / 25e-9)  # CW: N = Φ τ_eff, clicks = Φ T_gate
    assert tau_eff == pytest.approx(14.5e-9, rel=0.02)
    kw = dict(p_avg_w=1e-6 * 1e-2, duty=1e-2, period_s=1e-6, eps_aom=0.0, gating=True, n_g=dt.N_G, reflection_z_m=[0.0])
    rows = {}
    for start, mode_f in ((-12.5e-9, 0.035), (-7.5e-9, 0.431), (-2.5e-9, 0.471), (2.5e-9, 0.090)):
        out = gated_stray_photons(p, tg, t_offset_s=start, **kw)
        rows[start] = out["reflections_mode"] / cw["reflections_mode"]
        assert out["reflections_gate"] / cw["reflections_gate"] == pytest.approx(0.400, abs=1e-12)
        assert rows[start] == pytest.approx(mode_f, abs=0.005)  # review: bin-centre indicator, here exact edges
    assert max(rows.values()) / min(rows.values()) > 10
    # CW-like train (duty 1): mode factor 1 for lines (reflection) and distributed scattering
    tg2 = dict(tg, stray_rayleigh=True)
    full = gated_stray_photons(p, tg2, t_offset_s=0.0, **{**kw, "duty": 1.0, "p_avg_w": 1e-6})
    cwr = stray_photons(p, tg2, 1e-6)
    for col in ("reflections_mode", "rayleigh_mode", "photons_per_mode"):
        assert full[col] == pytest.approx(cwr[col], rel=1e-6)
