"""``probe_pulser`` (WP-B8): leak slots, toggles, compiled == reference, per-port gating of B10's stray table."""

import numpy as np
import pytest

import cranberry.blocks.dmx.testing as dt
from cranberry.blocks.fib.scattering import stray_photons
from cranberry.blocks.mod.probe_pulser import gate_factors, gated_stray_photons
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
