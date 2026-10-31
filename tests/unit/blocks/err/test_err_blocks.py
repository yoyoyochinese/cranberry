"""Unit tests of the ERR blocks (WP-B3): contract (D34, wiring, bands, toggles), static signal shapes,
compiled = reference (1e−12, VT-44 rule). Physics validation (T14 V1–V12) is in tests/vt/test_vt_43_*."""

import math

import numpy as np
import pytest

from cranberry.blocks.err import (
    amplitude_gains_from_table,
    amplitude_gains_opa,
    amplitude_gains_opo,
    gamma_coh_from_ledger,
    qnl_band_integrals,
)
from cranberry.blocks.err.testing import (
    DET,
    REF_A_DET,
    REF_B,
    bhd,
    dcf,
    dem,
    liu,
    pbd,
    qnl,
    qnl_flat,
    rfr,
    rig,
)
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.pipeline import iter_leaves
from cranberry.offline.band_ledger import BandLedger
from cranberry.offline.tables import Table

E = 1.602176634e-19
EXT_DCF = {"beam.lo.phi": "rad", "beam.pr.phi": "rad", "v.dp": "1", "v.voff": "V"}


def one_step(blocks, phases=None, f_sim=16e6, n=1):
    root, s, p = rig(blocks, f_sim)
    for k, v in (phases or {}).items():
        s[k] = v
    run_reference(root, s, p, n)
    return s, p, root


# ---------------------------------------------------------------- contract


def test_d34_missing_parameter_and_toggle_raise():
    obj = pbd()
    del obj["params"]["eta"]
    with pytest.raises(MissingParameterError):
        rig([obj], 1e6)
    obj = pbd()
    del obj["params"]["toggles"]["shot_noise"]
    with pytest.raises(MissingParameterError):
        rig([obj], 1e6)
    obj = pbd()
    del obj["params"]["mode"]
    with pytest.raises(MissingParameterError):
        rig([obj], 1e6)
    with pytest.raises(ValueError):
        rig([pbd(mode="bogus")], 1e6)


def test_null_input_with_toggle_on_raises():
    with pytest.raises(ValueError):
        rig([dcf(toggles={"rin": True})], 1e6)
    with pytest.raises(ValueError):
        rig([pbd(), dem(toggles={"rf_reference_phase_noise": True})], 1e6)


def test_demodulator_rejects_mode_mismatch():
    with pytest.raises(ValueError, match="must agree"):
        rig([pbd(rf=False), dem(rf=True)], 16e6)


def test_reads_writes_units_and_runtime_params():
    root, s, p = rig([pbd(noise=True), dem()], 16e6)
    det, d = list(iter_leaves(root))
    assert set(det.reads) == {"beam.pr.phi", "beam.p.phi"}
    assert "v.pbd" in det.writes and det.slot_units["v.pbd"] == "V"
    assert p["pbd.eta"] == 0.8 and s["v.pbd"].shape == (2,)
    assert d.reads == ("v.pbd",) and det.losses() == [] and d.bands() == []


def test_bands_sensor_and_rf_reference():
    root, _, _ = rig([pbd(noise=True), dem(), rfr()], 16e6)
    det, _, ref = list(iter_leaves(root))
    (b,) = det.bands()
    assert b.quantity == "sensor" and b.f_hi_hz == 8e6 and b.variance_outside == 0.0
    (r,) = ref.bands()
    assert r.quantity == "phase" and r.targets == ("beam.pr",)
    np.testing.assert_allclose(r.psd_outside(np.array([1e7])), 1e-12 + 1e-2 / 1e14)
    led = BandLedger()
    led.declare_all([b, r])
    led.define_metric("m", f_max_hz=1e9, weight=None)
    var = led.hf_variance("m", {"beam.pr": 1.0})
    assert var > 0 and gamma_coh_from_ledger(led, "m", {"beam.pr": 1.0}) == pytest.approx(math.exp(-var / 2))
    root, _, _ = rig([pbd(noise=False), rfr(on=False)], 16e6)
    assert all(leaf.bands() == [] for leaf in iter_leaves(root))


def test_amplitude_gain_helpers():
    ap, am = amplitude_gains_opa(1.0)
    assert ap * am == pytest.approx(1.0) and ap == pytest.approx(math.e)
    ap, am = amplitude_gains_opo(0.5)
    assert (ap, am) == pytest.approx((2.0, 2.0 / 3.0))
    with pytest.raises(ValueError):
        amplitude_gains_opo(1.0)
    p = np.array([0.0, 1.0])
    tab = Table("apm", {"p_w": p}, np.array([[1.0, 1.0], [3.0, 1.0 / 3.0]]), "1")
    assert amplitude_gains_from_table(tab, p_w=0.5) == pytest.approx((2.0, 2.0 / 3.0))


def test_qnl_band_integrals_flat_band():
    f = np.linspace(8e6, 300e6, 2001)
    one = np.ones_like(f)
    r = qnl_band_integrals(f, 2.0 * one, 0.3 * one, 10.0 * one, 0.1 * one, one)
    b = 292e6
    pm, pd = b * (2 * 5.15 + 0.1), b * 2 * 4.85
    assert r["p_mean"] == pytest.approx(pm) and r["p_delta"] == pytest.approx(pd)
    flat = {"q0": pm * pm / b, "q1": -2 * pm * pd / b, "q2": pd * pd / b}
    for k, v in flat.items():
        assert r[k] == pytest.approx(v)


# ---------------------------------------------------------------- static signal shapes (noise off)


@pytest.mark.parametrize("phi1", [0.0, 0.7, -2.1])
def test_parametric_heterodyne_baseband(phi1):
    s, p, _ = one_step([pbd(phi_opa_rad=0.2), dem(phi_d_rad=0.4)], {"beam.pr.phi": 0.5 * (phi1 + 0.2)})
    ap, am = REF_A_DET
    k = 1e3 * 1.0 * 0.8e-3 * 0.5 * (ap**2 - am**2)  # G ℜ η P_in sinh 2r_d
    np.testing.assert_allclose(s["v.pbd"], [k * math.cos(phi1), k * math.sin(phi1)], rtol=1e-12, atol=1e-14)
    assert float(s["dem.dem.i"]) == pytest.approx(k * math.cos(phi1 + 0.4), rel=1e-12, abs=1e-14)
    assert float(s["dem.dem.q"]) == pytest.approx(k * math.sin(phi1 + 0.4), rel=1e-12, abs=1e-14)


def test_bhd_dc_and_heterodyne_shapes():
    ap, am = REF_B
    gc, gs = (ap + am) / 2, (ap - am) / 2
    k2 = 1e3 * 2 * math.sqrt(1e-3 * 1e-9)
    u, v = 0.3, 1.1
    s, _, _ = one_step([bhd(mode="dc")], {"beam.pr.phi": u, "beam.lo.phi": v})
    assert float(s["v.bhd"]) == pytest.approx(k2 * (gc * math.cos(u - v) + gs * math.cos(u + v)), rel=1e-12)
    s, _, _ = one_step([bhd(a_plus=ap, a_minus=am), dem(inp="v.bhd", phi_d_rad=0.2)], {"beam.pr.phi": u, "beam.lo.phi": v})
    want = k2 * (gc * math.cos(u + 0.2 - v) + gs * math.cos(u + 0.2 + v))  # T14 §3.4, u′ = u + φ_d2
    assert float(s["dem.dem.i"]) == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("gating", [True, False])
def test_liu_homodyne_amplitude_and_gating(gating):
    s, _, _ = one_step([liu(gating=gating, v_off_v=0.01)], {"beam.pr.phi": 0.4, "beam.ref.phi": 0.1})
    p_eff = 1e-6 / 0.1 if gating else 0.1 * 1e-6
    k = 1e4 * 2 * math.sqrt(p_eff * 1e-3)
    assert float(s["v.liu"]) == pytest.approx(k * math.cos(0.3) + 0.01, rel=1e-12)


def test_liu_same_laser_cancels_different_laser_does_not():
    """D40: a laser phase common to probe and reference (same laser, matched delays) drops out; an
    independent one (different laser) enters the error 1:1."""
    ev = lambda pr, ref: float(one_step([liu()], {"beam.pr.phi": pr, "beam.ref.phi": ref})[0]["v.liu"])  # noqa: E731
    base = ev(0.4, 0.1)
    assert ev(0.4 + 2.3, 0.1 + 2.3) == pytest.approx(base, rel=1e-12)
    assert ev(0.4 + 0.2, 0.1) == pytest.approx(base * math.cos(0.5) / math.cos(0.3), rel=1e-12)


def test_dc_fringe_eq3_superset():
    """T15 eq. (1)–(3) with every term on: imbalance, coupler sign, μ, γ_coh, PER, offset."""
    kw = dict(kappa1=0.47, kappa2=0.53, t_arm1=0.9, t_arm2=0.97, eps_r=0.02, v_off_v=0.015, mu_abs=0.95,
              mu_phase_rad=0.05, gamma_coh=0.99, per_x=0.01, per_dchi_rad=1.3, coupler_sign=-1.0)
    phid = 0.8
    s, _, _ = one_step([dcf(**kw)], {"beam.lo.phi": phid})
    p1, p2 = 0.47 * 0.9 * 3.33e-3, 0.53 * 0.97 * 3.33e-3
    pd = p1 + p2
    eps = -(2 * 0.53 - 1) * (p1 - p2) / pd + 0.02
    vis = 4 * math.sqrt(0.53 * 0.47 * p1 * p2) / pd
    mu = 0.95 * np.exp(0.05j) * 0.99 * ((1 - 0.01) + 0.01 * np.exp(1.3j))
    want = 1e3 * pd * (eps - vis * abs(mu) * math.cos(phid + np.angle(mu))) + 0.015
    assert float(s["v.dcf"]) == pytest.approx(want, rel=1e-12)


# ---------------------------------------------------------------- toggles off => identity / zero source


def test_rf_reference_off_is_zero_source():
    s, _, _ = one_step([rfr(on=False)], n=100)
    assert float(s["rf.rfr.phi"]) == 0.0


@pytest.mark.parametrize("make", [
    lambda nz: pbd(noise=nz), lambda nz: pbd(mode="dither", noise=nz), lambda nz: bhd(noise=nz),
    lambda nz: bhd(mode="dither", noise=nz), lambda nz: liu(noise=nz), lambda nz: liu(mode="heterodyne", noise=nz),
    lambda nz: qnl(noise=False, toggles={"radiometer_noise": nz}), lambda nz: dcf(noise=nz)])
def test_noise_toggles_off_give_the_noiseless_signal(make):
    phases = {"beam.pr.phi": 0.3, "beam.lo.phi": 0.2, "beam.p.phi": 0.1}
    s_off, _, root = one_step([make(False)], phases)
    leaf = next(iter_leaves(root))
    assert leaf.draws is None
    s_on, _, _ = one_step([make(True)], phases)
    out = leaf.params["out"]
    assert not np.allclose(s_on[out], s_off[out], rtol=1e-6, atol=0)  # the noise is really there when on


@pytest.mark.parametrize("toggle,kw", [
    ("optical_imbalance", dict(eps_r=0.05, kappa1=0.4, kappa2=0.6)),
    ("per_crosscoupling", dict(per_x=0.02, per_dchi_rad=1.0)),
    ("visibility_hf_averaging", dict(gamma_coh=0.9)),
    ("saturation", dict(p_in_w=10e-3)),
])
def test_dc_fringe_feature_toggles_off_remove_the_feature(toggle, kw):
    phases = {"beam.lo.phi": 0.0}
    base, _, _ = one_step([dcf(kappa1=kw.get("kappa1", 0.5), kappa2=kw.get("kappa2", 0.5),
                               p_in_w=kw.get("p_in_w", 3.33e-3), toggles={"saturation": False})], phases)
    on, _, _ = one_step([dcf(**kw)], phases)
    off, _, _ = one_step([dcf(**kw, toggles={toggle: False, "saturation": toggle == "saturation" and False})], phases)
    if toggle == "optical_imbalance":  # ε → 0, 𝒱 keeps its (second-order) κ dependence
        assert float(off["v.dcf"]) != float(on["v.dcf"])
        assert float(on["v.dcf"]) - float(off["v.dcf"]) == pytest.approx(
            1e3 * 3.33e-3 * ((2 * 0.6 - 1) * (0.4 - 0.6) + 0.05), rel=1e-9)
    else:
        assert float(off["v.dcf"]) == pytest.approx(float(base["v.dcf"]), rel=1e-12)
        assert float(on["v.dcf"]) != pytest.approx(float(base["v.dcf"]), rel=1e-6)


def test_dc_fringe_rin_and_offset_inputs():
    obj = dcf(in_dp="v.dp", in_voff="v.voff", toggles={"rin": True, "offset_drift": True})
    root, s, p = rig([obj], 1e6, external=EXT_DCF)
    s["beam.lo.phi"], s["v.dp"], s["v.voff"] = 0.0, 0.1, 0.002
    run_reference(root, s, p, 1)
    assert float(s["v.dcf"]) == pytest.approx(3.33 * 1.1 + 0.002, rel=1e-12)
    obj = dcf(in_dp="v.dp", in_voff="v.voff")  # toggles off: inputs not even read
    root, s, p = rig([obj], 1e6, external=EXT_DCF)
    assert "v.dp" not in next(iter_leaves(root)).reads
    s["v.dp"], s["v.voff"] = 0.1, 0.002
    run_reference(root, s, p, 1)
    assert float(s["v.dcf"]) == pytest.approx(3.33, rel=1e-12)


def test_ram_and_leak_toggles_off_and_demod_offset():
    ph = {"beam.lo.phi": 0.5 * math.pi + 0.01}
    a, _, _ = one_step([qnl(ram_a=0.05)], ph)
    b, _, _ = one_step([qnl(ram_a=0.0)], ph)
    c, _, _ = one_step([qnl(ram_a=0.05, toggles={"ram_offset": True})], ph)
    np.testing.assert_array_equal(a["v.qnl"], b["v.qnl"])
    assert not np.allclose(c["v.qnl"], b["v.qnl"])
    s1, _, _ = one_step([pbd(rf=True, f_sh_leak_v=0.3), dem(rf=True)], n=3)
    s2, _, _ = one_step([pbd(rf=True, f_sh_leak_v=0.0), dem(rf=True)], n=3)
    assert float(s1["v.pbd"]) == float(s2["v.pbd"])
    s3, _, _ = one_step([pbd(), dem(v_off_v=0.1)])
    s4, _, _ = one_step([pbd(), dem(v_off_v=0.1, toggles={"demod_offset": True})])
    assert float(s4["dem.dem.i"]) - float(s3["dem.dem.i"]) == pytest.approx(0.1)


def test_demod_rf_reference_phase_toggle():
    objs = [rfr(), pbd(), dem(in_phase="rf.rfr.phi")]
    s_off, _, _ = one_step(objs, n=50)
    s_ref, _, _ = one_step([pbd(), dem()], n=50)
    assert float(s_off["dem.dem.i"]) == float(s_ref["dem.dem.i"])
    s_on, _, _ = one_step([rfr(), pbd(), dem(in_phase="rf.rfr.phi", toggles={"rf_reference_phase_noise": True})], n=50)
    k = np.hypot(*s_on["v.pbd"])
    phi = float(s_on["rf.rfr.phi"])
    a, b = s_on["v.pbd"]
    assert float(s_on["dem.dem.i"]) == pytest.approx(a * math.cos(phi) - b * math.sin(phi), rel=1e-12) and k > 0


# ---------------------------------------------------------------- compiled = reference (VT-44 rule)


@pytest.mark.numba
def test_compiled_equals_reference_all_blocks_with_noise():
    from cranberry.core.rates import PIECE_SAMPLES

    blocks = [
        rfr(), pbd(noise=True), dem(in_phase="rf.rfr.phi", order=2, toggles={"rf_reference_phase_noise": True,
                                                                          "demod_offset": True}, v_off_v=1e-3),
        pbd("pd2", "dither", noise=True, toggles={"ram_offset": True}, ram_a=1e-3), dem("d2", "v.pd2", order=1, rate_div=2),
        pbd("pd3", rf=True, noise=True, toggles={"f_sh_leak": True}, f_sh_leak_v=0.1), dem("d3", "v.pd3", rf=True, order=2),
        bhd(noise=True), bhd("b2", "dither", noise=True), bhd("b3", "dc", noise=True), bhd("b4", "dither", rf=True, noise=True),
        liu(noise=True), liu("l2", "heterodyne", noise=True, gating=False), liu("l3", "dither", rf=True, noise=True),
        qnl(noise=True), qnl("q2", "log", noise=True, toggles={"ram_offset": True}, ram_a=1e-2),
        qnl("q3", "linear", rf=True, noise=True, rate_div=3),
        dcf(noise=True, out_err="v.dcf_err", per_x=0.01, per_dchi_rad=0.3),
    ]
    n = PIECE_SAMPLES + 333
    r1, s1, p1 = rig(blocks, 16e6, seed=5)
    r2, s2, p2 = rig(blocks, 16e6, seed=5)
    for s in (s1, s2):
        s["beam.pr.phi"], s["beam.lo.phi"], s["beam.p.phi"] = 0.3, -0.2, 0.05
    run_reference(r1, s1, p1, n)
    run_compiled(r2, s2, p2, n)
    np.testing.assert_array_equal(s1.I, s2.I)
    scale = np.maximum(np.abs(s1.S), 1.0)
    assert np.max(np.abs(s1.S - s2.S) / scale) < 1e-12
    assert float(s1["dem.d3.i"]) != 0.0 and float(s1["v.dcf_err"]) != 0.0
