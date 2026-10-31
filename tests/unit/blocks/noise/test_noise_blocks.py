"""Unit tests of the WP-B2 NOISE blocks: contract (D34, toggles, bands, slots), toggles off ⇒ zero, source
independence (switching one source off leaves the others unchanged, G11/G13), compiled = reference (1e−12,
VT-44 rule), exact discretisation of the generators. Physics VTs: tests/vt/test_vt_34_laser.py, _37_."""

import copy
import math

import numpy as np
import pytest
import scipy.integrate

import cranberry.blocks.noise.testing as nt
from cranberry.blocks.err.testing import rig
from cranberry.blocks.noise._lgs import cascade, discretize, flicker, lowpass, resonance, section_psd, two_pole
from cranberry.blocks.noise.laser import laser_oob_variance, s_nu_model
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.pipeline import iter_leaves
from cranberry.offline.band_ledger import BandLedger

pytestmark = pytest.mark.numba
TOPO = nt.seg_topo({"a": (5.0, 1, "uniform"), "b": (9.0, 4, "per_sub"), "xtal": (0.045, 1, "uniform")})


def all_blocks():
    return [nt.laser(), nt.thermal(TOPO, ["a", "b"], {"a": 30.0, "b": 300.0}, rate_div=10),
            nt.crystal(), nt.pd_offset(),
            nt.mech(TOPO, ["a", "b"], spatial_coherence={"b": 0.3}),
            nt.thermo(TOPO, ["a", "b"]), nt.gawbs()]


def run(blocks, n, fs=1e6, compiled=True, seed=3, external=None):
    root, st, pa = rig(blocks, fs, seed=seed, external=external or {})
    (run_compiled if compiled else run_reference)(root, st, pa, n)
    return st, root


def series(blocks, n, slots, fs=1e6, seed=3):
    from cranberry.blocks.err.testing import record
    return record(blocks, fs, n, slots, seed=seed, external={})


# ------------------------------------------------------------------------------------------- contract


@pytest.mark.parametrize("obj, key", [
    (nt.laser(), "h0_hz2_per_hz"), (nt.laser(), "rin_psd"), (nt.laser(), "oob_pairs"),
    (nt.thermal(TOPO, ["a"], {"a": 1.0}), "r_T_k_per_s"), (nt.thermal(TOPO, ["a"], {"a": 1.0}), "tau_th_s"),
    (nt.crystal(), "dT_fwhm_k"), (nt.pd_offset(), "offset_v"), (nt.mech(TOPO, ["a"]), "k_eps_rad_per_m"),
    (nt.mech(TOPO, ["a"]), "rho_0"), (nt.mech(TOPO, ["b"], spatial_coherence={"b": 0.0}), "spatial_coherence"),
    (nt.thermo(TOPO, ["a"]), "accumulate"), (nt.gawbs(), "sigma2_rad2"),
])
def test_d34_missing_parameter_raises(obj, key):
    obj = copy.deepcopy(obj)
    del obj["params"][key]
    with pytest.raises(MissingParameterError):
        rig([obj], 1e5, external={})


@pytest.mark.parametrize("obj", all_blocks() + [{"block": "filter_cavity", "name": "fc", "params": {
    "segment": "fc", "beam": "lo", "dnu_fc_hz": 1e6, "t0": 0.9, "source": "test", "toggles": {}}}])
def test_d34_missing_toggle_raises(obj):
    obj = copy.deepcopy(obj)
    tg = obj["params"]["toggles"]
    if tg:
        tg.pop(next(iter(tg)))
    with pytest.raises(MissingParameterError):
        rig([obj], 1e5, external={})


def test_band_limits_enforced():
    with pytest.raises(ValueError, match="rin_psd"):
        rig([nt.laser(rin_psd=[[1e4, 1e-14]])], 1e5, external={})
    with pytest.raises(ValueError, match="dither"):
        rig([nt.laser(dither_f_hz=6e4)], 1e5, external={})
    with pytest.raises(ValueError, match="peak"):
        rig([nt.laser(peaks=[{"f_hz": 6e4, "q": 1, "area": 1, "quantity": "freq"}])], 1e5, external={})
    import cranberry.path.testing as pt
    topo = pt.topo_spec([pt.seg_spec("a", 1.0, drivers={"strain": "uniform"})],
                        [pt.beam_spec("b", [pt.step_spec("a", 1)], "o")], ["o"], pt.coef_spec("a", "strain", "*", 1.0))
    with pytest.raises(ValueError, match="no dT"):
        rig([nt.thermal(topo, ["a"], {"a": 1.0})], 1e5, external={})


def test_slots_units_and_shapes():
    st, root = run(all_blocks(), 1)
    lay = st.layout
    assert lay.slot("laser.dx1.phi").unit == "rad" and lay.slot("laser.dx1.power").unit == "W"
    assert lay.slot("drv.a.dT").shape == () and lay.slot("drv.b.dT").shape == (4,)
    assert lay.slot("drv.b.strain").shape == (4,) and lay.slot("drv.b.strain").unit == "1"
    assert lay.slot("opa.xt.dk").unit == "1/m" and lay.slot("v.voff").unit == "V"
    assert all(leaf.runtime_params == () for leaf in iter_leaves(root))


# ------------------------------------------------------------------------- compiled = reference


def test_compiled_equals_reference():
    n = 3000
    a, _ = run(all_blocks(), n, compiled=False)
    b, _ = run(all_blocks(), n, compiled=True)
    for name in a.layout.slot_names():
        x, y = np.asarray(a[name], dtype=float), np.asarray(b[name], dtype=float)
        np.testing.assert_allclose(y, x, rtol=1e-12, atol=1e-300, err_msg=name)
    assert np.std(np.asarray(b["drv.b.strain"])) > 0 and b["laser.dx1.phi"] != 0


# --------------------------------------------------------------------------- toggles off ⇒ zero


def test_toggles_off_give_zero_sources():
    blocks = [nt.laser(on=()), nt.thermal(TOPO, ["a", "b"], {"a": 3.0, "b": 3.0}, on=()), nt.crystal(on=()),
              nt.pd_offset(on=()), nt.mech(TOPO, ["a", "b"], on=(), spatial_coherence={"b": 0.3}),
              nt.thermo(TOPO, ["a", "b"], on=(), name="td2", accumulate=True)]
    rec = series(blocks, 500, ["laser.dx1.phi", "laser.dx1.dnu", "laser.dx1.dp", "laser.dx1.power", "drv.a.dT",
                               "drv.b.dT", "drv.xtal.dT", "v.voff", "drv.a.strain", "drv.b.strain", "opa.xt.dk"])
    for k, v in rec.items():
        if k == "laser.dx1.power":
            assert np.all(v == 3.33e-3)
        elif k == "opa.xt.dk":
            np.testing.assert_allclose(v, 2 * 2.784 * 0.01 / (0.045 * 2.7), rtol=1e-14)  # static detuning only
        else:
            assert np.all(v == 0.0), k


def test_gawbs_toggle_and_ledger_only():
    _, root = run([nt.gawbs(on=False)], 1)
    assert root.bands()[0].variance_outside == 0.0
    _, root = run([nt.gawbs(sigma2=2e-4)], 1)
    led = BandLedger(root.bands(), beams={"sqz": ()})
    led.define_metric("i", f_max_hz=1e9, weight=None)
    assert led.hf_variance("i", {"sqz": 1.0}) == pytest.approx(2e-4)


# ------------------------------------------------------- one source off leaves the others unchanged


def test_laser_sources_superpose_and_are_independent():
    n, s = 4000, ["laser.dx1.dnu", "laser.dx1.dp"]
    tot = series([nt.laser(on=("fm_white", "fm_flicker", "peaks", "rin"))], n, s)
    parts = [series([nt.laser(on=on)], n, s) for on in (("fm_white",), ("fm_flicker",), ("peaks",), ("rin",))]
    for k in s:
        np.testing.assert_allclose(sum(p[k] for p in parts), tot[k], rtol=0, atol=1e-9 * np.abs(tot[k]).max())
    assert np.std(parts[0]["laser.dx1.dnu"]) > 0 and np.std(parts[3]["laser.dx1.dp"]) > 0


def test_mech_parts_superpose_and_blocks_independent():
    n, s = 3000, ["drv.a.strain"]
    tot = series([nt.mech(TOPO, ["a", "b"], spatial_coherence={"b": 0.3})], n, s)
    m = series([nt.mech(TOPO, ["a", "b"], on=("mech", "correlation"), spatial_coherence={"b": 0.3})], n, s)
    p = series([nt.mech(TOPO, ["a", "b"], on=("peaks", "correlation"), spatial_coherence={"b": 0.3})], n, s)
    np.testing.assert_allclose(m["drv.a.strain"] + p["drv.a.strain"], tot["drv.a.strain"], rtol=1e-12, atol=1e-25)
    # another block switched off does not change this one
    alone = series([nt.laser()], n, ["laser.dx1.phi"])
    both = series([nt.laser(), nt.mech(TOPO, ["a"], on=())], n, ["laser.dx1.phi"])
    np.testing.assert_array_equal(alone["laser.dx1.phi"], both["laser.dx1.phi"])


def test_thermodynamic_accumulates_onto_mech():
    n, s = 2000, ["drv.a.strain"]
    base = series([nt.mech(TOPO, ["a"])], n, s)
    tot = series([nt.mech(TOPO, ["a"]), nt.thermo(TOPO, ["a"])], n, s)
    own = series([nt.thermo(TOPO, ["a"], accumulate=False)], n, s)
    np.testing.assert_allclose(tot["drv.a.strain"] - base["drv.a.strain"], own["drv.a.strain"], rtol=1e-6,
                               atol=1e-12 * np.abs(base["drv.a.strain"]).max())


# --------------------------------------------------------------------------------- deterministic parts


def test_thermal_ramp_lag_exact_and_per_sub_equal():
    fs, n = 100.0, 4000
    rec = series([nt.thermal(TOPO, ["a", "b"], {"a": 3.0, "b": 10.0}, on=("drift",))], n, ["drv.a.dT", "drv.b.dT"], fs=fs)
    t = np.arange(n) / fs
    for k, tau in (("drv.a.dT", 3.0), ("drv.b.dT", 10.0)):
        want = 1e-4 * (t - tau * (1 - np.exp(-t / tau)))
        np.testing.assert_allclose(rec[k].reshape(n, -1), want[:, None] * np.ones((1, rec[k].reshape(n, -1).shape[1])),
                                   rtol=1e-9, atol=1e-15)
    st, _ = run([nt.thermal(TOPO, ["b"], {"b": 10.0})], 777, fs=fs)
    assert np.all(st["drv.b.dT"] == st["drv.b.dT"][0])


def test_crystal_dk_and_pd_offset():
    fs, n = 100.0, 1000
    rec = series([nt.crystal(on=("drift",), tau_th_s=1.0), nt.pd_offset(on=("offset", "drift"))], n,
                 ["drv.xtal.dT", "opa.xt.dk", "v.voff"], fs=fs)
    np.testing.assert_allclose(rec["opa.xt.dk"] * 0.045 / 2, 2.784 * (0.01 + rec["drv.xtal.dT"]) / 2.7, rtol=1e-13)
    np.testing.assert_allclose(rec["v.voff"], 5e-3 + 1e-6 * np.arange(n) / fs, rtol=1e-13)


def test_bands_declared():
    root, _, _ = rig([nt.laser(oob_pairs=[{"target": "d1", "tau_s": 49e-9}, {"target": "d2", "tau_s": 1e-6}],
                               rin_targets=["lo"])] + all_blocks()[1:], 1e6, external={})
    decls = root.bands()
    led = BandLedger(decls)  # no double registration
    las = [d for d in decls if d.source == "las"]
    assert {(d.quantity, d.targets) for d in las} == {("freq", ("d1",)), ("freq", ("d2",)), ("rin", ("lo",))}
    s_nu = s_nu_model(nt.laser()["params"])
    f = np.array([1e6, 3e6])
    d1 = next(d for d in las if d.targets == ("d1",))
    np.testing.assert_allclose(d1.psd_outside(f), 4 * np.sin(np.pi * f * 49e-9) ** 2 * s_nu(f))
    assert d1.f_hi_hz == 5e5
    rin = next(d for d in las if d.quantity == "rin")
    np.testing.assert_allclose(rin.psd_outside(np.array([1e7, 5e8, 2e9])), [1e-14, 1e-15, 0.0],
                               rtol=2e-5, atol=1e-23)  # + the 177 kHz rin peak's tail
    assert {d.quantity for d in decls if d.source == "mn"} == {"strain"}
    assert {d.targets for d in decls if d.source == "th"} == {("a",), ("b",)}
    assert len(led.decls) == len(decls)


# ------------------------------------------------------------------------------------- generators


@pytest.mark.parametrize("sec", [lowpass(3.0, 2.0), two_pole(1e-3, 0.1, 1e3), resonance(40.0, 5.0, 1e-6),
                                 cascade(two_pole(1e-3, 0.1, 1e3), 10.0, 0.3, 1.0)])
def test_exact_discretisation_stationary(sec):
    phi, lq, p_half = discretize(sec, 1e-3)
    P = p_half @ p_half.T
    np.testing.assert_allclose(phi @ P @ phi.T + lq @ lq.T, P, rtol=1e-9, atol=1e-12 * np.abs(P).max())
    var = float(sec.C @ P @ sec.C)
    integ = scipy.integrate.quad(lambda f: section_psd(sec, f)[0], 0, np.inf, limit=500, points=None)[0]
    assert var == pytest.approx(integ, rel=1e-4)


def test_flicker_sum_approximates_one_over_f():
    f = np.geomspace(10.0, 1e4, 200)
    s = sum(section_psd(x, f) for x in flicker(2.0, 1.0, 1e5, 2))
    np.testing.assert_allclose(s * f / 2.0, 1.0, rtol=0.04)


def test_laser_oob_variance_white_closed_form():
    h0, tau, fa = 8e3, 49e-9, 5e5
    # ∫_a^∞ 4 sin²(πfτ) h0/f² df = 2h0[1/a − cos(2πaτ)/a + 2πτ(π/2 − Si(2πaτ))]
    from scipy.special import sici
    x = 2 * np.pi * fa * tau
    exact = 2 * h0 * ((1 - np.cos(x)) / fa + 2 * np.pi * tau * (np.pi / 2 - sici(x)[0]))
    got = laser_oob_variance(lambda f: h0 + 0 * f, tau, fa, 1e12)
    assert got == pytest.approx(exact, rel=1e-4)
    assert math.isclose(laser_oob_variance(lambda f: h0 + 0 * f, 0.0, 1.0, 2.0), 0.0)
