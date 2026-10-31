"""WP-B7 M4 POL blocks: segment_jones SMF refinements (static Haar U0, OU drift spectrum,
re-unitarisation), pdl_element, faraday_nonreciprocal, the beam_assembler J_nr_slot hook (F6 addition),
pmd_mode_overlap. Analytic forms, toggles off ⇒ identity, compiled == reference (1e−12), validation."""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.pol.pdl_element import pdl_matrix
from cranberry.blocks.pol.pmd_mode_overlap import make_pmd_mode_overlap, pulse_projection
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.offline.loss_registry import LossRegistry

ON = {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}


def R(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def rig(pipe, external=(), setup=None, fs=1e6, seed=5):
    lay = SlotLayout(fs)
    if setup:
        setup(lay)
    blk = build(spec_from_obj(pipe), lay, seed, external=list(external))
    return blk, State(lay), Params(lay)


# ------------------------------------------------------------------------------- segment_jones (SMF)


def smf_topo(n):
    segs = [pt.seg_spec(f"s{i}", 3.0, kind="smf", pol_model="jones", pol_params={"J_slot": f"J.s{i}"})
            for i in range(n)]
    return pt.topo_spec(segs, [pt.beam_spec("b", [pt.step_spec(s["id"]) for s in segs], "o", jones_out=(1, 0))],
                        ["o"], [])


def smf_cfg(u0="haar", drift="diffusion", D=0.0, fc=None):
    c = {"model": "smf", "axis_drivers": [], "axis_stokes": [1, 0, 0], "u0": u0, "drift": drift,
         "D_pol_rad2_per_s": D, "k_T_rad_per_m_K": 0.0}
    return c if fc is None else {**c, "f_c_hz": fc}


def sj_rig(segs, n=None, toggles=ON, fs=1e6, seed=5, name="sj"):
    topo = smf_topo(n or len(segs))
    return rig([{"block": "segment_jones", "name": name,
                 "params": {"topology": topo, "segments": segs, "toggles": dict(toggles)}}], fs=fs, seed=seed)


def test_static_haar_u0():
    segs = {f"s{i}": smf_cfg() for i in range(3)}
    blk, s, p = sj_rig(segs)
    run_reference(blk, s, p, 3)
    us = [np.array(s[f"J.s{i}"]) for i in range(3)]
    for u in us:
        np.testing.assert_allclose(u @ u.conj().T, np.eye(2), atol=1e-14)
        assert abs(np.linalg.det(u) - 1) < 1e-14
    assert not np.allclose(us[0], us[1])  # one fiber realisation per segment
    np.testing.assert_allclose(us[0], next(iter_leaves(blk)).u0["s0"], atol=0)
    # same seed ⇒ same fiber, independent of the other segments (order-independent sub-streams)
    blk2, s2, p2 = sj_rig({"s1": smf_cfg(), "s0": smf_cfg(u0="identity")}, n=2)
    run_reference(blk2, s2, p2, 1)
    np.testing.assert_array_equal(s2["J.s1"], s["J.s1"])
    np.testing.assert_array_equal(s2["J.s0"], np.eye(2))
    blk3, s3, _ = sj_rig(segs, seed=6)
    assert not np.allclose(next(iter_leaves(blk3)).u0["s0"], us[0])


def test_u0_does_not_shift_the_drift_draws():
    """The drift W is the same with and without a Haar U0 (U0 comes from its own sub-stream)."""
    (b1, s1, p1), (b2, s2, p2) = (sj_rig({"s0": smf_cfg(u0=u, D=1e3)}) for u in ("identity", "haar"))
    run_reference(b1, s1, p1, 50)
    run_reference(b2, s2, p2, 50)
    np.testing.assert_array_equal(s1["J.sj.s0.w"], s2["J.sj.s0.w"])
    np.testing.assert_allclose(s2["J.s0"], np.array(s1["J.s0"]) @ np.array(s2["J.sj.s0.u0"]), atol=1e-15)


def test_diffusion_stays_exactly_su2():
    blk, s, p = sj_rig({"s0": smf_cfg(D=1e5)})
    run_reference(blk, s, p, 2000)
    w = np.array(s["J.sj.s0.w"])
    assert w[1, 1] == np.conj(w[0, 0]) and w[0, 1] == -np.conj(w[1, 0])  # re-unitarised SU(2) form
    assert abs(abs(w[0, 0]) ** 2 + abs(w[1, 0]) ** 2 - 1) < 1e-15


def test_ou_drift_spectrum():
    """OU: stationary variance σ² = D/(4π f_c) per Stokes component, one-step increment variance
    2σ²(1 − r), r = e^{−2π f_c h}; W = exp(−iΩ·σ/2)."""
    D, fc, fs, n = 2.0, 50.0, 1e4, 40000
    blk, s, p = sj_rig({"s0": smf_cfg(u0="identity", drift="ou", D=D, fc=fc)}, fs=fs)
    om = np.empty((n, 3))
    for k in range(n):
        run_reference(blk, s, p, 1)
        om[k] = s["J.sj.s0.om"]
    sig2, r = D / (4 * math.pi * fc), math.exp(-2 * math.pi * fc / fs)
    assert np.var(om[n // 4:]) == pytest.approx(sig2, rel=0.08)  # ~400 correlation times
    assert np.var(np.diff(om, axis=0)) == pytest.approx(2 * sig2 * (1 - r), rel=0.03)
    assert 2 * sig2 * (1 - r) == pytest.approx(D / fs, rel=0.02)  # = diffusion with the same D_pol
    a, b, c = om[-1]
    th = math.sqrt(a * a + b * b + c * c)
    sig = (a * np.diag([1, -1]) + b * np.array([[0, 1], [1, 0]]) + c * np.array([[0, -1j], [1j, 0]])) / th
    np.testing.assert_allclose(s["J.s0"], math.cos(th / 2) * np.eye(2) - 1j * math.sin(th / 2) * sig, atol=1e-14)


def test_smf_toggle_off_identity_and_validation():
    blk, s, p = sj_rig({"s0": smf_cfg(), "s1": smf_cfg(drift="ou", D=1.0, fc=1.0)},
                       toggles={**ON, "pol_smf_drift": False})
    run_reference(blk, s, p, 5)
    for j in ("J.s0", "J.s1"):
        np.testing.assert_array_equal(s[j], np.eye(2))
    assert next(iter_leaves(blk)).draws is None
    with pytest.raises(MissingParameterError):
        sj_rig({"s0": {k: v for k, v in smf_cfg().items() if k != "u0"}})
    with pytest.raises(MissingParameterError):
        sj_rig({"s0": smf_cfg(drift="ou", D=1.0)})  # f_c_hz missing
    with pytest.raises(ValueError, match="drift"):
        sj_rig({"s0": smf_cfg(drift="brownian")})
    with pytest.raises(ValueError, match="f_c_hz"):
        sj_rig({"s0": smf_cfg(drift="ou", D=1.0, fc=0.0)})


@pytest.mark.numba
def test_smf_compiled_equals_reference():
    segs = {"s0": smf_cfg(D=50.0), "s1": smf_cfg(drift="ou", D=80.0, fc=300.0), "s2": smf_cfg(u0="identity", D=5.0)}
    (b1, s1, p1), (b2, s2, p2) = sj_rig(segs), sj_rig(segs)
    run_reference(b1, s1, p1, 400)
    run_compiled(b2, s2, p2, 400)
    for j in ("J.s0", "J.s1", "J.s2", "J.sj.s1.om"):
        np.testing.assert_allclose(s1[j], s2[j], rtol=0, atol=1e-12)


# --------------------------------------------------------------------------------------- pdl_element


def pdl_spec(on=True, pdl=0.5, ax=0.3, il=0.4):
    return {"block": "pdl_element", "name": "pdl", "params": {
        "out": "J.pdl", "pdl_db": pdl, "axis_rad": ax, "il_db": il, "segment": "pdl", "applies_to": ["gbs", "probe"],
        "source": "[AE]", "toggles": {"pol_pdl": on}}}


def test_pdl_matrix_and_average():
    blk, s, p = rig([pdl_spec()])
    run_reference(blk, s, p, 1)
    q = 10 ** (-0.05)
    M = R(0.3) @ np.diag([math.sqrt(2 / (1 + q)), math.sqrt(2 * q / (1 + q))]) @ R(-0.3)
    np.testing.assert_allclose(s["J.pdl"], M, atol=1e-15)
    np.testing.assert_allclose(M, M.T, atol=0)  # reciprocal
    # |M e|² averaged over the Poincaré sphere = 1; extremes T_max/T_min = PDL
    rng = np.random.default_rng(1)
    e = rng.normal(size=(200000, 2)) + 1j * rng.normal(size=(200000, 2))
    e /= np.linalg.norm(e, axis=1)[:, None]
    t = np.sum(np.abs(e @ M.T) ** 2, axis=1)
    assert np.mean(t) == pytest.approx(1.0, abs=3e-4)
    sv = np.linalg.svd(M, compute_uv=False) ** 2
    assert 10 * math.log10(sv[0] / sv[1]) == pytest.approx(0.5, abs=1e-12)
    leaf = next(iter_leaves(blk))
    np.testing.assert_allclose(leaf.jones(), math.sqrt(10 ** (-0.04)) * M, atol=1e-15)
    reg = LossRegistry(leaf.losses(), paths={"sq": {"segments": ["pdl"]}})
    assert reg.eta("sq", "gbs") == pytest.approx(10 ** (-0.04), rel=1e-14)
    p["pdl.axis_rad"] = 0.0  # runtime parameter
    run_reference(blk, s, p, 1)
    np.testing.assert_allclose(s["J.pdl"], pdl_matrix(0.5, 0.0), atol=1e-15)


def test_pdl_off_identity_and_validation():
    blk, s, p = rig([pdl_spec(on=False)])
    run_reference(blk, s, p, 1)
    np.testing.assert_array_equal(s["J.pdl"], np.eye(2))
    assert next(iter_leaves(blk)).losses()[0].enabled  # the averaged IL is not the PDL
    with pytest.raises(ValueError, match="T_max"):
        rig([pdl_spec(pdl=3.0, il=0.1)])
    with pytest.raises(MissingParameterError):
        rig([{**pdl_spec(), "params": {**pdl_spec()["params"], "toggles": {}}}])


# ----------------------------------------------------------------------------- faraday_nonreciprocal


def far_spec(on=True):
    return {"block": "faraday_nonreciprocal", "name": "far", "params": {
        "out": "J.far", "verdet_rad_per_T_m": 0.6, "B_par_T": 50e-6, "L_m": 63.0, "toggles": {"pol_faraday": on}}}


def test_faraday_rotation_and_toggle():
    blk, s, p = rig([far_spec()])
    run_reference(blk, s, p, 1)
    rho = 0.6 * 50e-6 * 63
    assert next(iter_leaves(blk)).rho_rad == pytest.approx(rho, rel=1e-15)
    assert rho == pytest.approx(1.89e-3, rel=1e-3)  # T09 §1.8: ≈ 2e−3 rad
    np.testing.assert_allclose(s["J.far"], R(rho), atol=1e-16)
    blk, s, p = rig([far_spec(on=False)])
    run_reference(blk, s, p, 1)
    np.testing.assert_array_equal(s["J.far"], np.eye(2))


# ------------------------------------------------------------------- beam_assembler J_nr_slot (F6 add)


def nr_rig(nr=True):
    U = np.array([[0.6, -0.8j], [-0.8j, 0.6]]) @ R(0.4)  # not symmetric: Uᵀ ≠ U
    pol = {"J": U.tolist(), **({"J_nr_slot": "J.far"} if nr else {})}
    seg = pt.seg_spec("f", 1.0, pol_model="jones", pol_params=pol)
    e_in = (0.8, 0.6j)
    beams = [pt.beam_spec("fw", [pt.step_spec("f", 1)], "o", jones_in=e_in, jones_out=(1, 0)),
             pt.beam_spec("bw", [pt.step_spec("f", -1)], "o", jones_in=e_in, jones_out=(1, 0))]
    topo = pt.topo_spec([seg], beams, ["o"], [])
    pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1), far_spec(),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.fw.phi", "beam.bw.phi"], 2)]
    return pipe, U, np.array(e_in)


@pytest.mark.parametrize("compiled", [False, pytest.param(True, marks=pytest.mark.numba)])
def test_nonreciprocal_factor_not_transposed(compiled):
    pipe, U, e = nr_rig()
    _, st, _ = pt.run_and_record(pipe, 1e6, 2, compiled=compiled)
    N = R(0.6 * 50e-6 * 63)
    np.testing.assert_allclose(st["beam.fw.jones"], N @ U @ e, atol=1e-15)
    np.testing.assert_allclose(st["beam.bw.jones"], U.T @ N @ e, atol=1e-15)
    pipe0, _, _ = nr_rig(nr=False)  # no J_nr_slot: unchanged F6 behaviour
    _, st0, _ = pt.run_and_record(pipe0, 1e6, 2, compiled=compiled)
    np.testing.assert_allclose(st0["beam.bw.jones"], U.T @ e, atol=1e-15)


def test_topology_rejects_bad_nr_slot():
    from cranberry.path.topology import Topology, TopologyError

    seg = pt.seg_spec("f", 1.0, pol_model="jones", pol_params={"J": np.eye(2).tolist(), "J_nr_slot": 3})
    with pytest.raises(TopologyError, match="J_nr_slot"):
        Topology.from_spec(pt.topo_spec([seg], [pt.beam_spec("b", [pt.step_spec("f")], "o", jones_out=(1, 0))],
                                        ["o"], []))


@pytest.mark.numba
def test_pdl_faraday_compiled_equals_reference():
    pipe = [pdl_spec(), far_spec()]
    (b1, s1, p1), (b2, s2, p2) = rig(pipe), rig(pipe)
    run_reference(b1, s1, p1, 3)
    run_compiled(b2, s2, p2, 3)
    for j in ("J.pdl", "J.far"):
        np.testing.assert_allclose(s1[j], s2[j], rtol=0, atol=1e-15)


# ---------------------------------------------------------------------------------- pmd_mode_overlap


def pmd(L, tp, ref="psp", D=0.1e-12 / math.sqrt(1e3), on=True):
    return make_pmd_mode_overlap({"D_PMD_s_per_sqrt_m": D, "L_m": L, "tau_p_fwhm_s": tp, "timing_reference": ref,
                                  "segment": "smf", "applies_to": ["gbs", "bhd"], "source": "[T09 §1.7 AE]",
                                  "toggles": {"pol_pmd": on}}, name="pmd")


def gauss(n, dt, tp, shift=0.0):
    t = (np.arange(n) - n // 2) * dt
    sig = tp / (2 * math.sqrt(2 * math.log(2)))
    u = np.exp(-(t - shift) ** 2 / (4 * sig ** 2)).astype(complex)
    return t, u / math.sqrt(np.sum(abs(u) ** 2) * dt)


@pytest.mark.parametrize("ref", ["psp", "centered"])
def test_pmd_closed_form_matches_pulse_projection(ref):
    """Equal PSP split with DGD τ, reference mode at one PSP (psp) or at the mean delay (centered)."""
    b = pmd(630.0, 100e-15, ref)
    tau, dt, n = b.tau_dgd_s, 2e-15, 4096
    t, u = gauss(n, dt, 100e-15)
    V = np.array([[1, 1], [1, -1]]) / math.sqrt(2)  # PSPs ±45°; input x = equal split

    def U(om):
        d = np.zeros((len(om), 2, 2), complex)
        d[:, 0, 0], d[:, 1, 1] = np.exp(0.5j * tau * om), np.exp(-0.5j * tau * om)
        return V @ d @ V.T

    f = pulse_projection(u, dt, U, (1, 0), (1, 0))
    ref_mode = gauss(n, dt, 100e-15, tau / 2 if ref == "psp" else 0.0)[1]  # psp: u(T − τ/2), the delayed PSP
    eta = abs(np.sum(np.conj(ref_mode) * f) * dt) ** 2
    assert eta == pytest.approx(b.eta, rel=1e-6)


def test_pmd_t09_numbers_and_toggle():
    assert 1 - pmd(63.0, 100e-15).eta == pytest.approx(0.042, abs=5e-4)  # T09 §4.2
    assert 1 - pmd(630.0, 100e-15).eta == pytest.approx(0.32, abs=5e-3)
    assert 1 - pmd(630.0, 1e-12).eta < 0.005
    assert 1 - pmd(630.0, 100e-15, D=0.04e-12 / math.sqrt(1e3)).eta == pytest.approx(0.066, abs=1e-3)
    assert pmd(63.0, 100e-15).dphi_cw_pulse_rad == pytest.approx(1.1e-2, rel=0.01)  # T09: ≤ 1e−2
    off = pmd(63.0, 100e-15, on=False)
    assert off.eta == 1.0 and off.dphi_cw_pulse_rad == 0.0
    reg = LossRegistry(pmd(63.0, 100e-15).losses(), paths={"sq": {"segments": ["smf"]}})
    assert reg.eta("sq", "bhd") == pytest.approx(pmd(63.0, 100e-15).eta, rel=1e-14)
    assert LossRegistry(off.losses(), paths={"sq": {"segments": ["smf"]}}).eta("sq", "gbs") == 1.0
    with pytest.raises(ValueError, match="timing_reference"):
        pmd(63.0, 100e-15, ref="peak")
