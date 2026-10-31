"""WP-B7 POL blocks: segment_jones, polarizer_projection, per_fringe. Analytic U, SMF diffusion
⟨Θ²⟩ = 2·D_pol·t, toggles off ⇒ identity, compiled == reference (1e−12), loss registration, validation."""

import math

import numpy as np
import pytest

import cranberry.core.layout as layout_mod
import cranberry.path.testing as pt
from cranberry.blocks.pol import mu_to_x_dchi
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.offline.loss_registry import LossRegistry

ON = {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}
OFF = dict.fromkeys(ON, False)


@pytest.fixture(autouse=True)
def pol_group(monkeypatch):
    monkeypatch.setattr(layout_mod, "SLOT_GROUPS", layout_mod.SLOT_GROUPS | {"pol"})


def topo(n_smf=1):
    segs = [pt.seg_spec("a", 8.4, n_sub=3, drivers={"dT": {"subs": [0, 1]}}, pol_model="jones",
                        pol_params={"J_slot": "J.a"})]
    segs += [pt.seg_spec(f"b{i}", 5.0, n_sub=2, drivers={"dT": "per_sub"}, spatial="distributed", pol_model="jones",
                         pol_params={"J_slot": f"J.b{i}"}, kind="smf") for i in range(n_smf)]
    path = [pt.step_spec(s["id"]) for s in segs]
    coef = sum((pt.coef_spec(s["id"], "dT", "*", 100.0) for s in segs), [])
    return pt.topo_spec(segs, [pt.beam_spec("s", path, "det", jones_out=[1, 0])], ["det"], coef)


PM = {"model": "pm", "axis_drivers": ["dT"], "dpsi_b0_rad": 0.3, "K_T_rad_per_m_K": -2.05,
      "theta_in_rad": 0.05, "theta_out_rad": -0.03}


def smf(D=10.0):
    return {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0.6, 0.0, 0.8], "D_pol_rad2_per_s": D,
            "k_T_rad_per_m_K": 1e-2}


def sj_spec(t, toggles=ON, D=10.0, n_smf=1):
    segs = {"a": PM, **{f"b{i}": smf(D) for i in range(n_smf)}}
    return {"block": "segment_jones", "name": "sj", "params": {"topology": t, "segments": segs, "toggles": toggles}}


def rig(pipe, external, setup=None, fs=1e6, seed=5):
    lay = SlotLayout(fs)
    if setup:
        setup(lay)
    blk = build(spec_from_obj(pipe), lay, seed, external=external)
    return blk, State(lay), Params(lay)


def sj_rig(toggles=ON, D=10.0, n_smf=1, fs=1e6):
    t = topo(n_smf)
    ext = ["drv.a.dT", *(f"drv.b{i}.dT" for i in range(n_smf))]

    def setup(lay):
        lay.add("drv.a.dT", (), np.float64, "K", init=0.2)
        for i in range(n_smf):
            lay.add(f"drv.b{i}.dT", (2,), np.float64, "K", init=0.4)

    return rig([sj_spec(t, toggles, D, n_smf)], ext, setup, fs)


def R(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def su2(v):
    v = np.asarray(v, float)
    th = np.linalg.norm(v)
    n = v / th
    sig = n[0] * np.diag([1, -1]) + n[1] * np.array([[0, 1], [1, 0]]) + n[2] * np.array([[0, -1j], [1j, 0]])
    return math.cos(th / 2) * np.eye(2) - 1j * math.sin(th / 2) * sig


def test_pm_matrix_is_the_analytic_form():
    blk, s, p = sj_rig(D=0.0)
    run_reference(blk, s, p, 2)
    T = 0.2 * 2 / 3  # placement {subs: [0, 1]} of N_sub = 3: weight 1/3 each
    want = R(-0.03) @ np.diag([1, np.exp(-1j * (0.3 - 2.05 * 8.4 * T))]) @ R(0.05)
    np.testing.assert_allclose(s["J.a"], want, atol=1e-14)
    # SMF with D = 0: pure temperature rotation δ = k_T·L·mean(drv) about the Stokes axis
    np.testing.assert_allclose(s["J.b0"], su2(1e-2 * 5.0 * 0.4 * np.array([0.6, 0, 0.8])), atol=1e-14)
    p["sj.a.theta_in_rad"] = 0.0  # runtime parameter: no rebuild
    run_reference(blk, s, p, 1)
    np.testing.assert_allclose(s["J.a"], R(-0.03) @ np.diag([1, np.exp(-1j * (0.3 - 2.05 * 8.4 * T))]), atol=1e-14)


@pytest.mark.parametrize("off", ["pol_birefringence_drift", "pol_crosstalk", "pol_smf_drift"])
def test_each_toggle_removes_its_effect(off):
    blk, s, p = sj_rig({**ON, off: False}, D=0.0)
    run_reference(blk, s, p, 1)
    T = 0.2 * 2 / 3
    bire = 0.0 if off == "pol_birefringence_drift" else 0.3 - 2.05 * 8.4 * T
    ti, to = (0.0, 0.0) if off == "pol_crosstalk" else (0.05, -0.03)
    np.testing.assert_allclose(s["J.a"], R(to) @ np.diag([1, np.exp(-1j * bire)]) @ R(ti), atol=1e-14)
    if off == "pol_smf_drift":
        np.testing.assert_array_equal(s["J.b0"], np.eye(2))


def test_all_toggles_off_is_identity():
    blk, s, p = sj_rig(OFF)
    run_reference(blk, s, p, 10)
    for j in ("J.a", "J.b0"):
        np.testing.assert_array_equal(s[j], np.eye(2))
    assert next(iter_leaves(blk)).draws is None


def test_smf_diffusion_rate():
    """Per step the output state moves by Θ² with mean 2·D_pol·h (2D diffusion on the Poincaré sphere)."""
    D, fs, n = 1.0, 1e4, 20000
    blk, s, p = sj_rig(D=D, fs=fs)
    p["sj.b0.k_T_rad_per_m_K"] = 0.0
    run_reference(blk, s, p, 1)
    th2 = np.empty(n)
    e_prev = np.array(s["J.b0"])[:, 0]
    for k in range(n):
        run_reference(blk, s, p, 1)
        e = np.array(s["J.b0"])[:, 0]
        th2[k] = 4 * (1 - abs(np.vdot(e_prev, e)) ** 2)  # Θ² ≈ 4 sin²(Θ/2), Θ ≪ 1
        e_prev = e
    assert np.mean(th2) == pytest.approx(2 * D / fs, rel=0.05)  # 3σ ≈ 2 %
    U = np.array(s["J.b0"])
    np.testing.assert_allclose(U @ U.conj().T, np.eye(2), atol=1e-12)
    assert abs(np.linalg.det(U) - 1) < 1e-12


@pytest.mark.numba
def test_segment_jones_compiled_equals_reference():
    (b1, s1, p1), (b2, s2, p2) = sj_rig(n_smf=2), sj_rig(n_smf=2)
    run_reference(b1, s1, p1, 300)
    run_compiled(b2, s2, p2, 300)
    for j in ("J.a", "J.b0", "J.b1"):
        np.testing.assert_allclose(s1[j], s2[j], rtol=0, atol=1e-12)
    assert not np.allclose(s1["J.b0"], s1["J.b1"])  # independent draws per segment


def test_segment_jones_validation():
    t = topo()
    with pytest.raises(MissingParameterError):
        build(spec_from_obj([{"block": "segment_jones", "name": "x", "params": {
            "topology": t, "segments": {"a": {k: v for k, v in PM.items() if k != "K_T_rad_per_m_K"}},
            "toggles": ON}}]), SlotLayout(1e6), 0)
    with pytest.raises(MissingParameterError):
        build(spec_from_obj([{"block": "segment_jones", "name": "x", "params": {
            "topology": t, "segments": {"a": PM}, "toggles": {"pol_crosstalk": True}}}]), SlotLayout(1e6), 0)
    t2 = pt.topo_spec([pt.seg_spec("c", 1.0, drivers={"dT": "uniform"})],
                      [pt.beam_spec("s", [pt.step_spec("c")], "det")], ["det"], pt.coef_spec("c", "dT", "*", 1.0))
    with pytest.raises(ValueError, match="J_slot"):
        build(spec_from_obj([{"block": "segment_jones", "name": "x", "params": {
            "topology": t2, "segments": {"c": PM}, "toggles": ON}}]), SlotLayout(1e6), 0)
    with pytest.raises(ValueError, match="model"):
        build(spec_from_obj([{"block": "segment_jones", "name": "x", "params": {
            "topology": t, "segments": {"a": {**PM, "model": "bowtie"}}, "toggles": ON}}]), SlotLayout(1e6), 0)


# ------------------------------------------------------------------------------- polarizer_projection


def polz_spec(drift=True, ext=True, er=28.0):
    return {"block": "polarizer_projection", "name": "polz", "params": {
        "in": "beam.mon.jones", "out_jones": "beam.polz.jones", "out_dp": "v.polz.dp", "er_db": er,
        "il_db": 0.3, "segment": "polz", "applies_to": ["probe"], "source": "[D51-8 Opneti]",
        "toggles": {"pol_drift_before_polarizer": drift, "pol_extinction": ext}}}


def polz_rig(**kw):
    def setup(lay):
        lay.add("beam.mon.jones", (2,), np.complex128, "1")

    return rig([polz_spec(**kw)], ["beam.mon.jones"], setup)


def test_polarizer_projection():
    blk, s, p = polz_rig()
    th = 0.3  # Poincaré angle Θ
    s["beam.mon.jones"] = [math.cos(th / 2), 1j * math.sin(th / 2)]
    run_reference(blk, s, p, 1)
    q = 10 ** (-28 / 10)
    assert float(s["v.polz.dp"]) == pytest.approx(math.cos(th / 2) ** 2 + q * math.sin(th / 2) ** 2 - 1, rel=1e-13)
    e = np.array(s["beam.polz.jones"])
    assert abs(np.linalg.norm(e) - 1) < 1e-15
    assert abs(e[1]) ** 2 == pytest.approx(q * math.sin(th / 2) ** 2 / (math.cos(th / 2) ** 2 + q * math.sin(th / 2) ** 2))
    reg = LossRegistry(next(iter_leaves(blk)).losses(), paths={"pr": {"segments": ["polz"]}})
    assert reg.eta("pr", "probe") == pytest.approx(10 ** (-0.03))


def test_polarizer_toggles_off_identity():
    blk, s, p = polz_rig(drift=False, ext=False)
    assert next(iter_leaves(blk)).reads == ()
    s["beam.mon.jones"] = [0.3, 0.9]
    run_reference(blk, s, p, 1)
    assert float(s["v.polz.dp"]) == 0.0
    np.testing.assert_array_equal(s["beam.polz.jones"], [1, 0])


@pytest.mark.numba
def test_polarizer_compiled_equals_reference():
    (b1, s1, p1), (b2, s2, p2) = polz_rig(), polz_rig()
    for s in (s1, s2):
        s["beam.mon.jones"] = [0.8, 0.6j]
    run_reference(b1, s1, p1, 3)
    run_compiled(b2, s2, p2, 3)
    assert abs(float(s1["v.polz.dp"]) - float(s2["v.polz.dp"])) < 1e-15
    np.testing.assert_allclose(s1["beam.polz.jones"], s2["beam.polz.jones"], atol=1e-15)


# ------------------------------------------------------------------------------------------ per_fringe


@pytest.mark.parametrize("mu", [0.99 + 0.01j, 0.5 - 0.3j, 1 - 2e-4 + 1e-6j, -0.9 + 0j, 0.2j])
def test_mu_inversion_is_exact(mu):
    x, d = mu_to_x_dchi(mu)
    assert 0 < x <= 1
    assert abs((1 - x) + x * np.exp(1j * d) - mu) < 1e-12
    assert mu_to_x_dchi(1.0) == (0.0, 0.0)


def pf_rig(on=True):
    def setup(lay):
        lay.add("beam.in.jones", (2,), np.complex128, "1", init=(0.995, 0.0998j))
        lay.add("J.p1", (2, 2), np.complex128, "1", init=tuple((R(0.02) @ np.diag([1, np.exp(-0.7j)])).ravel()))
        lay.add("J.p2", (2, 2), np.complex128, "1", init=tuple(np.diag([1, np.exp(1.9j)]).ravel()))
        lay.add("J.l1", (2, 2), np.complex128, "1", init=tuple((np.diag([1, np.exp(0.4j)]) @ R(-0.01)).ravel()))

    pf = {"block": "per_fringe", "name": "refD", "params": {
        "in_jones": "beam.in.jones", "j_pr": ["J.p1", "J.p2"], "j_lo": ["J.l1"], "beam_jones_in": [1, 0],
        "beam_jones_out": [1, 0], "out_x": "pol.refD.x", "out_dchi": "pol.refD.dchi",
        "toggles": {"per_crosscoupling": on}}}
    return rig([pf], ["beam.in.jones", "J.p1", "J.p2", "J.l1"], setup)


def test_per_fringe_matches_jones_overlap():
    blk, s, p = pf_rig()
    run_reference(blk, s, p, 1)
    e = np.array([0.995, 0.0998j])
    e = e / np.linalg.norm(e)
    Up = np.diag([1, np.exp(1.9j)]) @ R(0.02) @ np.diag([1, np.exp(-0.7j)])
    Ul = np.diag([1, np.exp(0.4j)]) @ R(-0.01)
    mu_full = np.vdot(Up @ e, Ul @ e)
    tp, tl = (Up @ [1, 0])[0], (Ul @ [1, 0])[0]
    mu = mu_full * tp * np.conj(tl) / abs(tp * tl)
    x, d = float(s["pol.refD.x"]), float(s["pol.refD.dchi"])
    assert abs((1 - x) + x * np.exp(1j * d) - mu) < 1e-14


def test_per_fringe_toggle_off_zero():
    blk, s, p = pf_rig(on=False)
    assert next(iter_leaves(blk)).reads == ()
    run_reference(blk, s, p, 1)
    assert float(s["pol.refD.x"]) == 0.0 and float(s["pol.refD.dchi"]) == 0.0


@pytest.mark.numba
def test_per_fringe_compiled_equals_reference():
    (b1, s1, p1), (b2, s2, p2) = pf_rig(), pf_rig()
    run_reference(b1, s1, p1, 2)
    run_compiled(b2, s2, p2, 2)
    for k in ("pol.refD.x", "pol.refD.dchi"):
        assert abs(float(s1[k]) - float(s2[k])) < 1e-12
