"""VT-38 (M1 part, WP-B7): polarization (T09 §5 tests 1, 3, 4, 5, 6, 8; G04; D51-2, D51-8).

* T09 test 1: PM differential drift. L = 63 m, L_b = 4.2 mm, k_T = −1.37e−3 /K, λ = 1545.32 nm, ΔT =
  +1 mK. With an orthogonal counter-propagating probe the locked combination moves by
  K_T·L·ΔT = −0.1291 rad (rel. 1e−6). With a same-polarization probe it moves by 0 (1e−12).
* T09 test 3: a 3° connector misalignment gives PER 25.61 dB (0.01 dB).
* T09 test 4: distributed crosstalk h·L = 1e−3 gives tanh(1e−3) → −30.00 dB (0.01 dB).
* T09 test 5: Ref D lock point with ε = 0.01, δ_b swept over 0…2π by the 915B coil temperature
  (K_T·L = −17.1 rad/K). The ``per_fringe`` slots give shift −arg μ = −arg(1 + ε e^{iδ_b}) (1e−9; in
  φ_D = φ_LO − φ_pr, which is T09's −arg(1 + ε e^{−iδ_b}) in φ_pr − φ_lo), maximum
  arcsin(0.01) = 10.0002 mrad, and 0 for ε = 0. Re(μ e^{iφ_D}) equals the fringe term summed over
  both polarizations of the full fields (1e−12, with crosstalk inside the arms).
* T09 test 6: polarization overlap η = cos²(Θ/2) (Θ = π/2 → 0.5, Θ = 0.2 → 0.990033; 1e−6) through
  a segment_jones rotation and the beam analyzer. A polarizer (ER = ∞) at Θ = 0.1 transmits
  cos²(0.05) = 0.997502.
* T09 test 8: all polarization toggles off ⇒ U = 1, and the prompt §7-11 shared-segment
  cancellation (VT-22) passes unchanged, bit-equal (1e−12) to the pol-off topology. It also passes
  with polarization on and identical beam states.

M4 rows (WP-B7):

* T09 test 2: 100 random SU(2) chains (100 ``segment_jones`` SMF segments, Haar U0 + strong diffusion,
  one new chain per sample): same-polarization backward probe phase = squeezed phase (1e−12),
  |t_s| = |t_p|; orthogonal probe error = 2·arg(a), a = xᵀUx (1e−12); PM diagonal ⇒ Σ Δψ_b.
  Faraday (``faraday_nonreciprocal``, ρ = 0.6·50e−6·63 = 1.89e−3 rad) in the chain breaks it by exactly
  arg(xᵀBNAx) − arg(xᵀAᵀNBᵀx) (median ≈ 2.8e−3 rad); circular worst case −2ρ = −3.78e−3 rad; off ⇒ 0.
* T09 test 7: CW probe phase amplitude κ1κ2 = 3.16e−3 rad (1 %); 100 fs pulse ≤ 1e−4·κ1κ2 and a
  satellite of power κ1²κ2² = 1e−5 at −1.36 ps; 1 ps pulse κ1κ2·0.277 (T09's "≤ 1e−4 at 1 ps" does not
  hold for a Gaussian pulse).
* PMD (``pmd_mode_overlap``, worst case, 0.1 ps/√km): 100 fs → 4.23 % at 63 m, 32.25 % at 630 m;
  1 ps → 0.44 % at 630 m; re-timed reference 19.6 % at 630 m.
* PDL (``pdl_element``): power range T_mean·(1 ± (T_max − T_min)/(T_max + T_min)) (T09 §1.4).
* pol_pdl / pol_faraday / pol_smf_drift (Haar U0, OU) off ⇒ every matrix is 1.
"""

import math

import numpy as np
import pytest

import cranberry.core.layout as layout_mod
import cranberry.path.testing as pt
from cranberry.blocks.loss.per_budget import per_db_from_angle
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State

OFF = {"pol_birefringence_drift": False, "pol_crosstalk": False, "pol_smf_drift": False}
ON = {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}


@pytest.fixture(autouse=True)
def pol_group(monkeypatch):
    """Slot group ``pol`` (proposed one-line addition to core.layout.SLOT_GROUPS, WP-B7 report)."""
    monkeypatch.setattr(layout_mod, "SLOT_GROUPS", layout_mod.SLOT_GROUPS | {"pol"})


def jseg(id, L, *, kind="pm", n_sub=1, drivers=None):
    return pt.seg_spec(id, L, n_sub=n_sub, drivers=drivers or {"dT": "uniform"}, pol_model="jones",
                       pol_params={"J_slot": f"J.{id}"}, kind=kind)


def pm(dpsi=0.0, K_T=-2.05, ti=0.0, to=0.0):
    return {"model": "pm", "axis_drivers": ["dT"], "dpsi_b0_rad": dpsi, "K_T_rad_per_m_K": K_T,
            "theta_in_rad": ti, "theta_out_rad": to}


def sj(topo, segments, toggles=ON, name="sj"):
    return {"block": "segment_jones", "name": name, "params": {"topology": topo, "segments": segments,
                                                               "toggles": dict(toggles)}}


def table(slot, values, unit="K"):
    v = np.asarray(values, dtype=float)
    return pt.pathtest("pathtest_table", f"tab.{slot}", slot, v.shape[1:], unit, n_rows=len(v)), {f"{slot}_table": v}


# ------------------------------------------------------------------------------------------- test 1


def test_t09_1_pm_differential_drift():
    L, lam, L_b, k_T = 63.0, 1545.32e-9, 4.2e-3, -1.37e-3
    K_T = k_T * 2 * math.pi / L_b  # T09 §1.2: K_T = k_T·2π/L_b
    seg = jseg("f", L, n_sub=4)
    back = dict(direction=-1)
    beams = [pt.beam_spec("sq", [pt.step_spec("f")], "o", jones_in=(1, 0), jones_out=(1, 0)),
             pt.beam_spec("orth", [pt.step_spec("f", **back)], "o", jones_in=(0, 1), jones_out=(0, 1)),
             pt.beam_spec("same", [pt.step_spec("f", **back)], "o", jones_in=(1, 0), jones_out=(1, 0))]
    topo = pt.topo_spec([seg], beams, ["o"], pt.coef_spec("f", "dT", "*", 44.0 * L))
    drv, tabs = table("drv.f.dT", [0.0] * 5 + [1e-3] * 15)
    src, tab2 = table("laser.las.phi", [0.0] * 20, "rad")
    pipe = [drv, src, sj(topo, {"f": pm(dpsi=0.7, K_T=K_T)}), *pt.path_pipeline(topo, {"las": "laser.las.phi"}),
            pt.record(["beam.sq.phi", "beam.orth.phi", "beam.same.phi"], 20)]
    rec, _, _ = pt.run_and_record(pipe, 1e6, 20, {**tabs, **tab2})
    d_orth = (rec["beam.sq.phi"] - rec["beam.orth.phi"])
    d_same = (rec["beam.sq.phi"] - rec["beam.same.phi"])
    want = K_T * L * 1e-3
    assert want == pytest.approx(-0.1291, abs=1e-4)
    assert d_orth[-1] - d_orth[2] == pytest.approx(want, rel=1e-6)
    assert abs(d_same[-1] - d_same[2]) < 1e-12 and np.max(np.abs(d_same)) < 1e-12


# ----------------------------------------------------------------------------------- tests 3 and 4


def per_budget(elements):
    from cranberry.blocks.loss.per_budget import make_per_budget

    return make_per_budget({"elements": elements, "loss": None, "toggles": {"pol_crosstalk": True}}, name="pb")


def test_t09_3_connector_misalignment():
    assert per_db_from_angle(math.radians(3.0)) == pytest.approx(25.61, abs=0.01)
    b = per_budget([{"id": "fc", "type": "connector", "theta_deg": 3.0, "source": "[AE]"}])
    assert b.per_total_db == pytest.approx(25.61, abs=0.01)


def test_t09_4_distributed_crosstalk():
    b = per_budget([{"id": "fib", "type": "fiber", "h_per_m": 1e-5, "L_m": 100.0, "source": "[Corning]"}])
    assert b.x_total == pytest.approx(math.tanh(1e-3), rel=1e-14)
    assert 10 * math.log10(b.x_total) == pytest.approx(-30.00, abs=0.01)


# ------------------------------------------------------------------------------------------- test 5


def refd_rig(eps, dT, *, ti=0.0, to=0.0, emission_time=True):
    """Ref D arms: pr = 915B coil (38 × 0.22 m, jones PM), lo = 0.5 m PM; per_fringe → pol.refD.*."""
    x = eps / (1 + eps)
    segs = [jseg("a1", 38 * 0.22), jseg("a2", 0.5)]
    beams = [pt.beam_spec("pr", [pt.step_spec("a1")], "pd", jones_out=(1, 0)),
             pt.beam_spec("lo", [pt.step_spec("a2")], "pd", jones_out=(1, 0))]
    coef = pt.coef_spec("a1", "dT", "*", 44.0 * 38 * 0.22) + pt.coef_spec("a2", "dT", "*", 22.0)
    topo = pt.topo_spec(segs, beams, ["pd"], coef)
    n = len(dT)
    d1, t1 = table("drv.a1.dT", dT)
    d2, t2 = table("drv.a2.dT", np.zeros(n))
    d3, t3 = table("laser.las.phi", np.zeros(n), "rad")
    pipe = [d1, d2, d3,
            sj(topo, {"a1": pm(dpsi=0.0, ti=ti, to=to), "a2": pm(dpsi=0.0, ti=-0.5 * ti, to=0.3 * to)}),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}, emission_time=emission_time),
            {"block": "per_fringe", "name": "refD", "params": {
                "in_jones": "beam.pol_in.jones", "j_pr": ["J.a1"], "j_lo": ["J.a2"], "beam_jones_in": [1, 0],
                "beam_jones_out": [1, 0], "out_x": "pol.refD.x", "out_dchi": "pol.refD.dchi",
                "toggles": {"per_crosscoupling": True}}},
            pt.record(["pol.refD.x", "pol.refD.dchi", "beam.pr.phi", "beam.lo.phi", "drv.a1.dT", "psi.a1.pr"], n)]

    def setup(lay):
        lay.add("beam.pol_in.jones", (2,), np.complex128, "1", init=(math.sqrt(1 - x), math.sqrt(x)))

    rec, _, _ = pt.run_and_record(pipe, 1e6, n, {**t1, **t2, **t3}, setup=setup, external=["beam.pol_in.jones"])
    return rec, x


def mu_of(rec):
    x, d = rec["pol.refD.x"], rec["pol.refD.dchi"]
    return (1 - x) + x * np.exp(1j * d)


def test_t09_5_refD_lock_point():
    eps, KL = 0.01, -2.05 * 38 * 0.22
    assert KL == pytest.approx(-17.138, abs=1e-3)  # T09 §4.1: −17 rad/K
    d_star = math.acos(-eps)  # extremum of arg(1 + ε e^{iδ})
    deltas = np.sort(np.concatenate([np.linspace(0, 2 * math.pi, 1001), [d_star, 2 * math.pi - d_star]]))
    rec, x = refd_rig(eps, deltas / KL)  # δ_b = Δψ_b,a1 − Δψ_b,a2 = K_T·L·ΔT
    moving = np.abs(np.exp(1j * deltas) - 1) > 1e-9  # μ = 1 at δ_b = 0 mod 2π: x is then reported as 0
    np.testing.assert_allclose(rec["pol.refD.x"][moving], x, rtol=1e-9)
    shift = -np.angle(mu_of(rec))
    np.testing.assert_allclose(shift, -np.angle(1 + eps * np.exp(1j * deltas)), rtol=0, atol=1e-9)
    assert np.max(np.abs(shift)) == pytest.approx(math.asin(eps), abs=1e-12)
    assert math.asin(eps) == pytest.approx(10.0002e-3, abs=1e-7)
    rec0, _ = refd_rig(0.0, deltas / KL)
    assert np.max(np.abs(rec0["pol.refD.x"])) == 0.0 and np.max(np.abs(np.angle(mu_of(rec0)))) == 0.0


def test_t09_5_fringe_equals_full_field_sum():
    """Re(μ e^{iφ_D}) (dc_fringe's eq. (3) term with the per_fringe slots) equals the PD fringe term
    summed over both polarizations of the full fields, here with crosstalk inside the arms."""
    eps, ti, to = 0.02, 0.04, -0.03
    dT = np.linspace(0, 0.5, 200)
    rec, x = refd_rig(eps, dT, ti=ti, to=to, emission_time=False)

    def R(t):
        return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])

    e = np.array([math.sqrt(1 - x), math.sqrt(x)])
    term = []
    for T, psi in zip(rec["drv.a1.dT"], rec["psi.a1.pr"]):
        U1 = R(to) @ np.diag([1, np.exp(-1j * (-2.05 * 38 * 0.22 * T))]) @ R(ti)
        U2 = R(0.3 * to) @ np.diag([1, 1]) @ R(-0.5 * ti)
        E1 = np.exp(1j * psi) * U1 @ e  # pr; ψ as segment_phase evaluates it (subsegment timing)
        E2 = U2 @ e  # lo (ψ = 0, its dT is 0)
        term.append(np.real(np.vdot(E1, E2)) / (np.linalg.norm(E1) * np.linalg.norm(E2)))
    phid = rec["beam.lo.phi"] - rec["beam.pr.phi"]
    np.testing.assert_allclose(np.real(mu_of(rec) * np.exp(1j * phid)), term, rtol=0, atol=1e-12)


# ------------------------------------------------------------------------------------------- test 6


@pytest.mark.parametrize("theta,eta", [(math.pi / 2, 0.5), (0.2, 0.990033)])
def test_t09_6_overlap(theta, eta):
    L, kT = 2.0, 0.5
    seg = jseg("s", L, kind="smf")
    topo = pt.topo_spec([seg], [pt.beam_spec("b", [pt.step_spec("s")], "o", jones_out=(1, 0))], ["o"],
                        pt.coef_spec("s", "dT", "*", 0.0))
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 0, 1], "u0": "identity", "drift": "diffusion", "D_pol_rad2_per_s": 0.0,
           "k_T_rad_per_m_K": kT}
    d1, t1 = table("drv.s.dT", [theta / (kT * L)] * 3)
    d2, t2 = table("laser.las.phi", [0.0] * 3, "rad")
    pipe = [d1, d2, sj(topo, {"s": smf}), *pt.path_pipeline(topo, {"las": "laser.las.phi"}),
            pt.record(["beam.b.power"], 3)]
    rec, _, _ = pt.run_and_record(pipe, 1e6, 3, {**t1, **t2})
    assert rec["beam.b.power"][-1] / 1e-3 == pytest.approx(math.cos(theta / 2) ** 2, abs=1e-12)
    assert math.cos(theta / 2) ** 2 == pytest.approx(eta, abs=1e-6)


def test_t09_6_polarizer_power():
    lay = SlotLayout(1e6)
    lay.add("beam.mon.jones", (2,), np.complex128, "1")
    pol = {"block": "polarizer_projection", "name": "polz", "params": {
        "in": "beam.mon.jones", "out_jones": "beam.polz.jones", "out_dp": "v.polz.dp", "er_db": 30.0,
        "il_db": 0.3, "segment": "polz", "applies_to": ["probe"], "source": "[D51-8]",
        "toggles": {"pol_drift_before_polarizer": True, "pol_extinction": False}}}
    blk = build(spec_from_obj([pol]), lay, 0, external=["beam.mon.jones"])
    s, p = State(lay), Params(lay)
    s["beam.mon.jones"] = [math.cos(0.05), math.sin(0.05)]  # Θ = 0.1 on the Poincaré sphere
    run_reference(blk, s, p, 1)
    assert 1 + float(s["v.polz.dp"]) == pytest.approx(0.997502, abs=1e-6)
    assert 1 + float(s["v.polz.dp"]) == pytest.approx(math.cos(0.05) ** 2, rel=1e-14)


# ------------------------------------------------------------------------------------------- test 8


def shared_topology(pol: bool):
    js = (lambda id, L, kind="pm", **kw: jseg(id, L, kind=kind, **kw)) if pol else \
        (lambda id, L, kind="pm", **kw: pt.seg_spec(id, L, kind=kind, drivers=kw.get("drivers") or {"dT": "uniform"},
                                                     n_sub=kw.get("n_sub", 1)))
    segs = [js("s1", 7.3, n_sub=3), js("s2", 2.1, kind="smf"),
            pt.seg_spec("t_sq", 1.7, n_sub=2, drivers={"dT": "uniform"}),
            pt.seg_spec("t_pr", 1.7, n_sub=2, drivers={"dT": "uniform"})]
    jout = {"jones_out": (1, 0)} if pol else {}
    path = [pt.step_spec("s1"), pt.step_spec("s2")]
    beams = [pt.beam_spec("sq", path + [pt.step_spec("t_sq")], "bhd", **jout),
             pt.beam_spec("pr", path + [pt.step_spec("t_pr")], "bhd", **jout)]
    coef = (pt.coef_spec("s1", "dT", "*", 320.0) + pt.coef_spec("s2", "dT", "*", 92.0)
            + pt.coef_spec("t_sq", "dT", "*", 75.0) + pt.coef_spec("t_pr", "dT", "*", 75.0))
    locks = [{"id": "co", "combo": {"sq": 1.0, "pr": -1.0}, "setpoint_rad": 0.3, "observer": "bhd"}]
    return pt.topo_spec(segs, beams, ["bhd"], coef, locks)


def run_shared(pol: bool, toggles, shared_sigma=1.0, n=300):
    topo = shared_topology(pol)
    pipe = [pt.pathtest("pathtest_white", f"n.{s['id']}", f"drv.{s['id']}.dT", (), "K",
                        sigma=1e-3 * (shared_sigma if s["id"] in ("s1", "s2") else 1.0)) for s in topo["segments"]]
    pipe.append(pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=1.0))
    if pol:
        smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "u0": "haar", "drift": "diffusion", "D_pol_rad2_per_s": 50.0,
               "k_T_rad_per_m_K": 3.0}
        pipe.append(sj(topo, {"s1": pm(dpsi=1.1, ti=0.05, to=-0.04), "s2": smf}, toggles))
    pipe += [*pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True), pt.record(["lock.co.err", "beam.sq.phi"], n)]
    rec, _, _ = pt.run_and_record(pipe, 1e8, n, seed=38)
    return rec


@pytest.mark.parametrize("toggles", [OFF, ON], ids=["all_off", "all_on"])
def test_t09_8_toggles_off_shared_segment_cancellation(toggles):
    on = run_shared(True, toggles)
    off = run_shared(True, toggles, shared_sigma=0.0)
    scale = np.max(np.abs(on["beam.sq.phi"]))
    assert np.max(np.abs(on["lock.co.err"] - off["lock.co.err"])) < 1e-12 * scale  # prompt §7-11
    if toggles is OFF:  # U = 1: identical to the pol_model-off topology
        plain = run_shared(False, toggles)
        np.testing.assert_allclose(on["lock.co.err"], plain["lock.co.err"], rtol=0, atol=1e-12 * scale)
        np.testing.assert_allclose(on["beam.sq.phi"], plain["beam.sq.phi"], rtol=0, atol=1e-12 * scale)


def test_t09_8_all_off_gives_identity():
    topo = shared_topology(True)
    lay = SlotLayout(1e6)
    for s in ("s1", "s2"):
        lay.add(f"drv.{s}.dT", (), np.float64, "K", init=0.37)
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "u0": "haar", "drift": "diffusion", "D_pol_rad2_per_s": 50.0,
           "k_T_rad_per_m_K": 3.0}
    blk = build(spec_from_obj([sj(topo, {"s1": pm(dpsi=1.1, ti=0.05, to=-0.04), "s2": smf}, OFF)]), lay, 1,
                external=["drv.s1.dT", "drv.s2.dT"])
    s, p = State(lay), Params(lay)
    run_reference(blk, s, p, 5)
    for j in ("J.s1", "J.s2"):
        np.testing.assert_array_equal(s[j], np.eye(2))
    assert not any(leaf.draws for leaf in iter_leaves(blk))


# ================================================================================= M4 rows (WP-B7)
# ------------------------------------------------------------------------------------------- test 2


def reciprocity_rig(n_seg=100, faraday=False, n_steps=100, pm=False, mid=None):
    """n_seg SMF segments (Haar U0 + strong diffusion: every sample is a new random SU(2) chain) or PM
    segments (diagonal U); optional Faraday segment `far` after segment `mid`. Beams: sq forward x → x,
    probe `same` backward x → x, probe `orth` backward y → y."""
    ids = [f"s{i}" for i in range(n_seg)]
    segs = [pt.seg_spec(s, 0.63, kind="pm" if pm else "smf", pol_model="jones", pol_params={"J_slot": f"J.{s}"})
            for s in ids]
    order = list(ids)
    if mid is not None:
        segs.append(pt.seg_spec("far", 0.0, pol_model="jones", pol_params={"J": np.eye(2).tolist(), "J_nr_slot": "J.far"}))
        order.insert(mid + 1, "far")
    fw = [pt.step_spec(s) for s in order]
    bw = [pt.step_spec(s, -1) for s in reversed(order)]
    beams = [pt.beam_spec("sq", fw, "o", jones_in=(1, 0), jones_out=(1, 0)),
             pt.beam_spec("same", bw, "o", jones_in=(1, 0), jones_out=(1, 0)),
             pt.beam_spec("orth", bw, "o", jones_in=(0, 1), jones_out=(0, 1))]
    topo = pt.topo_spec(segs, beams, ["o"], [])
    rng = np.random.default_rng(2)
    dpsi = rng.uniform(0, 2 * math.pi, n_seg)
    cfg = ({s: {"model": "pm", "axis_drivers": [], "dpsi_b0_rad": float(d), "K_T_rad_per_m_K": -2.05,
                "theta_in_rad": 0.0, "theta_out_rad": 0.0} for s, d in zip(ids, dpsi)} if pm else
           {s: {"model": "smf", "axis_drivers": [], "axis_stokes": [1, 0, 0], "u0": "haar", "drift": "diffusion",
                "D_pol_rad2_per_s": 1e7, "k_T_rad_per_m_K": 0.0} for s in ids})
    pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1), sj(topo, cfg)]
    if mid is not None:
        pipe.append(far_spec(faraday))
    pipe += [*pt.path_pipeline(topo, {"las": "laser.las.phi"}),
             pt.record([f"beam.{b}.{f}" for b in ("sq", "same", "orth") for f in ("phi", "power")], n_steps)]
    rec, st, _ = pt.run_and_record(pipe, 1e6, n_steps, seed=9)
    return rec, st, ids, dpsi


def far_spec(on=True):
    return {"block": "faraday_nonreciprocal", "name": "far", "params": {
        "out": "J.far", "verdet_rad_per_T_m": 0.6, "B_par_T": 50e-6, "L_m": 63.0, "toggles": {"pol_faraday": on}}}


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def chain(st, ids):
    U = np.eye(2, dtype=complex)
    for s in ids:
        U = np.array(st[f"J.{s}"]) @ U
    return U


def R(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def test_t09_2_reciprocity_random_chains():
    """100 random SU(2) chains of 100 segment_jones segments: same-polarization probe phase = squeezed
    phase (1e−12), |t_s| = |t_p|; orthogonal probe error = 2·arg(a), a = xᵀUx (1e−12)."""
    rec, st, ids, _ = reciprocity_rig()
    sq, same, orth = rec["beam.sq.phi"], rec["beam.same.phi"], rec["beam.orth.phi"]
    assert np.std(sq) > 0.5  # the chain really changes from sample to sample
    assert np.max(np.abs(wrap(sq - same))) < 1e-12
    np.testing.assert_allclose(rec["beam.sq.power"], rec["beam.same.power"], rtol=1e-12)
    np.testing.assert_allclose(wrap(sq - orth), wrap(2 * sq), rtol=0, atol=1e-12)  # t_s = a
    a = chain(st, ids)[0, 0]  # last sample, from the J slots directly
    assert abs(wrap(sq[-1] - np.angle(a))) < 1e-12
    assert abs(wrap(sq[-1] - orth[-1] - 2 * np.angle(a))) < 1e-12


def test_t09_2_pm_diagonal_error_is_dpsi_b():
    rec, _, _, dpsi = reciprocity_rig(pm=True, n_steps=3)
    assert abs(wrap(rec["beam.sq.phi"][-1] - rec["beam.orth.phi"][-1] - np.sum(dpsi))) < 1e-11
    assert abs(wrap(rec["beam.sq.phi"][-1] - rec["beam.same.phi"][-1])) < 1e-12


def test_t09_2_faraday_breaks_reciprocity():
    """Faraday segment (ρ = V_F·B·L = 1.89e−3 rad, T09 §1.8) in the middle of a random chain: the
    same-polarization error equals arg(xᵀB N A x) − arg(xᵀAᵀN Bᵀx) (1e−12), of order ρ; off ⇒ 0."""
    rho = 0.6 * 50e-6 * 63
    rec, st, ids, _ = reciprocity_rig(n_seg=20, faraday=True, mid=9)
    err = wrap(rec["beam.sq.phi"] - rec["beam.same.phi"])
    A, B, N = chain(st, ids[:10]), chain(st, ids[10:]), R(rho)
    want = np.angle((B @ N @ A)[0, 0]) - np.angle((A.T @ N @ B.T)[0, 0])
    assert abs(wrap(err[-1] - want)) < 1e-12
    assert 0.1 * rho < np.median(np.abs(err)) < 10 * rho  # ~2e−3 rad scale
    off, _, _, _ = reciprocity_rig(n_seg=20, faraday=False, mid=9)
    assert np.max(np.abs(wrap(off["beam.sq.phi"] - off["beam.same.phi"]))) < 1e-12


def test_t09_2_faraday_circular_worst_case():
    """x → QWP → σ+ → Faraday → QWP⁻¹ → x: forward e^{−iρ}, backward e^{+iρ}, error −2ρ = −3.78e−3 rad.
    The same rotation as a reciprocal element (static J) gives 0."""
    rho = 0.6 * 50e-6 * 63
    q1 = np.array([[1, 1j], [1j, 1]]) / math.sqrt(2)
    q2 = q1.conj().T

    def run(nr: bool):
        mid = {"J": np.eye(2).tolist(), "J_nr_slot": "J.far"} if nr else {"J": R(rho).tolist()}
        segs = [pt.seg_spec("q1", 0.0, pol_model="jones", pol_params={"J": q1.tolist()}),
                pt.seg_spec("far", 0.0, pol_model="jones", pol_params=mid),
                pt.seg_spec("q2", 0.0, pol_model="jones", pol_params={"J": q2.tolist()})]
        fw = [pt.step_spec(s) for s in ("q1", "far", "q2")]
        bw = [pt.step_spec(s, -1) for s in ("q2", "far", "q1")]
        topo = pt.topo_spec(segs, [pt.beam_spec("sq", fw, "o", jones_out=(1, 0)),
                                   pt.beam_spec("same", bw, "o", jones_out=(1, 0))], ["o"], [])
        pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1), far_spec(),
                *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "beam.same.phi"], 2)]
        rec, _, _ = pt.run_and_record(pipe, 1e6, 2)
        return rec["beam.sq.phi"][-1], rec["beam.same.phi"][-1]

    s, p = run(True)
    assert s == pytest.approx(-rho, abs=1e-15) and p == pytest.approx(rho, abs=1e-15)
    assert s - p == pytest.approx(-3.78e-3, abs=1e-5)
    s, p = run(False)
    assert abs(s - p) < 1e-15


# ------------------------------------------------------------------------------------------- test 7


def test_t09_7_cw_probe_vs_pulse():
    """Two PER 25 dB joints (κ = sin θ, tan²θ = 10^−2.5) ℓ = 1 m apart, DGD 1.36 ps/m. CW probe (same
    polarization, time loop): phase deviation amplitude κ1κ2 = 3.16e−3 rad (1 %) over a Δψ_b sweep.
    Pulse (``pulse_projection`` of the same J with U(Ω) = R(θo)·diag(1, e^{−iτΩ})·R(−θo)·J): the main
    pulse deviation is κ1κ2·O(τ), O = exp(−τ²/(8σ_t²)); 100 fs → ≤ 1e−4·κ1κ2, and a satellite of power
    κ1²κ2² = 1e−5 at −1.36 ps (it travelled on the fast axis, so it arrives early).
    T09 test 7 states ≤ 1e−4 also for τ_p = 1 ps; the Gaussian overlap gives O = 0.277 there (asserted
    below; 1e−4 needs τ_p ≤ 0.37 ps)."""
    from cranberry.blocks.loss.per_budget import angle_from_per_db
    from cranberry.blocks.pol.pmd_mode_overlap import pulse_projection

    th, tau, KL = angle_from_per_db(25.0), 1.36e-12, -2.05
    k2 = math.sin(th) ** 2
    assert k2 == pytest.approx(3.16e-3, rel=0.01)
    topo = pt.topo_spec([jseg("f", 1.0)], [pt.beam_spec("sq", [pt.step_spec("f")], "o", jones_out=(1, 0))], ["o"],
                        pt.coef_spec("f", "dT", "*", 0.0))
    dpsi = np.linspace(0, 2 * math.pi, 181)
    d1, t1 = table("drv.f.dT", dpsi / KL)
    d2, t2 = table("laser.las.phi", np.zeros(len(dpsi)), "rad")
    pipe = [d1, d2, sj(topo, {"f": pm(dpsi=0.0, K_T=KL, ti=th, to=th)}),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"})]
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(pipe), lay, 0)
    s, p = State(lay), Params(lay)
    for k, v in {**t1, **t2}.items():
        s[k] = v
    n, dt = 1 << 14, 5e-15
    t = (np.arange(n) - n // 2) * dt
    modes = {}
    for tp in (100e-15, 1e-12):
        sig = tp / (2 * math.sqrt(2 * math.log(2)))
        u = np.exp(-t ** 2 / (4 * sig ** 2))
        modes[tp] = (u / math.sqrt(np.sum(u ** 2) * dt), sig)
    Ro = R(th)
    cw, pulse, sat = [], {tp: [] for tp in modes}, None
    for _ in dpsi:
        run_reference(blk, s, p, 1)
        cw.append(float(s["beam.sq.phi"]))
        J = np.array(s["J.f"])

        def U(om, J=J):
            d = np.zeros((len(om), 2, 2), complex)
            d[:, 0, 0], d[:, 1, 1] = 1.0, np.exp(-1j * tau * om)
            return Ro @ d @ Ro.T @ J

        for tp, (u, sig) in modes.items():
            f = pulse_projection(u.astype(complex), dt, U, (1, 0), (1, 0))
            pulse[tp].append(np.angle(np.sum(u * f) * dt))
            if tp == 100e-15 and sat is None:
                sat = np.abs(f) ** 2, sig
    amp_cw = np.max(np.abs(cw))
    assert amp_cw == pytest.approx(k2, rel=0.01)  # T09: 3.16e−3 rad
    assert amp_cw == pytest.approx(math.asin(math.tan(th) ** 2), rel=1e-4)
    assert np.max(np.abs(pulse[100e-15])) <= 1e-4 * k2
    ov = math.exp(-tau ** 2 / (8 * modes[1e-12][1] ** 2))
    assert ov == pytest.approx(0.277, abs=1e-3)
    assert np.max(np.abs(pulse[1e-12])) == pytest.approx(math.asin(math.tan(th) ** 2 * ov), rel=1e-4)
    pw, sig = sat
    win = np.abs(t + tau) < 6 * sig
    assert np.sum(pw[win]) * dt == pytest.approx(math.sin(th) ** 4, rel=1e-6)  # κ1²κ2²
    assert math.sin(th) ** 4 == pytest.approx(1e-5, rel=0.01)
    assert t[win][np.argmax(pw[win])] == pytest.approx(-tau, abs=dt)
    assert np.sum(pw[np.abs(t - tau) < 6 * sig]) * dt < 1e-20  # nothing at +1.36 ps


# ------------------------------------------------------------------------------------- PMD, PDL, off


def test_pmd_mode_overlap_t09_numbers():
    """T09 §4.2 (worst case, D_PMD = 0.1 ps/√km): 100 fs → 4.2 % at 63 m, 32 % at 630 m; 1 ps < 0.5 %."""
    from cranberry.blocks.pol.pmd_mode_overlap import make_pmd_mode_overlap

    def loss(L, tp, ref="psp"):
        return 1 - make_pmd_mode_overlap({
            "D_PMD_s_per_sqrt_m": 0.1e-12 / math.sqrt(1e3), "L_m": L, "tau_p_fwhm_s": tp, "timing_reference": ref,
            "segment": "smf", "applies_to": ["gbs"], "source": "[T09]", "toggles": {"pol_pmd": True}}).eta

    assert loss(63.0, 100e-15) == pytest.approx(0.0423, abs=1e-4)
    assert loss(630.0, 100e-15) == pytest.approx(0.3225, abs=1e-4)
    assert loss(630.0, 1e-12) < 0.005 and loss(63.0, 1e-12) < 0.0005
    assert loss(630.0, 100e-15, "centered") == pytest.approx(0.196, abs=1e-3)  # re-timed reference


def test_pdl_power_drift_range():
    """T09 §1.4: through a PDL element the power moves in T_mean·(1 ± (T_max − T_min)/(T_max + T_min))
    as the input polarization turns; beam.power (beam_assembler, analyzer x) follows |M_xx e_x|²."""
    seg = pt.seg_spec("pdl", 0.0, pol_model="jones", pol_params={"J_slot": "J.pdl"})
    beams = [pt.beam_spec("b", [pt.step_spec("rot"), pt.step_spec("pdl")], "o", jones_out=(1, 0))]
    topo = pt.topo_spec([jseg("rot", 1.0, kind="smf"), seg], beams, ["o"], pt.coef_spec("rot", "dT", "*", 0.0))
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 0, 1], "u0": "identity", "drift": "diffusion",
           "D_pol_rad2_per_s": 0.0, "k_T_rad_per_m_K": 1.0}
    th = np.linspace(0, 2 * math.pi, 73)  # Stokes rotation angle of the input (linear) state
    d1, t1 = table("drv.rot.dT", th)
    d2, t2 = table("laser.las.phi", np.zeros(len(th)), "rad")
    pdl = {"block": "pdl_element", "name": "pdl", "params": {
        "out": "J.pdl", "pdl_db": 0.3, "axis_rad": 0.0, "il_db": 0.5, "segment": "pdl", "applies_to": ["gbs"],
        "source": "[AE]", "toggles": {"pol_pdl": True}}}
    pipe = [d1, d2, sj(topo, {"rot": smf}), pdl, *pt.path_pipeline(topo, {"las": "laser.las.phi"}),
            pt.record(["beam.b.power"], len(th))]
    rec, _, _ = pt.run_and_record(pipe, 1e6, len(th), {**t1, **t2})
    q = 10 ** (-0.03)
    np.testing.assert_allclose(rec["beam.b.power"] / 1e-3, 2 / (1 + q) * np.cos(th / 2) ** 2, atol=1e-14)
    tot = 2 / (1 + q) * np.cos(th / 2) ** 2 + 2 * q / (1 + q) * np.sin(th / 2) ** 2  # both analyzer ports
    assert tot.max() == pytest.approx(1 + (1 - q) / (1 + q), rel=1e-12)
    assert tot.min() == pytest.approx(1 - (1 - q) / (1 + q), rel=1e-12)


def test_m4_toggles_off_identity():
    """pol_pdl / pol_faraday / pol_smf_drift (Haar U0, OU) off ⇒ every M4 matrix is 1 (T09 test 8)."""
    topo = shared_topology(True)
    lay = SlotLayout(1e6)
    for s in ("s1", "s2"):
        lay.add(f"drv.{s}.dT", (), np.float64, "K", init=0.37)
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "u0": "haar", "drift": "ou",
           "D_pol_rad2_per_s": 50.0, "k_T_rad_per_m_K": 3.0, "f_c_hz": 10.0}
    pipe = [sj(topo, {"s1": pm(dpsi=1.1, ti=0.05, to=-0.04), "s2": smf}, OFF),
            {"block": "pdl_element", "name": "pdl", "params": {
                "out": "J.pdl", "pdl_db": 0.3, "axis_rad": 0.2, "il_db": 0.5, "segment": "pdl", "applies_to": ["gbs"],
                "source": "[AE]", "toggles": {"pol_pdl": False}}},
            far_spec(False)]
    blk = build(spec_from_obj(pipe), lay, 1, external=["drv.s1.dT", "drv.s2.dT"])
    s, p = State(lay), Params(lay)
    run_reference(blk, s, p, 5)
    for j in ("J.s1", "J.s2", "J.pdl", "J.far"):
        np.testing.assert_array_equal(s[j], np.eye(2))
