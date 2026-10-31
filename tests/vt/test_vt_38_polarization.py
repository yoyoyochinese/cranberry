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
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 0, 1], "D_pol_rad2_per_s": 0.0,
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
        smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "D_pol_rad2_per_s": 50.0,
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
    smf = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "D_pol_rad2_per_s": 50.0,
           "k_T_rad_per_m_K": 3.0}
    blk = build(spec_from_obj([sj(topo, {"s1": pm(dpsi=1.1, ti=0.05, to=-0.04), "s2": smf}, OFF)]), lay, 1,
                external=["drv.s1.dT", "drv.s2.dT"])
    s, p = State(lay), Params(lay)
    run_reference(blk, s, p, 5)
    for j in ("J.s1", "J.s2"):
        np.testing.assert_array_equal(s[j], np.eye(2))
    assert not any(leaf.draws for leaf in iter_leaves(blk))
