"""VT-22 (M0 part, WP-F6): co-propagation / shared-segment cancellation (prompt §7-11 restricted by
D51-4 to shared segments with equal segment × beam coefficients; 00_summary §4: T10 test 1, T05 V3,
T09 test 8). Later parts: WP-I4 (Ref B, M3), WP-B7 (polarization, M4).

* All segments shared, same direction and axis, equal coefficients, noise on every segment and on the
  laser: the lock combination φ_sq − φ_pr carries none of it (|e + setpoint| < 1e−12 relative to the
  beam-phase excursion; prompt tolerance 1e−6, laser common path 1e−12).
* Shared segments followed by one non-shared segment per beam (equal group delay): the lock error is
  bit-for-bit the same with the shared-segment noise on or off (same seed, < 1e−12 relative), i.e.
  shared noise never enters the lock's residual; the non-shared noise does.
* Unequal coefficients on a shared segment (1 %): the shared noise appears at 1 % (cancellation is
  computed, not assumed, D11; the full test is VT-23).
* Polarization off (segment_phase.pm_axis_drift and beam_assembler.polarization off, J = e^{iψ}·1, T09
  test 8) passes the same check with jones / pm_diagonal segments in the shared path; with polarization
  on and identical input / analyzer states it passes too.
"""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt

FS = 1.0e8
N = 400
PM = {"K_T_rad_per_m_K": -2.05, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.36e-12, "axis_drivers": ["dT"]}


def rand_su2(seed):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    q = q * (np.diag(r) / abs(np.diag(r)))
    return (q / np.sqrt(np.linalg.det(q))).tolist()


def topology(*, tails: bool, pr_scale: float = 1.0, pol: str = "off"):
    shared = [
        pt.seg_spec("s1", 7.3, n_sub=5, spatial="mixed", drivers={"dT": "uniform", "strain": "per_sub"}),
        pt.seg_spec("s2", 2.1, n_sub=2, drivers={"dT": "uniform"},
                    pol_model="pm_diagonal" if pol != "off" else "off", pol_params=PM if pol != "off" else {}),
        pt.seg_spec("s3", 0.9, drivers={"dT": "uniform"},
                    pol_model="jones" if pol != "off" else "off",
                    pol_params={"J": rand_su2(1)} if pol != "off" else {}),
    ]
    tail_segs = [pt.seg_spec("t_sq", 1.7, n_sub=3, drivers={"dT": "uniform"}),
                 pt.seg_spec("t_pr", 1.7, n_sub=3, drivers={"dT": "uniform"})] if tails else []
    ax = {"axis": "slow"} if pol != "off" else {}
    path = [pt.step_spec("s1"), pt.step_spec("s2", **ax), pt.step_spec("s3")]
    jout = [1, 0] if pol != "off" else None
    beams = [pt.beam_spec("sq", path + ([pt.step_spec("t_sq")] if tails else []), "bhd", jones_out=jout),
             pt.beam_spec("pr", path + ([pt.step_spec("t_pr")] if tails else []), "bhd", jones_out=jout)]
    coef = (pt.coef_spec("s1", "dT", "*", 320.0, 3e-13) + pt.coef_spec("s1", "strain", "*", 5e6)
            + pt.coef_spec("s2", "dT", "*", 92.0) + pt.coef_spec("s3", "dT", "*", 40.0)
            + pt.coef_spec("s1", "strain", "pr", 5e6 * pr_scale))
    if tails:
        coef += pt.coef_spec("t_sq", "dT", "*", 75.0) + pt.coef_spec("t_pr", "dT", "*", 75.0)
    locks = [{"id": "co", "combo": {"sq": 1.0, "pr": -1.0}, "setpoint_rad": 0.3, "observer": "bhd"}]
    return pt.topo_spec(shared + tail_segs, beams, ["bhd"], coef, locks)


def run(topo, sigmas, **toggles):
    drivers = []
    for s in topo["segments"]:
        for d, pl in s["drivers"].items():
            shape = (s["n_sub"],) if pl == "per_sub" else ()
            drivers.append(pt.pathtest("pathtest_white", f"n.{s['id']}.{d}", f"drv.{s['id']}.{d}", shape,
                                       "K" if d == "dT" else "1", sigma=sigmas.get(s["id"], 1.0) * (1e-3 if d == "dT" else 1e-7)))
    pipe = drivers + [pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=1.0),
                      *pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True, **toggles),
                      pt.record(["beam.sq.phi", "beam.pr.phi", "lock.co.err"], N)]
    rec, _, _ = pt.run_and_record(pipe, FS, N, seed=2022)
    return rec


def test_all_segments_shared_cancel_exactly():
    rec = run(topology(tails=False), {})
    excursion = np.max(np.abs(rec["beam.sq.phi"]))
    assert excursion > 1.0  # laser + path noise of order rad
    assert np.max(np.abs(rec["lock.co.err"] + 0.3)) < 1e-12 * excursion


def test_shared_noise_absent_from_lock_error_with_tails():
    topo = topology(tails=True)
    on = run(topo, {})
    off = run(topo, {"s1": 0.0, "s2": 0.0, "s3": 0.0})
    scale = np.max(np.abs(on["beam.sq.phi"]))
    assert np.max(np.abs(on["lock.co.err"] - off["lock.co.err"])) < 1e-12 * scale
    tails_off = run(topo, {"t_sq": 0.0, "t_pr": 0.0})
    assert np.std(on["lock.co.err"]) > 1e-3  # the non-shared tails do enter
    assert np.max(np.abs(tails_off["lock.co.err"] + 0.3)) < 1e-12 * scale


def test_unequal_coefficients_do_not_cancel():
    eq = run(topology(tails=False), {})
    neq = run(topology(tails=False, pr_scale=1.01), {})
    shared_strain = run(topology(tails=False, pr_scale=0.0), {})  # e = ψ_s1,strain(sq) − 0·… → strain part
    resid = neq["lock.co.err"] + 0.3
    assert np.max(np.abs(eq["lock.co.err"] + 0.3)) < 1e-12 * np.max(np.abs(eq["beam.sq.phi"]))
    np.testing.assert_allclose(resid, -0.01 * (shared_strain["lock.co.err"] + 0.3), rtol=1e-9, atol=1e-15)


@pytest.mark.parametrize("pol_on", [False, True])
def test_polarization_off_and_on(pol_on):
    topo = topology(tails=True, pol="jones")
    toggles = {"pm_axis_drift": pol_on, "polarization": pol_on}
    on = run(topo, {}, **toggles)
    off = run(topo, {"s1": 0.0, "s2": 0.0, "s3": 0.0}, **toggles)
    scale = np.max(np.abs(on["beam.sq.phi"]))
    assert np.max(np.abs(on["lock.co.err"] - off["lock.co.err"])) < 1e-12 * scale
    if not pol_on:  # J = e^{iψ}·1: identical to the pol_model-off topology
        plain = run(topology(tails=True), {})
        np.testing.assert_allclose(on["lock.co.err"], plain["lock.co.err"], rtol=0, atol=1e-12 * scale)
        assert math.isfinite(float(plain["lock.co.err"][-1]))


# ------------------------------------------------------------------------------------------------------------------
# WP-I4 (M3): Ref B part — the co-propagating probe in the refB preset, closed loops (refB.yaml header).
#
# The probe `pr_bhd` and the squeezed light (θ_sq bookkeeping beam `sqz`) share the OPO and the OPO → HD path SQ.
# `lock2` = φ_LO − φ_pr (lock 2 itself holds φ_LO − φ_prf, `lock2f`, M3B-05); metric (i) is Δ = φ_LO − θ_sq; the
# difference Δ − lock2 is `cop` = φ_pr − θ_sq at the HD,
# the part of the squeezing angle the locks cannot see. With table acoustics on the non-common paths (PUMP,
# PROBE, LO) running, the acoustics of the shared SQ (block env.shared) are switched on / off with the same seed:
#   * cop is identical to 1e−12 of the beam-phase excursion (shared noise cancels, D11: computed, not assumed);
#   * Δ and lock2 change (the shared noise is an in-loop disturbance that lock 2 sees one to one), but
#     Δ − lock2 = cop does not;
#   * the non-common noise does enter cop (lock-1 residual on the probe path), and the OPO detuning enters it
#     with the coefficient difference c_pr − c_sq only (the exception named by D51-4; VT-23). Since Phase 5
#     M3B-03 lock 1 senses Δ (its PD is after the OPO), so its actuator also responds to Δ; the path-only check
#     holds the lock-1 actuator (pid1 feedback off in both runs) and uses c_pr = the transit coefficient −1/κ.
# ------------------------------------------------------------------------------------------------------------------

import importlib.util  # noqa: E402
from pathlib import Path  # noqa: E402

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parents[1] / "integration" / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)


@pytest.mark.numba
@pytest.mark.integration
def test_refB_shared_path_noise_absent_from_cop():
    from cranberry.runner.run import run

    slots = ["lock.cop.err", "lock.lock2.err", "metric.mi.port0", "beam.sqz.phi", "drv.OPO.detuning"]

    def go(shared: bool, detuning: bool = False, lock1: bool = True):
        tg = {**rs.NOISE_OFF, "env.mech.mech": True, "env.shared.mech": shared, "opo_dn.length_noise": detuning,
              "pid1.feedback": lock1}
        return run(rs.load(run={"duration_s": 0.05}, toggles=tg, measure={"record": {"slots": slots, "decim": 1}})).series

    on, off = go(True), go(False)
    scale = np.max(np.abs(on["beam.sqz.phi"]))
    assert scale > 0.1                                                    # rad-scale beam-phase excursions
    np.testing.assert_allclose(on["lock.cop.err"], off["lock.cop.err"], rtol=0, atol=1e-12 * scale)
    d_on, d_off = on["metric.mi.port0"] - on["lock.lock2.err"], off["metric.mi.port0"] - off["lock.lock2.err"]
    np.testing.assert_allclose(d_on, on["lock.cop.err"], rtol=0, atol=1e-12 * scale)   # Δ − lock2 = cop
    np.testing.assert_allclose(d_on, d_off, rtol=0, atol=1e-12 * scale)
    assert np.std(on["metric.mi.port0"] - off["metric.mi.port0"]) > 1e-4  # shared noise is in-loop for lock 2
    assert np.std(on["lock.cop.err"]) > 1e-3                              # the non-common paths do enter cop
    # OPO detuning: cop changes by (c_pr − c_sq)·Δ_opo only (D51-4 exception; the full check is VT-23), lock-1
    # actuator held (M3B-03: a closed lock 1 also follows Δ)
    det, ref = go(True, detuning=True, lock1=False), go(True, lock1=False)
    p = rs.params_of(rs.load())
    c_res = -1.0 / p["kappa_rad_s"] + 1 / (p["kappa_rad_s"] * (1 + p["x"] ** 2))
    dcop, drv = det["lock.cop.err"] - ref["lock.cop.err"], det["drv.OPO.detuning"]
    assert np.std(drv) > 1e-3 * p["kappa_rad_s"]
    # read at the HD through ONE linear interpolation at (half the lumped OPO segment + SQ) group delay (F6, D68-3)
    mu = (0.5 * 0.5 * 1.0 + 1.0 * 1.0003) / 299792458.0 * 2.0e6
    want = c_res * ((1 - mu) * drv[1:] + mu * drv[:-1])
    np.testing.assert_allclose(dcop[1:], want, rtol=1e-9, atol=1e-12 * scale)
