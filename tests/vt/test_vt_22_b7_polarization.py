"""VT-22 (M4 part, WP-B7): shared-segment cancellation with every polarization block in the shared path
(prompt §7-11 / D51-4; T09 test 8; G04 "偏光をオフにしても通る").

Shared path: PM segment (``segment_jones`` pm: Δψ_b drift + end crosstalk), SMF segment (Haar U0 + OU
drift + temperature rotation), PDL element (``pdl_element``), Faraday element (``faraday_nonreciprocal``,
J_nr_slot); then one non-shared tail per beam. Squeezed and probe co-propagate with the same input /
analyzer states.

* All B7 toggles on: the lock error φ_sq − φ_pr is bit-for-bit the same with the shared-segment noise on
  or off (< 1e−12 relative to the beam-phase excursion); the polarization matrices really vary.
* All B7 toggles off (U = N = M = 1): the run equals the pol_model-off topology (1e−12), i.e. VT-22 passes
  unchanged.
"""

import numpy as np
import pytest

import cranberry.path.testing as pt

FS, N = 1e8, 300
ON = {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}
OFF = dict.fromkeys(ON, False)


def topology(pol: bool):
    def seg(id, L, kind="pm", n_sub=1, **pol_kw):
        return pt.seg_spec(id, L, kind=kind, n_sub=n_sub, drivers={"dT": "uniform"},
                           pol_model="jones" if pol else "off", pol_params=pol_kw if pol else {})

    segs = [seg("s1", 7.3, n_sub=3, J_slot="J.s1"), seg("s2", 2.1, kind="smf", J_slot="J.s2"),
            pt.seg_spec("pdl", 0.0, kind="component", pol_model="jones" if pol else "off",
                        pol_params={"J_slot": "J.pdl"} if pol else {}),
            pt.seg_spec("far", 0.0, kind="component", pol_model="jones" if pol else "off",
                        pol_params={"J": np.eye(2).tolist(), "J_nr_slot": "J.far"} if pol else {}),
            pt.seg_spec("t_sq", 1.7, n_sub=2, drivers={"dT": "uniform"}),
            pt.seg_spec("t_pr", 1.7, n_sub=2, drivers={"dT": "uniform"})]
    path = [pt.step_spec(s) for s in ("s1", "s2", "pdl", "far")]
    jout = {"jones_out": (1, 0)} if pol else {}
    beams = [pt.beam_spec("sq", path + [pt.step_spec("t_sq")], "bhd", **jout),
             pt.beam_spec("pr", path + [pt.step_spec("t_pr")], "bhd", **jout)]
    coef = (pt.coef_spec("s1", "dT", "*", 320.0) + pt.coef_spec("s2", "dT", "*", 92.0)
            + pt.coef_spec("t_sq", "dT", "*", 75.0) + pt.coef_spec("t_pr", "dT", "*", 75.0))
    locks = [{"id": "co", "combo": {"sq": 1.0, "pr": -1.0}, "setpoint_rad": 0.3, "observer": "bhd"}]
    return pt.topo_spec(segs, beams, ["bhd"], coef, locks)


def b7_blocks(topo, on: bool):
    sj = {"s1": {"model": "pm", "axis_drivers": ["dT"], "dpsi_b0_rad": 1.1, "K_T_rad_per_m_K": -2.05,
                 "theta_in_rad": 0.05, "theta_out_rad": -0.04},
          "s2": {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [0, 1, 0], "u0": "haar", "drift": "ou",
                 "D_pol_rad2_per_s": 1e5, "k_T_rad_per_m_K": 3.0, "f_c_hz": 1e3}}
    return [{"block": "segment_jones", "name": "sj", "params": {"topology": topo, "segments": sj,
                                                                "toggles": dict(ON if on else OFF)}},
            {"block": "pdl_element", "name": "pdl", "params": {
                "out": "J.pdl", "pdl_db": 0.5, "axis_rad": 0.4, "il_db": 0.5, "segment": "pdl",
                "applies_to": ["gbs", "probe"], "source": "[AE]", "toggles": {"pol_pdl": on}}},
            {"block": "faraday_nonreciprocal", "name": "far", "params": {
                "out": "J.far", "verdet_rad_per_T_m": 0.6, "B_par_T": 50e-6, "L_m": 63.0,
                "toggles": {"pol_faraday": on}}}]


def run(pol: bool, on: bool, shared_sigma=1.0):
    topo = topology(pol)
    pipe = [pt.pathtest("pathtest_white", f"n.{s['id']}", f"drv.{s['id']}.dT", (), "K",
                        sigma=1e-3 * (shared_sigma if s["id"] in ("s1", "s2") else 1.0))
            for s in topo["segments"] if s["drivers"]]
    pipe.append(pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=1.0))
    if pol:
        pipe += b7_blocks(topo, on)
    pipe += [*pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True),
             pt.record(["lock.co.err", "beam.sq.phi", "beam.sq.power"], N)]
    rec, st, _ = pt.run_and_record(pipe, FS, N, seed=22)
    return rec, st


def test_b7_on_shared_noise_cancels():
    on, st = run(True, True)
    off, _ = run(True, True, shared_sigma=0.0)
    scale = np.max(np.abs(on["beam.sq.phi"]))
    assert np.max(np.abs(on["lock.co.err"] - off["lock.co.err"])) < 1e-12 * scale  # prompt §7-11
    assert np.std(on["beam.sq.power"]) > 1e-6  # polarization really moves (OU drift, PDL)
    for j in ("J.s1", "J.s2", "J.pdl", "J.far"):
        assert not np.allclose(st[j], np.eye(2))


def test_b7_off_equals_pol_off_topology():
    on, _ = run(True, False)
    plain, _ = run(False, False)
    scale = np.max(np.abs(on["beam.sq.phi"]))
    np.testing.assert_allclose(on["lock.co.err"], plain["lock.co.err"], rtol=0, atol=1e-12 * scale)
    np.testing.assert_allclose(on["beam.sq.phi"], plain["beam.sq.phi"], rtol=0, atol=1e-12 * scale)
    np.testing.assert_allclose(on["beam.sq.power"], plain["beam.sq.power"], rtol=0, atol=0)


@pytest.mark.numba
def test_b7_compiled_equals_reference():
    topo = topology(True)
    pipe = [pt.pathtest("pathtest_white", f"n.{s['id']}", f"drv.{s['id']}.dT", (), "K", sigma=1e-3)
            for s in topo["segments"] if s["drivers"]]
    pipe += [pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=1.0), *b7_blocks(topo, True),
             *pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True),
             pt.record(["lock.co.err", "beam.sq.phi", "beam.pr.power"], N)]
    ref, _, _ = pt.run_and_record(pipe, FS, N, seed=22)
    cmp, _, _ = pt.run_and_record(pipe, FS, N, seed=22, compiled=True)
    for k in ref:
        np.testing.assert_allclose(cmp[k], ref[k], rtol=0, atol=1e-12 * max(1.0, np.max(np.abs(ref[k]))))
