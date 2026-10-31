"""WP-F6: segment_phase / beam_assembler / lock_error blocks — compiled == reference (1e−12, VT-44 rule
for kernels, PLAN §7.3), toggles off ⇒ identity, subsegment and emission timing, wiring errors."""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import PipelineCheckError, build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.path.locks import lock_error
from cranberry.path.topology import C_M_PER_S, Topology

PM = {"K_T_rad_per_m_K": -2.05, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.36e-12, "axis_drivers": ["dT", "tec"]}


def rand_unitary(rng):
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / abs(np.diag(r)))


def make(pipe, fs, *, external=(), setup=None, seed=3):
    lay = SlotLayout(fs)
    if setup:
        setup(lay)
    blk = build(spec_from_obj(pipe), lay, seed, external=external)
    return blk, State(lay), Params(lay)


def rich_topology():
    rng = np.random.default_rng(5)
    U1 = rand_unitary(rng)
    segs = [
        pt.seg_spec("c", 13.7, n_sub=7, spatial="mixed", actuator_sub=2,
                    drivers={"strain": "per_sub", "dT": "uniform", "stretcher": "actuator", "tec": {"subs": [4, 5]}}),
        pt.seg_spec("p", 3.3, n_sub=3, drivers={"dT": "uniform", "tec": "uniform"}, pol_model="pm_diagonal", pol_params=PM),
        pt.seg_spec("j1", 1.1, pol_model="jones", pol_params={"J": U1.tolist()}, drivers={"dT": "uniform"}),
        pt.seg_spec("j2", 0.7, pol_model="jones", pol_params={"J_slot": "J.j2"}),
        pt.seg_spec("f", 2.9, n_sub=2, kind="free_space", n_g=1.0),
    ]
    beams = [
        pt.beam_spec("sq", [pt.step_spec("c"), pt.step_spec("p", axis="slow"), pt.step_spec("j1"), pt.step_spec("j2"),
                            pt.step_spec("f")], "gbs", jones_out=[1, 0]),
        pt.beam_spec("pr", [pt.step_spec("f", -1), pt.step_spec("j2", -1), pt.step_spec("j1", -1),
                            pt.step_spec("p", -1, axis="fast"), pt.step_spec("c", -1, sub_in=5, sub_out=1)], "det",
                     jones_in=[0, 1], jones_out=[2 ** -0.5, 1j * 2 ** -0.5]),
        pt.beam_spec("lo", [pt.step_spec("p", axis="slow")], "det", source="las2", freq_offset_hz=80e6),
    ]
    coef = (pt.coef_spec("c", "strain", "*", 1.3e6, 2e-9) + pt.coef_spec("c", "dT", "*", 600.0, 5e-13)
            + pt.coef_spec("c", "dT", "pr", 610.0, 5.1e-13) + pt.coef_spec("c", "stretcher", "*", 6.49)
            + pt.coef_spec("c", "tec", "*", 600.0, 5e-13) + pt.coef_spec("p", "dT", "*", 145.0, 1.3e-13)
            + pt.coef_spec("p", "tec", "*", 145.0, 1.3e-13) + pt.coef_spec("j1", "dT", "*", 48.0))
    locks = [{"id": "lk", "combo": {"lo": 1.0, "pr": -1.0}, "setpoint_rad": math.pi / 2, "observer": "det"}]
    return pt.topo_spec(segs, beams, ["gbs", "det"], coef, locks), rng


def rich_pipeline(topo, **toggles):
    drivers = [pt.pathtest("pathtest_white", "w.strain", "drv.c.strain", (7,), "1", sigma=1e-6),
               pt.pathtest("pathtest_white", "w.cdT", "drv.c.dT", (), "K", sigma=1e-3),
               pt.pathtest("pathtest_white", "w.str", "drv.c.stretcher", (), "V", sigma=0.1),
               pt.pathtest("pathtest_white", "w.ctec", "drv.c.tec", (), "K", sigma=1e-3),
               pt.pathtest("pathtest_white", "w.pdT", "drv.p.dT", (), "K", sigma=1e-3),
               pt.pathtest("pathtest_white", "w.ptec", "drv.p.tec", (), "K", sigma=1e-3),
               pt.pathtest("pathtest_white", "w.j1", "drv.j1.dT", (), "K", sigma=1e-3),
               pt.pathtest("pathtest_white", "w.las", "laser.las.phi", (), "rad", sigma=0.3),
               pt.pathtest("pathtest_white", "w.las2", "laser.las2.phi", (), "rad", sigma=0.3)]
    return drivers + pt.path_pipeline(topo, {"las": "laser.las.phi", "las2": "laser.las2.phi"}, locks=True, **toggles)


@pytest.mark.numba
@pytest.mark.parametrize("toggles", [{}, {"phase": False, "group_delay": False, "pm_axis_drift": False,
                                           "emission_time": False, "polarization": False}])
def test_compiled_equals_reference(toggles):
    topo, rng = rich_topology()
    U2 = rand_unitary(rng)

    def setup(lay):
        lay.add("J.j2", (2, 2), np.complex128, "1")

    runs = []
    for runner in (run_reference, run_compiled):
        blk, st, pa = make(rich_pipeline(topo, **toggles), 2.0e8, external=["J.j2"], setup=setup)
        st["J.j2"] = U2
        runner(blk, st, pa, 1500)
        runs.append(st)
    a, b = runs
    scale = np.maximum(1.0, np.abs(a.S))
    assert np.max(np.abs(a.S - b.S) / scale) < 1e-12
    assert np.max(np.abs(a.Z - b.Z)) < 1e-12
    assert np.array_equal(a.I, b.I)
    assert abs(float(a["beam.sq.phi"])) > 0 and abs(float(a["lock.lk.err"])) > 0


def test_toggles_off_identity():
    """phase off → ψ = 0; group_delay off → nominal τ_g; polarization off → J = e^{iψ}·1; emission off →
    instantaneous; pm_axis_drift off → ψ_fast = ψ_slow for equal coefficients (T09 test 8)."""
    topo, rng = rich_topology()
    U2 = rand_unitary(rng)
    t = Topology.from_spec(topo)

    def run(**tg):
        blk, st, pa = make(rich_pipeline(topo, **tg), 2.0e8, external=["J.j2"],
                           setup=lambda lay: lay.add("J.j2", (2, 2), np.complex128, "1"))
        st["J.j2"] = U2
        hist = {"las": [], "las2": []}
        for _ in range(400):
            run_reference(blk, st, pa, 1)
            hist["las"].append(float(st["laser.las.phi"]))
            hist["las2"].append(float(st["laser.las2.phi"]))
        return st, hist

    st, hist = run(phase=False, group_delay=False, pm_axis_drift=False, polarization=False, emission_time=False)
    for name in st.layout.slot_names():
        if name.startswith("psi.") and not name.endswith(".sub"):
            assert float(st[name]) == 0.0
        if name.startswith("tau_g."):
            seg, beam = name.split(".")[1:]
            assert float(st[name]) == pytest.approx(next(s.tau_s for s in t.timing(beam).steps if s.segment == seg),
                                                    rel=1e-15)
    assert float(st["beam.sq.phi"]) == hist["las"][-1]  # instantaneous, no segment phases, no pol phase
    assert float(st["beam.lo.phi"]) == hist["las2"][-1]  # f_off phase vanishes without emission delay
    assert float(st["beam.pr.power"]) == 1e-3 and np.array_equal(st["beam.pr.jones"], [0, 1])

    # pm_axis_drift off: the fast-axis probe and slow-axis LO see the same p-segment phase
    st2, _ = run(pm_axis_drift=False)
    st3, _ = run(pm_axis_drift=True)
    assert float(st2["psi.p.pr"]) == pytest.approx(float(st2["psi.p.lo"]), rel=1e-12)
    assert float(st3["psi.p.pr"]) != pytest.approx(float(st3["psi.p.lo"]), rel=1e-6)


def test_point_actuator_timing_both_directions():
    """A stretcher step at subsegment a reaches the forward beam (exit x = L) after (N − a − ½)δ and the
    backward beam (exit x = 0) after (a + ½)δ (conv §9, T10 §3.2)."""
    N, a = 8, 2
    L = 2 * N * C_M_PER_S / 1.468 / 1e8  # δ = 2 samples at 100 MHz
    segs = [pt.seg_spec("c", L, n_sub=N, actuator_sub=a, drivers={"stretcher": "actuator"})]
    beams = [pt.beam_spec("f", [pt.step_spec("c", 1)], "x1"), pt.beam_spec("b", [pt.step_spec("c", -1)], "x0")]
    topo = pt.topo_spec(segs, beams, ["x1", "x0"], pt.coef_spec("c", "stretcher", "*", 1.0))
    pipe = [pt.pathtest("pathtest_table", "tab", "drv.c.stretcher", (), "V", n_rows=200),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1), *pt.path_pipeline(topo, {"las": "laser.las.phi"})]
    blk, st, pa = make(pipe, 1e8)
    st["drv.c.stretcher_table"][10:] = 1.0
    f, b = [], []
    for _ in range(60):
        run_reference(blk, st, pa, 1)
        f.append(float(st["beam.f.phi"]))
        b.append(float(st["beam.b.phi"]))
    # sample n written at step n: step starts at n = 10
    assert np.flatnonzero(np.array(f) > 0.5)[0] == 10 + 2 * (N - a) - 1
    assert np.flatnonzero(np.array(b) > 0.5)[0] == 10 + 2 * a + 1


def test_emission_delay_fractional_ramp_exact():
    """Linear interpolation is exact for a ramp: φ_b(t) = φ_L(t − τ_b) + 2π f_off τ_b to rounding."""
    L = 17.3
    segs = [pt.seg_spec("c", L, n_sub=3)]
    topo = pt.topo_spec(segs, [pt.beam_spec("b", [pt.step_spec("c")], "o", freq_offset_hz=1e6)], ["o"], [])
    pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=500),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"})]
    fs = 1e9
    blk, st, pa = make(pipe, fs)
    slope = 2.0e5  # rad/s
    st["laser.las.phi_table"][:] = slope * np.arange(500) / fs
    run_reference(blk, st, pa, 300)
    tau = 1.468 * L / C_M_PER_S
    t_last = 299 / fs
    assert float(st["beam.b.phi"]) == pytest.approx(slope * (t_last - tau) + 2 * math.pi * 1e6 * tau, rel=1e-12)
    assert float(st["beam.b.t_emit"]) == pytest.approx(t_last - tau, rel=1e-12)


def test_lock_error_function_and_block():
    topo, _ = rich_topology()
    t = Topology.from_spec(topo)
    lk = t.locks["lk"]
    assert lock_error(lk, {"lo": 1.0, "pr": 0.25}) == pytest.approx(0.75 - math.pi / 2)
    with pytest.raises(KeyError):
        lock_error(lk, {"lo": 1.0})
    blk, st, pa = make(rich_pipeline(topo), 2e8, external=["J.j2"],
                       setup=lambda lay: lay.add("J.j2", (2, 2), np.complex128, "1"))
    st["J.j2"] = np.eye(2)
    pa["le.lk.setpoint_rad"] = 0.1
    run_reference(blk, st, pa, 5)
    assert float(st["lock.lk.err"]) == pytest.approx(float(st["beam.lo.phi"]) - float(st["beam.pr.phi"]) - 0.1, abs=1e-12)


def test_wiring_errors():
    segs = [pt.seg_spec("c", 1.0, n_sub=4, drivers={"strain": "per_sub"}, spatial="distributed")]
    topo = pt.topo_spec(segs, [pt.beam_spec("b", [pt.step_spec("c")], "o")], ["o"], pt.coef_spec("c", "strain", "*", 1.0))
    las = pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1)
    bad_shape = [pt.pathtest("pathtest_table", "d", "drv.c.strain", (), "1", n_rows=1), las,
                 *pt.path_pipeline(topo, {"las": "laser.las.phi"})]
    with pytest.raises(TypeError, match="shape"):
        make(bad_shape, 1e8)
    with pytest.raises(ValueError, match="no slot given for source"):
        make([las, *pt.path_pipeline(topo, {})], 1e8)
    no_segment_phase = [las, pt.path_pipeline(topo, {"las": "laser.las.phi"})[1]]
    with pytest.raises(PipelineCheckError, match="psi.c.b"):
        make(no_segment_phase, 1e8)


REFD_YAML = """
preset: f6_refD_like
run: {duration_s: 2.0e-6, f_sim_hz: 1.0e8, chunk_samples: 64, seed: 9}
topology: &topo
  segments:
    - {id: A1, kind: pm, L_m: 10.5, n_g: 1.468, beta2_s2_per_m: 0.0, beta3_s3_per_m: 0.0, zone: lab, n_sub: 4,
       spatial: lumped, pol_model: "off", pol_params: {}, actuator_sub: 1, drivers: {dT: uniform, stretcher: actuator}}
    - {id: A2, kind: pm, L_m: 0.5, n_g: 1.468, beta2_s2_per_m: 0.0, beta3_s3_per_m: 0.0, zone: lab, n_sub: 1,
       spatial: lumped, pol_model: "off", pol_params: {}, drivers: {dT: uniform}}
  beams:
    - {id: pr, source: dx1, path: [{segment: A1, direction: 1}], freq_offset_hz: 0.0, jones_in: [1, 0],
       power_w: 1.6e-3, observer: pd}
    - {id: lo, source: dx1, path: [{segment: A2, direction: 1}], freq_offset_hz: 0.0, jones_in: [1, 0],
       power_w: 1.6e-3, observer: pd}
  observers: [{id: pd, kind: detector}]
  locks: [{id: lockD, combo: {lo: 1, pr: -1}, setpoint_rad: !pi_frac 0.5, observer: pd}]
  ports: []
  coef:
    - {segment: A1, driver: dT, beam: "*", effect: phase, value: !expr "44.0 * 10.5"}
    - {segment: A1, driver: dT, beam: "*", effect: group_delay, value: 0.0}
    - {segment: A1, driver: stretcher, beam: "*", effect: phase, value: 6.49}
    - {segment: A1, driver: stretcher, beam: "*", effect: group_delay, value: 0.0}
    - {segment: A2, driver: dT, beam: "*", effect: phase, value: !expr "44.0 * 0.5"}
    - {segment: A2, driver: dT, beam: "*", effect: group_delay, value: 0.0}
pipeline:
  - {block: pathtest_white, name: laser, params: {out: laser.dx1.phi, shape: [], unit: rad, sigma: 0.1}}
  - {block: pathtest_white, name: t1, params: {out: drv.A1.dT, shape: [], unit: K, sigma: 1.0e-4}}
  - {block: pathtest_white, name: t2, params: {out: drv.A2.dT, shape: [], unit: K, sigma: 1.0e-4}}
  - {block: pathtest_white, name: pzt, params: {out: drv.A1.stretcher, shape: [], unit: V, sigma: 1.0e-2}}
  - {block: segment_phase, name: sp, params: {topology: *topo, segments: all,
                                              toggles: {phase: true, group_delay: true, pm_axis_drift: true}}}
  - {block: beam_assembler, name: ba, params: {topology: *topo, beams: all, sources: {dx1: laser.dx1.phi},
                                               toggles: {emission_time: true, polarization: true}}}
  - {block: lock_error, name: le, params: {topology: *topo, locks: all}}
measure: {record: {slots: [beam.pr.phi, beam.lo.phi, lock.lockD.err], decim: 1}}
"""


@pytest.mark.numba
def test_yaml_preset_through_runner_and_round_trip():
    """A YAML config shares RunSpec.topology with the path blocks through a YAML anchor; it runs through
    the runner (reference == compiled) and round-trips (VT-F1 rule)."""
    from cranberry import config
    from cranberry.config.yaml_io import loads_yaml
    from cranberry.runner.run import run

    spec = config.load(loads_yaml(REFD_YAML))
    assert spec.pipeline.children[4].params["topology"] == spec.topology
    res = {m: run(spec, mode=m) for m in ("reference", "compiled")}
    s = res["reference"].series
    np.testing.assert_allclose(s["lock.lockD.err"], s["beam.lo.phi"] - s["beam.pr.phi"] - math.pi / 2, atol=1e-12)
    for k in s:
        np.testing.assert_allclose(res["compiled"].series[k], s[k], rtol=1e-12, atol=1e-12)
    assert config.load(loads_yaml(config.dump_yaml(spec))) == spec
