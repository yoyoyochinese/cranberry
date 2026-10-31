"""Smoke tests of ``runner.ablate``, ``runner.decompose`` (argument checks) and ``runner.convergence``
(WP-F4; PLAN §2.8, §6.4, G13)."""

import math

import pytest

import cranberry.core.testing  # noqa: F401
import cranberry.path.testing as pt
from cranberry.blocks.meas.testing import (
    TWO_PORT_GROUPS,
    metric_cfg,
    path_blocks,
    refa_source,
    run_spec,
    single_port_topology,
    toy,
    two_port_loop_spec,
)
from cranberry.runner.ablate import ablate, ablation_impact, list_toggles
from cranberry.runner.convergence import _double_topology, convergence, double_n_sub
from cranberry.runner.decompose import decompose_metric_ii

pytestmark = pytest.mark.numba
SMALL = dict(n_samples=1 << 13, record_decim=8)


def test_list_toggles_and_ablate_impact():
    spec = two_port_loop_spec(**SMALL)
    keys = list_toggles(spec)
    assert {"n_sq.enabled", "n_cp.enabled", "sp.phase", "ba.emission_time"} <= set(keys)
    assert "sp.pm_axis_drift" not in keys and "sp.pm_axis_drift" in list_toggles(spec, only_on=False)
    res = ablate(spec, ["n_sq.enabled", "n_cp.enabled"])
    assert list(res) == ["baseline", "n_sq.enabled", "n_cp.enabled"]
    imp = ablation_impact(res, ["ii.port1.std_rad", "i.port1.std_rad"])
    for t in ("n_sq.enabled", "n_cp.enabled"):
        row = imp[t]["ii.port1.std_rad"]
        assert row["delta"] < 0 and row["ablated"] == pytest.approx(row["baseline"] + row["delta"])
    full = ablation_impact(res)
    assert "lock.L1.err.std" in full["n_sq.enabled"]
    with pytest.raises(ValueError):
        ablate(spec, ["n_sq.enabled", "n_sq.enabled"])
    with pytest.raises(KeyError):
        ablation_impact({"x": res["baseline"]})


def test_decompose_argument_checks_and_portless():
    spec = two_port_loop_spec(**SMALL)
    with pytest.raises(ValueError, match="exactly the keys"):
        decompose_metric_ii(spec, groups={"common_points": []})
    with pytest.raises(ValueError, match="groups"):
        decompose_metric_ii(spec, groups={**TWO_PORT_GROUPS, "squeezed_only": ["n_cp.enabled"]})
    with pytest.raises(ValueError, match="slot, abs_max"):
        decompose_metric_ii(spec, validity=[{"slot": "x"}])
    portless = refa_metric_spec()
    with pytest.raises(KeyError, match="decompose"):
        decompose_metric_ii(portless)
    empty = {c: [] for c in TWO_PORT_GROUPS}
    assert decompose_metric_ii(portless, groups=empty, validity=[]) == {}
    # n_lock / n_las left out of the groups stay on in every run, so they are counted in both component
    # runs: the linearity check flags the incomplete grouping
    groups = {**empty, "squeezed_only": ["n_sq.enabled"], "common_points": ["n_cp.enabled"]}
    c = decompose_metric_ii(spec, groups=groups, validity=[])[1]
    assert c.variances_rad2["locked_path_residual"] == 0.0 and c.linearity_ratio > 1.05
    assert c.status == "decomposition_invalid" and "linearity" in c.reason


def refa_metric_spec(n=4000, f_sim=1.0e6, topo=None):
    pipe = [toy("toy_sine", "jit", out="drv.lo_seg.stretcher", unit="1", amp=0.01, freq_hz=1e3),
            toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
            *path_blocks({"las": "laser.las.phi"}),
            {"block": "metric_i", "name": "mi", "params": {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                                                          "setpoint_rad": 0.0}}]
    return run_spec(pipe, f_sim_hz=f_sim, n_samples=n, topology=topo or single_port_topology(),
                    offline={"blocks": [refa_source()]}, chunk_samples=1000,
                    measure={"record": {"slots": ["beam.lo.phi"], "decim": 1},
                             "metrics": {"i": metric_cfg("refA.sqz")},
                             "convergence": {"quantities": {"i.port0.std_rad": 1e-4, "i.port0.r_minus_db": 0.03},
                                             "seeds": [1, 2, 3]}})


def test_convergence_checks():
    from cranberry.config import FrozenExpressionsWarning

    with pytest.warns(FrozenExpressionsWarning):  # RunSpec(...) built in Python: no unresolved tree (WP-F2c)
        rep = convergence(refa_metric_spec())
    assert rep.checks["chunk"]["passed"] is True and rep.checks["chunk"]["values"]["chunk_samples"] == [1000, 334]
    f = rep.checks["f_sim"]
    assert f["passed"] is True and set(f["values"]) == {"x2", "x4"}
    assert f["values"]["x4"]["i.port0.std_rad"]["reference"] == pytest.approx(0.01 / math.sqrt(2), rel=1e-3)
    s = rep.checks["seed"]
    assert s["passed"] is False and s["values"]["seeds"] == [1, 2, 3]  # 10–30 seeds required to pass
    assert s["values"]["i.port0.std_rad"]["std"] == pytest.approx(0.0, abs=1e-15)  # deterministic jitter
    assert rep.checks["n_sub"]["passed"] is True
    assert rep.passed is False
    with pytest.raises(ValueError):
        convergence(refa_metric_spec(), ["mesh"])
    with pytest.raises(KeyError, match="seeds"):
        convergence(refa_metric_spec(), ["seed"], seeds=[])


def test_double_n_sub_remaps_indices():
    segs = [pt.seg_spec("c", 10.0, n_sub=4, drivers={"strain": {"point": 3}, "tec": {"subs": [0, 1]},
                                                     "stretcher": "actuator"}, actuator_sub=2, spatial="mixed")]
    beams = [pt.beam_spec("f", [pt.step_spec("c", 1, sub_in=1, sub_out=2)], "o1"),
             pt.beam_spec("b", [pt.step_spec("c", -1)], "o2")]
    obs = [{"id": "o1", "kind": "detector", "segment": "c", "sub": 2}, {"id": "o2", "kind": "detector", "segment": "c",
                                                                         "sub": 0}]
    coef = pt.coef_spec("c", "strain", "*", 1.0) + pt.coef_spec("c", "tec", "*", 1.0) + pt.coef_spec("c", "stretcher", "*", 1.0)
    topo = pt.topo_spec(segs, beams, obs, coef)
    d = _double_topology(topo)
    s = d["segments"][0]
    assert s["n_sub"] == 8 and s["actuator_sub"] == 4
    assert s["drivers"] == {"strain": {"point": 6}, "tec": {"subs": [0, 1, 2, 3]}, "stretcher": "actuator"}
    assert d["beams"][0]["path"][0]["sub_in"] == 2 and d["beams"][0]["path"][0]["sub_out"] == 5
    assert d["observers"][0]["sub"] == 5 and d["observers"][1]["sub"] == 0
    from cranberry.path.topology import Topology

    Topology.from_spec(d)  # still a valid topology
    assert topo["segments"][0]["n_sub"] == 4  # input untouched
    spec = double_n_sub(refa_metric_spec())
    assert spec.topology["segments"][0]["n_sub"] == 2
