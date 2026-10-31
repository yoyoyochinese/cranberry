"""Unit tests of the ``metric_i`` / ``metric_ii`` probes and ``MetricsResult`` (conv §10, G13, WP-F4)."""

import json
import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
import cranberry.path.testing as pt
from cranberry.blocks.meas.metrics import (
    MetricsResult,
    PhaseStats,
    make_metric_i,
    make_metric_ii,
    metric_blocks,
    metric_name,
    wrap_quadrature,
)
from cranberry.blocks.meas.testing import path_blocks, single_port_topology, toy, two_port_topology
from cranberry.core.block import MissingParameterError

pytestmark = pytest.mark.numba


def single_pipe(amp, setpoint):
    topo = single_port_topology()
    return topo, [toy("toy_sine", "jit", out="drv.lo_seg.stretcher", unit="1", amp=amp, freq_hz=1e4),
                  toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
                  *pt.path_pipeline(topo, {"las": "laser.las.phi"}),
                  {"block": "metric_i", "name": "mi", "params": {"topology": topo, "setpoint_rad": setpoint,
                                                                "ports": {"0": {"squeezed": "sqz", "lo": "lo"}}}},
                  pt.record(["metric.mi.port0", "beam.lo.phi", "beam.sqz.phi"], 400)]


@pytest.mark.parametrize("compiled", [False, True])
def test_metric_i_value_wrap_and_accumulators(compiled):
    topo, pipe = single_pipe(amp=3.0, setpoint=0.25)
    rec, st, pa = pt.run_and_record(pipe, 1e6, 400, compiled=compiled)
    x = rec["beam.lo.phi"] - rec["beam.sqz.phi"] - 0.25
    d = rec["metric.mi.port0"]
    np.testing.assert_allclose(d, wrap_quadrature(x), atol=1e-12)
    assert np.all(d > -math.pi / 2) and np.all(d <= math.pi / 2) and np.ptp(x) > math.pi  # wrapping exercised
    acc = st["metric.mi.port0.acc"]
    np.testing.assert_allclose(acc, [400, d.sum(), (d * d).sum(), np.cos(2 * d).sum(), np.sin(2 * d).sum()],
                               rtol=1e-12, atol=1e-12)


def test_metric_compiled_equals_reference_and_setpoint_is_runtime():
    _, pipe = single_pipe(amp=0.2, setpoint=0.0)
    a, sa, _ = pt.run_and_record(pipe, 1e6, 400, compiled=False)
    b, sb, _ = pt.run_and_record(pipe, 1e6, 400, compiled=True)
    np.testing.assert_allclose(a["metric.mi.port0"], b["metric.mi.port0"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(sa["metric.mi.port0.acc"], sb["metric.mi.port0.acc"], rtol=1e-12)
    blk = make_metric_i({"topology": single_port_topology(), "ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                         "setpoint_rad": 0.1}, name="mi")
    assert blk.setpoint_param(0) == "mi.port0.setpoint_rad" and blk.toggles == {}
    assert blk.reads == ("beam.lo.phi", "beam.sqz.phi")
    assert blk.writes == ("metric.mi.port0", "metric.mi.port0.acc")


def test_port_selection():
    topo = two_port_topology()
    mi = make_metric_i({"topology": topo, "ports": "all", "setpoint_rad": {1: 0.5}}, name="mi")
    assert [(p.port, p.plus, p.minus) for p in mi.pairs] == [(1, "lo1", "sqz1")] and mi.setpoints == {1: 0.5}
    mii = make_metric_ii({"topology": topo, "ports": "all", "setpoint_rad": 0.0}, name="mii")
    assert [(p.port, p.plus, p.minus) for p in mii.pairs] == [(1, "sqz1", "sqz0")]
    assert mii.beams_of(mii.pairs[0]) == {"sqz1": 1.0, "sqz0": -1.0} and mii.squeezed_beam(mii.pairs[0]) == "sqz1"
    explicit = make_metric_ii({"topology": topo, "ports": {"0": {"squeezed": "sqz0"}, "1": {"squeezed": "sqz1"}},
                               "setpoint_rad": 0.0}, name="m2")
    assert [(p.port, p.plus, p.minus) for p in explicit.pairs] == [(1, "sqz1", "sqz0")]
    with pytest.raises(MissingParameterError):
        make_metric_i({"topology": topo, "ports": "all", "setpoint_rad": {2: 0.0}}, name="x")
    with pytest.raises(ValueError, match="LO and a BHD"):
        make_metric_i({"topology": topo, "ports": [0], "setpoint_rad": 0.0}, name="x")
    with pytest.raises(ValueError, match="observer"):  # sqz1 is at gbs, pr1 at det1
        make_metric_i({"topology": topo, "ports": {1: {"squeezed": "sqz1", "lo": "pr1"}}, "setpoint_rad": 0.0}, name="x")
    with pytest.raises(ValueError, match="reference port 0"):
        make_metric_ii({"topology": topo, "ports": {1: {"squeezed": "sqz1"}}, "setpoint_rad": 0.0}, name="x")
    with pytest.raises(MissingParameterError):
        make_metric_i({"topology": topo, "ports": "all"}, name="x")
    portless = make_metric_ii({"topology": single_port_topology(), "ports": "all", "setpoint_rad": 0.0}, name="n")
    assert portless.pairs == [] and portless.writes == ()  # G13: not applicable, nothing computed


def test_duplicate_metric_port_is_an_error():
    topo = two_port_topology()
    a = make_metric_i({"topology": topo, "ports": "all", "setpoint_rad": 0.0}, name="a")
    b = make_metric_i({"topology": topo, "ports": [1], "setpoint_rad": 0.0}, name="b")
    assert metric_name("metric_i", 1) == "i.port1" and metric_name("ii", 3) == "ii.port3"
    with pytest.raises(ValueError, match="i.port1"):
        metric_blocks([a, b])


def test_metrics_result_json_and_summary():
    from cranberry.blocks.meas.metrics import Efficiency, MetricI, Overlap

    st = PhaseStats.from_series([0.01, -0.01, 0.02])
    r = MetricsResult(t_s=np.arange(3.0), metric_i={0: MetricI(0, "bhd", np.zeros(3), 0.0, st, {"lo": 1, "sqz": -1})},
                      efficiency={0: Efficiency(0, None, 0.9)}, overlap={1: Overlap(1, np.array([0.99, 0.98]))},
                      noise_levels={"i.port0": {"r_minus_db": -10.0, "v_prime_t": np.zeros(3), "f": np.ones(2)}})
    js = json.loads(json.dumps(r.to_json()))
    assert js["metric_i"]["0"]["std_rad"] == pytest.approx(float(st.std))
    assert "v_prime_t" not in js["noise_levels"]["i.port0"] and js["noise_levels"]["i.port0"]["f"] == [1.0, 1.0]
    assert js["overlap"]["1"]["min"] == pytest.approx(0.98)
    s = r.summary()
    assert s["i.port0.rms_rad"] == pytest.approx(float(st.rms)) and s["i.port0.r_minus_db"] == -10.0
    assert s["eff.port0.eta_bhd"] == 0.9 and "eff.port0.eta_gbs" not in s
    assert r.metric_i[0].rms_rad == pytest.approx(float(st.rms))
    with pytest.raises(ValueError):
        MetricsResult(t_s=np.zeros(0), metric_ii_status="maybe")


def test_path_blocks_helper_has_no_topology_param():
    for b in path_blocks({"las": "laser.las.phi"}, locks=True):
        assert "topology" not in b["params"]
