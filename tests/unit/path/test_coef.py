"""WP-F6: coefficient values — constant, callable, Table, config !expr / !fn (G16-a)."""

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import Expr, Fn, loads_yaml
from cranberry.core.registry import register_fn
from cranberry.offline.tables import Table
from cranberry.path.coef import evaluate_coef
from cranberry.path.testing import beam_spec, seg_spec, step_spec, topo_spec
from cranberry.path.topology import Topology, TopologyError


@register_fn("f6test_kT_L")
def _kt_l(kT, L_m):
    return kT * L_m


def test_constant_and_callable():
    assert evaluate_coef(3) == 3.0 and evaluate_coef(2.5, x=1.0) == 2.5
    assert evaluate_coef(lambda L_m: 44.0 * L_m, L_m=2.0, other=1) == 88.0
    assert evaluate_coef(lambda **kw: kw["a"] + kw["b"], a=1.0, b=2.0) == 3.0
    with pytest.raises(ValueError, match="needs coordinate"):
        evaluate_coef(lambda L_m: L_m)
    with pytest.raises(TypeError):
        evaluate_coef(True)
    with pytest.raises(ValueError, match="non-finite"):
        evaluate_coef(float("nan"))


def test_table_expr_fn():
    tab = Table("k", {"T_k": np.array([290.0, 300.0])}, np.array([40.0, 50.0]), "rad/K")
    assert evaluate_coef(tab, T_k=295.0, unused=1) == pytest.approx(45.0)
    with pytest.raises(ValueError, match="needs coordinate"):
        evaluate_coef(tab)
    assert evaluate_coef(Expr("k * L", {"k": 44.0}, ("L",)), L=2.0) == pytest.approx(88.0)
    assert evaluate_coef(Fn("f6test_kT_L", {"kT": 44.0}), L_m=3.0) == pytest.approx(132.0)


def test_yaml_topology_with_callable_coefficient():
    text = """
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 1}
topology:
  segments:
    - {id: a, kind: pm, L_m: 63.0, n_g: 1.468, beta2_s2_per_m: 0.0, beta3_s3_per_m: 0.0, zone: z, n_sub: 1,
       spatial: lumped, pol_model: "off", pol_params: {}, drivers: {dT: uniform}}
  beams:
    - {id: sq, source: las, path: [{segment: a, direction: 1}], freq_offset_hz: 0.0, jones_in: [1, 0],
       power_w: 1.0e-3, observer: o}
  observers: [{id: o, kind: gbs_input}]
  locks: []
  ports: []
  coef:
    - {segment: a, driver: dT, beam: "*", effect: phase, value: !expr {source: "kT * L", args: [L], where: {kT: 44.0}},
       coords: {L: 63.0}}
    - {segment: a, driver: dT, beam: "*", effect: group_delay, value: !expr "39e-15 * 63.0"}
pipeline: []
"""
    spec = config.load(loads_yaml(text))
    t = Topology.from_spec(spec.topology)
    assert t.coef_value(("a", "dT", "sq", "phase")) == pytest.approx(44.0 * 63.0)
    assert t.coef_value(("a", "dT", "sq", "group_delay")) == pytest.approx(2.457e-12)


def test_bad_coefficient_reported_at_parse():
    segs = [seg_spec("a", 1.0, drivers={"dT": "uniform"})]
    coef = [{"segment": "a", "driver": "dT", "beam": "*", "effect": "phase", "value": lambda L_m: L_m},
            {"segment": "a", "driver": "dT", "beam": "*", "effect": "group_delay", "value": 0.0}]
    with pytest.raises(TopologyError, match="needs coordinate"):
        Topology.from_spec(topo_spec(segs, [beam_spec("b", [step_spec("a")], "o")], ["o"], coef))
    coef[0]["coords"] = {"L_m": 2.0}
    t = Topology.from_spec(topo_spec(segs, [beam_spec("b", [step_spec("a")], "o")], ["o"], coef))
    assert t.coef_value(("a", "dT", "b", "phase")) == 2.0
