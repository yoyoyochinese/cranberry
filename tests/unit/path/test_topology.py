"""WP-F6: Topology parsing, validation, timing, registry helpers (PLAN §2.3.1, conv §9)."""

import copy

import numpy as np
import pytest

from cranberry.offline.band_ledger import BandLedger
from cranberry.offline.loss_registry import LossElement, LossRegistry
from cranberry.path.coef import CoefKey
from cranberry.path.testing import beam_spec, coef_spec, seg_spec, step_spec, topo_spec
from cranberry.path.topology import C_M_PER_S, Topology, TopologyError, as_topology


def base():
    segs = [seg_spec("a", 10.0, n_sub=4, drivers={"dT": "uniform", "stretcher": "actuator"}, actuator_sub=1),
            seg_spec("b", 2.0, n_sub=2, drivers={"strain": "per_sub"}, spatial="distributed")]
    beams = [beam_spec("sq", [step_spec("a", 1), step_spec("b", 1)], "gbs"),
             beam_spec("pr", [step_spec("b", -1), step_spec("a", -1)], "det")]
    coef = (coef_spec("a", "dT", "*", 440.0, 4e-13) + coef_spec("a", "stretcher", "*", 6.49)
            + coef_spec("b", "strain", "*", 1e7))
    locks = [{"id": "l1", "combo": {"pr": 1.0}, "setpoint_rad": 0.5, "observer": "det"}]
    return topo_spec(segs, beams, ["gbs", "det"], coef, locks)


def test_parse_and_queries():
    t = Topology.from_spec(base())
    assert set(t.segments) == {"a", "b"} and t.beams_on("a") == ("sq", "pr")
    assert t.segments["a"].drivers == {"dT": ("uniform",), "stretcher": ("point", 1)}
    assert t.step("pr", "a").subs == (3, 2, 1, 0)
    assert t.coef_value(CoefKey("a", "dT", "pr", "phase")) == 440.0
    assert not t.has_ports and t.locks["l1"].setpoint_rad == 0.5
    np.testing.assert_allclose(t.segments["a"].weights("dT"), [0.25] * 4)
    np.testing.assert_allclose(t.segments["a"].weights("stretcher"), [0, 1, 0, 0])
    assert t.segments["b"].weights("strain") is None
    assert as_topology(t) is t


def test_timing():
    t = Topology.from_spec(base())
    tm = t.timing("sq")
    ta, tb = 1.468 * 10 / C_M_PER_S, 1.468 * 2 / C_M_PER_S
    assert tm.tau_s == pytest.approx(ta + tb, rel=1e-15)
    assert tm.steps[0].after_s == pytest.approx(tb, rel=1e-15) and tm.steps[1].after_s == 0.0
    d = ta / 4
    np.testing.assert_allclose(tm.steps[0].to_exit_s, [3.5 * d, 2.5 * d, 1.5 * d, 0.5 * d], rtol=1e-14)
    # counter-propagating: same subsegments, traversed in the opposite order
    tp = t.timing("pr")
    assert tp.steps[1].subs == (3, 2, 1, 0)
    np.testing.assert_allclose(tp.steps[1].to_exit_s, tm.steps[0].to_exit_s, rtol=1e-14)


def test_partial_traversal():
    spec = base()
    spec["beams"][0]["path"][0] = step_spec("a", 1, sub_in=1, sub_out=2)
    t = Topology.from_spec(spec)
    st = t.timing("sq").steps[0]
    assert st.subs == (1, 2) and st.tau_s == pytest.approx(0.5 * 1.468 * 10 / C_M_PER_S)


def test_wildcard_override_and_python_mapping():
    spec = base()
    spec["coef"] += [{"segment": "a", "driver": "dT", "beam": "pr", "effect": "phase", "value": 400.0}]
    t = Topology.from_spec(spec)
    assert t.coef_value(CoefKey("a", "dT", "pr", "phase")) == 400.0
    assert t.coef_value(CoefKey("a", "dT", "sq", "phase")) == 440.0
    spec2 = base()
    spec2["coef"] = {CoefKey(k["segment"], k["driver"], k["beam"], k["effect"]): k["value"] for k in spec["coef"][:-1]}
    t2 = Topology.from_spec(spec2)
    assert t2.coef_value(CoefKey("b", "strain", "sq", "phase")) == 1e7


@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda s: s["coef"].pop(0), r"missing coefficients.*\(a, dT, sq, phase\)"),
        (lambda s: s["coef"].append({"segment": "a", "driver": "dT", "beam": "zz", "effect": "phase", "value": 1}),
         "does not traverse"),
        (lambda s: s["coef"].append({"segment": "a", "driver": "eom", "beam": "sq", "effect": "phase", "value": 1}),
         "not declared"),
        (lambda s: s["beams"][0]["path"].append(step_spec("a", 1)), "twice"),
        (lambda s: s["beams"][0]["path"].__setitem__(0, step_spec("a", 0)), r"direction must be \+1 or -1"),
        (lambda s: s["beams"][0]["path"].__setitem__(0, step_spec("a", -1, sub_in=0, sub_out=2)), "sub_in >= sub_out"),
        (lambda s: s["segments"][0].__setitem__("pol_model", "pm_diagonal"), "pm_diagonal pol_params"),
        (lambda s: s["segments"][1].__setitem__("spatial", "lumped"), "lumped forbids per_sub"),
        (lambda s: s["segments"][1]["drivers"].__setitem__("strain", "uniform"), "distributed needs per_sub"),
        (lambda s: s["locks"][0].__setitem__("combo", {"sq": 1.0}), "not evaluated at the lock observer"),
        (lambda s: s["ports"].append({"index": 1, "squeezed": "sq", "lo": None, "probe": "pr",
                                      "gbs_observer": "gbs", "bhd_observer": None}), "reference port 0"),
        (lambda s: s["segments"][0].pop("zone"), r"missing required key\(s\) \['zone'\]"),
        (lambda s: s["segments"][0]["drivers"].__setitem__("stretcher", {"point": 9}), "point must be"),
        (lambda s: s["segments"][0].__setitem__("bogus", 1), "unknown key"),
        (lambda s: s.pop("ports"), "missing required key"),
    ],
)
def test_validation_errors(mutate, match):
    spec = copy.deepcopy(base())
    mutate(spec)
    with pytest.raises(TopologyError, match=match):
        Topology.from_spec(spec)


def test_pm_axis_and_jones_requirements():
    pm = {"K_T_rad_per_m_K": -2.05, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.36e-12, "axis_drivers": ["dT"]}
    segs = [seg_spec("p", 1.0, drivers={"dT": "uniform"}, pol_model="pm_diagonal", pol_params=pm),
            seg_spec("j", 1.0, pol_model="jones", pol_params={"J": [[1, 0], [0, 1]]})]
    beams = [beam_spec("b", [step_spec("p", 1, axis="slow"), step_spec("j", 1)], "o", jones_out=[1, 0])]
    good = topo_spec(segs, beams, ["o"], coef_spec("p", "dT", "*", 1.0))
    t = Topology.from_spec(good)
    assert t.step("b", "p").axis == "slow" and t.segments["j"].pol_params["J"].dtype == np.complex128
    for mut, match in [
        (lambda s: s["beams"][0]["path"][0].pop("axis"), "needs axis"),
        (lambda s: s["beams"][0]["path"][1].__setitem__("axis", "slow"), "only used on pm_diagonal"),
        (lambda s: s["beams"][0].pop("jones_out"), "jones_out"),
        (lambda s: s["beams"][0].__setitem__("jones_in", [1, 1]), "unit norm"),
        (lambda s: s["segments"][1].__setitem__("pol_params", {}), "exactly one of J"),
    ]:
        spec = copy.deepcopy(good)
        mut(spec)
        with pytest.raises(TopologyError, match=match):
            Topology.from_spec(spec)


def test_pm_axis_ok():
    pm = {"K_T_rad_per_m_K": -2.05, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.36e-12, "axis_drivers": ["dT"]}
    segs = [seg_spec("p", 1.0, drivers={"dT": "uniform"}, pol_model="pm_diagonal", pol_params=pm)]
    t = Topology.from_spec(topo_spec(segs, [beam_spec("b", [step_spec("p", 1, axis="fast")], "o")], ["o"],
                                     coef_spec("p", "dT", "*", 1.0)))
    assert t.step("b", "p").axis == "fast" and t.segments["p"].pol_params["axis_drivers"] == ("dT",)


def test_registry_helpers():
    t = Topology.from_spec(base())
    reg = LossRegistry()
    reg.add(LossElement("a.il", 0.9, "component", frozenset({"gbs", "probe"}), "UNVERIFIED", True, segment="a"))
    reg.add(LossElement("b.il", 0.8, "component", frozenset({"gbs", "probe"}), "UNVERIFIED", True, segment="b"))
    paths = t.register_loss_paths(reg)
    assert [p.id for p in paths] == ["sq", "pr"] and paths[0].observers == ("gbs",)
    assert reg.eta("sq", "gbs") == pytest.approx(0.72)
    led = BandLedger()
    t.register_band_beams(led)
    assert led.beams == {"sq": ("a", "b"), "pr": ("b", "a")}


def test_setpoint_alias_and_ports():
    spec = base()
    spec["locks"][0] = {"id": "l1", "combo": {"pr": 1.0}, "setpoint": 0.25, "observer": "det"}
    spec["ports"] = [{"index": 0, "squeezed": "sq", "lo": None, "probe": "pr", "gbs_observer": "gbs",
                      "bhd_observer": None}]
    t = Topology.from_spec(spec)
    assert t.locks["l1"].setpoint_rad == 0.25 and t.has_ports and t.ports[0].probe == "pr"
