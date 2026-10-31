"""Block contract, registry and pipeline combinators (PLAN §2.4–§2.5, G09)."""

from collections.abc import Mapping
from typing import Any

import pytest

from cranberry.core import registry
from cranberry.core.block import BaseBlock, Block, BlockSpec, MissingParameterError
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import (
    Parallel,
    ParallelSpec,
    Seq,
    SeqSpec,
    Stack,
    StackSpec,
    from_spec,
    iter_leaves,
    spec_from_obj,
    spec_to_obj,
)
from cranberry.core.rng import stream
from cranberry.core.state import Params, State


class Gain(BaseBlock):
    """Test block: out = gain * in (runtime param), with an 'enabled' toggle."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        super().__init__("test_gain", name, params, reads=(params["in"],), writes=(params["out"],),
                         runtime_params=("gain",))

    def init(self, layout, rng):
        super().init(layout, rng)
        if self.params["out"] not in layout:
            layout.add(self.params["out"])

    def __call__(self, state, params):
        g = float(params[self.param_name("gain")]) if self.toggles["enabled"] else 1.0
        state[self.params["out"]] = g * state[self.params["in"]]
        return state


@registry.register_block("test_gain")
def make_test_gain(params, *, name=None):
    return Gain(params, name)


def gain(name, inp, out, g=2.0, enabled=True):
    return make_test_gain({"in": inp, "out": out, "gain": g, "toggles": {"enabled": enabled}}, name=name)


def setup(block, extra_slots=("v.in",)):
    lay = SlotLayout(1.0e3)
    for s in extra_slots:
        lay.add(s)
    for leaf in iter_leaves(block):
        leaf.init(lay, stream(1, leaf.name))
    lay.freeze()
    block.bind(lay)
    return State(lay), Params(lay)


def test_blockspec_round_trip_and_validation():
    s = BlockSpec("pid", "lock1.pid", {"kp": 1.0}, rate_div=2)
    assert s.to_obj() == {"block": "pid", "name": "lock1.pid", "params": {"kp": 1.0}, "rate_div": 2}
    assert BlockSpec.from_obj(s.to_obj()) == s
    assert BlockSpec.from_obj({"block": "pid"}) == BlockSpec("pid", "pid", {})
    for bad in [dict(kind="Pid", name="x"), dict(kind="pid", name=""), dict(kind="pid", name="x", rate_div=0)]:
        with pytest.raises(ValueError):
            BlockSpec(**bad)
    with pytest.raises(ValueError):
        BlockSpec.from_obj({"block": "pid", "typo": 1})


def test_factory_and_registry():
    b = registry.make_block(BlockSpec("test_gain", "g1", {"in": "v.in", "out": "ctrl.g1.u", "gain": 3.0,
                                                          "toggles": {"enabled": True}}, rate_div=4))
    assert isinstance(b, Block)
    assert b.rate_div == 4 and b.spec.rate_div == 4 and b.name == "g1" and b.kind == "test_gain"
    assert b.reads == ("v.in",) and b.writes == ("ctrl.g1.u",) and b.runtime_params == ("gain",)
    assert b.losses() == [] and b.bands() == [] and b.kernel is None
    with pytest.raises(KeyError):
        registry.get_factory("nope")
    with pytest.raises(ValueError):
        registry.register_block("seq")
    with pytest.raises(ValueError):
        registry.register_block("test_gain")(lambda params, name=None: None)


def test_d34_no_defaults():
    with pytest.raises(MissingParameterError):
        make_test_gain({"in": "v.in", "out": "ctrl.x.u", "toggles": {"enabled": True}})  # gain missing
    with pytest.raises(MissingParameterError):
        make_test_gain({"in": "v.in", "out": "ctrl.x.u", "gain": 1.0})  # toggle missing
    with pytest.raises(ValueError):
        make_test_gain({"in": "v.in", "out": "ctrl.x.u", "gain": 1.0, "toggles": {"enabled": True, "x": False}})


def test_reference_step_in_place_and_toggle_off_is_identity():
    seq = gain("a", "v.in", "ctrl.a.u", 2.0) >> gain("b", "ctrl.a.u", "ctrl.b.u", 5.0, enabled=False)
    st, p = setup(seq)
    st["v.in"] = 1.5
    out = seq(st, p)
    assert out is st
    assert float(st["ctrl.a.u"]) == 3.0 and float(st["ctrl.b.u"]) == 3.0
    p["a.gain"] = 4.0  # runtime parameter change without rebuilding
    seq(st, p)
    assert float(st["ctrl.b.u"]) == 6.0


def test_rate_div_hold():
    blk = registry.make_block(BlockSpec("test_gain", "slow", {"in": "v.in", "out": "ctrl.slow.u", "gain": 1.0,
                                                              "toggles": {"enabled": True}}, rate_div=2))
    seq = Seq([blk])
    st, p = setup(seq)
    got = []
    for n in range(4):
        st["n"] = n
        st["v.in"] = float(n + 1)
        seq(st, p)
        got.append(float(st["ctrl.slow.u"]))
    assert got == [1.0, 1.0, 3.0, 3.0]


def test_seq_parallel_operators_flatten_and_union():
    a, b, c = gain("a", "v.in", "ctrl.a.u"), gain("b", "v.in", "ctrl.b.u"), gain("c", "v.in", "ctrl.c.u")
    s = a >> b >> c
    assert isinstance(s, Seq) and [x.name for x in s.children] == ["a", "b", "c"]
    par = a | b | c
    assert isinstance(par, Parallel) and len(par.children) == 3
    assert par.writes == ("ctrl.a.u", "ctrl.b.u", "ctrl.c.u") and par.reads == ("v.in",)
    assert set(par.toggles) == {"a.enabled", "b.enabled", "c.enabled"}


def test_stack_expansion_matches_t16():
    c = [gain(f"c{k}", "v.in" if k == 1 else f"ctrl.c{k-1}.u", f"ctrl.c{k}.u") for k in (1, 2, 3)]
    a = [gain(f"a{k}", f"ctrl.c{k}.u", f"act.a{k}.x") for k in (1, 2, 3)]
    st = Stack(c, a, name="refC")
    e = st.expand()
    # c1 >> (a1 | (c2 >> (a2 | (c3 >> a3))))
    assert isinstance(e, Seq) and e.children[0] is c[0]
    p1 = e.children[1]
    assert isinstance(p1, Parallel) and p1.children[0] is a[0]
    s2 = p1.children[1]
    assert s2.children[0] is c[1] and s2.children[1].children[0] is a[1]
    s3 = s2.children[1].children[1]
    assert [x.name for x in s3.children] == ["c3", "a3"]
    assert [x.name for x in iter_leaves(st)] == ["c1", "a1", "c2", "a2", "c3", "a3"]
    state, params = setup(st)
    state["v.in"] = 1.0
    st(state, params)
    assert float(state["act.a3.x"]) == 16.0  # 2^4 through c1, c2, c3, a3


def test_pipeline_spec_obj_round_trip_and_from_spec():
    obj = [
        {"block": "test_gain", "name": "g", "params": {"in": "v.in", "out": "ctrl.g.u", "gain": 1.0,
                                                        "toggles": {"enabled": True}}},
        {"parallel": [{"block": "test_gain", "name": "h", "params": {"in": "ctrl.g.u", "out": "ctrl.h.u",
                                                                       "gain": 2.0, "toggles": {"enabled": True}}}],
         "name": "par"},
        {"stack": {"stages": [{"block": "test_gain", "name": "s1", "params": {"in": "v.in", "out": "ctrl.s1.u",
                                                                             "gain": 1.0, "toggles": {"enabled": True}}}],
                   "actuators": [{"block": "test_gain", "name": "x1", "rate_div": 2,
                                  "params": {"in": "ctrl.s1.u", "out": "act.x1.y", "gain": 1.0,
                                             "toggles": {"enabled": True}}}]}},
        {"seq": [], "name": "empty_named"},
    ]
    spec = spec_from_obj(obj)
    assert isinstance(spec, SeqSpec) and isinstance(spec.children[1], ParallelSpec)
    assert isinstance(spec.children[2], StackSpec)
    assert spec_to_obj(spec) == obj
    assert spec_from_obj(spec_to_obj(spec)) == spec
    built = from_spec(SeqSpec(spec.children[:3]))
    assert built.spec == SeqSpec(spec.children[:3])
    with pytest.raises(ValueError):
        spec_from_obj({"stack": [obj[0]]})
    with pytest.raises(ValueError):
        spec_from_obj({"what": 1})


def test_composite_rejects_non_blocks_and_empty():
    with pytest.raises(TypeError):
        Seq([object()])
    with pytest.raises(ValueError):
        Parallel([])
    with pytest.raises(ValueError):
        Stack([gain("a", "v.in", "ctrl.a.u")], [])


def test_framework_stubs_name_their_wp():
    registry.ensure_loaded()
    assert {"delay_buffer", "segment_phase", "beam_assembler", "metric_i", "metric_ii"} <= set(registry.registered_kinds())
    with pytest.raises(MissingParameterError):  # delay_buffer is implemented (WP-F1); D34: params required
        registry.make_block(BlockSpec("delay_buffer", "d", {}))
    with pytest.raises(NotImplementedError, match="WP-F6"):
        registry.make_block(BlockSpec("segment_phase", "sp", {}))
