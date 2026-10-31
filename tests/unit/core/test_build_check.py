"""core.pipeline.build / check_pipeline: the build-time checks of PLAN §2.5 (G09), toggle overrides,
RNG stream per leaf (G11)."""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401  (registers toy_* kinds)
from cranberry.core.block import BlockSpec
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import (
    PipelineCheckError,
    SeqSpec,
    apply_toggle_overrides,
    build,
    iter_leaves,
    spec_from_obj,
)
from cranberry.core.rng import stream
from cranberry.core.state import Params, State
from cranberry.core.testing import toy


def delay(name, inp, out, tau_s=0.0, method="linear"):
    return {"block": "delay_buffer", "name": name, "params": {"in": inp, "out": out, "tau_s": tau_s, "method": method}}


def gain(name, inp, out, unit="V", g=1.0):
    return toy("toy_gain", name, in_=inp, out=out, unit=unit, gain=g)


def loop_obj(with_delay=True):
    """plant(integrator) -> sensor -> pid -> [delay] -> act -> (back edge) plant."""
    act_in = "delay.d.out" if with_delay else "ctrl.pid.u"
    obj = [
        toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="rad", k=1.0),
        gain("sensor", "psi.plant.x", "v.y", unit="V"),
        toy("toy_pid", "pid", in_="v.y", out="ctrl.pid.u", unit="V", kp=0.1, ki=0.0, kd=0.0, setpoint=0.0),
    ]
    if with_delay:
        obj.append(delay("d", "ctrl.pid.u", "delay.d.out", 2e-6))
    obj.append(gain("act", act_in, "act.u"))
    return obj


def try_build(obj, external=(), pre=()):
    lay = SlotLayout(1.0e6)
    for s in pre:
        lay.add(s)
    return build(spec_from_obj(obj), lay, 1, external=external), lay


def violations(obj, **kw):
    with pytest.raises(PipelineCheckError) as ei:
        try_build(obj, **kw)
    return ei.value.violations


def test_feedback_through_delay_buffer_builds():
    blk, lay = try_build(loop_obj(True))
    assert lay.frozen and blk.layout is lay
    assert [leaf.name for leaf in iter_leaves(blk)] == ["plant", "sensor", "pid", "d", "act"]


def test_zero_delay_loop_rejected():
    v = violations(loop_obj(False))
    assert len(v) == 1 and "zero-delay feedback loop" in v[0] and "'plant'" in v[0]


def test_read_before_write_rejected():
    # 'late' reads v.b which is written later and does not depend on 'late': an ordering error
    obj = [toy("toy_sine", "src", out="v.a", unit="V", amp=1.0, freq_hz=1.0),
           gain("late", "v.b", "v.c"), gain("b", "v.a", "v.b")]
    v = violations(obj)
    assert len(v) == 1 and v[0].startswith("read-before-write") and "'late'" in v[0]
    try_build([obj[0], obj[2], obj[1]])  # correct order builds


def test_missing_writer_and_external():
    obj = [gain("g", "v.in", "v.out")]
    v = violations(obj, pre=("v.in",))
    assert "no block writes" in v[0]
    try_build(obj, external=("v.in",), pre=("v.in",))
    v = violations(obj)  # not even registered
    assert any("unregistered slot 'v.in'" in m for m in v)


def test_parallel_disjoint_writes_and_no_cross_reads():
    src = toy("toy_sine", "src", out="v.a", unit="V", amp=1.0, freq_hz=1.0)
    v = violations([src, {"parallel": [gain("p", "v.a", "v.x"), gain("q", "v.x", "v.y")], "name": "par"}])
    assert any("reads ['v.x'] written by sibling" in m for m in v)
    try_build([src, {"parallel": [gain("p", "v.a", "v.x"), gain("q", "v.a", "v.y")]}])


def test_parallel_write_overlap_detected_by_check():
    # two writers of the same slot need one block that does not register it; use a toy_sum with out
    # pre-registered externally so the layout does not reject the duplicate add
    from cranberry.core.pipeline import Parallel, check_pipeline
    from cranberry.core.testing import make_toy_gain

    a = make_toy_gain({"in": "v.a", "out": "v.x", "unit": "V", "gain": 1.0, "toggles": {"enabled": True}}, name="a")
    b = make_toy_gain({"in": "v.a", "out": "v.y", "unit": "V", "gain": 1.0, "toggles": {"enabled": True}}, name="b")
    b.writes = ("v.x",)  # declare an overlapping write
    lay = SlotLayout(1.0)
    lay.add("v.a")
    par = Parallel([a, b], name="par")
    for leaf in iter_leaves(par):
        leaf.init(lay, stream(0, leaf.name))
    with pytest.raises(PipelineCheckError, match="both write"):
        check_pipeline(par, lay, external=("v.a",))


def test_unit_tag_mismatch():
    obj = [toy("toy_sine", "src", out="v.a", unit="V", amp=1.0, freq_hz=1.0),
           {"block": "delay_buffer", "name": "d", "params": {"in": "v.a", "out": "delay.d.out", "tau_s": 0.0,
                                                             "method": "linear", "unit": "rad"}}]
    v = violations(obj)
    assert any("expects unit 'rad' for slot 'v.a'" in m for m in v)


def test_duplicate_leaf_names():
    obj = [toy("toy_sine", "s", out="v.a", unit="V", amp=1.0, freq_hz=1.0),
           toy("toy_sine", "s", out="v.b", unit="V", amp=1.0, freq_hz=1.0)]
    assert "duplicate leaf name 's'" in violations(obj)[0]


def test_all_violations_reported_together():
    obj = [gain("g1", "v.nope", "v.o1"), gain("late", "v.o3", "v.o2"), gain("w", "v.o1", "v.o3")]
    v = violations(obj, pre=("v.nope",))
    assert len(v) == 2


def test_each_leaf_gets_its_own_stream():
    obj = [toy("toy_noise_ar1", "n1", out="v.n1", unit="V", a=0.0, b=1.0),
           toy("toy_noise_ar1", "n2", out="v.n2", unit="V", a=0.0, b=1.0)]
    blk, lay = try_build(obj)
    st, p = State(lay), Params(lay)
    leaves = list(iter_leaves(blk))
    for leaf in leaves:
        leaf.prepare(st, p, 0, 16)
    np.testing.assert_array_equal(st["rng.n1.buf"][:16], stream(1, "n1").standard_normal(16))
    np.testing.assert_array_equal(st["rng.n2.buf"][:16], stream(1, "n2").standard_normal(16))


def test_apply_toggle_overrides():
    spec = spec_from_obj([toy("toy_sine", "a.src", out="v.a", unit="V", amp=1.0, freq_hz=1.0),
                          {"stack": {"stages": [gain("c1", "v.a", "ctrl.c1.u")],
                                     "actuators": [gain("x1", "ctrl.c1.u", "act.x1")]}}])
    out = apply_toggle_overrides(spec, {"a.src.enabled": False, "x1.enabled": False})
    leaves = {s.name: s for s in _leaf_specs(out)}
    assert leaves["a.src"].params["toggles"] == {"enabled": False}
    assert leaves["x1"].params["toggles"] == {"enabled": False}
    assert leaves["c1"].params["toggles"] == {"enabled": True}
    assert apply_toggle_overrides(spec, {}) is spec
    with pytest.raises(KeyError):
        apply_toggle_overrides(spec, {"nope.enabled": False})


def _leaf_specs(s):
    if isinstance(s, BlockSpec):
        yield s
    elif isinstance(s, SeqSpec) or hasattr(s, "children"):
        for c in s.children:
            yield from _leaf_specs(c)
    else:
        for c in (*s.stages, *s.actuators):
            yield from _leaf_specs(c)
