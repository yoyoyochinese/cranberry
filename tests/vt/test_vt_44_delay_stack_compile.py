"""VT-44 (framework part, WP-F1): compiled loop = pure-Python reference (relative 1e−12, PLAN §6.3),
Stack = its Seq/Parallel expansion, Stack = parallel F_k = C_k⋯C_1 in the linear range (T16 §4.4),
kernel reuse across builds with the same structure, multirate hold in the compiled loop.

Delay-buffer accuracy is in tests/unit/core/test_delay.py; actuator blocks are WP-B4's part of VT-44;
f_sim convergence is in test_vt_f2_chunk_invariance.py.
"""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.block import BaseBlock
from cranberry.core.compile import (
    NotCompilableError,
    compile_loop,
    loop_cache_size,
    run_compiled,
    run_reference,
)
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import Seq, build, iter_leaves, spec_from_obj
from cranberry.core.rates import PIECE_SAMPLES
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba
FS = 1.0e6


def pid(name, inp, kp, ki, kd=0.0, rate_div=1):
    return toy("toy_pid", name, rate_div, in_=inp, out=f"ctrl.{name}.u", unit="V", kp=kp, ki=ki, kd=kd, setpoint=0.0)


def gain(name, inp, out, g):
    return toy("toy_gain", name, in_=inp, out=out, unit="V", gain=g)


def nontrivial_obj():
    """Two noise sources, a sine, a sum, a decimated PID, a fractional Thiran delay, a saturating
    actuator, an integrator plant (feedback via back edge), plus a Stack and a Parallel."""
    return [
        toy("toy_noise_ar1", "n1", out="v.n1", unit="V", a=0.999, b=0.05),
        toy("toy_noise_ar1", "n2", 3, out="v.n2", unit="V", a=0.9, b=0.2),
        toy("toy_sine", "dist", out="v.d", unit="V", amp=0.3, freq_hz=2.1e3),
        toy("toy_sum", "sens", ins=["psi.plant.x", "v.n1", "v.n2", "v.d"], gains=[1.0, 1.0, 0.5, 1.0],
            out="v.y", unit="V"),
        pid("pid", "v.y", kp=0.4, ki=3e3, kd=2e-7, rate_div=2),
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out",
                                                           "tau_s": 3.37e-6, "method": "thiran1"}},
        toy("toy_sat", "act", in_="delay.dl.out", out="act.u", unit="V", gain=1.0, limit=0.8),
        toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="V", k=2e4),
        {"stack": {"stages": [pid("c1", "v.y", 0.5, 1e2), pid("c2", "ctrl.c1.u", 0.2, 1e1, rate_div=4)],
                   "actuators": [gain("a1", "ctrl.c1.u", "act.a1", 2.0), gain("a2", "ctrl.c2.u", "act.a2", 3.0)]}},
        {"parallel": [gain("p1", "v.y", "v.p1", 1.5),
                      {"seq": [{"block": "delay_buffer", "name": "dl2", "params": {
                          "in": "v.y", "out": "delay.dl2.out", "tau_s": 2.5e-6, "method": "linear"}},
                               gain("p2", "delay.dl2.out", "v.p2", -1.0)], "name": "branch2"}]},
    ]


def make(obj, seed=11, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def assert_states_close(a, b, rtol=1e-12):
    np.testing.assert_array_equal(a.I, b.I)
    for x, y in ((a.S, b.S), (a.Z, b.Z)):
        np.testing.assert_allclose(x, y, rtol=rtol, atol=rtol * max(1.0, float(np.max(np.abs(y), initial=0.0))))


def test_compiled_equals_reference_nontrivial_pipeline():
    n = 2 * PIECE_SAMPLES + 777  # crosses piece boundaries
    b1, s1, p1 = make(nontrivial_obj())
    b2, s2, p2 = make(nontrivial_obj())
    run_reference(b1, s1, p1, n)
    run_compiled(b2, s2, p2, n)
    assert s1.n == s2.n == n and s1.t == s2.t == n / FS
    assert all(abs(float(s1[k])) > 0 for k in ("psi.plant.x", "act.u", "act.a2", "v.p2"))  # pipeline is alive
    assert_states_close(s1, s2)


def test_runtime_param_change_mid_run_matches():
    b1, s1, p1 = make(nontrivial_obj())
    b2, s2, p2 = make(nontrivial_obj())
    run_reference(b1, s1, p1, 500)
    run_compiled(b2, s2, p2, 500)
    for p in (p1, p2):
        p.set({"pid.kp": 0.1, "act.limit": 0.05})
    run_reference(b1, s1, p1, 500)
    run_compiled(b2, s2, p2, 500)
    assert_states_close(s1, s2)


def test_kernel_reuse_same_structure():
    obj = nontrivial_obj()
    b1, s1, p1 = make(obj, seed=1)
    loop1 = compile_loop(b1)
    size = loop_cache_size()
    obj2 = nontrivial_obj()
    obj2[4]["params"]["kp"] = 0.05  # runtime parameter value differs
    b2, s2, p2 = make(obj2, seed=2)  # and the seed
    assert compile_loop(b2) is loop1 and loop_cache_size() == size
    b3, s3, p3 = make(obj2, seed=2)
    run_compiled(b2, s2, p2, 1000, loop1)
    run_reference(b3, s3, p3, 1000)
    assert_states_close(s3, s2)
    obj3 = nontrivial_obj()
    obj3[5]["params"]["tau_s"] = 5e-6  # structural parameter differs -> new loop
    b4, _, _ = make(obj3)
    assert compile_loop(b4) is not loop1


def test_cache_key_by_content_not_repr_or_id():
    """Phase 5 F1-4: two 5000-point Tables that differ in one middle element (numpy's repr abbreviates
    them identically) get different structure keys and loops; an unregistered callable in a structural
    param has no stable key (ids are reused), so that pipeline bypasses the cache; a registered function
    is keyed by its registry name."""
    from cranberry.core.compile import structure_key
    from cranberry.core.registry import register_fn
    from cranberry.offline.tables import Table

    v1 = np.linspace(0.0, 1.0, 5000)
    v2 = v1.copy()
    v2[2500] += 1e-3
    assert repr(v1) == repr(v2)

    def with_param(value, g=1.0):
        obj = [gain("g", "v.x", "v.y", g)]
        obj[0]["params"]["structural"] = value
        lay = SlotLayout(FS)
        lay.add("v.x", (), np.float64, "V")
        blk = build(spec_from_obj(obj), lay, 0, external=["v.x"])
        return blk, lay

    t1 = Table("tab", {"f_hz": np.arange(5000.0)}, v1, "1")
    t2 = Table("tab", {"f_hz": np.arange(5000.0)}, v2, "1")
    (b1, l1), (b2, l2), (b3, l3) = with_param(t1), with_param(t2), with_param(t1, g=2.0)
    assert structure_key(b1, l1) != structure_key(b2, l2)
    assert structure_key(b1, l1) == structure_key(b3, l3)  # gain is runtime: same structure
    assert compile_loop(b1) is not compile_loop(b2) and compile_loop(b3) is compile_loop(b1)
    assert structure_key(*with_param(v1)) != structure_key(*with_param(v2))  # bare arrays too

    keys = {structure_key(*with_param(lambda x, c=c: c * x)) for c in range(3)}
    assert keys == {None}  # unregistered callables: no cache key
    n0 = loop_cache_size()
    b4, _ = with_param(lambda x: x)
    assert compile_loop(b4) is not compile_loop(b4) and loop_cache_size() == n0

    @register_fn("vt44_f1_4_square")
    def square(x):
        return x * x

    k1 = structure_key(*with_param(square))
    assert k1 is not None and k1 == structure_key(*with_param(square))


def stack_vs_expansion_objs():
    stages = [pid("c1", "v.e", 0.5, 2e3), pid("c2", "ctrl.c1.u", 0.3, 5e2), pid("c3", "ctrl.c2.u", 0.2, 1e2)]
    acts = [gain("a1", "ctrl.c1.u", "act.a1", 1.0), gain("a2", "ctrl.c2.u", "act.a2", 2.0),
            gain("a3", "ctrl.c3.u", "act.a3", 3.0)]
    src = toy("toy_noise_ar1", "e", out="v.e", unit="V", a=0.95, b=1.0)
    stacked = [src, {"stack": {"stages": stages, "actuators": acts}}]
    expanded = [src, stages[0], {"parallel": [acts[0], {"seq": [stages[1], {"parallel": [acts[1], {"seq": [
        stages[2], acts[2]]}]}]}]}]
    return stacked, expanded


def test_stack_equals_its_expansion_bitwise():
    stacked, expanded = stack_vs_expansion_objs()
    b1, s1, p1 = make(stacked)
    b2, s2, p2 = make(expanded)
    assert [x.name for x in iter_leaves(b1)] == [x.name for x in iter_leaves(b2)]
    run_compiled(b1, s1, p1, 3000)
    run_compiled(b2, s2, p2, 3000)
    np.testing.assert_array_equal(s1.S, s2.S)
    np.testing.assert_array_equal(s1.I, s2.I)


def test_stack_equals_parallel_in_linear_range():
    """Stacked u_k = C_k u_{k−1} ≡ parallel branches F_k = C_k⋯C_1 (T16 §4.4), linear blocks only."""
    stacked, _ = stack_vs_expansion_objs()
    src = stacked[0]

    def chain(k):  # branch k: copies of C_1..C_k then a_k, slot names suffixed by branch
        blocks = []
        inp = "v.e"
        for j, (kp, ki) in enumerate([(0.5, 2e3), (0.3, 5e2), (0.2, 1e2)][:k], start=1):
            nm = f"b{k}c{j}"
            blocks.append(pid(nm, inp, kp, ki))
            inp = f"ctrl.{nm}.u"
        blocks.append(gain(f"b{k}a", inp, f"act.b{k}", float(k)))
        return {"seq": blocks, "name": f"branch{k}"}

    parallel = [src, {"parallel": [chain(1), chain(2), chain(3)]}]
    b1, s1, p1 = make(stacked)
    b2, s2, p2 = make(parallel)
    a_stack, a_par = [], []
    for _ in range(200):
        run_compiled(b1, s1, p1, 7)
        run_compiled(b2, s2, p2, 7)
        a_stack.append([float(s1[f"act.a{k}"]) for k in (1, 2, 3)])
        a_par.append([float(s2[f"act.b{k}"]) for k in (1, 2, 3)])
    a_stack, a_par = np.array(a_stack), np.array(a_par)
    assert np.max(np.abs(a_stack)) > 0.1
    np.testing.assert_allclose(a_par, a_stack, rtol=1e-12, atol=1e-12 * np.max(np.abs(a_stack)))


def test_multirate_hold_compiled():
    obj = [toy("toy_sine", "s", out="v.s", unit="V", amp=1.0, freq_hz=1e4),
           toy("toy_gain", "slow", 5, in_="v.s", out="v.slow", unit="V", gain=1.0)]
    b1, s1, p1 = make(obj)
    got, src = [], []
    loop = compile_loop(b1)
    for _ in range(23):
        run_compiled(b1, s1, p1, 1, loop)
        got.append(float(s1["v.slow"]))
        src.append(float(s1["v.s"]))
    expect = [src[(n // 5) * 5] for n in range(23)]  # updated at n % 5 == 0, held in between
    np.testing.assert_array_equal(got, expect)


class PythonOnly(BaseBlock):
    def __init__(self):
        super().__init__("py_only", "py", {}, writes=("v.py",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add("v.py")

    def __call__(self, state, params):
        state["v.py"] = float(state["v.py"]) + 1.0
        return state


def test_python_only_leaf_not_compilable_but_reference_runs():
    lay = SlotLayout(FS)
    blk = build(Seq([PythonOnly()]), lay, 0)
    with pytest.raises(NotCompilableError, match="py"):
        compile_loop(blk)
    st, p = State(lay), Params(lay)
    run_reference(blk, st, p, 10)
    assert float(st["v.py"]) == 10.0 and st.n == 10
