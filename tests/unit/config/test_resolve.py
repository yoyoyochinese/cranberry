"""Reference / callable resolution with cycle detection (G12, prompt §5.5)."""

import math

import numpy as np
import pytest

from cranberry.config import ConfigError
from cranberry.config.expr import ExprError, compile_expr
from cranberry.config.resolve import as_callable, find_markers, normalize, parse_path, path_str, resolve
from cranberry.config.yaml_io import Expr, Fn, PiFrac, Sweep, loads_yaml


def test_plan_example_where_and_sibling_reference():
    tree = loads_yaml("""
topology:
  segments:
    - {id: arm1_915b, L_m: 8.4, N_sub: 8}
    - {id: arm2, L_m: !expr "arm1_915b.L_m - dL", where: {dL: 9.0}}
    - {id: arm3, L_m: !expr "2 * N_sub_local", N_sub_local: 3}
""")
    out = resolve(tree)
    segs = out["topology"]["segments"]
    assert segs[1] == {"id": "arm2", "L_m": 8.4 - 9.0}  # where: removed once used
    assert segs[2]["L_m"] == 6 and type(segs[2]["L_m"]) is int


def test_absolute_paths_dotted_block_names_and_where_references():
    tree = loads_yaml("""
run: {f_sim_hz: 2.0e6}
pump: {x: 0.5}
pipeline:
  - {block: pid, name: lock1.pid, params: {kp: 2.0}}
  - {block: pid, name: lock2.pid, params: {kp: !expr "lock1.pid.params.kp * 3"}}
  - {block: x, name: n, params: {dt: !expr "1 / run.f_sim_hz"}}
  - {block: x, name: m, params: {y: !expr "2 * px", where: {px: pump.x}}}
""")
    out = resolve(tree)
    assert out["pipeline"][1]["params"]["kp"] == 6.0
    assert out["pipeline"][2]["params"]["dt"] == 0.5e-6
    assert out["pipeline"][3]["params"] == {"y": 1.0}
    assert "where" not in out["pipeline"][3]


def test_chained_expressions_resolve_in_dependency_order():
    tree = {"c": Expr("b + 1"), "b": Expr("a * 2"), "a": Expr("1 + 1")}
    assert resolve(tree) == {"c": 5, "b": 4, "a": 2}


def test_cycles_are_errors_naming_the_cycle():
    with pytest.raises(ConfigError, match="circular"):
        resolve({"a": Expr("b + 1"), "b": Expr("a + 1")})
    with pytest.raises(ConfigError, match="circular"):
        resolve({"a": Expr("a + 1")})
    with pytest.raises(ConfigError, match=r"seg\.s1\.L_m"):
        resolve({"seg": [{"id": "s1", "L_m": Expr("s2.L_m")}, {"id": "s2", "L_m": Expr("s1.L_m")}]})


def test_unknown_names_and_bad_references():
    with pytest.raises(ConfigError, match="cannot resolve name 'nope'"):
        resolve({"a": Expr("nope + 1")})
    with pytest.raises(ConfigError, match="mapping/list"):
        resolve({"m": {"k": 1}, "a": Expr("m * 2")})
    with pytest.raises(ConfigError, match="!fn"):
        resolve({"f": Fn("x"), "a": Expr("f * 2")})
    with pytest.raises(ConfigError, match="where-binding"):
        resolve({"a": Expr("q", {"q": "no.such.path"})})
    with pytest.raises(ConfigError, match="where: must be a mapping"):
        resolve({"a": Expr("1"), "where": 3})
    with pytest.raises(ExprError, match=r"a: !expr"):  # evaluation errors name the path
        resolve({"a": Expr("1 / z"), "z": 0})


def test_sweep_dependence_is_deferred_and_where_kept():
    tree = loads_yaml("""
segs:
  - {id: a, L_m: 8.4}
  - {id: b, L_m: !expr "a.L_m - dL", where: {dL: !sweep [0.0, 15.0, 4]}, k: !expr "2 * 3"}
x: !expr "segs.b.L_m * 2"
""")
    out = resolve(tree)
    b = out["segs"][1]
    assert b["L_m"] == Expr("a.L_m - dL") and b["where"] == {"dL": Sweep(0.0, 15.0, 4)} and b["k"] == 6
    assert out["x"] == Expr("segs.b.L_m * 2")  # depends on a deferred expression → deferred too
    assert [path_str(out, p) for p, _ in find_markers(out)] == ["segs.b.where.dL"]
    # fixing the axis value resolves both
    out["segs"][1]["where"]["dL"] = 4.4
    again = resolve(out)
    assert again["segs"][1] == {"id": "b", "L_m": 4.0, "k": 6} and again["x"] == 8.0
    assert resolve(again) == again  # idempotent


def test_pi_frac_becomes_float_everywhere():
    out = resolve({"s": PiFrac(0.5), "l": [PiFrac(-1)], "f": Fn("g", {"phi": PiFrac(2)}), "e": Expr("2 * s")})
    assert out == {"s": math.pi / 2, "l": [-math.pi], "f": Fn("g", {"phi": 2 * math.pi}), "e": math.pi}


def test_callable_expr_captures_resolved_bindings():
    tree = loads_yaml("""
opo: {L0: 0.00249}
loss: !expr {source: "L0 + k * x", where: {L0: opo.L0, k: 0.00222}, args: [x]}
""")
    out = resolve(tree)
    assert out["loss"] == Expr("L0 + k * x", {"L0": 0.00249, "k": 0.00222}, ("x",))
    f = as_callable(out["loss"], ("x",))
    assert f(2.0) == pytest.approx(0.00249 + 0.00222 * 2.0)
    assert resolve(out) == out
    with pytest.raises(ConfigError, match="callable"):
        resolve({**tree, "y": Expr("loss * 2")})
    with pytest.raises(ConfigError, match="args"):
        as_callable(out["loss"], ("P",))


def test_as_callable_forms(blocks):
    assert as_callable(Fn("cfgtest_affine", {"a": 2.0, "b": 1.0}), ("x",))(3.0) == 7.0
    assert as_callable(0.5, ("x",))(123.0) == 0.5
    assert as_callable(blocks.affine, ("x",))(3.0) == 3.0
    with pytest.raises(ConfigError):
        as_callable("text", ("x",))


def test_normalize_python_values(blocks):
    import functools

    v = normalize({"t": (1, 2), "a": np.array([1.0, 2.0]), "s": np.float64(0.5), "i": np.int64(3),
                   "f": blocks.affine, "p": functools.partial(blocks.affine, a=2.0),
                   "e": compile_expr("2 * x", ("x",))})
    assert v == {"t": [1, 2], "a": [1.0, 2.0], "s": 0.5, "i": 3, "f": Fn("cfgtest_affine"),
                 "p": Fn("cfgtest_affine", {"a": 2.0}), "e": Expr("2 * x", {}, ("x",))}
    assert type(v["s"]) is float and type(v["i"]) is int
    with pytest.raises(ConfigError, match="unregistered callable"):
        normalize({"f": lambda x: x}, "pipeline[0].params")
    with pytest.raises(ConfigError, match="positional"):
        normalize(functools.partial(blocks.affine, 1.0))
    with pytest.raises(ConfigError, match="set"):
        normalize({1, 2})
    with pytest.raises(ConfigError, match="key"):
        normalize({1: 2})


def test_parse_path_round_trips_path_str():
    tree = {"pipeline": [{"block": "a", "name": "lock1.pid", "params": {"kp": 1}},
                         {"stack": {"stages": [{"block": "b", "name": "s", "params": {"x": 2}}], "actuators": []}}]}
    for p in [("pipeline", 0, "params", "kp"), ("pipeline", 1, "stack", "stages", 0, "params", "x")]:
        assert parse_path(tree, path_str(tree, p)) == p
    assert parse_path(tree, "pipeline[1].stack.stages[0].params.x") == ("pipeline", 1, "stack", "stages", 0, "params", "x")
    with pytest.raises(ConfigError):
        parse_path(tree, "pipeline.nope.params")
