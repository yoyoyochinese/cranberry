"""YAML tags, safe loading and dumping (G12, prompt §5.2)."""

import pytest
import yaml

from cranberry.config.yaml_io import (
    YAML_TAGS,
    Choice,
    Expr,
    Fn,
    PiFrac,
    Sweep,
    SweepLog,
    canonical_json,
    dump_yaml_obj,
    load_yaml,
    loads_yaml,
)


def test_all_tags_load_and_round_trip(tmp_path):
    text = """
a: !sweep [1e-8, 1e-2]
b: !sweep [0, 15, 4]
c: !sweep_log {lo: 1.0e-3, hi: 10, n: 5}
d: !choice [+1, -1, !pi_frac 0.5]
e: !expr "arm1.L_m - dL"
f: !expr {source: "a*x + b", where: {a: 2, b: pump.x}, args: [x]}
g: !fn {name: refB_loss_vs_x, args: {k: 1.5}}
h: !fn refB_loss_vs_x
i: !pi_frac -0.25
j: [1, 2.5, "s", true, null]
"""
    obj = loads_yaml(text)
    assert obj["a"] == Sweep(1e-8, 1e-2)
    assert obj["b"] == Sweep(0, 15, 4)
    assert obj["c"] == SweepLog(1e-3, 10, 5)
    assert obj["d"] == Choice((1, -1, PiFrac(0.5)))
    assert obj["e"] == Expr("arm1.L_m - dL")
    assert obj["f"] == Expr("a*x + b", {"a": 2, "b": "pump.x"}, ("x",))
    assert obj["g"] == Fn("refB_loss_vs_x", {"k": 1.5})
    assert obj["h"] == Fn("refB_loss_vs_x")
    assert obj["i"] == PiFrac(-0.25)
    text2 = dump_yaml_obj(obj)
    for tag in YAML_TAGS:
        assert tag in text2
    assert loads_yaml(text2) == obj
    p = tmp_path / "x.yaml"
    p.write_text(text2)
    assert load_yaml(p) == obj


def test_floats_without_dot_and_strings_that_look_like_floats():
    obj = loads_yaml("a: 1e-8\nb: 1E+3\nc: -2.5e-3\nd: '1e-8'\ne: .5\nf: 7\n")
    assert obj == {"a": 1e-8, "b": 1e3, "c": -2.5e-3, "d": "1e-8", "e": 0.5, "f": 7}
    assert type(obj["f"]) is int
    back = loads_yaml(dump_yaml_obj(obj))
    assert back == obj and type(back["d"]) is str and type(back["a"]) is float
    # dumped floats are readable by plain YAML-1.1 parsers as floats too
    assert yaml.safe_load(dump_yaml_obj({"x": 1e-8}))["x"] == 1e-8


def test_timestamps_stay_strings():
    assert loads_yaml("d: 2026-10-31\n") == {"d": "2026-10-31"}


def test_float_round_trip_exact():
    vals = [0.1, 1 / 3, 2.0**-1074, 1.7976931348623157e308, -0.0, 6.02214076e23]
    assert loads_yaml(dump_yaml_obj({"v": vals}))["v"] == vals


@pytest.mark.parametrize(
    "text",
    [
        "a: !!python/object/apply:os.system ['true']",
        "a: !!python/name:os.system",
        "a: !!python/object:collections.OrderedDict {}",
        "a: !unknown_tag 1",
    ],
)
def test_unsafe_and_unknown_tags_rejected(text):
    with pytest.raises(yaml.YAMLError):
        loads_yaml(text)


@pytest.mark.parametrize(
    "text",
    [
        "a: 1\na: 2\n",  # duplicate keys are not silently overwritten
        "a: !sweep [1]",
        "a: !sweep [1, 2, 3, 4]",
        "a: !sweep [x, 2]",
        "a: !sweep [1, 2, 0]",
        "a: !sweep [1, 2, 2.5]",
        "a: !sweep_log [0, 1]",
        "a: !sweep_log [-1, 1]",
        "a: !sweep {lo: 1}",
        "a: !choice []",
        "a: !choice 1",
        "a: !expr ''",
        "a: !expr {where: {}}",
        "a: !expr {source: x, args: [1x]}",
        "a: !fn {args: {}}",
        "a: !pi_frac abc",
        "a: !pi_frac [1]",
    ],
)
def test_malformed_tags_and_duplicates(text):
    with pytest.raises(yaml.YAMLError):
        loads_yaml(text)


def test_alias_bomb_bounded():
    lines = ["a0: &a0 [x, x, x, x, x, x, x, x, x, x]"]
    for i in range(1, 9):
        prev = f"*a{i - 1}"
        lines.append(f"a{i}: &a{i} [{', '.join([prev] * 10)}]")
    with pytest.raises(ValueError, match="too large"):
        loads_yaml("\n".join(lines))


def test_single_document_only():
    with pytest.raises(ValueError):
        loads_yaml("a: 1\n---\nb: 2\n")


def test_canonical_json_key_order_independent():
    assert canonical_json({"a": 1, "b": Sweep(0, 1)}) == canonical_json({"b": Sweep(0, 1), "a": 1})
    with pytest.raises(TypeError):
        canonical_json({"a": object()})
