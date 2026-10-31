"""dB → linear on load (G16-i, ``config.db``): names, values, scope, round trip."""

import math

import pytest

from cranberry import config
from cranberry.config import ConfigError
from cranberry.config.db import db_value_to_linear, linear_name, linearize
from cranberry.config.expr import compile_expr
from cranberry.config.resolve import as_callable, resolve
from cranberry.config.yaml_io import Expr, Fn, loads_yaml
from cranberry.core import units

TEXT = """
preset: dbtest
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 1}
topology:
  segments:
    - {id: s1, L_m: 2.0, il_db: 3.0}
offline:
  laser: {pn_dbc_per_hz: -100.0, rin_db_per_hz: -140.0, il_db: [0.0, 10.0, [20.0]]}
  opo: {loss_vs_x_db: !expr {source: "a * x", where: {a: 2.0}, args: [x]}}
  derived: {g_db: !expr "topology.segments.s1.il_db * 2"}
pipeline:
  - {block: cfgtest_shaper, name: sh, params: {in: t, out: v.sh, curve: !fn {name: cfgtest_affine, args: {a_db: 10.0}}}}
toggles: {}
measure: {esa: {rbw_db: 0.0}}
"""


def test_names():
    assert linear_name("il_db") == "il_lin"
    assert linear_name("IL_dB") == "IL_lin"
    assert linear_name("rin_db_per_hz") == "rin_per_hz"
    assert linear_name("pn_dbc_per_hz") == "pn_rad2_per_hz"
    assert linear_name("f_sim_hz") is None and linear_name("_db") is None and linear_name("dbx") is None
    for name in ("il_db", "rin_db_per_hz", "pn_dbc_per_hz"):  # same notion of "dB name" as core.units
        assert units.is_db_name(name)


def test_values_follow_core_units():
    assert db_value_to_linear("il_db", 3.0) == float(units.db_to_lin(3.0))
    assert db_value_to_linear("rin_db_per_hz", -140.0) == pytest.approx(1e-14)
    assert db_value_to_linear("pn_dbc_per_hz", -100.0) == pytest.approx(2e-10)  # S_phi = 2 L (conv §7)
    assert db_value_to_linear("sq_db", -10) == pytest.approx(0.1) and type(db_value_to_linear("sq_db", 0)) is float
    for bad in ("3 dB", True, None, {"a": 1}, Fn("f")):
        with pytest.raises(ConfigError, match="il_lin"):
            db_value_to_linear("il_db", bad, "x.il_db")


def test_load_converts_and_renames_everywhere_but_where_and_run():
    spec = config.load(loads_yaml(TEXT))
    seg = spec.topology["segments"][0]
    assert seg == {"id": "s1", "L_m": 2.0, "il_lin": pytest.approx(10 ** 0.3)}
    laser = spec.offline["laser"]
    assert set(laser) == {"pn_rad2_per_hz", "rin_per_hz", "il_lin"}
    assert laser["il_lin"] == [1.0, 10.0, [100.0]]
    assert spec.offline["derived"] == {"g_lin": pytest.approx(10 ** 0.6)}  # the expr saw il_db = 3 dB
    f = spec.offline["opo"]["loss_vs_x_lin"]  # a callable is wrapped
    assert as_callable(f, ("x",))(1.5) == pytest.approx(10 ** (3.0 / 10))
    assert spec.pipeline.children[0].params["curve"] == Fn("cfgtest_affine", {"a_lin": 10.0})
    assert spec.measure == {"esa": {"rbw_lin": 1.0}}


def test_round_trip_is_lossless_and_idempotent(tmp_path):
    spec = config.load(loads_yaml(TEXT))
    text = config.dump_yaml(spec)
    assert "_db" not in text and "il_lin: 1.9952623149688795" in text
    again = config.load(loads_yaml(text))
    assert again == spec and config.dump_yaml(again) == text
    assert linearize(config.spec.spec_to_tree(spec)) == config.spec.spec_to_tree(spec)


def test_python_to_spec_converts_the_same_way(blocks):
    from cranberry.config import RunSettings

    sh = blocks.shaper({"in": "t", "out": "v.sh", "curve": compile_expr("x", ("x",)), "gain_db": 20.0}, name="sh")
    spec = config.to_spec([sh], run=RunSettings(1e-3, 1e6, 64, 1), offline={"tap": {"il_db": 10.0}})
    assert spec.pipeline.children[0].params["gain_lin"] == pytest.approx(100.0)
    assert spec.offline == {"tap": {"il_lin": 10.0}}


def test_errors_name_the_path():
    with pytest.raises(ConfigError, match=r"offline\.tap\.il_db: both 'il_db' \(dB\) and 'il_lin'"):
        config.load({**loads_yaml(TEXT), "offline": {"tap": {"il_db": 1.0, "il_lin": 2.0}}})
    with pytest.raises(ConfigError, match=r"offline\.tap\.il_db: a dB parameter needs a number"):
        config.load({**loads_yaml(TEXT), "offline": {"tap": {"il_db": "high"}}})


def test_where_bindings_and_expr_names_are_not_renamed():
    out = resolve(
        {"offline": {"a": Expr("x_db + 1", {"x_db": 2.0}, ("y",)), "where": {"k_db": 1.0}, "b": Expr("k_db")}}
    )
    assert out["offline"]["a"] == Expr("x_db + 1", {"x_db": 2.0}, ("y",))
    assert out["offline"]["b"] == 1.0 and "where" not in out["offline"]
    assert math.isclose(db_value_to_linear("x_db", 0.0), 1.0)
