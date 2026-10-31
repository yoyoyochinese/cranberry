"""``runner.sweep``: axis expansion (cartesian) and deterministic per-point seeds (PLAN §2.8, D54)."""

import math

import pytest

from cranberry import config
from cranberry.config import ConfigError
from cranberry.config.yaml_io import loads_yaml
from cranberry.config.yaml_io import Choice, Sweep, SweepLog
from cranberry.runner.sweep import axis_values, expand, point_seed, sweep

TEXT = """
preset: sweeptest
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 11, seed_mode: per_point}
topology:
  segments:
    - {id: a, L_m: 8.4}
    - {id: b, L_m: !expr "a.L_m - dL", where: {dL: !sweep [0.0, 4.0, 3]}}
  locks:
    - {id: l1, setpoint: !choice [!pi_frac 0.5, !pi_frac -0.5]}
pipeline:
  - {block: cfgtest_noise, name: n1, params: {out: drv.n1.x, sigma: !sweep_log [0.01, 1.0, 3], toggles: {enabled: true}}}
"""


def test_axis_values():
    assert axis_values(Sweep(0, 1, 3)) == [0.0, 0.5, 1.0]
    assert axis_values(SweepLog(1e-3, 1e-1, 3)) == pytest.approx([1e-3, 1e-2, 1e-1])
    assert axis_values(Choice((1, "x"))) == [1, "x"]
    with pytest.raises(ConfigError, match="no point count"):
        axis_values(Sweep(0, 1), "p.q")


def test_cartesian_expansion_resolves_each_point():
    spec = config.load(loads_yaml(TEXT))
    pts = expand(spec)
    assert len(pts) == 3 * 2 * 3
    assert list(pts[0].values) == ["topology.segments.b.where.dL", "topology.locks.l1.setpoint", "pipeline.n1.params.sigma"]
    p = pts[1 * 6 + 1 * 3 + 2]  # dL=2.0, setpoint=-π/2, sigma=1.0
    assert p.values == {"topology.segments.b.where.dL": 2.0, "topology.locks.l1.setpoint": -math.pi / 2,
                        "pipeline.n1.params.sigma": 1.0}
    assert p.spec.topology["segments"][1] == {"id": "b", "L_m": 8.4 - 2.0}
    assert p.spec.topology["locks"][0]["setpoint"] == -math.pi / 2
    assert p.spec.pipeline.children[0].params["sigma"] == 1.0
    for q in pts:
        config.validate(q.spec)  # every point is a single run


def test_point_seeds_deterministic_distinct_and_grid_independent():
    spec = config.load(loads_yaml(TEXT))
    pts = expand(spec)
    seeds = [p.spec.run.seed for p in pts]
    assert len(set(seeds)) == len(seeds)
    assert seeds == [p.spec.run.seed for p in expand(config.load(loads_yaml(TEXT)))]
    assert all(0 <= s < 2**63 for s in seeds)
    for p in pts:
        assert p.spec.run.seed == point_seed(11, p.values)
    # a point's seed depends on its coordinates, not on the rest of the grid
    finer = expand(spec, {"topology.segments.b.where.dL": [2.0, 3.0]})
    same = [q for q in finer if q.values == pts[1 * 6 + 1 * 3 + 2].values][0]
    assert same.spec.run.seed == pts[1 * 6 + 1 * 3 + 2].spec.run.seed
    # a different root seed changes every point
    other = expand(config.load({**loads_yaml(TEXT), "run": {**loads_yaml(TEXT)["run"], "seed": 12}}))
    assert all(a.spec.run.seed != b.spec.run.seed for a, b in zip(pts, other))


def test_explicit_axes_and_seed_axis():
    spec = config.load(loads_yaml(TEXT))
    before = config.dump_yaml(spec)
    pts = expand(spec, {"run.f_sim_hz": [1e6, 2e6], "run.seed": [5], "pipeline.n1.params.sigma": [0.3]})
    assert config.dump_yaml(spec) == before  # expansion never mutates the spec
    assert len(pts) == 3 * 2 * 1 * 2
    assert {p.spec.run.f_sim_hz for p in pts} == {1e6, 2e6}
    assert {p.spec.run.seed for p in pts} == {5}  # verbatim
    with pytest.raises(ConfigError, match="no config value"):
        expand(spec, {"pipeline.nope.params.x": [1]})
    with pytest.raises(ConfigError, match="no values"):
        expand(spec, {"run.seed": []})


def test_range_without_n_needs_explicit_values():
    spec = config.load(loads_yaml(TEXT.replace("!sweep [0.0, 4.0, 3]", "!sweep [0.0, 4.0]")))
    with pytest.raises(ConfigError, match=r"topology\.segments\.b\.where\.dL.*no point count"):
        expand(spec)
    assert len(expand(spec, {"topology.segments.b.where.dL": [1.0]})) == 6


def test_no_axes_gives_single_point_with_derived_seed():
    spec = config.load({"run": {"duration_s": 1e-3, "f_sim_hz": 1e6, "chunk_samples": 8, "seed": 1}, "pipeline": []})
    (p,) = expand(spec)
    assert p.values == {} and p.spec.run.seed == point_seed(1, {})


def test_sweep_runs_each_point_via_runner_run():
    spec = config.load(loads_yaml(TEXT))
    try:
        results = sweep(spec)
    except NotImplementedError as e:
        assert "WP-F1" in str(e)
    else:  # once WP-F1 lands
        assert len(results) == 18


# ------------------------------------------------------------------------- WP-F2b: seed mode, nesting


def test_seed_mode_is_required_with_axes_and_validated():
    no_mode = TEXT.replace(", seed_mode: per_point", "")
    spec = config.load(loads_yaml(no_mode))
    assert spec.run.seed_mode is None and "seed_mode" not in config.dump_yaml(spec)
    with pytest.raises(ConfigError, match=r"seed_mode is required.*per_point \| common"):
        expand(spec)
    single = config.load(loads_yaml(no_mode.replace("!sweep [0.0, 4.0, 3]", "1.0")
                                    .replace("!choice [!pi_frac 0.5, !pi_frac -0.5]", "0.0")
                                    .replace("!sweep_log [0.01, 1.0, 3]", "0.1")))
    assert len(expand(single)) == 1  # no axis: no mode needed
    with pytest.raises(ConfigError, match="seed_mode is required"):  # explicit axes count too
        expand(single, {"run.f_sim_hz": [1e6]})
    assert len(expand(spec, {"run.seed": [1, 2]})) == 36  # a seed axis is verbatim: no mode needed
    with pytest.raises(ConfigError, match="seed_mode must be one of"):
        config.load(loads_yaml(TEXT.replace("per_point", "random")))


def test_seed_mode_common_gives_every_point_the_root_seed():
    spec = config.load(loads_yaml(TEXT.replace("per_point", "common")))
    pts = expand(spec)
    assert {p.spec.run.seed for p in pts} == {11}
    assert all(p.spec.run.seed_mode == "common" for p in pts)
    text = config.dump_yaml(spec)
    assert "seed_mode: common" in text and config.load(loads_yaml(text)) == spec


def test_sweep_marker_inside_expr_where_and_fn_args_are_axes():
    text = """
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 1, seed_mode: common}
pipeline:
  - {block: cfgtest_noise, name: n1, params: {out: drv.n1.x, toggles: {enabled: true},
     sigma: !expr {source: "2 * s", where: {s: !sweep [1.0, 2.0, 2]}}}}
  - {block: cfgtest_shaper, name: sh, params: {in: drv.n1.x, out: v.sh,
     curve: !fn {name: cfgtest_affine, args: {a: !choice [1.0, 3.0], b: !expr "offline.opo.L0"}}}}
offline: {opo: {L0: 0.5}}
"""
    spec = config.load(loads_yaml(text))
    assert sorted(config.sweep_markers(spec)) == ["pipeline.n1.params.sigma", "pipeline.n1.params.sigma.where.s",
                                                  "pipeline.sh.params.curve.args.a"]
    assert spec.pipeline.children[1].params["curve"].args["b"] == 0.5  # resolved although a sibling is an axis
    pts = expand(spec)
    assert [p.values for p in pts] == [
        {"pipeline.n1.params.sigma.where.s": s, "pipeline.sh.params.curve.args.a": a}
        for s in (1.0, 2.0) for a in (1.0, 3.0)
    ]
    last = pts[-1].spec.pipeline.children
    assert last[0].params["sigma"] == 4.0
    assert last[1].params["curve"] == config.yaml_io.Fn("cfgtest_affine", {"a": 3.0, "b": 0.5})
    for p in pts:
        config.validate(p.spec)


def test_db_parameters_stay_db_axes_and_each_point_is_linearised():
    text = """
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 1, seed_mode: per_point}
offline: {tap: {il_db: !sweep [0.0, 10.0, 2], k: !expr "il_db * 2"}}
pipeline: []
"""
    spec = config.load(loads_yaml(text))
    assert "il_db" in spec.offline["tap"]  # unresolved axis: kept as written
    pts = expand(spec)
    assert [p.values for p in pts] == [{"offline.tap.il_db": 0.0}, {"offline.tap.il_db": 10.0}]
    assert pts[1].spec.offline["tap"] == {"il_lin": 10.0, "k": 20.0}  # the expr saw the dB value
    # a single-run spec holds il_lin; an explicit axis may still name il_db (values in dB)
    single = config.load(loads_yaml(text.replace("!sweep [0.0, 10.0, 2]", "3.0")))
    assert single.offline["tap"]["il_lin"] == pytest.approx(10 ** 0.3)
    q = expand(single, {"offline.tap.il_db": [20.0]})
    assert q[0].values == {"offline.tap.il_db": 20.0}
    assert q[0].spec.offline["tap"]["il_lin"] == pytest.approx(100.0)
    assert q[0].spec.offline["tap"]["k"] == 6.0  # resolved at load, before the sweep


def test_recorded_sweep_writes_one_record_per_point(tmp_path):
    import json

    from cranberry.runner import record

    spec = config.load(loads_yaml(TEXT.replace("per_point", "common")))
    results = sweep(spec, {"topology.segments.b.where.dL": [1.0], "topology.locks.l1.setpoint": [0.0],
                           "pipeline.n1.params.sigma": [0.1, 0.2]}, record_root=tmp_path)
    dirs = [r.run_dir for r in results]
    assert len(set(dirs)) == 2 and all(d.parent == tmp_path for d in dirs)
    for r in results:
        inv = json.loads((r.run_dir / "invocation.json").read_text())
        assert inv["seed_mode"] == "common" and inv["n_points"] == 2
        assert inv["sweep_point"]["pipeline.n1.params.sigma"] == r.spec.pipeline.children[0].params["sigma"]
        assert record.load_record(r.run_dir) == r.spec
        assert config.load(r.run_dir / "config.input.yaml") == spec  # the sweep's input document
