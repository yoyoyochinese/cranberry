"""RunSpec from config trees, ``to_spec``, and D34 enforcement at config level (no silent defaults)."""

import pytest

from cranberry import config
from cranberry.config import ConfigError, RunSettings
from cranberry.config.spec import apply_toggles, iter_block_specs, spec_from_tree
from cranberry.config.yaml_io import Fn, Sweep, loads_yaml
from cranberry.core.pipeline import SeqSpec

RUN = {"duration_s": 1.0e-3, "f_sim_hz": 1.0e6, "chunk_samples": 64, "seed": 3}


def noise(name="n1", **params):
    p = {"out": f"drv.{name}.x", "sigma": 1.0, "toggles": {"enabled": True}}
    p.update(params)
    return {"block": "cfgtest_noise", "name": name, "params": p}


def tree(**over):
    t = {"preset": "t", "run": dict(RUN), "pipeline": [noise()]}
    t.update(over)
    return t


def test_minimal_tree_and_types_normalised():
    spec = config.load(tree(run={"duration_s": 1, "f_sim_hz": 1000000, "chunk_samples": 64.0, "seed": 3}))
    assert spec.run == RunSettings(1.0e-3 * 1000, 1.0e6, 64, 3)
    assert type(spec.run.duration_s) is float and type(spec.run.chunk_samples) is int
    assert isinstance(spec.pipeline, SeqSpec) and spec.topology == {} and spec.overlays == ()


@pytest.mark.parametrize(
    ("t", "match"),
    [
        ({"run": RUN}, r"missing required section\(s\) \['pipeline'\]"),
        ({"pipeline": []}, r"\['run'\]"),
        ({"run": {k: v for k, v in RUN.items() if k != "seed"}, "pipeline": []}, r"missing run\.seed"),
        ({"run": {**RUN, "dt": 1}, "pipeline": []}, "unknown key"),
        ({"run": {**RUN, "seed": Sweep(0, 1, 2)}, "pipeline": []}, "sweep markers are not allowed"),
        ({"run": {**RUN, "seed": 1.5}, "pipeline": []}, "integers"),
        ({"run": {**RUN, "f_sim_hz": -1.0}, "pipeline": []}, "f_sim_hz"),
        ({"run": RUN, "pipeline": [], "pipline": []}, "unknown top-level"),
        ({"run": RUN, "pipeline": [{"blok": "x"}]}, "pipeline"),
        ({"run": RUN, "pipeline": [], "toggles": {"n1": True}}, "<block name>.<toggle>"),
        ({"run": RUN, "pipeline": [], "toggles": {"n1.enabled": "yes"}}, "true/false"),
        ({"run": RUN, "pipeline": [], "topology": [1]}, "topology: must be a mapping"),
        ({"run": RUN, "pipeline": [], "preset": 3}, "preset"),
    ],
)
def test_tree_errors_name_the_key(t, match):
    with pytest.raises(ConfigError, match=match):
        config.load(t)


def test_spec_from_tree_rejects_unapplied_overlays():
    with pytest.raises(ConfigError, match="overlays"):
        spec_from_tree({**tree(), "overlays": ["x.yaml"]})


def test_validate_d34_missing_toggle_and_parameter_name_block_and_key():
    config.validate(config.load(tree()))
    p = noise()
    del p["params"]["toggles"]
    with pytest.raises(ConfigError, match=r"pipeline\[0\]: block 'n1' \(cfgtest_noise\).*enabled.*D34"):
        config.validate(config.load(tree(pipeline=[p])))
    p = noise()
    del p["params"]["sigma"]
    with pytest.raises(ConfigError, match=r"block 'n1'.*sigma.*D34"):
        config.validate(config.load(tree(pipeline=[p])))
    with pytest.raises(ConfigError, match="unknown block kind 'nope'"):
        config.validate(config.load(tree(pipeline=[{"block": "nope", "name": "x"}])))


def test_validate_names_toggles_and_markers():
    with pytest.raises(ConfigError, match="duplicate block name 'n1'"):
        config.validate(config.load(tree(pipeline=[noise(), {"parallel": [noise()]}])))
    with pytest.raises(ConfigError, match="no block named 'zz'"):
        config.validate(config.load(tree(toggles={"zz.enabled": False})))
    with pytest.raises(ConfigError, match="no toggle 'nope'"):
        config.validate(config.load(tree(toggles={"n1.nope": False})))
    with pytest.raises(ConfigError, match=r"unresolved .*pipeline\.n1\.params\.sigma"):
        config.validate(config.load(tree(pipeline=[noise(sigma=Sweep(0.1, 1.0, 3))])))


def test_apply_toggles_and_iter_block_specs():
    spec = config.load(tree(pipeline=[noise("a.b"), {"stack": {"stages": [noise("s")], "actuators": [noise("x")]}}],
                            toggles={"a.b.enabled": False}))
    assert [(p, b.name) for p, b in iter_block_specs(spec.pipeline)] == [
        ("pipeline[0]", "a.b"), ("pipeline[1].stack.stages[0]", "s"), ("pipeline[1].stack.actuators[0]", "x")]
    out = apply_toggles(spec.pipeline, spec.toggles)
    assert out.children[0].params["toggles"] == {"enabled": False}
    assert spec.pipeline.children[0].params["toggles"] == {"enabled": True}  # spec untouched


def test_to_spec_matches_load_and_rejects_lambdas(blocks):
    run = RunSettings(**RUN)
    py = [blocks.noise({"out": "drv.n1.x", "sigma": 1.0, "toggles": {"enabled": True}}, name="n1")]
    assert config.to_spec(py, run=run, preset="t") == config.load(tree())
    seq = py[0] >> blocks.shaper({"in": "drv.n1.x", "out": "v.s", "curve": blocks.affine}, name="s")
    spec = config.to_spec(seq, run=run, preset="t")
    assert spec.pipeline.children[1].params["curve"] == Fn("cfgtest_affine")
    assert config.load(loads_yaml(config.dump_yaml(spec))) == spec
    bad = py[0] >> blocks.shaper({"in": "drv.n1.x", "out": "v.s", "curve": lambda x: 2 * x}, name="s")
    with pytest.raises(ConfigError, match=r"block 's'.*curve.*unregistered callable.*lambda"):
        config.to_spec(bad, run=run)
    with pytest.raises(ConfigError, match="RunSettings"):
        config.to_spec(py, run=RUN)
