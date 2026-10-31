"""VT-F1 (PLAN §6.2, prompt §5.2; protected list "YAML config loading and its equivalence with
Python-list pipelines"): a pipeline built from Python objects and the same pipeline written in YAML
give identical ``RunSpec``s and identical built pipelines; ``spec == load(dump_yaml(spec))``; and the
same seed gives bit-identical results (tolerance: exact, ``numpy.array_equal`` on float64).

Test blocks (kinds ``vtf1_*``) stand in for physics blocks; the reference loop is the pure-Python
step of ``core.block`` (``runner.run`` from WP-F1 is used additionally once it exists).
"""

import math

import numpy as np
import pytest

from cranberry import config
from cranberry.config import RunSettings
from cranberry.config.expr import compile_expr
from cranberry.config.resolve import as_callable
from cranberry.config.spec import apply_toggles
from cranberry.config.yaml_io import loads_yaml
from cranberry.core import registry
from cranberry.core.block import BaseBlock, require
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import Parallel, Seq, Stack, from_spec, iter_leaves
from cranberry.core.rng import stream
from cranberry.core.state import Params, State


class _Noise(BaseBlock):
    toggle_names = ("enabled",)

    def __init__(self, params, name):
        require(params, "out", block=name or "vtf1_noise")
        super().__init__("vtf1_noise", name, params, writes=(params["out"],), runtime_params=("sigma",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"])

    def __call__(self, state, params):
        x = self.rng.standard_normal()
        state[self.params["out"]] = float(params[self.param_name("sigma")]) * x if self.toggles["enabled"] else 0.0
        return state


class _Shaper(BaseBlock):
    def __init__(self, params, name):
        require(params, "in", "out", "curve", "offset_rad", block=name or "vtf1_shaper")
        super().__init__("vtf1_shaper", name, params, reads=(params["in"],), writes=(params["out"],))
        self.curve = as_callable(params["curve"], ("x",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"])

    def __call__(self, state, params):
        state[self.params["out"]] = self.curve(float(state[self.params["in"]])) + self.params["offset_rad"]
        return state


@registry.register_block("vtf1_noise")
def make_vtf1_noise(params, *, name=None):
    return _Noise(params, name)


@registry.register_block("vtf1_shaper")
def make_vtf1_shaper(params, *, name=None):
    return _Shaper(params, name)


@registry.register_fn("vtf1_affine")
def vtf1_affine(x, a, b):
    return a * x + b


def _run_block(block, run: RunSettings):
    lay = SlotLayout(run.f_sim_hz)
    for leaf in iter_leaves(block):
        leaf.init(lay, stream(run.seed, leaf.name))
    lay.freeze()
    block.bind(lay)
    st, p = State(lay), Params(lay)
    out = {s: np.empty(run.n_samples) for s in block.writes}
    for n in range(run.n_samples):
        st["n"] = n
        st["t"] = n * lay.dt_s
        block(st, p)
        for s in out:
            out[s][n] = float(st[s])
    return out


def _run_spec(spec):
    return _run_block(from_spec(apply_toggles(spec.pipeline, spec.toggles)), spec.run)


RUN = RunSettings(duration_s=2.0e-4, f_sim_hz=1.0e5, chunk_samples=8, seed=20261031)

YAML = """
preset: vtf1
run: {duration_s: 2.0e-4, f_sim_hz: 1.0e5, chunk_samples: 8, seed: 20261031}
topology:
  segments:
    - {id: arm1, kind: pm, L_m: 8.4}
    - {id: arm2, kind: pm, L_m: !expr "arm1.L_m - dL", where: {dL: 9.0}}
  locks:
    - {id: lockD, combo: {lo: 1, pr: -1}, setpoint: !pi_frac 0.5}
pipeline:
  - {block: vtf1_noise, name: laser.noise, params: {out: laser.dx1.phi, sigma: 0.3, toggles: {enabled: true}}}
  - {block: vtf1_noise, name: env.noise, params: {out: drv.env.dT, sigma: 2.0, toggles: {enabled: true}}}
  - block: vtf1_shaper
    name: s1
    params: {in: laser.dx1.phi, out: v.s1, offset_rad: !pi_frac 0.25,
             curve: !expr {source: "k * x + c", where: {k: 2.0, c: 0.5}, args: [x]}}
  - parallel:
      - {block: vtf1_shaper, name: p1, params: {in: v.s1, out: ctrl.p1.u, offset_rad: 0.0,
                                                curve: !fn {name: vtf1_affine, args: {a: -1.5, b: 0.1}}}}
      - {block: vtf1_shaper, name: p2, rate_div: 3, params: {in: v.s1, out: ctrl.p2.u, offset_rad: 0.0, curve: 0.25}}
    name: branches
  - stack:
      stages:
        - {block: vtf1_shaper, name: c1, params: {in: ctrl.p1.u, out: ctrl.c1.u, offset_rad: 0.0, curve: !fn {name: vtf1_affine, args: {a: 0.5, b: 0.0}}}}
        - {block: vtf1_shaper, name: c2, params: {in: ctrl.c1.u, out: ctrl.c2.u, offset_rad: 0.0, curve: !fn {name: vtf1_affine, args: {a: 0.5, b: 0.0}}}}
      actuators:
        - {block: vtf1_shaper, name: a1, params: {in: ctrl.c1.u, out: act.a1.x, offset_rad: 0.0, curve: !expr {source: "sin(x)", args: [x]}}}
        - seq:
            - {block: vtf1_shaper, name: a2, params: {in: ctrl.c2.u, out: act.a2.x, offset_rad: 0.0, curve: !expr {source: "x * x", args: [x]}}}
          name: a2chain
toggles: {env.noise.enabled: false}
measure: {esa_note: {f_s_hz: 3.0e6}}
"""


def _python_pipeline():
    from functools import partial

    shaper = make_vtf1_shaper
    noise_l = make_vtf1_noise({"out": "laser.dx1.phi", "sigma": 0.3, "toggles": {"enabled": True}}, name="laser.noise")
    noise_e = make_vtf1_noise({"out": "drv.env.dT", "sigma": 2.0, "toggles": {"enabled": True}}, name="env.noise")
    s1 = shaper({"in": "laser.dx1.phi", "out": "v.s1", "offset_rad": math.pi / 4,
                 "curve": compile_expr("k * x + c", ("x",), {"k": 2.0, "c": 0.5})}, name="s1")
    p1 = shaper({"in": "v.s1", "out": "ctrl.p1.u", "offset_rad": 0.0,
                 "curve": partial(vtf1_affine, a=-1.5, b=0.1)}, name="p1")
    p2 = shaper({"in": "v.s1", "out": "ctrl.p2.u", "offset_rad": 0.0, "curve": 0.25}, name="p2")
    p2.set_rate_div(3)
    half = partial(vtf1_affine, a=0.5, b=0.0)
    c1 = shaper({"in": "ctrl.p1.u", "out": "ctrl.c1.u", "offset_rad": 0.0, "curve": half}, name="c1")
    c2 = shaper({"in": "ctrl.c1.u", "out": "ctrl.c2.u", "offset_rad": 0.0, "curve": half}, name="c2")
    a1 = shaper({"in": "ctrl.c1.u", "out": "act.a1.x", "offset_rad": 0.0,
                 "curve": compile_expr("sin(x)", ("x",))}, name="a1")
    a2 = shaper({"in": "ctrl.c2.u", "out": "act.a2.x", "offset_rad": 0.0,
                 "curve": compile_expr("x * x", ("x",))}, name="a2")
    return [noise_l, noise_e, s1, Parallel([p1, p2], name="branches"),
            Stack([c1, c2], [a1, Seq([a2], name="a2chain")])]


def _python_spec(pipeline):
    return config.to_spec(
        pipeline,
        run=RUN,
        topology={"segments": [{"id": "arm1", "kind": "pm", "L_m": 8.4}, {"id": "arm2", "kind": "pm", "L_m": 8.4 - 9.0}],
                  "locks": [{"id": "lockD", "combo": {"lo": 1, "pr": -1}, "setpoint": math.pi / 2}]},
        toggles={"env.noise.enabled": False},
        measure={"esa_note": {"f_s_hz": 3.0e6}},
        preset="vtf1",
    )


def _assert_identical(a, b):
    assert a.keys() == b.keys()
    for k in a:
        assert a[k].dtype == b[k].dtype == np.float64
        assert np.array_equal(a[k], b[k]), k  # bit-identical (no NaNs here)
        assert a[k].tobytes() == b[k].tobytes(), k


def test_yaml_and_python_give_identical_specs(tmp_path):
    py_spec = _python_spec(_python_pipeline())
    yaml_spec = config.load(loads_yaml(YAML))
    assert py_spec == yaml_spec
    assert config.spec_hash(py_spec) == config.spec_hash(yaml_spec)
    f = tmp_path / "vtf1.yaml"
    f.write_text(YAML)
    assert config.load(f) == yaml_spec


def test_round_trip_spec_equals_load_dump_spec(tmp_path):
    spec = config.load(loads_yaml(YAML))
    text = config.dump_yaml(spec)
    path = tmp_path / "dumped.yaml"
    path.write_text(text)
    again = config.load(path)
    assert again == spec
    assert config.dump_yaml(again) == text  # dump is a fixed point
    py_spec = _python_spec(_python_pipeline())
    assert config.load(loads_yaml(config.dump_yaml(py_spec))) == py_spec


def test_identical_built_pipelines():
    py_spec = _python_spec(_python_pipeline())
    yaml_spec = config.load(loads_yaml(YAML))
    a = from_spec(apply_toggles(py_spec.pipeline, py_spec.toggles))
    b = from_spec(apply_toggles(yaml_spec.pipeline, yaml_spec.toggles))
    assert a.spec == b.spec
    la, lb = list(iter_leaves(a)), list(iter_leaves(b))
    assert [(x.kind, x.name, x.reads, x.writes, x.rate_div, x.runtime_params, dict(x.toggles)) for x in la] == [
        (x.kind, x.name, x.reads, x.writes, x.rate_div, x.runtime_params, dict(x.toggles)) for x in lb
    ]
    assert dict(a.toggles) == dict(b.toggles) and dict(a.toggles)["env.noise.enabled"] is False


def test_same_seed_bit_identical_results():
    pipeline = _python_pipeline()
    py_spec = _python_spec(pipeline)
    yaml_spec = config.load(loads_yaml(YAML))
    # the Python objects themselves (real callables, toggle override applied by rebuilding env.noise)
    direct = Seq(pipeline[:1] + [make_vtf1_noise({"out": "drv.env.dT", "sigma": 2.0, "toggles": {"enabled": False}},
                                                 name="env.noise")] + pipeline[2:])
    r_direct = _run_block(direct, RUN)
    r_py = _run_spec(py_spec)
    r_yaml = _run_spec(yaml_spec)
    r_rt = _run_spec(config.load(loads_yaml(config.dump_yaml(yaml_spec))))
    _assert_identical(r_direct, r_py)
    _assert_identical(r_py, r_yaml)
    _assert_identical(r_yaml, r_rt)
    assert np.all(r_yaml["drv.env.dT"] == 0.0) and np.std(r_yaml["laser.dx1.phi"]) > 0  # toggle took effect
    # a different seed changes the result (the check is not vacuous)
    other = _run_spec(config.load({**loads_yaml(YAML), "run": {**loads_yaml(YAML)["run"], "seed": 1}}))
    assert not np.array_equal(other["laser.dx1.phi"], r_yaml["laser.dx1.phi"])


def test_extends_composition_equals_the_flat_file(tmp_path):
    """WP-F2b: a config that ``extends`` a base and edits it (merge, ``!delete``, ``!replace``) gives
    the same spec as the flat file; its dump carries no ``extends`` and reloads to the same spec."""
    flat = config.load(loads_yaml(YAML))
    base = YAML.replace("toggles: {env.noise.enabled: false}", "toggles: {laser.noise.enabled: false}")
    base = base.replace("measure: {esa_note: {f_s_hz: 3.0e6}}", "measure: {esa_note: {f_s_hz: 1.0}, old: {x: 1}}")
    (tmp_path / "base.yaml").write_text(base)
    child = tmp_path / "child.yaml"
    child.write_text("extends: ./base.yaml\n"
                     "toggles: !replace {env.noise.enabled: false}\n"
                     "measure:\n  esa_note: {f_s_hz: 3.0e6}\n  old: !delete\n")
    spec = config.load(child)
    assert spec == flat and config.spec_hash(spec) == config.spec_hash(flat)
    text = config.dump_yaml(spec)
    assert "\nextends:" not in text and config.load(loads_yaml(text)) == flat
    assert config.load(config.read_input(child)) == flat
    _assert_identical(_run_spec(spec), _run_spec(flat))


def test_unnamed_lambda_cannot_round_trip():
    pipeline = _python_pipeline()
    pipeline[2] = make_vtf1_shaper({"in": "laser.dx1.phi", "out": "v.s1", "offset_rad": 0.0,
                                    "curve": lambda x: 2 * x}, name="s1")
    with pytest.raises(config.ConfigError, match=r"block 's1'.*curve.*unnamed / unregistered callable"):
        _python_spec(pipeline)


def test_with_runner_run_when_available():
    from cranberry.runner.run import run

    spec = config.load(loads_yaml(YAML))
    try:
        r1 = run(spec)
    except NotImplementedError as e:
        assert "WP-F1" in str(e)
        pytest.skip("runner.run not implemented yet (WP-F1)")
    r2 = run(config.load(loads_yaml(config.dump_yaml(spec))))
    for k in r1.series:
        assert np.array_equal(r1.series[k], r2.series[k])
