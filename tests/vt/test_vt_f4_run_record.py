"""VT-F4 (PLAN §6.2, §2.11, prompt §6; protected list "logging of config, seed, and commit hash"):
a run record holds the input and resolved config, the seeds, the git commit (+ dirty flag and diff),
the environment and the results, and the run can be reproduced from the record alone (tolerance:
bit-identical results).

Test blocks (kind ``vtf4_noise``) stand in for physics blocks; the loop is the pure-Python reference
step, and ``runner.run(spec, record_root=...)`` (WP-F1 + WP-F2b wiring) is checked end to end: it
writes the record itself and re-running the recorded config reproduces the stored results.
"""

import json
import subprocess

import numpy as np
import pytest

from cranberry import config
from cranberry.config.spec import apply_toggles
from cranberry.core import registry
from cranberry.core.block import BaseBlock, require
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import from_spec, iter_leaves
from cranberry.core.rng import stream, stream_key
from cranberry.core.state import Params, State
from cranberry.runner import record


class _Noise(BaseBlock):
    toggle_names = ("enabled",)

    def __init__(self, params, name):
        require(params, "out", block=name or "vtf4_noise")
        super().__init__("vtf4_noise", name, params, writes=(params["out"],), runtime_params=("sigma",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"])

    def __call__(self, state, params):
        x = self.rng.standard_normal()
        state[self.params["out"]] = float(params[self.param_name("sigma")]) * x if self.toggles["enabled"] else 0.0
        return state


@registry.register_block("vtf4_noise")
def make_vtf4_noise(params, *, name=None):
    return _Noise(params, name)


def _simulate(spec):
    block = from_spec(apply_toggles(spec.pipeline, spec.toggles))
    lay = SlotLayout(spec.run.f_sim_hz)
    for leaf in iter_leaves(block):
        leaf.init(lay, stream(spec.run.seed, leaf.name))
    lay.freeze()
    block.bind(lay)
    st, p = State(lay), Params(lay)
    out = {s: np.empty(spec.run.n_samples) for s in block.writes}
    for n in range(spec.run.n_samples):
        st["n"] = n
        block(st, p)
        for s in out:
            out[s][n] = float(st[s])
    return out


CONFIG = """
preset: vtf4
run: {duration_s: 1.0e-4, f_sim_hz: 1.0e6, chunk_samples: 16, seed: 424242}
topology:
  segments:
    - {id: arm1, L_m: 8.4}
    - {id: arm2, L_m: !expr "arm1.L_m - dL", where: {dL: 9.0}}
pipeline:
  - {block: vtf4_noise, name: laser.pn, params: {out: laser.dx1.phi, sigma: 0.1, toggles: {enabled: true}}}
  - {block: vtf4_noise, name: env.pn, params: {out: drv.arm1.dT, sigma: 1.0, toggles: {enabled: true}}}
"""
OVERLAY = "pipeline: [{name: env.pn, params: {sigma: 0.5}}]\nrun: {seed: 7}\n"


def _files(tmp_path):
    cfg = tmp_path / "cfg" / "vtf4.yaml"
    cfg.parent.mkdir()
    cfg.write_text(CONFIG)
    ov = tmp_path / "cfg" / "half_env.yaml"
    ov.write_text(OVERLAY)
    return cfg, ov


def test_record_is_complete(tmp_path):
    cfg, ov = _files(tmp_path)
    spec = config.load(cfg, overlays=[ov])
    results = _simulate(spec)
    run_dir = record.record_run(spec, tmp_path / "runs", config_input=config.read_input(cfg, [ov]),
                                results=results, metrics={"std_phi": float(np.std(results["laser.dx1.phi"]))})
    assert run_dir.parent == tmp_path / "runs" and "_vtf4_" in run_dir.name
    for key in ("config_input", "config_resolved", "seeds", "git", "env", "results", "metrics"):
        assert (run_dir / record.RECORD_FILES[key]).is_file(), key

    # seeds: root + every block's stream key
    seeds = json.loads((run_dir / "seeds.json").read_text())
    assert seeds == {"root": 7, "streams": {"laser.pn": list(stream_key("laser.pn")), "env.pn": list(stream_key("env.pn"))}}
    # git: commit hash of the checkout holding the package, dirty flag, diff when dirty
    git = json.loads((run_dir / "git.json").read_text())
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=record._PKG_DIR, capture_output=True, text=True).stdout.strip()
    assert git["commit"] == head and len(head) == 40
    assert isinstance(git["dirty"], bool) and (git["diff"] is not None) == git["dirty"]
    # env: interpreter and dependency versions
    env = json.loads((run_dir / "env.json").read_text())
    assert env["versions"]["numpy"] == np.__version__
    for name in ("scipy", "numba", "llvmlite", "pyyaml"):
        assert env["versions"][name], name
    assert env["python_version"]
    # configs: resolved = the spec; input = the config as written + overlays (location independent)
    assert config.load(run_dir / "config.resolved.yaml") == spec
    assert config.load(run_dir / "config.input.yaml") == spec
    text_in = (run_dir / "config.input.yaml").read_text()
    assert "!expr" in text_in and "overlays" in text_in  # input keeps expressions and the overlay list


def test_record_reproduces_the_result(tmp_path):
    cfg, ov = _files(tmp_path)
    spec = config.load(cfg, overlays=[ov])
    results = _simulate(spec)
    run_dir = record.record_run(spec, tmp_path / "runs", config_input=config.read_input(cfg, [ov]), results=results)

    replay_spec = record.load_record(run_dir)
    assert replay_spec == spec
    replay = _simulate(replay_spec)
    stored = np.load(run_dir / "results.npz")
    assert set(stored.files) == set(replay)
    for k in replay:
        assert stored[k].tobytes() == replay[k].tobytes(), k  # bit-identical
    # the record is tamper-evident w.r.t. the seed
    (run_dir / "seeds.json").write_text(json.dumps({"root": 1, "streams": {}}))
    with pytest.raises(config.ConfigError, match="seeds.json"):
        record.load_record(run_dir)


def test_python_built_run_records_spec_as_input(tmp_path):
    spec = config.load(config.yaml_io.loads_yaml(CONFIG))
    run_dir = record.record_run(spec, tmp_path / "runs")
    assert config.load(run_dir / "config.input.yaml") == spec == record.load_record(run_dir)
    with pytest.raises(FileExistsError):  # records are never overwritten
        record.write_record(run_dir, git={})


def test_runner_run_without_record_root_writes_nothing(tmp_path, monkeypatch):
    from cranberry.runner.run import run

    monkeypatch.chdir(tmp_path)
    result = run(config.load(config.yaml_io.loads_yaml(CONFIG)))
    assert result.run_dir is None and list(tmp_path.iterdir()) == []  # D34: no silent default directory


REC_OVERLAY = "measure: {record: {slots: [laser.dx1.phi, drv.arm1.dT], decim: 1}}\n"


def test_recorded_runner_run_end_to_end(tmp_path):
    """WP-F2b: ``run(spec, record_root=...)`` writes the complete record, and re-running the recorded
    config reproduces the stored results bit for bit (also equal to the reference loop)."""
    from cranberry.runner.run import run

    cfg, ov = _files(tmp_path)
    rec = tmp_path / "cfg" / "rec.yaml"
    rec.write_text(REC_OVERLAY)
    spec = config.load(cfg, overlays=[ov, rec])
    result = run(spec, record_root=tmp_path / "runs", config_input=config.read_input(cfg, [ov, rec]))
    run_dir = result.run_dir
    assert run_dir is not None and run_dir.parent == tmp_path / "runs" and "_vtf4_" in run_dir.name
    for key in ("config_input", "config_resolved", "seeds", "git", "env", "results"):
        assert (run_dir / record.RECORD_FILES[key]).is_file(), key
    seeds = json.loads((run_dir / "seeds.json").read_text())
    assert seeds["root"] == 7 == spec.run.seed
    assert json.loads((run_dir / "git.json").read_text())["commit"]
    assert config.load(run_dir / "config.input.yaml") == spec == record.load_record(run_dir)

    stored = np.load(run_dir / "results.npz")
    assert set(stored.files) == {"laser.dx1.phi", "drv.arm1.dT", "t_s"}
    replay = run(record.load_record(run_dir))  # from the record alone
    reference = _simulate(spec)
    for k in ("laser.dx1.phi", "drv.arm1.dT"):
        assert stored[k].tobytes() == replay.series[k].tobytes() == result.series[k].tobytes(), k
        assert stored[k].tobytes() == reference[k].tobytes(), k  # runner == reference loop
    assert stored["t_s"].tobytes() == replay.t_s.tobytes()
    assert np.std(stored["laser.dx1.phi"]) > 0


def test_runner_run_from_a_path_records_the_file_as_input(tmp_path):
    from cranberry.runner.run import run

    cfg, _ = _files(tmp_path)
    result = run(cfg, record_root=tmp_path / "runs")
    text_in = (result.run_dir / "config.input.yaml").read_text()
    assert "!expr" in text_in  # the input as written, not the resolved spec
    assert config.load(result.run_dir / "config.input.yaml") == config.load(cfg) == record.load_record(result.run_dir)
