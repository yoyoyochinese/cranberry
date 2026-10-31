"""``runner.record`` (unit level; the full record/reproduce check is VT-F4) and CLI wiring of ``load``."""

import json
import subprocess

import numpy as np
import pytest

from cranberry import __main__ as cli
from cranberry import config
from cranberry.runner import record

CFG = """
preset: clitest
run: {duration_s: 1.0e-5, f_sim_hz: 1.0e6, chunk_samples: 4, seed: 2}
pipeline:
  - {block: cfgtest_noise, name: n1, params: {out: drv.n1.x, sigma: 1.0, toggles: {enabled: true}}}
"""


def _spec():
    return config.load(config.yaml_io.loads_yaml(CFG))


def test_run_dir_name_and_uniqueness(tmp_path):
    spec = _spec()
    a, b = record.make_run_dir(spec, tmp_path), record.make_run_dir(spec, tmp_path)
    assert a != b and a.is_dir() and b.is_dir()
    stamp, preset, h = a.name.split("_")[:3]
    assert stamp.endswith("Z") and preset == "clitest" and h[:8] == config.spec_hash(spec)[:8]


def test_git_and_env_info():
    g = record.git_info()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=record._PKG_DIR, capture_output=True, text=True).stdout.strip()
    assert g["commit"] == head and isinstance(g["dirty"], bool)
    assert (g["diff"] is not None) == g["dirty"] and (g["diff_sha256"] is not None) == g["dirty"]
    e = record.env_info()
    assert e["versions"]["numpy"] == np.__version__ and e["versions"]["pyyaml"]
    assert e["python_version"] and "numpy" in e["packages"]


def test_git_info_outside_a_repo_records_the_error(tmp_path, monkeypatch):
    monkeypatch.setattr(record, "_PKG_DIR", tmp_path)
    g = record.git_info()
    assert g["commit"] is None and "error" in g


def test_write_record_kinds_and_no_overwrite(tmp_path):
    spec = _spec()
    record.write_record(tmp_path, config_resolved=spec, metrics={"a": np.float64(1.5), "v": np.arange(2)},
                        results={"x": np.arange(3.0)}, config_input="raw: text\n")
    assert config.load(tmp_path / "config.resolved.yaml") == spec
    assert json.loads((tmp_path / "metrics.json").read_text()) == {"a": 1.5, "v": [0, 1]}
    assert np.load(tmp_path / "results.npz")["x"].tolist() == [0.0, 1.0, 2.0]
    assert (tmp_path / "config.input.yaml").read_text() == "raw: text\n"
    with pytest.raises(FileExistsError):
        record.write_record(tmp_path, metrics={})
    with pytest.raises(KeyError, match="unknown record part"):
        record.write_record(tmp_path, typo={})


def test_seeds_info_lists_every_block_stream():
    from cranberry.core.rng import stream_key

    s = record.seeds_info(_spec())
    assert s == {"root": 2, "streams": {"n1": list(stream_key("n1"))}}


def test_cli_loads_validates_and_delegates(tmp_path, capsys):
    good = tmp_path / "c.yaml"
    good.write_text(CFG)
    # WP-F1 has landed: run and sweep now execute
    assert cli.main(["run", str(good), "--seed", "5"]) == 0
    assert cli.main(["sweep", str(good)]) == 0
    ov = tmp_path / "ov.yaml"
    ov.write_text("pipeline: [{name: n1, params: {toggles: null}}]\n")
    assert cli.main(["run", str(good), "--overlay", str(ov)]) == 2
    err = capsys.readouterr().err
    assert "block 'n1'" in err and "enabled" in err and "D34" in err
    bad = tmp_path / "bad.yaml"
    bad.write_text("preset: x\n")
    assert cli.main(["run", str(bad)]) == 2
    assert "missing required section" in capsys.readouterr().err
