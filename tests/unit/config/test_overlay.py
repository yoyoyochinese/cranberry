"""Preset + overlay deep merge (G12) and overlay application in ``config.load``."""

import pytest

from cranberry import config
from cranberry.config import ConfigError
from cranberry.config.overlay import REPO_ROOT, display_ref, find_overlay, merge
from cranberry.config.yaml_io import Sweep


def test_mappings_merge_recursively_and_null_deletes():
    base = {"a": {"x": 1, "y": {"p": 1, "q": 2}}, "b": 1, "c": 3}
    ov = {"a": {"y": {"q": 20, "r": 30}, "z": 5}, "b": None, "d": {"k": None, "m": 1}}
    out = merge(base, ov)
    assert out == {"a": {"x": 1, "y": {"p": 1, "q": 20, "r": 30}, "z": 5}, "c": 3, "d": {"m": 1}}
    assert base == {"a": {"x": 1, "y": {"p": 1, "q": 2}}, "b": 1, "c": 3}  # inputs untouched


def test_plain_lists_are_replaced():
    assert merge({"l": [1, 2, 3]}, {"l": [9]}) == {"l": [9]}
    assert merge({"l": [{"a": 1}]}, {"l": [{"a": 2}]}) == {"l": [{"a": 2}]}  # no key → replace


def test_keyed_lists_merge_by_id_or_name():
    base = {"segments": [{"id": "s1", "L_m": 1.0, "zone": "lab"}, {"id": "s2", "L_m": 2.0}]}
    ov = {"segments": [{"id": "s2", "L_m": 5.0}, {"id": "s3", "L_m": 3.0, "tmp": None}, {"id": "s1", "zone": None}]}
    assert merge(base, ov) == {
        "segments": [{"id": "s1", "L_m": 1.0}, {"id": "s2", "L_m": 5.0}, {"id": "s3", "L_m": 3.0}]
    }
    pipe = [{"block": "pid", "name": "lock1.pid", "params": {"kp": 1.0, "ki": 2.0}}, {"block": "dac", "name": "dac"}]
    out = merge({"pipeline": pipe}, {"pipeline": [{"name": "lock1.pid", "params": {"kp": 3.0}}]})
    assert out["pipeline"][0] == {"block": "pid", "name": "lock1.pid", "params": {"kp": 3.0, "ki": 2.0}}
    assert out["pipeline"][1] == {"block": "dac", "name": "dac"}


def test_mixed_key_list_is_replaced_and_duplicate_keys_error():
    base = {"l": [{"name": "a"}, {"stack": {}}]}
    assert merge(base, {"l": [{"name": "a", "x": 1}]}) == {"l": [{"name": "a", "x": 1}]}
    with pytest.raises(ConfigError, match="duplicate"):
        merge({"l": [{"id": "a"}, {"id": "a"}]}, {"l": [{"id": "a", "x": 1}]})


def test_markers_are_leaves():
    assert merge({"a": {"lo": 1}}, {"a": Sweep(0, 1)}) == {"a": Sweep(0, 1)}
    assert merge({"a": Sweep(0, 1)}, {"a": {"lo": 1}}) == {"a": {"lo": 1}}


def _write(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


BASE = """
preset: testbase
run: {duration_s: 1.0e-3, f_sim_hz: 1.0e6, chunk_samples: 64, seed: 1}
topology:
  segments:
    - {id: s1, L_m: 1.0}
pipeline:
  - {block: cfgtest_noise, name: n1, params: {out: drv.n1.x, sigma: 1.0, toggles: {enabled: true}}}
"""


def test_load_applies_overlays_in_order_and_records_them(tmp_path):
    base = _write(tmp_path / "base.yaml", BASE + "overlays: [ov1.yaml]\n")
    _write(tmp_path / "ov1.yaml", "topology: {segments: [{id: s1, L_m: 2.0}]}\nrun: {seed: 7}\n")
    ov2 = _write(tmp_path / "sub" / "ov2.yaml", "topology: {segments: [{id: s2, L_m: 3.0}]}\nrun: {seed: 8}\n")
    spec = config.load(base, overlays=[ov2])
    assert spec.run.seed == 8  # CLI overlay after in-file overlays
    assert spec.topology["segments"] == [{"id": "s1", "L_m": 2.0}, {"id": "s2", "L_m": 3.0}]
    assert spec.overlays == (str((tmp_path / "ov1.yaml").resolve()), str(ov2.resolve()))
    # dump records them as applied_overlays; loading the dump does not re-apply anything
    text = config.dump_yaml(spec)
    assert "applied_overlays" in text and "\noverlays:" not in text
    assert config.load(_write(tmp_path / "resolved.yaml", text)) == spec
    # read_input gives a self-contained input document equivalent to (file, overlays)
    assert config.load(config.read_input(base, [ov2])) == spec


def test_overlay_errors(tmp_path):
    base = _write(tmp_path / "base.yaml", BASE)
    with pytest.raises(ConfigError, match="not found"):
        config.load(base, overlays=["nope/missing"])
    bad = _write(tmp_path / "bad.yaml", "overlays: [x.yaml]\n")
    with pytest.raises(ConfigError, match="cannot list overlays"):
        config.load(base, overlays=[bad])
    lst = _write(tmp_path / "lst.yaml", "- 1\n")
    with pytest.raises(ConfigError, match="mapping"):
        config.load(base, overlays=[lst])


def test_overlay_search_path(tmp_path, monkeypatch):
    axes = REPO_ROOT / "presets" / "axes"
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "local.yaml", "a: 1\n")
    assert find_overlay("local.yaml") == (tmp_path / "local.yaml").resolve()
    assert find_overlay("presets/axes/.gitkeep") == (axes / ".gitkeep").resolve()
    assert display_ref(axes / ".gitkeep") == "presets/axes/.gitkeep"
    assert display_ref(tmp_path / "local.yaml") == str((tmp_path / "local.yaml").resolve())
