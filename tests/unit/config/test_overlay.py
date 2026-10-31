"""Preset + overlay deep merge (G12) and overlay application in ``config.load``."""

import pytest

from cranberry import config
from cranberry.config import ConfigError
from cranberry.config.overlay import REPO_ROOT, display_ref, find_overlay, merge
from cranberry.config.yaml_io import Sweep, loads_yaml


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


# ------------------------------------------------------------------- WP-F2b: !delete / !replace


def test_delete_mapping_key_and_keyed_element():
    base = {"a": {"x": 1, "y": 2}, "segs": [{"id": "s1", "L_m": 1.0}, {"id": "s2"}, {"id": "s3"}]}
    # in flow style a bare tag needs a space before ',' / '}' (YAML syntax): {x: !delete }
    ov = loads_yaml("a: {x: !delete }\nsegs: [!delete s2, !delete {id: s3}, {id: s4}]\n")
    assert merge(base, ov) == {"a": {"y": 2}, "segs": [{"id": "s1", "L_m": 1.0}, {"id": "s4"}]}
    assert merge({"a": 1, "b": 2}, loads_yaml("a: !delete\n")) == {"b": 2}
    pipe = [{"block": "pid", "name": "lock1.pid"}, {"block": "dac", "name": "dac"}]
    assert merge({"pipeline": pipe}, loads_yaml("pipeline: [!delete lock1.pid]\n")) == {"pipeline": [pipe[1]]}


def test_delete_errors_when_the_target_does_not_exist():
    with pytest.raises(ConfigError, match="nope: !delete of a key the base does not have"):
        merge({"a": 1}, loads_yaml("nope: !delete\n"))
    with pytest.raises(ConfigError, match=r"segs: !delete 's9': no element with id 's9'"):
        merge({"segs": [{"id": "s1"}]}, loads_yaml("segs: [!delete s9]\n"))
    with pytest.raises(ConfigError, match="keyed list"):
        merge({"l": [1, 2]}, loads_yaml("l: [!delete x]\n"))
    with pytest.raises(ConfigError, match="keyed-list elements"):
        merge({"a": 1}, loads_yaml("a: !delete x\n"))
    with pytest.raises(ConfigError, match="nothing to act on"):  # inside a newly added subtree
        merge({"a": 1}, loads_yaml("b:\n  c: !delete\n"))


def test_replace_mapping_list_and_keyed_element():
    base = {"p": {"kp": 1.0, "ki": 2.0}, "segs": [{"id": "s1", "L_m": 1.0, "zone": "lab"}, {"id": "s2"}]}
    ov = loads_yaml("p: !replace {kd: 3.0}\nsegs: [!replace {id: s1, L_m: 5.0}, {id: s2, L_m: 2.0}]\n")
    assert merge(base, ov) == {"p": {"kd": 3.0}, "segs": [{"id": "s1", "L_m": 5.0}, {"id": "s2", "L_m": 2.0}]}
    # a whole keyed list replaced instead of merged by key
    assert merge(base, loads_yaml("segs: !replace [{id: s9}]\n"))["segs"] == [{"id": "s9"}]


def test_replace_errors_when_the_target_does_not_exist():
    with pytest.raises(ConfigError, match="q: !replace of a key the base does not have"):
        merge({"p": {}}, loads_yaml("q: !replace {a: 1}\n"))
    with pytest.raises(ConfigError, match=r"!replace id='s9': no such element"):
        merge({"segs": [{"id": "s1"}]}, loads_yaml("segs: [!replace {id: s9, L_m: 1.0}]\n"))
    with pytest.raises(Exception, match="mapping or a list"):
        loads_yaml("a: !replace 3\n")


def test_directives_in_overlay_files_and_errors_without_a_base(tmp_path):
    base = _write(tmp_path / "base.yaml", BASE + "offline: {tap: {x: 1}, old: 2}\n")
    ov = _write(tmp_path / "ov.yaml", "offline:\n  old: !delete\n  tap: !replace {y: 3}\n"
                                      "topology: {segments: [!delete s1]}\n")
    spec = config.load(base, overlays=[ov])
    assert spec.offline == {"tap": {"y": 3}} and spec.topology == {"segments": []}
    typo = _write(tmp_path / "typo.yaml", "offline:\n  olld: !delete\n")
    with pytest.raises(ConfigError, match=r"overlay .*typo\.yaml: offline: olld: !delete"):
        config.load(base, overlays=[typo])
    lone = _write(tmp_path / "lone.yaml", BASE + "offline:\n  x: !delete\n")
    with pytest.raises(ConfigError, match=r"!delete at offline\.x .*nothing to act on"):
        config.load(lone)
    # the input document keeps the directives and round-trips through YAML
    from cranberry.config.yaml_io import dump_yaml_obj

    assert loads_yaml(dump_yaml_obj(loads_yaml(ov.read_text()))) == loads_yaml(ov.read_text())


# ------------------------------------------------------------------------------- WP-F2b: extends


@pytest.fixture
def presets(tmp_path, monkeypatch):
    """A private presets/ tree (find_base searches overlay.PRESETS_DIR)."""
    from cranberry.config import overlay

    root = tmp_path / "presets"
    monkeypatch.setattr(overlay, "PRESETS_DIR", root)
    _write(root / "reference" / "refT.yaml", BASE + "offline: {tap: {il_db: 3.0, k: 1}}\n")
    _write(root / "axes" / "seed" / "s9.yaml", "run: {seed: 9}\n")
    return root


def test_extends_by_name_merges_own_keys_then_overlays(tmp_path, presets):
    child = _write(tmp_path / "cfg" / "child.yaml", """
extends: refT
preset: refT_variant
offline: {tap: {k: 2}}
pipeline: [{name: n1, params: {sigma: 0.5}}]
overlays: [ov.yaml]
""")
    _write(tmp_path / "cfg" / "ov.yaml", "offline: {tap: {k: 3}}\n")
    spec = config.load(child)
    assert spec.preset == "refT_variant"
    assert spec.offline == {"tap": {"il_lin": pytest.approx(10 ** 0.3), "k": 3}}  # base, own keys, overlay
    assert spec.pipeline.children[0].params["sigma"] == 0.5 and spec.run.seed == 1
    text = config.dump_yaml(spec)
    assert "\nextends:" not in text and config.load(_write(tmp_path / "dump.yaml", text)) == spec
    # read_input pins the base to its file (location independent) and reproduces the spec
    inp = config.read_input(child)
    assert inp["extends"] == str((presets / "reference" / "refT.yaml").resolve())
    assert config.load(inp) == spec
    # preset: is inherited when the child does not set it
    bare = _write(tmp_path / "cfg" / "bare.yaml", "extends: refT\n")
    assert config.load(bare).preset == "testbase"


def test_extends_by_path_chain_and_base_overlays(tmp_path, presets):
    _write(presets / "reference" / "refU.yaml", "extends: refT.yaml\noverlays: [seed/s9]\n"
                                                "offline: {tap: {k: 5}}\n")
    child = _write(tmp_path / "c.yaml", f"extends: {presets / 'reference' / 'refU.yaml'}\n"
                                        "offline: {tap: {il_db: !delete }}\n")
    spec = config.load(child)
    assert spec.run.seed == 9 and spec.offline == {"tap": {"k": 5}}
    assert spec.overlays == (str((presets / "axes" / "seed" / "s9.yaml").resolve()),)  # base's overlay recorded


def test_extends_errors(tmp_path, presets):
    _write(presets / "a" / "loopA.yaml", "extends: loopB\n")
    _write(presets / "b" / "loopB.yaml", "extends: loopA\n")
    with pytest.raises(ConfigError, match=r"circular chain .*loopA\.yaml -> .*loopB\.yaml -> .*loopA\.yaml"):
        config.load(presets / "a" / "loopA.yaml")
    with pytest.raises(ConfigError, match="no preset named 'nope'"):
        config.load(_write(tmp_path / "x.yaml", "extends: nope\n"))
    _write(presets / "future" / "refT.yaml", "run: {}\n")
    with pytest.raises(ConfigError, match="ambiguous"):
        config.load(_write(tmp_path / "y.yaml", "extends: refT\n"))
    with pytest.raises(ConfigError, match="not found"):
        config.load(_write(tmp_path / "z.yaml", "extends: ../missing.yaml\n"))
    with pytest.raises(ConfigError, match="cannot list overlays or extends"):
        config.load(presets / "reference" / "refT.yaml", overlays=[_write(tmp_path / "o.yaml", "extends: refT\n")])
