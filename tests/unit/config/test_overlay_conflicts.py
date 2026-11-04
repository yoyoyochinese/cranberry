"""Phase 6 fix (reviews/phase6_wp6g_ab.md §5, G-3..G-6): config.load detects conflicts between overlays (a later
discard of what an earlier overlay wrote; block keys merged into a seq/stack) and names both; ``overrides:`` accepts an
intended discard; scalar point overrides stay allowed. Also the structured alias ``!expr "<path to a list>"``."""

from __future__ import annotations

import pytest

from cranberry import config
from cranberry.config import ConfigError
from cranberry.config.resolve import resolve
from cranberry.config.yaml_io import loads_yaml

BASE = {
    "run": {"duration_s": 1e-5, "f_sim_hz": 1e6, "chunk_samples": 1000, "seed": 1},
    "offline": {"params": {"a": 1.0, "lst": [1.0, 2.0]}},
    "pipeline": [{"block": "ramp", "name": "x", "params": {"out": "v.x", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                           "period_s": 1.0, "offset_v": 1.0, "toggles": {"enabled": True}}}],
}


def ov(tmp_path, name, text):
    p = tmp_path / f"{name}.yaml"
    p.write_text(text)
    return str(p)


SET_X = "pipeline:\n  - name: x\n    params: {offset_v: 2.0}\n"
REPLACE_X = ("pipeline:\n  - !replace\n    name: x\n    seq:\n      - block: ramp\n        name: x\n        params: {out: v.x, "
             "waveform: sawtooth, amplitude_v: 0.0, period_s: 1.0, offset_v: 3.0, toggles: {enabled: true}}\n")


def test_replace_after_modify_names_both(tmp_path):
    a, b = ov(tmp_path, "a", SET_X), ov(tmp_path, "b", REPLACE_X)
    with pytest.raises(ConfigError) as e:
        config.load(BASE, overlays=[a, b])
    assert "a.yaml" in str(e.value) and "b.yaml" in str(e.value) and "pipeline.x" in str(e.value)


def test_params_into_replaced_seq_names_both(tmp_path):
    a, b = ov(tmp_path, "a", REPLACE_X), ov(tmp_path, "b", SET_X)
    with pytest.raises(ConfigError) as e:
        config.load(BASE, overlays=[a, b])
    assert "a.yaml" in str(e.value) and "b.yaml" in str(e.value) and "seq / stack" in str(e.value)


def test_list_replaced_by_later_overlay(tmp_path):
    a, b = ov(tmp_path, "a", "offline:\n  params: {lst: [3.0]}\n"), ov(tmp_path, "b", "offline:\n  params: {lst: [4.0]}\n")
    with pytest.raises(ConfigError, match="a.yaml"):
        config.load(BASE, overlays=[a, b])


def test_overrides_declaration_accepts(tmp_path):
    a = ov(tmp_path, "a", SET_X)
    b = ov(tmp_path, "b", f"overrides: [{a}]\n" + REPLACE_X)
    spec = config.load(BASE, overlays=[a, b])
    assert spec.pipeline.children[0].children[0].params["offset_v"] == 3.0
    star = ov(tmp_path, "c", "overrides: ['*']\n" + REPLACE_X)
    config.load(BASE, overlays=[a, star])


def test_scalar_point_override_and_base_changes_are_not_conflicts(tmp_path):
    a = ov(tmp_path, "a", "offline:\n  params: {a: !sweep [1.0, 2.0, 2]}\n")
    b = ov(tmp_path, "b", "offline:\n  params: {a: 5.0}\n")
    assert config.load(BASE, overlays=[a, b]).offline["params"]["a"] == 5.0
    # one overlay replacing a node of the BASE config, and a refinement inside an earlier overlay's seq
    r = ov(tmp_path, "r", REPLACE_X)
    inner = ov(tmp_path, "i", "pipeline:\n  - name: x\n    seq:\n      - name: x\n        params: {offset_v: 7.0}\n")
    spec = config.load(BASE, overlays=[r, inner])
    assert spec.pipeline.children[0].children[0].params["offset_v"] == 7.0


def test_structured_alias():
    t = loads_yaml("""
offline:
  params:
    f: 3.0
    q: !sweep [1, 2, 2]
    modes: [{f_hz: !expr "f", q: 4.0, r: 1.0}]
    deferred: [{f_hz: !expr "f", q: !expr "offline.params.q", r: 1.0}]
a: {modes: !expr "offline.params.modes", d: !expr "offline.params.deferred", x: !expr "2 * offline.params.f"}
""")
    r = resolve(t)
    assert r["a"]["modes"] == [{"f_hz": 3.0, "q": 4.0, "r": 1.0}] and r["a"]["x"] == 6.0
    assert type(r["a"]["d"]).__name__ == "Expr"  # deferred until the sweep fixes q
    with pytest.raises(ConfigError, match="contains expressions"):
        resolve(loads_yaml("p: {m: [{f: !expr 'g'}], g: 1.0}\nb: !expr '2 * p.m'\n"))
