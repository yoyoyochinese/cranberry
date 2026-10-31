"""offline.tables: Table interpolation, spec_hash, npz + JSON cache with hash invalidation (G16-g)."""

import json

import numpy as np
import pytest

from cranberry.offline.tables import CACHE_DIR, Table, TableCache, canonical_json, code_version, spec_hash


def bilinear_table(value_shape=()):
    x = np.array([0.0, 1.0, 2.5, 4.0])
    y = np.array([-1.0, 0.0, 3.0])
    X, Y = np.meshgrid(x, y, indexing="ij")
    base = 1.0 + 2.0 * X - 0.5 * Y + 0.25 * X * Y
    values = base.reshape(base.shape + (1,) * len(value_shape)) * np.ones(value_shape) if value_shape else base
    return Table("bilin", {"x_hz": x, "y_s": y}, values, "1", {"sources": ["[test]"]})


def test_table_validation_and_read_only():
    t = bilinear_table()
    assert t.axes == ("x_hz", "y_s") and t.value_shape == ()
    assert not t.values.flags.writeable and not t.grid["x_hz"].flags.writeable
    with pytest.raises(ValueError):
        Table("bad", {"x": [0.0, 0.0, 1.0]}, np.zeros(3), "1")
    with pytest.raises(ValueError):
        Table("bad", {"x": [0.0, 1.0]}, np.zeros(3), "1")
    with pytest.raises(ValueError):
        Table("bad", {}, np.zeros(3), "1")
    src = np.zeros(2)
    t2 = Table("copy", {"x": [0.0, 1.0]}, src, "1")
    src[0] = 5.0
    assert t2.values[0] == 0.0  # copied


def test_interp_multilinear_exact_and_broadcast():
    t = bilinear_table()
    f = lambda x, y: 1.0 + 2.0 * x - 0.5 * y + 0.25 * x * y  # bilinear: reproduced exactly
    xs = np.array([0.0, 0.3, 1.7, 4.0])
    ys = np.array([-1.0, 2.2, 0.5, 3.0])
    np.testing.assert_allclose(t.interp(x_hz=xs, y_s=ys), f(xs, ys), rtol=1e-14)
    out = t.interp(x_hz=xs[:, None], y_s=ys[None, :])
    assert out.shape == (4, 4)
    np.testing.assert_allclose(out, f(xs[:, None], ys[None, :]), rtol=1e-14)
    assert t.interp(x_hz=2.5, y_s=0.0).shape == ()
    with pytest.raises(ValueError):
        t.interp(x_hz=4.1, y_s=0.0)  # no extrapolation
    with pytest.raises(ValueError):
        t.interp(x_hz=1.0)
    with pytest.raises(ValueError):
        t.interp(x_hz=1.0, y_s=0.0, z=1.0)


def test_interp_value_shape_and_single_point_axis():
    t = bilinear_table(value_shape=(2, 2))
    v = t.interp(x_hz=np.array([0.5, 1.5]), y_s=1.0)
    assert v.shape == (2, 2, 2)
    np.testing.assert_allclose(v[:, 0, 1], 1.0 + 2.0 * np.array([0.5, 1.5]) - 0.5 + 0.25 * np.array([0.5, 1.5]))
    # single-point axis (e.g. CW limit table with one Δt point, T03)
    s = Table("one", {"dt_s": [0.0], "f_hz": [1.0, 3.0]}, np.array([[10.0, 30.0]]), "1")
    assert float(s.interp(dt_s=0.0, f_hz=2.0)) == pytest.approx(20.0)
    with pytest.raises(ValueError):
        s.interp(dt_s=1e-9, f_hz=2.0)


def test_spec_hash_properties():
    a = {"source": {"L": 0.08, "alpha_per_w": 9.06, "grid": [1, 2, 3]}, "name": "vpm"}
    b = {"name": "vpm", "source": {"grid": (1, 2, 3), "alpha_per_w": np.float64(9.06), "L": 0.08}}
    h = spec_hash(a, "0.1.0|b/1")
    assert len(h) == 64 and h == spec_hash(b, "0.1.0|b/1")  # order / tuple / numpy-scalar independent
    assert h != spec_hash({**a, "source": {**a["source"], "L": 0.0800000001}}, "0.1.0|b/1")
    assert h != spec_hash(a, "0.1.0|b/2")  # code version is part of the key
    assert spec_hash({"x": np.arange(3.0)}, "v") != spec_hash({"x": [0.0, 1.0, 2.0]}, "v")  # arrays tagged
    assert spec_hash({"x": float("nan")}, "v") == spec_hash({"x": float("nan")}, "v")
    with pytest.raises(TypeError):
        spec_hash({"f": lambda x: x}, "v")
    with pytest.raises(TypeError):
        spec_hash({1: 2}, "v")
    with pytest.raises(ValueError):
        spec_hash(a, "")
    assert canonical_json({"b": 1, "a": [1.5, None, True]}) == '{"a":[1.5,null,true],"b":1}'


def test_code_version():
    v = code_version("opo_cw_spectrum", 3)
    assert v.endswith("|opo_cw_spectrum/3")
    with pytest.raises(ValueError):
        code_version("", 1)


def test_cache_get_or_build_and_invalidation(tmp_path):
    cache = TableCache(tmp_path / "_cache")
    calls = []

    def builder(scale):
        def build():
            calls.append(scale)
            x = np.linspace(0, 1, 5)
            return Table("lin", {"x_hz": x}, scale * x, "rad", {"sources": ["[test]"], "code_version": "v1"})

        return build

    sub = {"table": "lin", "scale": 2.0}
    key = spec_hash(sub, "v1")
    t1 = cache.get_or_build(key, builder(2.0))
    assert calls == [2.0] and t1.provenance["spec_hash"] == key
    npz, js = cache.paths(key)
    assert npz.is_file() and js.is_file()
    meta = json.loads(js.read_text())
    assert meta["spec_hash"] == key and meta["unit"] == "rad" and meta["axes"] == ["x_hz"]
    t2 = cache.get_or_build(key, builder(2.0))
    assert calls == [2.0]  # served from cache
    np.testing.assert_array_equal(t2.values, t1.values)
    np.testing.assert_array_equal(t2.grid["x_hz"], t1.grid["x_hz"])
    assert t2.provenance["sources"] == ["[test]"] and t2.name == "lin"

    # changed config subtree → new key → rebuilt (G16-g invalidation on hash change)
    key2 = spec_hash({**sub, "scale": 3.0}, "v1")
    t3 = cache.get_or_build(key2, builder(3.0))
    assert calls == [2.0, 3.0] and float(t3.interp(x_hz=1.0)) == pytest.approx(3.0)
    # changed code version → new key too
    key3 = spec_hash(sub, "v2")
    assert key3 != key
    cache.get_or_build(key3, builder(2.0))
    assert calls == [2.0, 3.0, 2.0]

    # corrupt / mismatching entries are rebuilt
    js.write_text("{not json")
    cache.get_or_build(key, builder(2.0))
    assert calls[-1] == 2.0 and len(calls) == 4
    meta = json.loads(js.read_text())
    meta["spec_hash"] = key2
    js.write_text(json.dumps(meta))
    assert cache.load(key) is None

    # no temporary files left behind; invalidate / clear
    assert not [p for p in cache.root.iterdir() if p.suffix == ".tmp"]
    assert cache.invalidate(key) and not cache.has(key)
    assert cache.clear() == 2 and list(cache.root.iterdir()) == []


def test_cache_rejects_bad_input(tmp_path):
    cache = TableCache(tmp_path)
    with pytest.raises(ValueError):
        cache.paths("not-a-hash")
    key = spec_hash({"a": 1}, "v")
    with pytest.raises(TypeError):
        cache.get_or_build(key, lambda: np.zeros(3))
    wrong = Table("t", {"x": [0.0, 1.0]}, np.zeros(2), "1", {"spec_hash": "0" * 64})
    with pytest.raises(ValueError):
        cache.get_or_build(key, lambda: wrong)
    assert CACHE_DIR.parts == ("runs", "_cache")
