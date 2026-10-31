"""runner.run minimal path (WP-F1): RunSpec → build → chunks → recorded series."""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.pipeline import spec_from_obj
from cranberry.core.testing import toy
from cranberry.runner.run import run

OBJ = [toy("toy_sine", "s", out="v.s", unit="V", amp=2.0, freq_hz=1e3),
       toy("toy_gain", "g", in_="v.s", out="v.g", unit="V", gain=0.5)]


def spec(measure=None, chunk=64):
    return RunSpec(run=RunSettings(duration_s=1e-3, f_sim_hz=1e5, chunk_samples=chunk, seed=3),
                   pipeline=spec_from_obj(OBJ), measure=measure or {})


@pytest.mark.numba
def test_run_records_and_returns_final_state():
    res = run(spec({"record": {"slots": ["v.g", "t"], "decim": 4}}))
    assert res.mode == "compiled" and res.run_dir is None
    np.testing.assert_allclose(res.t_s, np.arange(0, 100, 4) / 1e5)
    np.testing.assert_allclose(res.series["v.g"], np.sin(2 * np.pi * 1e3 * res.t_s), rtol=1e-12, atol=1e-15)
    assert res.final_state.n == 100
    assert set(res.seeds) == {"root", "s", "g"}


def test_run_without_record_and_reference_mode():
    res = run(spec(), mode="reference", keep_state=False)
    assert res.series == {} and len(res.t_s) == 0 and res.final_state is None and res.mode == "reference"


def test_run_rejects_bad_inputs():
    with pytest.raises(ValueError):
        run(spec(), mode="fast")
    with pytest.raises(KeyError, match="decim"):
        run(spec({"record": {"slots": ["v.g"]}}))
    with pytest.raises(ValueError):
        run(spec({"record": {"slots": ["v.g"], "decim": 0}}))


def test_run_validates_the_spec_first():
    import dataclasses

    from cranberry.config import ConfigError
    from cranberry.config.yaml_io import Sweep

    with pytest.raises(ConfigError, match="toggles.g.nope"):
        run(dataclasses.replace(spec(), toggles={"g.nope": True}))
    swept = spec_from_obj([OBJ[0], toy("toy_gain", "g", in_="v.s", out="v.g", unit="V", gain=Sweep(0.0, 1.0, 3))])
    with pytest.raises(ConfigError, match="runner.sweep"):
        run(dataclasses.replace(spec(), pipeline=swept))


def test_run_writes_a_record_only_when_asked(tmp_path):
    import json

    from cranberry import config
    from cranberry.runner import record

    res = run(spec({"record": {"slots": ["v.g"], "decim": 4}}), mode="reference", record_root=tmp_path / "runs",
              invocation={"argv": ["run", "x.yaml"]})
    assert res.run_dir is not None and res.run_dir.parent == tmp_path / "runs"
    names = {p.name for p in res.run_dir.iterdir()}
    assert names == {"config.input.yaml", "config.resolved.yaml", "seeds.json", "git.json", "env.json",
                     "results.npz", "invocation.json"}
    stored = np.load(res.run_dir / "results.npz")
    assert stored["v.g"].tobytes() == res.series["v.g"].tobytes() and stored["t_s"].tobytes() == res.t_s.tobytes()
    assert record.load_record(res.run_dir) == res.spec == config.load(res.run_dir / "config.input.yaml")
    assert json.loads((res.run_dir / "invocation.json").read_text()) == {"argv": ["run", "x.yaml"]}
    assert run(spec(), mode="reference").run_dir is None
    assert len(list((tmp_path / "runs").iterdir())) == 1  # nothing written without record_root
