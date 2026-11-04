"""D34 guard (Phase 6 fix of reviews/phase6_wp6g_ab.md): the structural refactors of the controller / actuator params
(G-3/G-4, D91) and the config-layer conflict detection must not change the numbers of any existing preset.

Golden digests in ``data/preset_identity.json`` were computed at dev 30aef9f, before the refactor:

* ``spec``: SHA-256 of the canonical resolved spec tree. Asserted for the presets that do not extend F_common
  (Ref A/B/C/D): their config must load to the very same tree.
* ``run``: SHA-256 of every recorded series (raw float bytes), the metric summary (repr of every float) and the CHK
  event codes of a short compiled run (metric settle time + 2000 samples). Asserted for every preset without sweep
  axes: the F_* nominals (whose specs gain explicit controller keys equal to the old device defaults, so the spec
  digest changes but the run must be bit-identical), Ref A, Ref A tap, refB/C/D_nominal; the bad presets in the slow
  suite.

Regenerate (only after a deliberate, documented numerical change):
``uv run python tests/integration/test_preset_identity.py --regen``.
"""

from __future__ import annotations

import dataclasses
import glob
import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).with_name("data") / "preset_identity.json"
F_COMMON_USERS = ("presets/future/", "presets/bad/")


def presets() -> list[str]:
    pats = ("presets/reference/*.yaml", "presets/future/*_nominal.yaml", "presets/bad/*.yaml")
    return sorted(str(Path(p).relative_to(ROOT)) for pat in pats for p in glob.glob(str(ROOT / pat)))


def _settle(spec) -> float:
    from cranberry.config.spec import iter_block_specs

    return max((float(b.params.get("t_settle_s") or 0.0) for _, b in iter_block_specs(spec.pipeline)), default=0.0)


def spec_digest(path: str) -> str:
    from cranberry import config
    from cranberry.config.spec import spec_to_tree
    from cranberry.config.yaml_io import canonical_json

    return hashlib.sha256(canonical_json(spec_to_tree(config.load(ROOT / path))).encode()).hexdigest()


def run_digest(path: str) -> str | None:
    from cranberry import config
    from cranberry.runner.run import run

    spec = config.load(ROOT / path)
    if config.sweep_markers(spec):
        return None
    dur = _settle(spec) + 2000.0 / spec.run.f_sim_hz
    measure = {k: v for k, v in (spec.measure or {}).items() if k not in ("esa", "pulse_sample")}  # sweep-time hooks
    spec = dataclasses.replace(spec, run=dataclasses.replace(spec.run, duration_s=dur), measure=measure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run(spec, mode="compiled")
    h = hashlib.sha256()
    for k in sorted(r.series):
        h.update(k.encode())
        h.update(np.ascontiguousarray(np.asarray(r.series[k])).tobytes())
    summ = r.metrics.summary() if r.metrics is not None else {}
    h.update(json.dumps({k: repr(v) for k, v in sorted(summ.items())}).encode())
    h.update(json.dumps(sorted(e["code"] for e in r.events)).encode())
    return h.hexdigest()


def _golden() -> dict:
    return json.loads(DATA.read_text()) if DATA.exists() else {}


def _is_f(path: str) -> bool:
    return path.startswith(F_COMMON_USERS)


@pytest.mark.parametrize("path", [p for p in presets() if not _is_f(p)])
def test_reference_spec_unchanged(path):
    assert spec_digest(path) == _golden()[path]["spec"]


def _run_cases(bad: bool):
    g = _golden()
    return [p for p in presets() if p in g and g[p]["run"] is not None and p.startswith("presets/bad/") == bad]


@pytest.mark.numba
@pytest.mark.integration
@pytest.mark.parametrize("path", _run_cases(bad=False))
def test_preset_run_bit_identical(path):
    assert run_digest(path) == _golden()[path]["run"]


@pytest.mark.slow
@pytest.mark.numba
@pytest.mark.parametrize("path", _run_cases(bad=True))
def test_bad_preset_run_bit_identical(path):
    assert run_digest(path) == _golden()[path]["run"]


if __name__ == "__main__":
    if "--regen" not in sys.argv:
        sys.exit("usage: test_preset_identity.py --regen")
    out = {}
    for p in presets():
        out[p] = {"spec": spec_digest(p), "run": run_digest(p)}
        print(p, out[p], flush=True)
    DATA.parent.mkdir(exist_ok=True)
    DATA.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
