"""Run records: config, seeds, git commit, environment, results (protected list "logging of config,
seed, and commit hash"; PLAN §2.11, prompt §6). Owner: WP-F2.

Layout of ``runs/<UTC time>_<preset>_<spechash8>/`` (``runs/`` is gitignored), file names fixed here.

Use by the runner (``runner.run.run(spec, record_root=...)``, wired in WP-F2b)::

    run_dir = record_run(spec, root, config_input=read_input(path, overlays))   # before the loop
    ...
    write_record(run_dir, results={...arrays...}, metrics={...}, tables={...}, ...)  # after

``config.input.yaml`` is the input as written (``extends``/overlays as references, dB values as
written); ``config.resolved.yaml`` is the spec actually run. They load to the same spec unless the
run was started with an override (CLI ``--seed``) or is a sweep point, which ``invocation.json``
then states.

Files are written once (exclusive create): a record is never silently overwritten. A run is
reproduced from its directory with :func:`load_record` (VT-F4).
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import re
import subprocess
import sys
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec

__all__ = [
    "RUNS_DIR",
    "RECORD_FILES",
    "make_run_dir",
    "git_info",
    "env_info",
    "seeds_info",
    "write_record",
    "record_run",
    "load_record",
]

RUNS_DIR = Path("runs")

#: Files of a run directory (PLAN §2.11).
RECORD_FILES: dict[str, str] = {
    "config_input": "config.input.yaml",
    "config_resolved": "config.resolved.yaml",
    "seeds": "seeds.json",
    "git": "git.json",
    "env": "env.json",
    "tables": "tables.json",
    "results": "results.npz",
    "metrics": "metrics.json",
    "loss_registry": "loss_registry.json",
    "band_ledger": "band_ledger.json",
    # WP-F2b: how the run was started when the resolved config alone does not say it (CLI argv and
    # overrides such as --seed, the coordinates of a sweep point, the seed mode).
    "invocation": "invocation.json",
}

#: Distributions whose versions are listed first in env.json (all installed ones follow).
_KEY_PACKAGES = ("cranberry", "numpy", "scipy", "numba", "llvmlite", "pyyaml", "matplotlib")
_PKG_DIR = Path(__file__).resolve().parents[1]


def make_run_dir(spec: RunSpec, root: Path = RUNS_DIR) -> Path:
    """Create and return ``root/<UTC time>_<preset>_<spechash8>`` (a numeric suffix if it exists)."""
    from cranberry.config import spec_hash

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S.%fZ")
    preset = re.sub(r"[^A-Za-z0-9_.-]", "_", spec.preset or "nopreset")
    base = Path(root) / f"{stamp}_{preset}_{spec_hash(spec)[:8]}"
    Path(root).mkdir(parents=True, exist_ok=True)
    path, k = base, 0
    while True:
        try:
            path.mkdir()
            return path
        except FileExistsError:
            k += 1
            path = base.with_name(f"{base.name}-{k}")


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=_PKG_DIR, check=True, capture_output=True, text=True, errors="replace", timeout=30
    ).stdout


def git_info() -> dict[str, Any]:
    """{"commit": sha, "dirty": bool, "diff_sha256": ..., "diff": ..., "untracked": [...]} of the
    working tree holding the ``cranberry`` package. Outside a git checkout: commit/dirty are None and
    "error" says why (the record is still written, never silently without this file)."""
    try:
        commit = _git("rev-parse", "HEAD").strip()
        status = _git("status", "--porcelain")
        diff = _git("diff", "HEAD", "--binary")
        branch = _git("rev-parse", "--abbrev-ref", "HEAD").strip()
    except (OSError, subprocess.SubprocessError) as e:
        return {"commit": None, "dirty": None, "diff_sha256": None, "diff": None, "untracked": None,
                "error": f"{type(e).__name__}: {e}"}
    untracked = [line[3:] for line in status.splitlines() if line.startswith("??")]
    dirty = bool(status.strip())
    return {
        "commit": commit,
        "branch": branch,
        "dirty": dirty,
        "diff_sha256": hashlib.sha256(diff.encode("utf-8")).hexdigest() if dirty else None,
        "diff": diff if dirty else None,
        "untracked": untracked,
    }


def env_info() -> dict[str, Any]:
    """Versions of Python, numpy, scipy, numba, llvmlite, pyyaml and the platform (+ all packages)."""
    packages: dict[str, str] = {}
    for dist in importlib.metadata.distributions():
        name = dist.metadata["Name"]
        if name:
            packages[name.lower()] = dist.version
    key = {name: packages.get(name) for name in _KEY_PACKAGES}
    return {
        "python": sys.version,
        "python_version": platform.python_version(),
        "implementation": platform.python_implementation(),
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "versions": key,
        "packages": dict(sorted(packages.items())),
    }


def seeds_info(spec: RunSpec) -> dict[str, Any]:
    """Root seed and the stream key of every block (``core.rng.stream_key(block name)``, G11)."""
    from cranberry.config.spec import iter_block_specs
    from cranberry.core.rng import stream_key

    return {
        "root": int(spec.run.seed),
        "streams": {b.name: list(stream_key(b.name)) for _, b in iter_block_specs(spec.pipeline)},
    }


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "to_json") and callable(obj.to_json):  # LossRegistry / BandLedger (WP-F3)
        return _jsonable(obj.to_json())
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if hasattr(obj, "__dataclass_fields__"):
        from dataclasses import asdict

        return _jsonable(asdict(obj))
    raise TypeError(f"cannot write {type(obj).__name__} to a JSON record")


def _yaml_text(value: Any) -> str:
    from cranberry.config import dump_yaml
    from cranberry.config.spec import RunSpec
    from cranberry.config.yaml_io import dump_yaml_obj

    if isinstance(value, str):
        return value
    if isinstance(value, RunSpec):
        return dump_yaml(value)
    return dump_yaml_obj(value)


def write_record(run_dir: Path, **parts: Any) -> None:
    """Write the parts named by :data:`RECORD_FILES` keys into ``run_dir``.

    ``*.yaml`` parts take text, a ``RunSpec`` or a config tree; ``*.json`` parts take JSON-able data
    (numpy values converted; objects with ``to_json()`` are asked for it); ``results`` takes a mapping
    of arrays (``numpy.savez``). Existing files are never overwritten (``FileExistsError``).
    """
    unknown = sorted(set(parts) - set(RECORD_FILES))
    if unknown:
        raise KeyError(f"unknown record part(s) {unknown}; known: {sorted(RECORD_FILES)}")
    run_dir = Path(run_dir)
    for key, value in parts.items():
        path = run_dir / RECORD_FILES[key]
        if path.suffix == ".yaml":
            with path.open("x", encoding="utf-8") as f:
                f.write(_yaml_text(value))
        elif path.suffix == ".json":
            text = json.dumps(_jsonable(value), indent=2, sort_keys=True, allow_nan=True)
            with path.open("x", encoding="utf-8") as f:
                f.write(text + "\n")
        elif path.suffix == ".npz":
            if not isinstance(value, Mapping):
                raise TypeError(f"record part {key!r} must be a mapping of arrays")
            with path.open("xb") as f:
                np.savez(f, **{str(k): np.asarray(v) for k, v in value.items()})
        else:  # pragma: no cover - RECORD_FILES only has the three kinds
            raise ValueError(f"unsupported record file {path.name}")


def record_run(spec: RunSpec, root: Path = RUNS_DIR, *, config_input: Any = None, **parts: Any) -> Path:
    """Create the run directory and write config.input/resolved, seeds, git and env (+ ``parts``).

    Args:
        spec: the resolved spec of this run.
        root: parent of run directories.
        config_input: the input document (``config.read_input(...)``); None → the spec itself
            (Python-built pipelines have no separate input).
        parts: further :data:`RECORD_FILES` parts available now (results, metrics, ...).
    """
    run_dir = make_run_dir(spec, root)
    write_record(
        run_dir,
        config_input=spec if config_input is None else config_input,
        config_resolved=spec,
        seeds=seeds_info(spec),
        git=git_info(),
        env=env_info(),
        **parts,
    )
    return run_dir


def load_record(run_dir: Path) -> RunSpec:
    """The ``RunSpec`` of a recorded run (from ``config.resolved.yaml``; checked against seeds.json)."""
    from cranberry.config import ConfigError, load

    run_dir = Path(run_dir)
    spec = load(run_dir / RECORD_FILES["config_resolved"])
    seeds_path = run_dir / RECORD_FILES["seeds"]
    if seeds_path.is_file():
        root = json.loads(seeds_path.read_text(encoding="utf-8"))["root"]
        if root != spec.run.seed:
            raise ConfigError(f"{run_dir}: seeds.json root {root} != run.seed {spec.run.seed}")
    return spec
