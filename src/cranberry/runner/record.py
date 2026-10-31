"""Run records: config, seeds, git commit, environment, results (protected list "logging of config,
seed, and commit hash"; PLAN §2.11, prompt §6). Owner: WP-F2.

Layout of ``runs/<UTC time>_<preset>_<spechash8>/`` (``runs/`` is gitignored), file names fixed here.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec

__all__ = [
    "RUNS_DIR",
    "RECORD_FILES",
    "make_run_dir",
    "git_info",
    "env_info",
    "write_record",
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
}


def make_run_dir(spec: RunSpec, root: Path = RUNS_DIR) -> Path:
    """Create and return ``root/<UTC time>_<preset>_<spechash8>``. WP-F2."""
    raise NotImplementedError("runner.record.make_run_dir is implemented by WP-F2")


def git_info() -> dict[str, Any]:
    """{"commit": sha, "dirty": bool, "diff_sha256": ..., "diff": ...} of the working tree. WP-F2."""
    raise NotImplementedError("runner.record.git_info is implemented by WP-F2")


def env_info() -> dict[str, Any]:
    """Versions of Python, numpy, scipy, numba, llvmlite, pyyaml and the platform. WP-F2."""
    raise NotImplementedError("runner.record.env_info is implemented by WP-F2")


def write_record(run_dir: Path, **parts: Any) -> None:
    """Write the parts named by :data:`RECORD_FILES` keys into ``run_dir``. WP-F2."""
    raise NotImplementedError("runner.record.write_record is implemented by WP-F2")
