"""``run(spec) -> RunResult``: offline stage → build → time loop (chunks) → measurement → record
(protected list "experiment runner"; PLAN §2.8). Owner: WP-F1 (records written via WP-F2's
``runner.record``).

The runner is a pure function of the spec: same spec (incl. seed) → bit-identical results,
independent of ``chunk_samples`` (VT-F2).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from cranberry.blocks.meas.metrics import MetricsResult
    from cranberry.config.spec import RunSpec

__all__ = ["RunResult", "run"]


@dataclass
class RunResult:
    """Result of one run.

    Attributes:
        spec: the resolved spec that was run.
        run_dir: ``runs/<UTC time>_<preset>_<spechash8>/`` or None when recording is off.
        t_s: time axis of the recorded (decimated) series [s].
        series: {slot name: decimated time series} (true path phases, errors, u, δψ, ...).
        psd: {name: (f [Hz], one-sided PSD)}.
        metrics: metric (i)/(ii), overlap, stray light, η (``blocks.meas.metrics.MetricsResult``).
        events: lock / CHK events as JSON-serialisable dicts.
        table_hashes: {table name: hash} of offline tables used.
        seeds: {"root": entropy, "<stream name>": key, ...} (G11).
    """

    spec: RunSpec
    run_dir: Path | None
    t_s: np.ndarray
    series: Mapping[str, np.ndarray] = field(default_factory=dict)
    psd: Mapping[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    metrics: MetricsResult | None = None
    events: list[dict[str, Any]] = field(default_factory=list)
    table_hashes: Mapping[str, str] = field(default_factory=dict)
    seeds: Mapping[str, Any] = field(default_factory=dict)


def run(spec: RunSpec) -> RunResult:
    """Run one fully resolved configuration (no sweep markers left). Implemented by WP-F1."""
    raise NotImplementedError("runner.run is implemented by WP-F1")
