"""``RunSpec``: the pure-data configuration tree (PLAN §2.1, §2.9). Shape fixed by WP-F0; loading,
overlays, tags and serialisation are WP-F2.

YAML and Python-list pipelines both become a ``RunSpec``; ``build`` turns it into blocks, and
``to_spec`` turns Python-built blocks back into one. Round trip: ``spec == load(dump_yaml(spec))``
and the same seed gives bit-identical results (VT-F1).

Top-level YAML keys ↔ fields::

    preset:   <str>                       # name of the preset this config started from (or null)
    overlays: [<path or axis/value>, ...] # applied in order by WP-F2 (G12); kept for the record
    run:      {duration_s, f_sim_hz, chunk_samples, seed}
    topology: {segments: [...], beams: [...], locks: [...], ports: [...], observers: [...], coef: ...}
    offline:  {...}                       # offline-stage inputs (sources, fibers, detectors, ...)
    pipeline: [...]                       # core.pipeline dict form → SeqSpec
    toggles:  {"<block name>.<toggle>": bool, ...}   # overrides of block toggles (ablate)
    measure:  {...}                       # measurement / metric settings

D34: there are no code defaults for any modelling choice; every field of :class:`RunSettings` is
required and presets give every value. ``topology``/``offline``/``measure`` are kept as plain mappings
here (parsed by WP-F6 ``path.topology.Topology.from_spec``, WP-F3 and WP-F4 respectively).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from cranberry.core.pipeline import SeqSpec

__all__ = ["RunSettings", "RunSpec"]


@dataclass(frozen=True)
class RunSettings:
    """``run:`` section.

    Attributes:
        duration_s: simulated time [s].
        f_sim_hz: base sample rate of the time loop [Hz] (G10, per preset).
        chunk_samples: samples per outer-loop chunk (results must not depend on it, VT-F2).
        seed: root seed of all RNG streams (G11), >= 0.
    """

    duration_s: float
    f_sim_hz: float
    chunk_samples: int
    seed: int

    def __post_init__(self) -> None:
        if not self.duration_s > 0:
            raise ValueError("run.duration_s must be > 0")
        if not self.f_sim_hz > 0:
            raise ValueError("run.f_sim_hz must be > 0")
        if isinstance(self.chunk_samples, bool) or int(self.chunk_samples) != self.chunk_samples or self.chunk_samples < 1:
            raise ValueError("run.chunk_samples must be an int >= 1")
        if isinstance(self.seed, bool) or int(self.seed) != self.seed or self.seed < 0:
            raise ValueError("run.seed must be an int >= 0")

    @property
    def n_samples(self) -> int:
        """Total number of base samples, round(duration_s · f_sim_hz)."""
        return int(round(self.duration_s * self.f_sim_hz))


@dataclass(frozen=True)
class RunSpec:
    """Complete, resolved configuration of one run (after overlays and callables, WP-F2)."""

    run: RunSettings
    pipeline: SeqSpec
    topology: Mapping[str, Any] = field(default_factory=dict)
    offline: Mapping[str, Any] = field(default_factory=dict)
    toggles: Mapping[str, bool] = field(default_factory=dict)
    measure: Mapping[str, Any] = field(default_factory=dict)
    preset: str | None = None
    overlays: tuple[str, ...] = ()
