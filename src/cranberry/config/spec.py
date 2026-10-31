"""``RunSpec``: the pure-data configuration tree (PLAN §2.1, §2.9). Shape fixed by WP-F0; loading,
overlays, tags and serialisation are WP-F2.

YAML and Python-list pipelines both become a ``RunSpec``; ``build`` turns it into blocks, and
``to_spec`` turns Python-built blocks back into one. Round trip: ``spec == load(dump_yaml(spec))``
and the same seed gives bit-identical results (VT-F1).

Top-level YAML keys ↔ fields::

    preset:   <str>                       # name of the preset this config started from (or null)
    overlays: [<path or axis/value>, ...] # applied in order by WP-F2 (G12); kept for the record
    applied_overlays: [...]               # overlays already merged (written by dump_yaml; not re-applied)
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

import copy
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field, fields
from typing import Any

from cranberry.core.block import BlockSpec
from cranberry.core.pipeline import ParallelSpec, PipelineSpec, SeqSpec, StackSpec, spec_from_obj, spec_to_obj

__all__ = [
    "ConfigError",
    "RunSettings",
    "RunSpec",
    "TOP_LEVEL_KEYS",
    "RUN_KEYS",
    "spec_from_tree",
    "spec_to_tree",
    "iter_block_specs",
    "apply_toggles",
]


class ConfigError(ValueError):
    """A configuration is invalid (missing modelling choice, unknown key, bad reference, ...).

    Messages name the offending path / block and key (D34: no silent defaults)."""


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


#: Allowed top-level keys of a config tree. ``overlays`` are applied by ``config.load`` (G12);
#: ``applied_overlays`` records overlays already merged (written by ``dump_yaml``, never re-applied).
TOP_LEVEL_KEYS: tuple[str, ...] = (
    "preset", "overlays", "applied_overlays", "run", "topology", "offline", "pipeline", "toggles", "measure",
)
#: Keys of ``run:``; all required (D34).
RUN_KEYS: tuple[str, ...] = tuple(f.name for f in fields(RunSettings))


def _is_plain_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def spec_from_tree(tree: Mapping[str, Any]) -> RunSpec:
    """Build a :class:`RunSpec` from a merged and resolved config tree (overlays already applied).

    Raises :class:`ConfigError` naming the key for unknown top-level keys, missing ``run`` /
    ``pipeline`` sections or ``run`` keys (D34), and malformed sections.
    """
    if not isinstance(tree, Mapping):
        raise ConfigError(f"config must be a mapping, got {type(tree).__name__}")
    unknown = [k for k in tree if k not in TOP_LEVEL_KEYS]
    if unknown:
        raise ConfigError(f"unknown top-level key(s) {unknown}; allowed: {list(TOP_LEVEL_KEYS)}")
    if tree.get("overlays"):
        raise ConfigError("'overlays' must be applied before building a RunSpec (use config.load)")
    missing = [k for k in ("run", "pipeline") if k not in tree]
    if missing:
        raise ConfigError(f"missing required section(s) {missing} (D34: presets give every value)")
    run = tree["run"]
    if not isinstance(run, Mapping):
        raise ConfigError("run: must be a mapping {duration_s, f_sim_hz, chunk_samples, seed}")
    missing = [k for k in RUN_KEYS if k not in run]
    if missing:
        raise ConfigError(f"missing run.{', run.'.join(missing)} (D34: no code defaults; presets give every value)")
    unknown = [k for k in run if k not in RUN_KEYS]
    if unknown:
        raise ConfigError(f"unknown key(s) in run: {unknown}; allowed: {list(RUN_KEYS)}")
    for k in RUN_KEYS:
        if not _is_plain_number(run[k]):
            raise ConfigError(
                f"run.{k} must be a number, got {run[k]!r} (sweep markers are not allowed in run:; "
                "use runner.sweep(spec, axes={'run.%s': [...]}))" % k
            )
    if any(isinstance(run[k], float) and not float(run[k]).is_integer() for k in ("chunk_samples", "seed")):
        raise ConfigError("run.chunk_samples and run.seed must be integers")
    try:
        settings = RunSettings(
            duration_s=float(run["duration_s"]),
            f_sim_hz=float(run["f_sim_hz"]),
            chunk_samples=int(run["chunk_samples"]),
            seed=int(run["seed"]),
        )
    except ValueError as e:
        raise ConfigError(str(e)) from None
    pipe = tree["pipeline"]
    if not isinstance(pipe, (list, Mapping)):
        raise ConfigError("pipeline: must be a list of blocks / composites")
    try:
        pipeline = spec_from_obj(pipe)
    except (ValueError, TypeError) as e:
        raise ConfigError(f"pipeline: {e}") from None
    if not isinstance(pipeline, SeqSpec):
        pipeline = SeqSpec((pipeline,))
    sections: dict[str, Any] = {}
    for key in ("topology", "offline", "measure", "toggles"):
        value = tree.get(key)
        if value is None:
            value = {}
        if not isinstance(value, Mapping):
            raise ConfigError(f"{key}: must be a mapping, got {type(value).__name__}")
        sections[key] = dict(value)
    for k, v in sections["toggles"].items():
        if not isinstance(k, str) or "." not in k:
            raise ConfigError(f"toggles: key {k!r} must be '<block name>.<toggle>'")
        if not isinstance(v, bool):
            raise ConfigError(f"toggles.{k} must be true/false, got {v!r}")
    preset = tree.get("preset")
    if preset is not None and not isinstance(preset, str):
        raise ConfigError(f"preset: must be a name (string), got {preset!r}")
    applied = tree.get("applied_overlays") or []
    if not isinstance(applied, (list, tuple)) or not all(isinstance(a, str) for a in applied):
        raise ConfigError("applied_overlays: must be a list of strings")
    return RunSpec(
        run=settings,
        pipeline=pipeline,
        topology=sections["topology"],
        offline=sections["offline"],
        toggles=sections["toggles"],
        measure=sections["measure"],
        preset=preset,
        overlays=tuple(applied),
    )


def spec_to_tree(spec: RunSpec) -> dict[str, Any]:
    """Inverse of :func:`spec_from_tree` (``overlays`` are written as ``applied_overlays``).
    Returns a deep copy: editing the tree never changes ``spec``."""
    tree: dict[str, Any] = {"preset": spec.preset}
    if spec.overlays:
        tree["applied_overlays"] = list(spec.overlays)
    tree["run"] = {k: getattr(spec.run, k) for k in RUN_KEYS}
    tree["topology"] = dict(spec.topology)
    tree["offline"] = dict(spec.offline)
    tree["pipeline"] = spec_to_obj(spec.pipeline)
    tree["toggles"] = dict(spec.toggles)
    tree["measure"] = dict(spec.measure)
    return copy.deepcopy(tree)


def iter_block_specs(spec: PipelineSpec, path: str = "pipeline") -> Iterator[tuple[str, BlockSpec]]:
    """Depth-first ``(path, BlockSpec)`` of every leaf in execution order (Stack: stage, actuator, ...)."""
    if isinstance(spec, BlockSpec):
        yield path, spec
    elif isinstance(spec, (SeqSpec, ParallelSpec)):
        tag = "seq" if isinstance(spec, SeqSpec) else "parallel"
        for i, c in enumerate(spec.children):
            sub = f"{path}[{i}]" if (tag == "seq" and path == "pipeline") else f"{path}.{tag}[{i}]"
            yield from iter_block_specs(c, sub)
    elif isinstance(spec, StackSpec):
        for i, (s, a) in enumerate(zip(spec.stages, spec.actuators)):
            yield from iter_block_specs(s, f"{path}.stack.stages[{i}]")
            yield from iter_block_specs(a, f"{path}.stack.actuators[{i}]")
    else:
        raise TypeError(f"not a pipeline spec: {spec!r}")


def _map_leaves(spec: PipelineSpec, fn) -> PipelineSpec:
    if isinstance(spec, BlockSpec):
        return fn(spec)
    if isinstance(spec, SeqSpec):
        return SeqSpec(tuple(_map_leaves(c, fn) for c in spec.children), spec.name)
    if isinstance(spec, ParallelSpec):
        return ParallelSpec(tuple(_map_leaves(c, fn) for c in spec.children), spec.name)
    if isinstance(spec, StackSpec):
        return StackSpec(
            tuple(_map_leaves(c, fn) for c in spec.stages), tuple(_map_leaves(c, fn) for c in spec.actuators), spec.name
        )
    raise TypeError(f"not a pipeline spec: {spec!r}")


def apply_toggles(pipeline: SeqSpec, toggles: Mapping[str, bool]) -> SeqSpec:
    """Return ``pipeline`` with ``RunSpec.toggles`` overrides ``{"<block name>.<toggle>": bool}``
    written into the leaves' ``params["toggles"]`` (used by the builder and ``ablate``).

    Block names may contain dots; the toggle name is the part after the last dot. An override for an
    unknown block or for a toggle the block's spec does not list is a :class:`ConfigError`.
    """
    names = {b.name for _, b in iter_block_specs(pipeline)}
    per_block: dict[str, dict[str, bool]] = {}
    for key, value in toggles.items():
        block, _, toggle = key.rpartition(".")
        if block not in names:
            raise ConfigError(f"toggles.{key}: no block named {block!r} in the pipeline; blocks: {sorted(names)}")
        per_block.setdefault(block, {})[toggle] = bool(value)

    def override(b: BlockSpec) -> BlockSpec:
        if b.name not in per_block:
            return b
        given = dict(b.params.get("toggles") or {})
        for t in per_block[b.name]:
            if t not in given:
                raise ConfigError(
                    f"toggles.{b.name}.{t}: block {b.name!r} ({b.kind}) has no toggle {t!r} in params.toggles "
                    f"(known: {sorted(given)})"
                )
        given.update(per_block[b.name])
        return BlockSpec(b.kind, b.name, {**b.params, "toggles": given}, b.rate_div)

    return _map_leaves(pipeline, override)  # type: ignore[return-value]
