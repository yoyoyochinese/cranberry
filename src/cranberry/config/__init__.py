"""Configuration layer (PLAN §2.9; protected list "YAML config loading and its equivalence with
Python-list pipelines"). Owner of the bodies: WP-F2.

Public interface (PLAN §4.2):

    load(path_or_dict, overlays=()) -> RunSpec
    dump_yaml(spec) -> str
    to_spec(pipeline, *, run, topology=..., offline=..., toggles=..., measure=..., preset=None) -> RunSpec

Helpers: ``read_input`` (the input document recorded as ``config.input.yaml``), ``validate`` (D34
checks of a single-run spec by instantiating every block), ``spec_hash``, ``sweep_markers``,
``ConfigError``; ``resolve_input`` / ``input_tree`` (the unresolved tree a spec carries, so sweeps
and f_sim rescaling re-evaluate dependent expressions, WP-F2c / D63-9); tags in ``config.yaml_io``; evaluator in ``config.expr``; merge, ``extends:`` and
the ``!delete`` / ``!replace`` directives in ``config.overlay``; references in ``config.resolve``;
dB → linear on load in ``config.db``.

Pipeline from YAML → ``load`` → RunSpec; pipeline from Python blocks → ``to_spec`` → RunSpec. Both
paths go through the same normalisation and :func:`~cranberry.config.resolve.resolve`, so equal
inputs give equal specs, and ``load(dump_yaml(spec)) == spec`` (VT-F1).
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from cranberry.config.overlay import (
    apply_overlays,
    check_no_directives,
    display_ref,
    expand_extends,
    find_base,
    find_overlay,
)
from cranberry.config.resolve import find_markers, normalize, path_str, resolve
from cranberry.config.spec import (
    ConfigError,
    RunSettings,
    RunSpec,
    apply_toggles,
    iter_block_specs,
    spec_from_tree,
    spec_to_tree,
)
from cranberry.config.yaml_io import Expr, canonical_json, dump_yaml_obj, load_yaml

if TYPE_CHECKING:
    from cranberry.core.block import Block

__all__ = [
    "ConfigError",
    "RunSettings",
    "RunSpec",
    "load",
    "dump_yaml",
    "to_spec",
    "read_input",
    "validate",
    "spec_hash",
    "sweep_markers",
    "resolve_input",
    "input_tree",
    "FrozenExpressionsWarning",
]


class FrozenExpressionsWarning(UserWarning):
    """A sweep / f_sim rescaling ran on resolved values only: ``!expr`` / ``!fn`` values that depend on
    a swept quantity keep the value they had when the spec was resolved (``config.input_tree``)."""


def _read(path_or_dict: str | Path | Mapping[str, Any]) -> tuple[dict[str, Any], Path | None]:
    if isinstance(path_or_dict, Mapping):
        return copy.deepcopy(dict(path_or_dict)), None
    path = Path(path_or_dict)
    try:
        tree = load_yaml(path)
    except FileNotFoundError:
        raise ConfigError(f"config file {str(path)!r} not found") from None
    except ConfigError:
        raise
    except Exception as e:  # yaml.YAMLError / ValueError carry position and reason
        raise ConfigError(f"{path}: {e}") from None
    if not isinstance(tree, Mapping):
        raise ConfigError(f"{path}: a config must be a YAML mapping, got {type(tree).__name__}")
    return dict(tree), path.resolve().parent


def _overlay_refs(tree: Mapping[str, Any], overlays: Sequence[str | Path]) -> list[str | Path]:
    own = tree.get("overlays") or []
    if not isinstance(own, (list, tuple)) or not all(isinstance(r, str) for r in own):
        raise ConfigError("overlays: must be a list of overlay file references")
    return [*own, *overlays]


def read_input(path_or_dict: str | Path | Mapping[str, Any], overlays: Sequence[str | Path] = ()) -> dict[str, Any]:
    """The input document of a run (``config.input.yaml``): the config as written, with ``overlays``
    extended by the extra ``overlays`` and every reference (``overlays``, ``extends``) made location
    independent, so that ``load(read_input(x, o)) == load(x, o)``. dB values stay as written."""
    tree, base_dir = _read(path_or_dict)
    refs = [display_ref(find_overlay(r, base_dir)) for r in _overlay_refs(tree, overlays)]
    if refs:
        tree["overlays"] = refs
    else:
        tree.pop("overlays", None)
    if "extends" in tree:
        tree["extends"] = display_ref(find_base(tree["extends"], base_dir))
    return tree


def load(path_or_dict: str | Path | Mapping[str, Any], overlays: Sequence[str | Path] = ()) -> RunSpec:
    """Load a YAML file (or an already-parsed mapping), apply ``extends`` and ``overlays`` (G12),
    resolve callables, convert dB parameters to linear (G16-i), and return the ``RunSpec``.

    Order (``config.overlay``): the ``extends:`` base (recursively, its own overlays applied), this
    file's keys merged on top (``!delete`` / ``!replace`` / ``null`` act on the base), then the
    config's own ``overlays:`` list, then ``overlays``. Already merged overlays (``applied_overlays:``,
    written by :func:`dump_yaml`) are only recorded. ``preset:`` is the display name (D34: no default
    file). Sweep markers stay in the spec (``runner.sweep`` expands them; dB names stay dB until then,
    ``config.db``); :func:`validate` checks a single-run spec.
    """
    tree, base_dir = _read(path_or_dict)
    origin = () if base_dir is None or isinstance(path_or_dict, Mapping) else (Path(path_or_dict).resolve(),)
    tree = expand_extends(tree, base_dir, _stack=origin)
    check_no_directives(tree, "a config without a base (extends:)")
    refs = _overlay_refs(tree, overlays)
    tree.pop("overlays", None)
    merged, applied = apply_overlays(tree, refs, base_dir)
    previous = merged.get("applied_overlays") or []
    if not isinstance(previous, list):
        raise ConfigError("applied_overlays: must be a list of strings")
    merged["applied_overlays"] = [*previous, *applied]
    return resolve_input(normalize(merged))


def dump_yaml(spec: RunSpec) -> str:
    """Serialise a ``RunSpec`` to YAML text such that ``load(dump_yaml(s)) == s``."""
    return dump_yaml_obj(spec_to_tree(spec))


def to_spec(
    pipeline: Block | Sequence[Block],
    *,
    run: RunSettings,
    topology: Mapping[str, Any] | None = None,
    offline: Mapping[str, Any] | None = None,
    toggles: Mapping[str, bool] | None = None,
    measure: Mapping[str, Any] | None = None,
    preset: str | None = None,
) -> RunSpec:
    """Turn a Python-built pipeline (a Block or a list of Blocks = Seq) into a ``RunSpec``.
    Callables must be registered (``@register_fn``) or built with ``config.expr.compile_expr`` to
    round-trip (G12); an unnamed lambda raises :class:`ConfigError` naming the block and key."""
    from cranberry.core.pipeline import Seq, SeqSpec, spec_to_obj

    if not isinstance(run, RunSettings):
        raise ConfigError(f"to_spec: run must be a RunSettings, got {type(run).__name__}")
    if isinstance(pipeline, Sequence):
        pipe_spec = SeqSpec(tuple(b.spec for b in pipeline))
    elif isinstance(pipeline, Seq) and pipeline._name is None:
        pipe_spec = pipeline.spec
    else:
        pipe_spec = SeqSpec((pipeline.spec,))
    pipe_obj = spec_to_obj(pipe_spec)
    for path, b in iter_block_specs(pipe_spec):  # name the block in callable errors
        normalize(b.params, f"{path} (block {b.name!r}).params")
    tree = {
        "preset": preset,
        "run": {"duration_s": run.duration_s, "f_sim_hz": run.f_sim_hz, "chunk_samples": run.chunk_samples,
                "seed": run.seed, **({"seed_mode": run.seed_mode} if run.seed_mode is not None else {})},
        "topology": dict(topology or {}),
        "offline": dict(offline or {}),
        "pipeline": pipe_obj,
        "toggles": dict(toggles or {}),
        "measure": dict(measure or {}),
    }
    return resolve_input(normalize(tree))


def resolve_input(tree: Mapping[str, Any]) -> RunSpec:
    """Resolve an unresolved config tree (``extends`` / overlays merged, expressions symbolic) into a
    ``RunSpec`` that carries the tree as ``RunSpec.unresolved`` (WP-F2c)."""
    tree = copy.deepcopy(dict(tree))
    return dataclasses.replace(spec_from_tree(resolve(tree)), unresolved=tree)


def input_tree(spec: RunSpec) -> dict[str, Any] | None:
    """The tree to re-resolve ``spec`` from after changing values in it (WP-F2c, D63-9), or None.

    Returns a deep copy of ``spec.unresolved`` with its ``run:`` section replaced by ``spec.run``, so
    edits of the run settings (``dataclasses.replace(spec, run=...)``: CLI ``--seed``, a seed mode,
    f_sim) are kept and expressions using ``run.*`` see them. Every expression / ``!fn`` / dB value of
    the tree is evaluated again when the tree is resolved.

    None (expressions are frozen at their resolved values; callers warn with
    :class:`FrozenExpressionsWarning`) when

    * the spec carries no tree (built directly with ``RunSpec(...)`` / ``spec_from_tree``), or
    * the tree no longer matches the spec outside ``run:`` (a section was edited with
      ``dataclasses.replace``; the edit would otherwise be lost).

    A spec loaded from an already *resolved* document (``dump_yaml`` output, a record's
    ``config.resolved.yaml``) carries that document: it holds values, not expressions, so dependent
    values stay frozen without a warning. Sweep from the input document (``config.read_input`` /
    a record's ``config.input.yaml``) to re-evaluate them.
    """
    if spec.unresolved is None:
        return None
    if dataclasses.replace(spec_from_tree(resolve(spec.unresolved)), run=spec.run) != spec:
        return None
    tree = copy.deepcopy(dict(spec.unresolved))
    tree["run"] = spec_to_tree(spec)["run"]
    return tree


def spec_hash(spec: RunSpec) -> str:
    """SHA-256 (hex) of the canonical JSON of the spec tree (independent of key order)."""
    return hashlib.sha256(canonical_json(spec_to_tree(spec)).encode("utf-8")).hexdigest()


def sweep_markers(spec: RunSpec) -> dict[str, Any]:
    """{path: marker} of everything still unresolved in ``spec``: sweep axes and expressions that
    depend on them (empty for a single-run spec)."""
    tree = spec_to_tree(spec)
    out = {path_str(tree, p): m for p, m in find_markers(tree)}
    out.update({path_str(tree, p): m for p, m in find_markers(tree, (Expr,)) if not m.args})
    return out


def validate(spec: RunSpec) -> None:
    """D34 checks of a single-run spec; raises :class:`ConfigError` naming block and key.

    * no sweep markers / deferred expressions left (use ``runner.sweep``);
    * block names unique (they key RNG streams and runtime parameters, G11);
    * ``toggles`` overrides name existing blocks and toggles;
    * every block instantiates from its spec: unknown kinds, missing modelling choices and toggles
      (``MissingParameterError``) and invalid values are reported with the pipeline path.
    """
    from cranberry.core.block import MissingParameterError
    from cranberry.core.registry import ensure_loaded, make_block

    left = sweep_markers(spec)
    if left:
        raise ConfigError(f"unresolved sweep axes / expressions {sorted(left)}: expand them with runner.sweep")
    seen: dict[str, str] = {}
    for path, b in iter_block_specs(spec.pipeline):
        if b.name in seen:
            raise ConfigError(f"duplicate block name {b.name!r} at {seen[b.name]} and {path}")
        seen[b.name] = path
    pipeline = apply_toggles(spec.pipeline, spec.toggles)
    ensure_loaded()
    for path, b in iter_block_specs(pipeline):
        try:
            make_block(b)
        except MissingParameterError as e:
            raise ConfigError(f"{path}: block {b.name!r} ({b.kind}): {e.args[0] if e.args else e}") from None
        except KeyError as e:
            raise ConfigError(f"{path}: block {b.name!r}: {e.args[0] if e.args else e}") from None
        except (ValueError, TypeError) as e:
            raise ConfigError(f"{path}: block {b.name!r} ({b.kind}): {e}") from None
