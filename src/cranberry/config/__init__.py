"""Configuration layer (PLAN §2.9; protected list "YAML config loading and its equivalence with
Python-list pipelines"). Owner of the bodies: WP-F2.

Public interface (PLAN §4.2):

    load(path_or_dict, overlays=()) -> RunSpec
    dump_yaml(spec) -> str
    to_spec(pipeline, *, run, topology=..., offline=..., toggles=..., measure=..., preset=None) -> RunSpec
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from cranberry.config.spec import RunSettings, RunSpec

if TYPE_CHECKING:
    from cranberry.core.block import Block

__all__ = ["RunSettings", "RunSpec", "load", "dump_yaml", "to_spec"]


def load(path_or_dict: str | Path | Mapping[str, Any], overlays: Sequence[str | Path] = ()) -> RunSpec:
    """Load a YAML file (or an already-parsed mapping), apply ``overlays`` in order (G12), resolve
    callables, and return the ``RunSpec``. Implemented by WP-F2."""
    raise NotImplementedError("config.load (YAML → RunSpec) is implemented by WP-F2")


def dump_yaml(spec: RunSpec) -> str:
    """Serialise a ``RunSpec`` to YAML text such that ``load(dump_yaml(s)) == s``. WP-F2."""
    raise NotImplementedError("config.dump_yaml is implemented by WP-F2")


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
    Callables must be registered (``@register_fn``) to round-trip (G12). WP-F2."""
    raise NotImplementedError("config.to_spec is implemented by WP-F2")
