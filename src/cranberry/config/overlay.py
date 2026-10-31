"""Preset + overlay deep merge (G12). Owner: WP-F2.

Rules: start from ``preset``, apply ``overlays`` in order; mappings merge recursively; lists are
replaced, except lists of mappings carrying ``name`` or ``id`` (segments, blocks), which merge by that
key; a ``null`` value deletes the key. The merged result is saved as ``config.resolved.yaml``.
Overlays live in ``presets/axes/<axis>/<value>.yaml``; no preset file is a "default" (D34).

Details:

* A list merges by key only if every element of both lists is a mapping with the key (``id`` is
  tried first, then ``name``); matched elements merge recursively (so ``null`` inside an element
  deletes that element's key), unmatched overlay elements are appended in overlay order. Otherwise
  the overlay list replaces the base list.
* Tag markers (``!sweep``, ``!expr``, ...) are leaf values: they replace, never merge.
* Overlay references (``overlays:`` in a config, ``--overlay`` on the CLI) are file paths, tried in
  order: absolute; relative to the including file; relative to the repository root; ``presets/axes/
  <ref>[.yaml]`` (so ``fiber/smf`` names ``presets/axes/fiber/smf.yaml``); relative to the CWD.
* Overlay files are plain partial trees; they may not themselves list ``overlays``.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from cranberry.config.spec import ConfigError
from cranberry.config.yaml_io import load_yaml

__all__ = ["merge", "REPO_ROOT", "PRESETS_DIR", "find_overlay", "display_ref", "apply_overlays"]

#: Repository root (``src/cranberry/config/overlay.py`` → three levels up from ``src``).
REPO_ROOT = Path(__file__).resolve().parents[3]
PRESETS_DIR = REPO_ROOT / "presets"

_LIST_KEYS = ("id", "name")


def _list_key(base: list[Any], overlay: list[Any]) -> str | None:
    items = base + overlay
    if not items or not all(isinstance(x, Mapping) for x in items):
        return None
    for key in _LIST_KEYS:
        if all(key in x for x in items):
            return key
    return None


def _merge_value(base: Any, overlay: Any) -> Any:
    if isinstance(base, Mapping) and isinstance(overlay, Mapping):
        return merge(base, overlay)
    if isinstance(base, list) and isinstance(overlay, list):
        key = _list_key(base, overlay)
        if key is None:
            return copy.deepcopy(overlay)
        out = [copy.deepcopy(x) for x in base]
        index = {x[key]: i for i, x in enumerate(out)}
        if len(index) != len(out):
            raise ConfigError(f"duplicate {key!r} values in a keyed list: {[x[key] for x in out]}")
        for item in overlay:
            if item[key] in index:
                out[index[item[key]]] = merge(out[index[item[key]]], item)
            else:
                index[item[key]] = len(out)
                out.append(_strip_nulls(item))
        return out
    return _strip_nulls(overlay)


def _strip_nulls(value: Any) -> Any:
    """Copy of an overlay value added where the base has nothing (``null`` = delete → just omit)."""
    if isinstance(value, Mapping):
        return {k: _strip_nulls(v) for k, v in value.items() if v is not None}
    return copy.deepcopy(value)


def merge(base: Mapping[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Deep-merge ``overlay`` onto ``base`` (rules above) and return a new tree."""
    if not isinstance(base, Mapping) or not isinstance(overlay, Mapping):
        raise ConfigError("merge needs two mappings")
    out = {k: copy.deepcopy(v) for k, v in base.items()}
    for key, value in overlay.items():
        if value is None:
            out.pop(key, None)
        elif key in out:
            out[key] = _merge_value(out[key], value)
        else:
            out[key] = _strip_nulls(value)
    return out


def find_overlay(ref: str | Path, base_dir: Path | None = None) -> Path:
    """Locate overlay file ``ref`` (search order in the module docstring)."""
    p = Path(ref)
    candidates: list[Path] = []
    if p.is_absolute():
        candidates.append(p)
    else:
        if base_dir is not None:
            candidates.append(Path(base_dir) / p)
        candidates.append(REPO_ROOT / p)
        axes = PRESETS_DIR / "axes" / p
        candidates.append(axes if axes.suffix in (".yaml", ".yml") else axes.with_name(axes.name + ".yaml"))
        candidates.append(Path.cwd() / p)
    for c in candidates:
        if c.is_file():
            return c.resolve()
    raise ConfigError(f"overlay {str(ref)!r} not found; tried: {[str(c) for c in candidates]}")


def display_ref(path: Path) -> str:
    """Path relative to the repository root when inside it (machine independent), else absolute."""
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def apply_overlays(
    tree: Mapping[str, Any], refs: Sequence[str | Path], base_dir: Path | None = None
) -> tuple[dict[str, Any], list[str]]:
    """Merge the overlay files ``refs`` onto ``tree`` in order; returns (merged tree, display refs)."""
    out = dict(tree)
    applied: list[str] = []
    for ref in refs:
        path = find_overlay(ref, base_dir)
        try:
            ov = load_yaml(path)
        except Exception as e:  # yaml errors carry the position
            raise ConfigError(f"overlay {display_ref(path)}: {e}") from None
        if ov is None:
            ov = {}
        if not isinstance(ov, Mapping):
            raise ConfigError(f"overlay {display_ref(path)}: must be a mapping, got {type(ov).__name__}")
        if "overlays" in ov or "applied_overlays" in ov:
            raise ConfigError(f"overlay {display_ref(path)}: overlay files cannot list overlays")
        out = merge(out, ov)
        applied.append(display_ref(path))
    return out, applied
