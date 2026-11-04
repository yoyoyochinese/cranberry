"""Preset + overlay deep merge (G12), ``extends:`` and merge directives. Owner: WP-F2 (F2b).

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
* Overlay files are plain partial trees; they may not themselves list ``overlays`` or ``extends``.

Overlay conflicts (Phase 6, reviews/phase6_wp6g_ab.md §5 G-3..G-6; :func:`apply_overlays`). The overlays applied in
one ``load`` (a config's own ``overlays:`` then the extra ones; a base reached through ``extends:`` counts as the
base, its overlays included) are checked against each other, and a conflict is a :class:`ConfigError` naming both
overlays:

* **discard**: a later overlay discards wholesale what an earlier overlay wrote — a ``!replace`` (mapping, list or
  keyed-list element), a ``!delete`` / ``null``, a plain list replacing a list, or a value of another kind (scalar ↔
  mapping/list) — at or above a path the earlier overlay set. Example: ``controller/x`` sets ``adc0.params.device``,
  then a gate overlay ``!replace``-s ``adc0`` (the device would be silently lost).
* **structure**: a later overlay merges block keys (``block``, ``params``, ``rate_div``) into a pipeline node that
  is a ``seq`` / ``stack`` (or ``seq`` / ``stack`` into a leaf block), typically because an earlier overlay
  ``!replace``-d it: e.g. ``probe/detect_heterodyne`` makes ``liu0`` a seq, a gate overlay then sets ``liu0.params``.
  Never allowed (the merged node would not be a valid pipeline element).

A later overlay writing a *scalar* over an earlier overlay's scalar (a point of an axis: ``gate/borealis_si`` after
``gate/borealis_window``) is not a conflict, nor is any overlay changing the base config. An intentional discard is
declared in the later overlay with the top-level key ``overrides: [<overlay ref>, ...]`` (refs as in ``overlays:``,
e.g. ``gate/borealis_window``; ``"*"`` = any earlier overlay): its discards of what those overlays wrote are
accepted. ``overrides`` is consumed while loading (never part of the merged tree).

Directives (WP-F2b; tags in ``config.yaml_io``), valid wherever a partial tree is merged onto a
base (overlays, and the own keys of a file that ``extends:`` a base):

* ``key: !delete`` removes ``key`` from the base mapping; unlike ``null`` it is an error when the
  base has no such key (catches typos).
* ``- !delete <label>`` / ``- !delete {id|name: <label>}`` inside a keyed list removes that element
  (error if absent, and in a list that is not keyed).
* ``key: !replace {...}`` / ``key: !replace [...]`` replaces the base value wholesale instead of
  merging (error if the base has no ``key``); ``- !replace {name: x, ...}`` in a keyed list replaces
  element ``x`` wholesale, keeping its position (error if absent).
* A directive with nothing to act on (inside a newly added subtree, or left in a config that has no
  base) is an error (:func:`check_no_directives`).

``extends:`` (WP-F2b; :func:`find_base`, :func:`expand_extends`)::

    extends: refD                         # a preset name: the unique presets/**/refD.yaml
    extends: ../reference/refD.yaml       # a path: relative to this file, then the repository root

Resolution order of one config: (1) the base is loaded recursively the same way (its own
``extends`` chain first, then its own keys, then its own ``overlays``); (2) this file's keys are
merged on top (:func:`merge`: ``null`` and the directives work); (3) this file's ``overlays`` and the
extra overlays of ``config.load`` are applied. ``preset:`` stays the display name; a file that does
not set it inherits the base's. A chain that returns to a file already being loaded is an error
naming the cycle. ``extends`` is consumed while loading: a resolved spec (and ``dump_yaml``) never
carries it, so reloading a dump never merges a base again. ``config.read_input`` keeps ``extends``
but rewrites it to the base's repository-relative path (location independent, pinned to one file).
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from cranberry.config.spec import ConfigError
from cranberry.config.yaml_io import OVERLAY_DIRECTIVES, Delete, Expr, Fn, Replace, load_yaml

__all__ = [
    "merge",
    "REPO_ROOT",
    "PRESETS_DIR",
    "find_overlay",
    "display_ref",
    "apply_overlays",
    "find_base",
    "expand_extends",
    "check_no_directives",
    "footprint",
]

#: Repository root (``src/cranberry/config/overlay.py`` → three levels up from ``src``).
REPO_ROOT = Path(__file__).resolve().parents[3]
PRESETS_DIR = REPO_ROOT / "presets"

_LIST_KEYS = ("id", "name")


# ---------------------------------------------------------------------------------------------- merge


def _keyed(x: Any) -> Any:
    """The mapping a keyed-list key is read from (directives unwrapped)."""
    if isinstance(x, Replace):
        return x.value
    if isinstance(x, Delete) and x.field is not None:
        return {x.field: x.label}
    return x


def _list_key(base: list[Any], overlay: list[Any]) -> str | None:
    # a bare ``!delete <label>`` matches whichever key the rest of the list is keyed by
    items = [_keyed(x) for x in base + overlay if not (isinstance(x, Delete) and x.field is None)]
    if not items or not all(isinstance(x, Mapping) for x in items):
        return None
    for key in _LIST_KEYS:
        if all(key in x for x in items):
            return key
    return None


def _find(items: list[Any], key: str, label: Any) -> int | None:
    for i, x in enumerate(items):
        if x.get(key) == label:
            return i
    return None


def _merge_keyed(base: list[Any], overlay: list[Any], key: str) -> list[Any]:
    out = [copy.deepcopy(x) for x in base]
    labels = [x[key] for x in out]
    if len({repr(x) for x in labels}) != len(labels):
        raise ConfigError(f"duplicate {key!r} values in a keyed list: {labels}")
    for item in overlay:
        if isinstance(item, Delete):
            if item.label is None:
                raise ConfigError(f"a bare !delete in a keyed list must name the element (!delete <{key}>)")
            i = _find(out, key, item.label)
            if i is None:
                raise ConfigError(
                    f"!delete {item.label!r}: no element with {key} {item.label!r}; present: {[x[key] for x in out]}"
                )
            del out[i]
        elif isinstance(item, Replace):
            label = item.value[key]
            i = _find(out, key, label)
            if i is None:
                raise ConfigError(
                    f"!replace {key}={label!r}: no such element to replace; present: {[x[key] for x in out]}"
                )
            out[i] = _strip_nulls(item.value)
        else:
            i = _find(out, key, item[key])
            if i is None:
                out.append(_strip_nulls(item))
            else:
                out[i] = merge(out[i], item)
    return out


def _merge_value(base: Any, overlay: Any) -> Any:
    if isinstance(overlay, Replace):
        return _strip_nulls(overlay.value)
    if isinstance(base, Mapping) and isinstance(overlay, Mapping):
        return merge(base, overlay)
    if isinstance(base, list) and isinstance(overlay, list):
        key = _list_key(base, overlay)
        if key is not None:
            return _merge_keyed(base, overlay, key)
        if any(isinstance(x, OVERLAY_DIRECTIVES) for x in overlay):
            raise ConfigError("!delete / !replace elements need a keyed list (every element with 'id' or 'name')")
        return copy.deepcopy(overlay)
    return _strip_nulls(overlay)


def _strip_nulls(value: Any) -> Any:
    """Copy of an overlay value added where the base has nothing (``null`` = delete → just omit).
    Directives inside it have nothing to act on → error."""
    if isinstance(value, Replace):
        value = value.value
    check_no_directives(value, "a value added where the base has nothing")
    if isinstance(value, Mapping):
        return {k: _strip_nulls(v) for k, v in value.items() if v is not None}
    return copy.deepcopy(value)


def check_no_directives(tree: Any, where: str = "the config") -> None:
    """Raise :class:`ConfigError` if a ``!delete`` / ``!replace`` is left in ``tree`` (no base to act on)."""
    stack: list[tuple[str, Any]] = [("", tree)]
    while stack:
        path, x = stack.pop()
        if isinstance(x, OVERLAY_DIRECTIVES):
            tag = "!delete" if isinstance(x, Delete) else "!replace"
            raise ConfigError(
                f"{tag} at {path or '<root>'} in {where}: nothing to act on (directives apply only in overlays "
                "or in a file with extends:, onto an existing key / element)"
            )
        if isinstance(x, Mapping):
            stack.extend((f"{path}.{k}" if path else str(k), v) for k, v in x.items())
        elif isinstance(x, (list, tuple)):
            stack.extend((f"{path}[{i}]", v) for i, v in enumerate(x))
        elif isinstance(x, Expr):
            stack.extend((f"{path}.where.{k}", v) for k, v in x.where.items())
        elif isinstance(x, Fn):
            stack.extend((f"{path}.args.{k}", v) for k, v in x.args.items())


def merge(base: Mapping[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Deep-merge ``overlay`` onto ``base`` (rules above) and return a new tree."""
    if not isinstance(base, Mapping) or not isinstance(overlay, Mapping):
        raise ConfigError("merge needs two mappings")
    out = {k: copy.deepcopy(v) for k, v in base.items()}
    for key, value in overlay.items():
        if isinstance(value, Delete):
            if value.label is not None:
                raise ConfigError(f"{key}: !delete <label> is for keyed-list elements; write '{key}: !delete'")
            if key not in out:
                raise ConfigError(f"{key}: !delete of a key the base does not have; present: {sorted(out)}")
            del out[key]
        elif isinstance(value, Replace) and key not in out:
            raise ConfigError(
                f"{key}: !replace of a key the base does not have (write the value without !replace to add it); "
                f"present: {sorted(out)}"
            )
        elif value is None:
            out.pop(key, None)
        elif key in out:
            try:
                out[key] = _merge_value(out[key], value)
            except ConfigError as e:
                raise ConfigError(f"{key}: {e}") from None
        else:
            out[key] = _strip_nulls(value)
    return out


# ------------------------------------------------------------------------------------------- overlays


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


_BLOCK_KEYS = ("block", "params", "rate_div")
_STRUCT_KEYS = ("seq", "stack")


def _fmt(path: tuple[Any, ...]) -> str:
    return ".".join(str(p[1]) if isinstance(p, tuple) else str(p) for p in path) or "<root>"


def footprint(base: Any, ov: Mapping[str, Any]) -> tuple[list[tuple[Any, ...]], list[tuple[Any, ...]]]:
    """(writes, discards) of merging ``ov`` onto ``base``: the leaf paths it sets and the paths whose previous content
    it discards wholesale (module docstring, "Overlay conflicts"). Keyed-list elements appear as (key, label)."""
    writes: list[tuple[Any, ...]] = []
    discards: list[tuple[Any, ...]] = []

    def walk(b: Any, o: Mapping[str, Any], p: tuple[Any, ...]) -> None:
        for k, v in o.items():
            q = p + (k,)
            bk = b.get(k) if isinstance(b, Mapping) else None
            if isinstance(v, (Replace, Delete)) or v is None:
                discards.append(q)
                if isinstance(v, Replace):
                    writes.append(q)
            elif isinstance(v, Mapping) and isinstance(bk, Mapping):
                walk(bk, v, q)
            elif isinstance(v, list) and isinstance(bk, list) and (key := _list_key(bk, v)) is not None:
                for item in v:
                    if isinstance(item, Delete):
                        discards.append(q + ((key, item.label),))
                    elif isinstance(item, Replace):
                        lab = q + ((key, item.value[key]),)
                        discards.append(lab)
                        writes.append(lab)
                    else:
                        i = _find(bk, key, item[key])
                        if i is None:
                            writes.append(q + ((key, item[key]),))
                        else:
                            walk(bk[i], item, q + ((key, item[key]),))
            else:
                if bk is not None and (isinstance(v, (list, Mapping)) or isinstance(bk, (list, Mapping))):
                    discards.append(q)  # a list replacing a list, or a change of kind
                writes.append(q)

    walk(base, ov, ())
    return writes, discards


def _mixed_nodes(tree: Any, path: tuple[Any, ...] = ()) -> dict[tuple[Any, ...], list[str]]:
    """{path: block keys} of every mapping holding both seq/stack and block keys (an invalid pipeline element)."""
    out: dict[tuple[Any, ...], list[str]] = {}
    if isinstance(tree, Mapping):
        if any(k in tree for k in _STRUCT_KEYS) and any(k in tree for k in _BLOCK_KEYS):
            out[path] = [k for k in _BLOCK_KEYS if k in tree]
        for k, v in tree.items():
            out.update(_mixed_nodes(v, path + (k,)))
    elif isinstance(tree, list):
        for i, v in enumerate(tree):
            lab = _keyed(v).get("name", _keyed(v).get("id", i)) if isinstance(_keyed(v), Mapping) else i
            out.update(_mixed_nodes(v, path + (("name", lab),)))
    return out


def apply_overlays(
    tree: Mapping[str, Any], refs: Sequence[str | Path], base_dir: Path | None = None
) -> tuple[dict[str, Any], list[str]]:
    """Merge the overlay files ``refs`` onto ``tree`` in order; returns (merged tree, display refs). Conflicts
    between the overlays raise :class:`ConfigError` (module docstring, "Overlay conflicts")."""
    out = dict(tree)
    applied: list[str] = []
    written: dict[tuple[Any, ...], str] = {}   # leaf path → the overlay that last wrote it
    replaced: dict[tuple[Any, ...], str] = {}  # path → the overlay that last !replace-d it
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
        if "overlays" in ov or "applied_overlays" in ov or "extends" in ov:
            raise ConfigError(f"overlay {display_ref(path)}: overlay files cannot list overlays or extends")
        me = display_ref(path)
        ov = dict(ov)
        allowed = ov.pop("overrides", None) or []
        if not isinstance(allowed, list) or not all(isinstance(r, str) for r in allowed):
            raise ConfigError(f"overlay {me}: overrides: must be a list of overlay references")
        allow_all = "*" in allowed
        allowed_refs = {display_ref(find_overlay(r, path.parent)) for r in allowed if r != "*"}
        writes, discards = footprint(out, ov)
        for d in discards:
            for q, who in written.items():
                if who != me and q[: len(d)] == d and not allow_all and who not in allowed_refs:
                    raise ConfigError(
                        f"overlay conflict: {me} discards {_fmt(d)} wholesale (!replace / !delete / list or kind "
                        f"change), but the earlier overlay {who} set {_fmt(q)} there; its value would be silently lost. "
                        f"Apply them in the other order if that composes, write a combined overlay, or declare "
                        f"'overrides: [{who}]' in {me} if the discard is intended")
        before = _mixed_nodes(out)
        try:
            out = merge(out, ov)
        except ConfigError as e:
            raise ConfigError(f"overlay {me}: {e}") from None
        for node, keys in _mixed_nodes(out).items():
            if node not in before:
                who = next((w for q, w in replaced.items() if q == node), None)
                cause = f"the earlier overlay {who} replaced it (!replace)" if who else "it is so in the base config"
                raise ConfigError(
                    f"overlay conflict: {me} merges block keys {keys} into {_fmt(node)}, which is a seq / stack "
                    f"({cause}); the overlays cannot be combined in this order (write a combined overlay)")
        for d in discards:
            for q in [q for q in written if q[: len(d)] == d]:
                del written[q]
        for w in writes:
            written[w] = me
        for d in discards:
            if d in writes:
                replaced[d] = me
        applied.append(me)
    return out, applied


# -------------------------------------------------------------------------------------------- extends


def _is_path_ref(ref: str) -> bool:
    return "/" in ref or "\\" in ref or ref.endswith((".yaml", ".yml"))


def find_base(ref: str, base_dir: Path | None = None) -> Path:
    """Locate the base file of ``extends: <ref>``.

    * A path (contains ``/`` or ends in ``.yaml``/``.yml``): absolute; else relative to the including
      file's directory, then to the repository root.
    * A preset name: the unique ``presets/**/<ref>.yaml`` (or ``.yml``); none or several → error.
    """
    if not isinstance(ref, str) or not ref.strip():
        raise ConfigError(f"extends: must be a preset name or a relative path, got {ref!r}")
    if _is_path_ref(ref):
        p = Path(ref)
        if p.is_absolute():
            candidates = [p]
        else:
            candidates = [*([Path(base_dir) / p] if base_dir is not None else []), REPO_ROOT / p]
        for c in candidates:
            if c.is_file():
                return c.resolve()
        raise ConfigError(f"extends: {ref!r} not found; tried: {[str(c) for c in candidates]}")
    hits = sorted({*PRESETS_DIR.glob(f"**/{ref}.yaml"), *PRESETS_DIR.glob(f"**/{ref}.yml")})
    if not hits:
        raise ConfigError(f"extends: no preset named {ref!r} (looked for presets/**/{ref}.yaml)")
    if len(hits) > 1:
        raise ConfigError(
            f"extends: preset name {ref!r} is ambiguous: {[display_ref(h) for h in hits]}; write the relative path"
        )
    return hits[0].resolve()


def _load_base(path: Path, stack: tuple[Path, ...]) -> dict[str, Any]:
    """Fully composed tree of a base file: its extends chain, its own keys, then its overlays."""
    try:
        tree = load_yaml(path)
    except Exception as e:  # yaml errors carry the position
        raise ConfigError(f"extends: {display_ref(path)}: {e}") from None
    if not isinstance(tree, Mapping):
        raise ConfigError(f"extends: {display_ref(path)}: a config must be a YAML mapping")
    tree = expand_extends(tree, path.parent, _stack=stack)
    refs = tree.pop("overlays", None) or []
    if not isinstance(refs, list) or not all(isinstance(r, str) for r in refs):
        raise ConfigError(f"extends: {display_ref(path)}: overlays: must be a list of overlay file references")
    merged, applied = apply_overlays(tree, refs, path.parent)
    previous = merged.get("applied_overlays") or []
    if applied or previous:
        merged["applied_overlays"] = [*previous, *applied]
    return merged


def expand_extends(
    tree: Mapping[str, Any], base_dir: Path | None = None, *, _stack: tuple[Path, ...] = ()
) -> dict[str, Any]:
    """Resolve ``extends:`` of ``tree`` (module docstring): the base (fully composed, its overlays
    applied and recorded in ``applied_overlays``) with ``tree``'s own keys merged on top. ``tree``'s
    own ``overlays`` are kept for the caller to apply. Without ``extends``, ``tree`` is returned."""
    tree = dict(tree)
    if "extends" not in tree:
        return tree
    path = find_base(tree.pop("extends"), base_dir)
    if path in _stack:
        cycle = " -> ".join(display_ref(p) for p in (*_stack[_stack.index(path):], path))
        raise ConfigError(f"extends: circular chain {cycle}")
    base = _load_base(path, (*_stack, path))
    own_overlays = tree.pop("overlays", None)
    try:
        out = merge(base, tree)
    except ConfigError as e:
        raise ConfigError(f"extends {display_ref(path)}: {e}") from None
    if own_overlays is not None:
        out["overlays"] = own_overlays
    return out
