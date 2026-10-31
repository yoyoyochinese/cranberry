"""Device presets ``devices/<name>.yaml`` (protected list "device presets"; prompt §5.3, PLAN §2.9–§2.10).

File format (every value carries its source, PLAN §7.3)::

    meta: {device: <name>, title: <str>, document: <str>}     # free-form metadata, no sources
    <section>:                                               # adc, dac, pid, ramp, loop, ...
      <key>: {value: <scalar | list | null>, source: "<manual + page | Dnn | Tnn | UNVERIFIED ...>",
              sweep: [lo, hi]?, note: <str>?}

* Values are SI and linear (G16-i); key names carry unit suffixes (``_hz``, ``_s``, ``_v``, ``_ohm``).
* ``value: null`` means "not a device fact / not documented"; consumers that need it must get it from
  the preset (D34) and raise :class:`~cranberry.core.block.MissingParameterError` otherwise.
* ``sweep`` is the documented uncertainty range of a value (e.g. DSC1 latency 5.7–8.3 µs, D52).

:func:`load_device` returns the full tree (with sources) for the run record; :func:`device_values`
strips it to ``{section: {key: value}}``, the form taken by ``hardware.mapping`` and by blocks
configured with ``device: <name>`` (PLAN §2.9).
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml

__all__ = ["DEVICES_DIR", "LEAF_KEYS", "DeviceFormatError", "available_devices", "load_device", "device_values", "sweep_of"]

#: Repository ``devices/`` directory (src/cranberry/hardware/devices.py → repo root).
DEVICES_DIR: Path = Path(__file__).resolve().parents[3] / "devices"

#: Allowed keys of a value leaf.
LEAF_KEYS: frozenset[str] = frozenset({"value", "source", "sweep", "note"})


class DeviceFormatError(ValueError):
    """A device preset violates the file format (e.g. a value without ``source``)."""


def available_devices(root: Path | None = None) -> tuple[str, ...]:
    """Names of the presets in ``root`` (default :data:`DEVICES_DIR`)."""
    return tuple(sorted(p.stem for p in (root or DEVICES_DIR).glob("*.yaml")))


def _is_number(s: str) -> bool:
    try:
        float(s)
    except ValueError:
        return False
    return True


def _check(node: Any, path: str) -> None:
    if not isinstance(node, Mapping):
        raise DeviceFormatError(f"{path}: expected a section or a {{value, source}} leaf, got {node!r}")
    if "value" in node:
        extra = set(node) - LEAF_KEYS
        if extra:
            raise DeviceFormatError(f"{path}: unknown leaf keys {sorted(extra)}")
        src = node.get("source")
        if not isinstance(src, str) or not src.strip():
            raise DeviceFormatError(f"{path}: value without a non-empty 'source' (PLAN §7.3)")
        sweep = node.get("sweep")
        if sweep is not None and (not isinstance(sweep, list) or len(sweep) != 2):
            raise DeviceFormatError(f"{path}: 'sweep' must be [lo, hi], got {sweep!r}")
        value = node["value"]
        for v in (value if isinstance(value, list) else [value]) + list(sweep or []):
            if isinstance(v, str) and _is_number(v):
                raise DeviceFormatError(f"{path}: number {v!r} parsed as a string (YAML 1.1: write 1.0e+6, not 1.0e6)")
        return
    if not node:
        raise DeviceFormatError(f"{path}: empty section")
    for key, child in node.items():
        _check(child, f"{path}.{key}")


def load_device(name: str, root: Path | None = None) -> dict[str, Any]:
    """Load and validate ``<root>/<name>.yaml``; returns the full tree including sources."""
    path = (root or DEVICES_DIR) / f"{name}.yaml"
    if not path.is_file():
        raise FileNotFoundError(f"no device preset {name!r} at {path}; available {available_devices(root)}")
    with path.open(encoding="utf-8") as fh:
        tree = yaml.safe_load(fh)
    if not isinstance(tree, dict) or not isinstance(tree.get("meta"), dict):
        raise DeviceFormatError(f"{path}: top level must be a mapping with a 'meta' section")
    if tree["meta"].get("device") != name:
        raise DeviceFormatError(f"{path}: meta.device {tree['meta'].get('device')!r} != file name {name!r}")
    for key, section in tree.items():
        if key != "meta":
            _check(section, key)
    return tree


def device_values(tree: Mapping[str, Any]) -> dict[str, Any]:
    """Strip a loaded tree to plain values: ``{section: {key: value}}`` (``meta`` kept as is)."""

    def strip(node: Any) -> Any:
        if isinstance(node, Mapping) and "value" in node:
            return node["value"]
        return {k: strip(v) for k, v in node.items()}

    return {k: (dict(v) if k == "meta" else strip(v)) for k, v in tree.items()}


def sweep_of(tree: Mapping[str, Any], dotted: str) -> list[Any] | None:
    """The documented ``sweep`` range of ``"<section>.<key>"`` in a loaded tree (None if absent)."""
    node: Any = tree
    for part in dotted.split("."):
        node = node[part]
    return node.get("sweep")
