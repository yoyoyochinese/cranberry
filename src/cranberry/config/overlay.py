"""Preset + overlay deep merge (G12). Owner: WP-F2.

Rules: start from ``preset``, apply ``overlays`` in order; mappings merge recursively; lists are
replaced, except lists of mappings carrying ``name`` or ``id`` (segments, blocks), which merge by that
key; a ``null`` value deletes the key. The merged result is saved as ``config.resolved.yaml``.
Overlays live in ``presets/axes/<axis>/<value>.yaml``; no preset file is a "default" (D34).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

__all__ = ["merge"]


def merge(base: Mapping[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Deep-merge ``overlay`` onto ``base`` (rules above) and return a new tree. WP-F2."""
    raise NotImplementedError("config.overlay.merge is implemented by WP-F2")
