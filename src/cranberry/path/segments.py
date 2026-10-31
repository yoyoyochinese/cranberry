"""``segment_phase`` block: drivers × coefficients → ψ (segment × beam, per subsegment) and τ_g
(PLAN §2.3.2, §4.3 PATH, conv §9, T08, T10). Owner: WP-F6.

Reads ``drv.<seg>.<driver>`` (per subsegment, with history rings), writes ``psi.<seg>.<beam>`` [rad]
and ``tau_g.<seg>`` [s]. Spatial structure (lumped / distributed / mixed) and zones come from the
Topology. TEC writes temperature, not phase (G08); "write phase directly" is reproduced by a
coefficient preset (equal axis coefficients, zero group-delay effect), not by a separate code path.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from cranberry.core.block import Block
from cranberry.core.registry import register_block

__all__ = ["make_segment_phase"]


@register_block("segment_phase")
def make_segment_phase(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``segment_phase`` block. Implemented by WP-F6."""
    raise NotImplementedError("block 'segment_phase' is implemented by WP-F6")
