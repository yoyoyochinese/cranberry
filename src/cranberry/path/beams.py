"""``beam_assembler`` block: beam phase at its observer with emission-time evaluation (PLAN §2.3.2,
conv §7, §9). Owner: WP-F6.

    φ_b(t) = φ_L(t − τ_b) + Σ_{(i,m) ∈ path b} ψ_{i,m,b}(t − τ_{b,after(i,m)})

τ_b: group delay source → observer; τ_{b,after(i,m)}: from the exit of subsegment (i, m) to the
observer. Counter-propagating beams (d = −1) read the same subsegments at other times. Histories are
read from the rings of ``laser.*`` and ``drv.*`` with fractional interpolation (``core.delay``). The
frequency dependence of the propagation phase is NOT added on top (no double counting, conv §7).
Jones matrices (Jᵀ for d = −1) multiply the input Jones vector (G04). Writes ``beam.<id>.*``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from cranberry.core.block import Block
from cranberry.core.registry import register_block

__all__ = ["make_beam_assembler"]


@register_block("beam_assembler")
def make_beam_assembler(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``beam_assembler`` block. Implemented by WP-F6."""
    raise NotImplementedError("block 'beam_assembler' is implemented by WP-F6")
