"""Locks as linear combinations of beam phases (conv §9). Owner: WP-F6.

e(t) = Σ_b c_b φ_b(t) − setpoint. No wrapping here (controllers see the continuous phase error;
metrics wrap to (−π, π] themselves).
"""

from __future__ import annotations

from collections.abc import Mapping

from cranberry.path.topology import Lock

__all__ = ["lock_error"]


def lock_error(lock: Lock, beam_phases_rad: Mapping[str, float]) -> float:
    """e = Σ_b combo[b]·φ_b − setpoint [rad] for one sample. Implemented by WP-F6."""
    raise NotImplementedError("path.locks.lock_error is implemented by WP-F6")
