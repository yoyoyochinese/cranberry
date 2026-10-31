"""Offline stage (prompt §5.1, PLAN §2.2): tables + cache, SourceOutput (G07), LossRegistry (G06),
BandLedger (G05), Gaussian covariance operations. Owner: WP-F3.

The package ``__init__`` stays import-light (no submodule imports). It only holds the error types
shared by the registries, so that build-time bookkeeping errors (VT-F3) can be caught uniformly.
"""

from __future__ import annotations

__all__ = ["DoubleRegistrationError", "MissingRegistrationError"]


class DoubleRegistrationError(ValueError):
    """Something was registered twice where it must be counted once (VT-F3).

    Raised by ``LossRegistry`` (duplicate element id; second mode-overlap element of the same kind on
    one path pair, G06) and by ``BandLedger`` (second band declaration of the same source and quantity
    on a shared target, G05). Raised at build time, never during the time loop.
    """


class MissingRegistrationError(LookupError):
    """A mandatory registration is missing (e.g. a modulator on the squeezed path without a
    ``LossElement``, prompt §5.5, D10, G06)."""
