"""Band ledger: bookkeeping of which noise bandwidth the time loop represents (G05). Owner: WP-F3.

Every noise-source block declares, through ``Block.bands()``, the band [f_lo, f_hi] that its time
series represents and how much noise lies outside it (one-sided PSD, conv §7, or its integrated
variance). The out-of-band variance is added in exactly one place, ``phase_to_quadrature``
(``blocks/meas``), per metric and per beam combination (χ_HF; upper limit ~1/τ_p for pulses, ~κ_c
for OPO modes, conv §6). Blocks never add out-of-band noise themselves (PLAN §7.3). Registering the
same (source, quantity, target) band twice is a build-time error (VT-F3). The DC-fringe visibility
correction γ_coh (T15) is one ledger row referenced by ``dc_fringe``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = ["BandDecl", "BandLedger"]


@dataclass(frozen=True)
class BandDecl:
    """One band declaration.

    Attributes:
        source: declaring block name (noise source).
        quantity: what the noise is in: ``"phase"`` [rad], ``"freq"`` [Hz], ``"rin"`` [1],
            ``"strain"`` [1], ``"dT"`` [K], ... (``unit`` gives the unit explicitly).
        targets: beam ids or segment ids the noise acts on.
        f_lo_hz: lower edge of the band represented in the time loop [Hz] (≈ 1 / duration).
        f_hi_hz: upper edge [Hz] (≤ Nyquist of the block's rate, f_sim / (2·rate_div)).
        unit: unit of the quantity (PSD unit is unit²/Hz).
        psd_outside: one-sided PSD S(f) [unit²/Hz] valid outside [f_lo, f_hi], or None.
        variance_outside: integrated out-of-band variance [unit²], or None. Exactly one of
            ``psd_outside`` / ``variance_outside`` must be given (use 0.0 for "none").
    """

    source: str
    quantity: str
    targets: tuple[str, ...]
    f_lo_hz: float
    f_hi_hz: float
    unit: str
    psd_outside: Callable[[np.ndarray], np.ndarray] | None = None
    variance_outside: float | None = None

    def __post_init__(self) -> None:
        if not (0.0 <= self.f_lo_hz < self.f_hi_hz):
            raise ValueError(f"band of {self.source!r}: need 0 <= f_lo_hz < f_hi_hz")
        if (self.psd_outside is None) == (self.variance_outside is None):
            raise ValueError(f"band of {self.source!r}: give exactly one of psd_outside / variance_outside")
        object.__setattr__(self, "targets", tuple(self.targets))


class BandLedger:
    """Collects ``BandDecl`` rows and yields out-of-band variance per metric and beam set."""

    def __init__(self, decls: Iterable[BandDecl] = ()):
        self.decls: list[BandDecl] = []
        for d in decls:
            self.declare(d)

    def declare(self, decl: BandDecl) -> None:
        """Add a declaration (duplicate (source, quantity, target) is an error, VT-F3)."""
        raise NotImplementedError("offline.band_ledger.BandLedger.declare is implemented by WP-F3")

    def hf_variance(self, metric: str, beams: Iterable[str]) -> float:
        """Out-of-band phase variance [rad²] to add for ``metric`` and the beam combination ``beams``."""
        raise NotImplementedError("offline.band_ledger.BandLedger.hf_variance is implemented by WP-F3")

    def to_json(self) -> dict[str, Any]:
        """JSON-serialisable record (``band_ledger.json``)."""
        raise NotImplementedError("WP-F3")
