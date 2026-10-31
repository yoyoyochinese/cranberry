"""Band ledger: bookkeeping of which noise bandwidth the time loop represents (G05). Owner: WP-F3.

Every noise-source block declares, through ``Block.bands()``, the band [f_lo, f_hi] that its time
series represents and how much noise lies outside it (one-sided PSD, conv §7, or its integrated
variance). The out-of-band variance is added in exactly one place, ``phase_to_quadrature``
(``blocks/meas``), per metric and per beam combination (χ_HF; upper limit ~1/τ_p for pulses, ~κ_c
for OPO modes, conv §6). Blocks never add out-of-band noise themselves (PLAN §7.3). Registering the
same (source, quantity, target) band twice is a build-time error (VT-F3). The DC-fringe visibility
correction γ_coh (T15) is one ledger row referenced by ``dc_fringe``.

Semantics fixed by WP-F3
    * **Which side is added.** Only the part *above* ``f_hi_hz`` is added, up to the metric's upper
      limit ``f_max_hz`` (conv §6: set by the flat range of the source spectrum, ~1/τ_p for pulses,
      ~κ_c for OPO modes). Noise below ``f_lo_hz`` (≈ 1/duration) is slower than the simulated record
      and appears as a static offset, not as jitter (conv §6 lower limit 1/T, T04 eq. (4)).
    * **Out-of-band variance of one declaration** for metric m:
      ``variance_outside`` as given (independent of m), or ∫_{f_hi}^{f_max(m)} S(f)·w_m(f) df for
      ``psd_outside`` (w_m = sideband weighting of the metric, T04 §3.2, 1 for a flat source).
      Quantity ``"freq"`` [Hz] is converted with S_φ = S_ν/f² (conv §7); other quantities (rin,
      strain, dT, ...) are not phase and are ignored by :meth:`BandLedger.hf_variance` (use
      :meth:`BandLedger.oob_variance` and the block's coefficients for them).
    * **Beam combination.** A metric phase is a signed sum of beam phases, e.g. metric (i)
      δ = φ_LO,k − θ_sq,k → ``{"lo1": +1, "sqz1": -1}``. A declaration affects beam b if one of its
      targets is b or a segment of b (:meth:`BandLedger.define_beam`). It enters with
      c = Σ_{b affected} coef_b and contributes c²·v (common-mode noise on both beams of a pair
      cancels at high frequency; independent declarations add in variance). Sources whose high
      frequency part decorrelates between beams (laser noise with delay τ: 4 sin²(πfτ), T05) declare
      it with the decorrelation already folded into ``psd_outside`` on a dedicated target.
    * **Duplicates.** Two declarations with the same ``source`` and ``quantity`` that share a target
      raise :class:`~cranberry.offline.DoubleRegistrationError`, whether or not their bands overlap
      (each declaration's "outside" would cover the other's band, so even disjoint bands of one source
      would double count).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from cranberry.offline import DoubleRegistrationError

__all__ = ["BandDecl", "BandLedger", "MetricBand", "PHASE_QUANTITIES", "integrate_psd"]

#: Quantities that :meth:`BandLedger.hf_variance` converts to phase variance [rad²].
PHASE_QUANTITIES: tuple[str, ...] = ("phase", "freq")


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
        if not self.targets:
            raise ValueError(f"band of {self.source!r}: targets must be non-empty")
        if self.variance_outside is not None and not (float(self.variance_outside) >= 0.0):
            raise ValueError(f"band of {self.source!r}: variance_outside must be >= 0")
        if self.quantity == "freq" and self.variance_outside is not None and self.variance_outside != 0.0:
            raise ValueError(
                f"band of {self.source!r}: a 'freq' band needs psd_outside (S_ν cannot be turned into phase "
                "variance from an integrated value)"
            )


@dataclass(frozen=True)
class MetricBand:
    """Upper limit and sideband weighting of the band that enters one metric's jitter (conv §6).

    Attributes:
        name: metric name (e.g. ``"i.port1"``, ``"ii.port3"``).
        f_max_hz: upper integration limit [Hz] (flat range of the source spectrum; ~1/τ_p for pulses;
            ~κ_c for OPO modes; conv §6, T04 §3.4).
        weight: w(f) sideband weighting (T04 §3.2) or None for a flat source spectrum (w = 1).
    """

    name: str
    f_max_hz: float
    weight: Callable[[np.ndarray], np.ndarray] | None

    def __post_init__(self) -> None:
        if not (np.isfinite(self.f_max_hz) and self.f_max_hz > 0):
            raise ValueError(f"metric {self.name!r}: f_max_hz must be finite and > 0")


def integrate_psd(
    psd: Callable[[np.ndarray], np.ndarray],
    f_a_hz: float,
    f_b_hz: float,
    *,
    weight: Callable[[np.ndarray], np.ndarray] | None = None,
    divide_f2: bool = False,
    n_per_decade: int = 400,
) -> float:
    """∫_{f_a}^{f_b} S(f)·w(f) [/f²] df with the trapezoid rule on a log-spaced grid (f_a > 0).

    Integrated in ln f (∫ S f d ln f), so power laws over many decades are handled; the relative
    error for smooth spectra is ~(ln10/n_per_decade)² ≈ 3e−5 at the default resolution.
    """
    if f_b_hz <= f_a_hz:
        return 0.0
    if f_a_hz <= 0:
        raise ValueError("integrate_psd: lower limit must be > 0")
    decades = np.log10(f_b_hz / f_a_hz)
    n = max(64, int(np.ceil(decades * n_per_decade)) + 1)
    f = np.geomspace(f_a_hz, f_b_hz, n)
    s = np.broadcast_to(np.asarray(psd(f), dtype=float), f.shape)
    if weight is not None:
        s = s * np.broadcast_to(np.asarray(weight(f), dtype=float), f.shape)
    if divide_f2:
        s = s / (f * f)
    return float(np.trapezoid(s * f, np.log(f)))


class BandLedger:
    """Collects ``BandDecl`` rows and yields out-of-band variance per metric and beam set."""

    def __init__(
        self,
        decls: Iterable[BandDecl] = (),
        *,
        metrics: Mapping[str, MetricBand] | None = None,
        beams: Mapping[str, Iterable[str]] | None = None,
        n_per_decade: int = 400,
    ):
        self.decls: list[BandDecl] = []
        self.metrics: dict[str, MetricBand] = {}
        self.beams: dict[str, tuple[str, ...]] = {}
        self.n_per_decade = int(n_per_decade)
        for name, m in (metrics or {}).items():
            if m.name != name:
                raise ValueError(f"metrics[{name!r}] has name {m.name!r}")
            self.metrics[name] = m
        for b, segs in (beams or {}).items():
            self.define_beam(b, segs)
        for d in decls:
            self.declare(d)

    # ---------------------------------------------------------------------------------------------
    # registration

    def declare(self, decl: BandDecl) -> None:
        """Add a declaration (duplicate (source, quantity, target) is an error, VT-F3)."""
        if not isinstance(decl, BandDecl):
            raise TypeError(f"expected BandDecl, got {type(decl).__name__}")
        for d in self.decls:
            shared = set(d.targets) & set(decl.targets)
            if d.source == decl.source and d.quantity == decl.quantity and shared:
                raise DoubleRegistrationError(
                    f"band of {decl.source!r} ({decl.quantity}) declared twice on {sorted(shared)}: "
                    f"[{d.f_lo_hz}, {d.f_hi_hz}] Hz and [{decl.f_lo_hz}, {decl.f_hi_hz}] Hz (G05)"
                )
        self.decls.append(decl)

    def declare_all(self, decls: Iterable[BandDecl]) -> None:
        """Add several declarations (e.g. ``block.bands()``)."""
        for d in decls:
            self.declare(d)

    def define_metric(
        self, name: str, *, f_max_hz: float, weight: Callable[[np.ndarray], np.ndarray] | None
    ) -> MetricBand:
        """Declare metric ``name`` with its upper limit and sideband weighting (``None`` = flat)."""
        if name in self.metrics:
            raise DoubleRegistrationError(f"metric {name!r} defined twice")
        m = MetricBand(name, float(f_max_hz), weight)
        self.metrics[name] = m
        return m

    def define_beam(self, beam: str, segments: Iterable[str]) -> None:
        """Declare which segments beam ``beam`` traverses (so segment targets reach it)."""
        if beam in self.beams:
            raise DoubleRegistrationError(f"beam {beam!r} defined twice")
        self.beams[beam] = tuple(segments)

    # ---------------------------------------------------------------------------------------------
    # queries

    def find(self, *, source: str | None = None, quantity: str | None = None, target: str | None = None) -> list[BandDecl]:
        """Declarations matching all given filters (e.g. the γ_coh row read by ``dc_fringe``)."""
        return [
            d
            for d in self.decls
            if (source is None or d.source == source)
            and (quantity is None or d.quantity == quantity)
            and (target is None or target in d.targets)
        ]

    def affects(self, decl: BandDecl, beam: str) -> bool:
        """True if ``decl`` acts on ``beam`` (a target is the beam or one of its segments)."""
        segs = self.beams.get(beam, ())
        return any(t == beam or t in segs for t in decl.targets)

    def _metric(self, metric: str) -> MetricBand:
        try:
            return self.metrics[metric]
        except KeyError:
            raise KeyError(f"metric {metric!r} not defined (known: {sorted(self.metrics)})") from None

    def oob_variance(self, decl: BandDecl, metric: str) -> float:
        """Out-of-band variance of one declaration for ``metric`` in the declaration's own unit²
        (``"freq"`` → rad² via S_φ = S_ν/f²), above f_hi up to the metric's f_max."""
        m = self._metric(metric)
        if decl.variance_outside is not None:
            return float(decl.variance_outside)
        return integrate_psd(
            decl.psd_outside,
            decl.f_hi_hz,
            m.f_max_hz,
            weight=m.weight,
            divide_f2=decl.quantity == "freq",
            n_per_decade=self.n_per_decade,
        )

    def _coefficients(self, beams: Iterable[str] | Mapping[str, float]) -> tuple[dict[str, float], bool]:
        if isinstance(beams, Mapping):
            return {str(b): float(c) for b, c in beams.items()}, True
        if isinstance(beams, str):
            beams = [beams]
        return {str(b): 1.0 for b in beams}, False

    def hf_variance(self, metric: str, beams: Iterable[str] | Mapping[str, float]) -> float:
        """Out-of-band phase variance [rad²] to add for ``metric`` and the beam combination ``beams``.

        ``beams``: {beam: signed coefficient} (e.g. ``{"lo1": 1, "sqz1": -1}``), or an iterable of
        beam ids with coefficient +1 each — then a declaration acting on more than one listed beam is
        ambiguous (common mode cancels or adds depending on the signs) and raises ``ValueError``.
        """
        coef, signed = self._coefficients(beams)
        total = 0.0
        for d in self.decls:
            if d.quantity not in PHASE_QUANTITIES:
                continue
            hit = [b for b in coef if self.affects(d, b)]
            if not hit:
                continue
            if len(hit) > 1 and not signed:
                raise ValueError(
                    f"band of {d.source!r} acts on several beams {hit}: give signed coefficients {{beam: ±1}}"
                )
            c = sum(coef[b] for b in hit)
            if c != 0.0:
                total += c * c * self.oob_variance(d, metric)
        return total

    def to_json(self) -> dict[str, Any]:
        """JSON-serialisable record (``band_ledger.json``)."""
        rows = []
        for d in self.decls:
            row: dict[str, Any] = {
                "source": d.source,
                "quantity": d.quantity,
                "unit": d.unit,
                "targets": list(d.targets),
                "f_lo_hz": float(d.f_lo_hz),
                "f_hi_hz": float(d.f_hi_hz),
                "outside": "variance" if d.variance_outside is not None else "psd",
                "psd_outside": None if d.psd_outside is None else getattr(d.psd_outside, "__qualname__", repr(d.psd_outside)),
                "variance_outside": None if d.variance_outside is None else float(d.variance_outside),
                "oob_variance_by_metric": {m: self.oob_variance(d, m) for m in self.metrics},
            }
            rows.append(row)
        return {
            "decls": rows,
            "metrics": {
                n: {"f_max_hz": m.f_max_hz, "weight": None if m.weight is None else getattr(m.weight, "__qualname__", repr(m.weight))}
                for n, m in self.metrics.items()
            },
            "beams": {b: list(s) for b, s in self.beams.items()},
        }
