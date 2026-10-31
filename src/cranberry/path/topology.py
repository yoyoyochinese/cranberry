"""Static ``Topology``: segments, beams, coefficients, locks, ports, observers (protected list "per-path
phase model"; PLAN §2.3.1, conv §9, D11, D47). Owner: WP-F6.

The Topology is fixed at build time; the dynamic quantities live in State slots (``drv.*``, ``psi.*``,
``tau_g.*``, ``J.*``, ``beam.*``; ``core.layout``). Cancellation on shared segments is computed, never
assumed (D11): beam phases are sums over their paths evaluated at the times each beam passes each
subsegment (``path.beams``), and segment × beam coefficients may differ (VT-23).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cranberry.path.coef import CoefKey, CoefValue

__all__ = [
    "SEGMENT_KINDS",
    "SPATIAL_STRUCTURES",
    "POL_MODELS",
    "Segment",
    "PathStep",
    "Beam",
    "Lock",
    "Port",
    "Observer",
    "Topology",
]

SEGMENT_KINDS: tuple[str, ...] = ("pm", "smf", "free_space", "opo", "crystal", "demux", "component")
#: Spatial structure of in-segment noise (conv §9 一様 / 独立 / 混合 = lumped / distributed / mixed).
SPATIAL_STRUCTURES: tuple[str, ...] = ("lumped", "distributed", "mixed")
#: Polarization model per segment (G04): full Jones, PM diagonal fast path, or off (J = e^{iψ}·1).
POL_MODELS: tuple[str, ...] = ("jones", "pm_diagonal", "off")


@dataclass(frozen=True)
class Segment:
    """One path segment (conv §9).

    Attributes:
        id: segment id (used in slot names ``drv.<id>.*``, ``psi.<id>.<beam>``).
        kind: one of :data:`SEGMENT_KINDS`.
        length_m: length L [m].
        n_g: group index.
        beta2_s2_per_m, beta3_s3_per_m: dispersion [s²/m], [s³/m].
        zone: thermal zone id (segments in one zone share ambient temperature).
        n_sub: number of subsegments N_sub (timing of counter-propagation, T08/T10).
        spatial: one of :data:`SPATIAL_STRUCTURES`.
        pol_model: one of :data:`POL_MODELS`; ``pol_params`` holds its parameters.
        actuator_sub: subsegment index x_a of an actuator on this segment, or None.
    """

    id: str
    kind: str
    length_m: float
    n_g: float
    beta2_s2_per_m: float
    beta3_s3_per_m: float
    zone: str
    n_sub: int
    spatial: str
    pol_model: str
    pol_params: Mapping[str, Any] = field(default_factory=dict)
    actuator_sub: int | None = None

    @property
    def tau_g0_s(self) -> float:
        """Nominal group delay n_g·L/c [s]."""
        return self.n_g * self.length_m / 299_792_458.0


@dataclass(frozen=True)
class PathStep:
    """One traversal of a segment by a beam: direction d = +1 (source → detector/GBS) or −1, entering
    and leaving subsegment indices (conv §9). A counter-propagating beam samples the same subsegments
    at different times; it sees Jᵀ (reciprocity, G04)."""

    segment: str
    direction: int
    sub_in: int
    sub_out: int


@dataclass(frozen=True)
class Beam:
    """A beam: source, ordered path, frequency offset, input polarization, nominal power.

    Attributes:
        id: beam id (``beam.<id>.*`` slots).
        source: laser / pump / source block name it originates from.
        path: ordered :class:`PathStep` tuple.
        freq_offset_hz: optical frequency offset from the source carrier [Hz] (conv §4 sign).
        jones_in: input Jones vector (complex, shape (2,)).
        power_w: nominal power [W].
        observer: observer id where the beam phase is evaluated.
    """

    id: str
    source: str
    path: tuple[PathStep, ...]
    freq_offset_hz: float
    jones_in: np.ndarray
    power_w: float
    observer: str


@dataclass(frozen=True)
class Lock:
    """Lock = linear combination of beam phases held at a setpoint (conv §9).

    error e(t) = Σ_b combo[b]·φ_b(t) − setpoint_rad, evaluated at ``observer``. Device-specific
    error signs are converted inside the device's block (conv §9).
    """

    id: str
    combo: Mapping[str, float]
    setpoint_rad: float
    observer: str


@dataclass(frozen=True)
class Port:
    """Port k (D12, D31, G13). Port 0 is the reference port and always exists in port-based configs.

    ``squeezed`` / ``lo`` / ``probe`` are beam ids (None if absent); ``gbs_observer`` / ``bhd_observer``
    are observer ids (None if absent)."""

    index: int
    squeezed: str
    lo: str | None
    probe: str | None
    gbs_observer: str | None
    bhd_observer: str | None


@dataclass(frozen=True)
class Observer:
    """Point where beam phases are evaluated (detector, GBS input, interaction point).

    ``kind`` ∈ {"detector", "gbs_input", "interaction"}; ``segment``/``sub`` locate it."""

    id: str
    kind: str
    segment: str | None = None
    sub: int | None = None


@dataclass(frozen=True)
class Topology:
    """Static description of the optical layout (PLAN §2.3.1)."""

    segments: Mapping[str, Segment]
    beams: Mapping[str, Beam]
    coef: Mapping[CoefKey, CoefValue]
    locks: Mapping[str, Lock]
    ports: Mapping[int, Port]
    observers: Mapping[str, Observer]

    @classmethod
    def from_spec(cls, tree: Mapping[str, Any]) -> Topology:
        """Parse ``RunSpec.topology``. Implemented by WP-F6."""
        raise NotImplementedError("path.topology.Topology.from_spec is implemented by WP-F6")

    @property
    def has_ports(self) -> bool:
        """False for port-less configs (Ref A–D): metric (ii) is then "not applicable" (G13)."""
        return bool(self.ports)
