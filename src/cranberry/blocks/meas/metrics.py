"""Metric output containers and the metric blocks (conv §10, D12, D31, D46, D47, G13). Owner: WP-F4.

Metrics (conv §10):

* (i)  per port k at its BHD observer: δ_i,k(t) = wrap_π(φ_LO,k − θ_sq,k − Δ_set) → :class:`MetricI`.
* (ii) per port k ≥ 1 at the GBS inputs, pulses of the same GBS shot:
  Θ_k(t) = wrap_π(θ_sq,k − θ_sq,0 − Θ_k,set) → :class:`MetricII`, with the four components
  (:data:`METRIC_II_COMPONENTS`) obtained by reruns + linearity check (G13, ``runner.decompose``):

  1. ``locked_path_residual``     residual of the locked path noise
  2. ``common_points``            motion of the common points (probe split, reference delivery)
  3. ``laser_decorrelation_oob``  laser phase decorrelation φ_L(t − τ_j) − φ_L(t − τ_k) outside the
                                  lock band (+ in-band 16 sin²(πf(T_a − T_b)) sin²(πf(τ_k − τ_0)) S_φ
                                  term, D51-5)
  4. ``squeezed_only``            segments only the squeezed light traverses (incl. demux phase)

* alongside (ii): wave-packet overlap between ports (:class:`Overlap`), stray light in both units of
  D46 (:class:`StrayLight`), η_k at GBS input and at BHD separately (:class:`Efficiency`), and each
  port's own σ (D46).

Port convention (G13): port-based configs always include reference port 0. Port-less configs
(Ref A–D) report metric (i) under port key 0 when they have a squeezed beam and a BHD, and
``metric_ii_status = "not_applicable"`` with an empty ``metric_ii``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cranberry.core.block import Block
from cranberry.core.registry import register_block

__all__ = [
    "METRIC_II_COMPONENTS",
    "METRIC_II_STATUS",
    "wrap_pi",
    "MetricI",
    "MetricIIComponents",
    "MetricII",
    "Overlap",
    "StrayLight",
    "Efficiency",
    "MetricsResult",
    "make_metric_i",
    "make_metric_ii",
]

METRIC_II_COMPONENTS: tuple[str, ...] = (
    "locked_path_residual",
    "common_points",
    "laser_decorrelation_oob",
    "squeezed_only",
)
#: ``ok``: computed; ``not_applicable``: no ports (G13); ``decomposition_invalid``: saturation / lock
#: loss / linearity check failed; ``not_run``: decomposition not requested.
METRIC_II_STATUS: tuple[str, ...] = ("ok", "not_applicable", "decomposition_invalid", "not_run")


def wrap_pi(x):
    """Wrap angles to (−π, π] (conv §10 wrap_π)."""
    return np.pi - np.mod(np.pi - np.asarray(x, dtype=float), 2.0 * np.pi)


@dataclass
class MetricI:
    """Metric (i) of port k: δ_i,k(t) [rad] at its BHD observer (time axis = ``MetricsResult.t_s``)."""

    port: int
    observer: str
    delta_rad: np.ndarray
    setpoint_rad: float

    @property
    def rms_rad(self) -> float:
        return float(np.sqrt(np.mean(np.square(self.delta_rad))))


@dataclass
class MetricIIComponents:
    """Variance split of Θ_k [rad²] (G13).

    Attributes:
        variances_rad2: {component name (METRIC_II_COMPONENTS): variance}.
        total_var_rad2: variance of the all-on run.
        linearity_ratio: Σ components / total (valid if within 5 %, T04 V12).
        status: "ok" or "decomposition_invalid".
        reason: why invalid (saturation, lock loss, linearity), or None.
    """

    variances_rad2: Mapping[str, float]
    total_var_rad2: float
    linearity_ratio: float
    status: str
    reason: str | None = None


@dataclass
class MetricII:
    """Metric (ii) of port k ≥ 1: Θ_k(t) [rad] at the GBS input, plus its decomposition."""

    port: int
    observer: str
    theta_rad: np.ndarray
    setpoint_rad: float
    components: MetricIIComponents | None = None

    @property
    def rms_rad(self) -> float:
        return float(np.sqrt(np.mean(np.square(self.theta_rad))))


@dataclass
class Overlap:
    """Wave-packet overlap of port k with port 0 at the GBS input (∈ [0, 1]; series or constant)."""

    port: int
    overlap: np.ndarray


@dataclass
class StrayLight:
    """Stray (probe-derived) light in port k (D46): clicks per detection gate (both polarizations) and
    photons in the squeezed mode (N of conv §2)."""

    port: int
    clicks_per_gate: float
    photons_per_mode: float


@dataclass
class Efficiency:
    """Total transmission of port k's squeezed light to the GBS input and to its BHD (D46, G06)."""

    port: int
    eta_gbs: float | None
    eta_bhd: float | None


@dataclass
class MetricsResult:
    """All metrics of one run."""

    t_s: np.ndarray
    metric_i: dict[int, MetricI] = field(default_factory=dict)
    metric_ii: dict[int, MetricII] = field(default_factory=dict)
    metric_ii_status: str = "not_run"
    sigma_rad: dict[int, float] = field(default_factory=dict)
    overlap: dict[int, Overlap] = field(default_factory=dict)
    stray_light: dict[int, StrayLight] = field(default_factory=dict)
    efficiency: dict[int, Efficiency] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.metric_ii_status not in METRIC_II_STATUS:
            raise ValueError(f"metric_ii_status {self.metric_ii_status!r} not in {METRIC_II_STATUS}")

    def to_json(self) -> dict[str, Any]:
        """Summary for ``metrics.json`` (series go to results.npz). WP-F4."""
        raise NotImplementedError("MetricsResult.to_json is implemented by WP-F4")


@register_block("metric_i")
def make_metric_i(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the metric (i) probe block. Implemented by WP-F4."""
    raise NotImplementedError("block 'metric_i' is implemented by WP-F4")


@register_block("metric_ii")
def make_metric_ii(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the metric (ii) probe block. Implemented by WP-F4."""
    raise NotImplementedError("block 'metric_ii' is implemented by WP-F4")
