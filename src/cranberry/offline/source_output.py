"""Common output of every source model (G07, PLAN §2.2, §4.4). Owner: WP-F3.

The time loop and the measurement layer never branch on the source type: every SRC offline block
(single-mode reference, waveguide OPA, supermodes, CW OPO, pulsed OPO) returns a ``SourceOutput``.

Quadrature units (conv §2): ħ = 1, vacuum covariance = (1/2)·𝟙. Covariance ordering for N modes is
``xxpp``: (q_1, ..., q_N, p_1, ..., p_N). Angles follow conv §3–§4 (θ_sq = squeezing angle, mod π).

Basis kinds (:data:`SOURCE_BASES`):

* ``"cw_sideband"``: ``grid`` = sideband frequencies f_s [Hz]; ``cov`` shape (n_f, 2, 2).
* ``"temporal"``  : ``grid`` = time grid T [s]; ``modes`` (N, n_T) orthonormal (∫|u|² dT = 1);
  ``cov`` shape (2N, 2N).
* ``"single_mode"``: ``grid`` = None; ``cov`` shape (2, 2) (Ref A single-mode reference, N = 1).

``opo.pulsed = off`` must reproduce the CW OPO result (T17 V-P1; test of WP-B11).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cranberry.offline.tables import Table

__all__ = ["SOURCE_BASES", "SourceOutput"]

SOURCE_BASES: tuple[str, ...] = ("cw_sideband", "temporal", "single_mode")


@dataclass(frozen=True)
class SourceOutput:
    """Source model output.

    Attributes:
        basis: one of :data:`SOURCE_BASES`.
        grid: f_s grid [Hz] (cw_sideband), T grid [s] (temporal), None (single_mode).
        cov: covariance in the basis (module docstring for shapes), vacuum = 1/2·𝟙.
        phi_opa0_rad: per-OPA/OPO phase offset φ_OPA,0 (conv §4): θ_amp = (φ_p + φ_OPA,0)/2.
        modes: temporal mode functions (temporal only), else None.
        pump_path_phase_rad: {segment id: phase [rad] of that pump-path segment entering θ_sq as ½ψ}
            (crystal ψ_p,xtal; resonant pump cavity −arctan(Δ_p/κ_p); conv §4).
        gain_table: A±(P) amplitude-gain table for error signals, or None.
        coef_tables: segment × beam coefficient tables (e.g. OPO detuning, G16-a), by name.
        tail_into_next_slot: fraction of the mode in the next slot (N11), or None (CW).
        provenance: source tags / spec hash.
    """

    basis: str
    grid: np.ndarray | None
    cov: np.ndarray
    phi_opa0_rad: float  # required: a modelling choice per OPA/OPO (conv §4), no default (D34)
    modes: np.ndarray | None = None
    pump_path_phase_rad: Mapping[str, float] = field(default_factory=dict)
    gain_table: Table | None = None
    coef_tables: Mapping[str, Table] = field(default_factory=dict)
    tail_into_next_slot: float | None = None
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.basis not in SOURCE_BASES:
            raise ValueError(f"basis {self.basis!r} not in {SOURCE_BASES}")

    def project(self, lo_mode: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Project onto the LO mode → (V_−, V_+, θ_off) [variances in conv §2 units, θ_off rad].

        For ``cw_sideband`` each output is an array over ``grid``; otherwise 0-d. ``lo_mode`` is a
        temporal mode function on ``grid`` (temporal) or None (cw_sideband / single_mode). WP-F3.
        """
        raise NotImplementedError("offline.source_output.SourceOutput.project is implemented by WP-F3")
