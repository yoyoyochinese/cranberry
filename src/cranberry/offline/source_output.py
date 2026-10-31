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

Frame of ``cov`` and meaning of θ_off (WP-F3; T03 §2 "θ_off = angle of the minimum − θ_sq,0")
    ``cov`` is expressed in the lab quadrature frame of conv §2 **at pump phase φ_p = 0**. The
    nominal squeezing angle there is θ_sq,0 = φ_OPA,0/2 + π/2 (conv §4). The time loop adds the
    dynamic part (φ_p(t)/2, pump-path segments ½ψ, path phases Σψ; PLAN §2.3) to θ_sq,0, so
    :meth:`SourceOutput.project` returns θ_off = θ_min − θ_sq,0 wrapped to (−π/2, π/2]: the static
    rotation of the LO-mode ellipse that the lock does not see (OPO detuning, temporal-mode
    mismatch, GDD; T03, T17). For an ideal single-mode source θ_off = 0.

Temporal projection (T17 §2.2, T03 §2): the LO mode u(T) (complex allowed, normalised here) defines
â_LO = ∫u*(T) â(T) dT = Σ_n c_n â_n + √(1 − Σ|c_n|²) v̂ with c_n = ∫u* f_n dT (Riemann sum with
weights ``np.gradient(grid)``, exact for the uniform bins of T17) and v̂ vacuum (the part of u outside
span{f_n} sees vacuum — that is the mode-mismatch pure loss, T12 §1.3).

``opo.pulsed = off`` must reproduce the CW OPO result (T17 V-P1; test of WP-B11).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cranberry.offline import gaussian
from cranberry.offline.tables import Table

__all__ = ["SOURCE_BASES", "SourceOutput", "wrap_half_pi"]

SOURCE_BASES: tuple[str, ...] = ("cw_sideband", "temporal", "single_mode")


def wrap_half_pi(angle_rad):
    """Wrap a quadrature angle (defined mod π) to (−π/2, π/2]."""
    a = np.asarray(angle_rad, dtype=float)
    return np.pi / 2 - np.mod(np.pi / 2 - a, np.pi)


@dataclass(frozen=True)
class SourceOutput:
    """Source model output.

    Attributes:
        basis: one of :data:`SOURCE_BASES`.
        grid: f_s grid [Hz] (cw_sideband), T grid [s] (temporal), None (single_mode).
        cov: covariance in the basis (module docstring for shapes and frame), vacuum = 1/2·𝟙.
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
        cov = np.array(self.cov, dtype=float)
        if not np.isfinite(float(self.phi_opa0_rad)):
            raise ValueError("phi_opa0_rad must be finite")
        if self.basis == "single_mode":
            if self.grid is not None or self.modes is not None:
                raise ValueError("single_mode: grid and modes must be None")
            if cov.shape != (2, 2):
                raise ValueError(f"single_mode: cov must be (2, 2), got {cov.shape}")
            grid = None
            modes = None
        elif self.basis == "cw_sideband":
            grid = _ascending(self.grid, "cw_sideband grid (f_s [Hz])")
            if self.modes is not None:
                raise ValueError("cw_sideband: modes must be None")
            if cov.shape != (grid.size, 2, 2):
                raise ValueError(f"cw_sideband: cov must be (n_f, 2, 2) = ({grid.size}, 2, 2), got {cov.shape}")
            modes = None
        else:  # temporal
            grid = _ascending(self.grid, "temporal grid (T [s])")
            if grid.size < 2:
                raise ValueError("temporal: the T grid needs at least 2 points")
            if self.modes is None:
                raise ValueError("temporal: modes (N, n_T) are required")
            modes = np.array(self.modes)
            if modes.ndim != 2 or modes.shape[1] != grid.size:
                raise ValueError(f"temporal: modes must be (N, {grid.size}), got {modes.shape}")
            n = modes.shape[0]
            if cov.shape != (2 * n, 2 * n):
                raise ValueError(f"temporal: cov must be (2N, 2N) = ({2 * n}, {2 * n}), got {cov.shape}")
            gram = (modes.conj() * np.gradient(grid)) @ modes.T
            if not np.allclose(gram, np.eye(n), atol=1e-6):
                raise ValueError("temporal: modes must be orthonormal (∫ f_m* f_n dT = δ_mn within 1e-6)")
            modes.flags.writeable = False
        if not np.allclose(cov, np.swapaxes(cov, -1, -2), rtol=0.0, atol=1e-12 * max(1.0, np.max(np.abs(cov)))):
            raise ValueError("cov must be symmetric")
        if self.tail_into_next_slot is not None and not (0.0 <= self.tail_into_next_slot <= 1.0):
            raise ValueError("tail_into_next_slot must be in [0, 1] or None")
        cov.flags.writeable = False
        object.__setattr__(self, "cov", cov)
        object.__setattr__(self, "grid", grid)
        object.__setattr__(self, "modes", modes)
        object.__setattr__(self, "phi_opa0_rad", float(self.phi_opa0_rad))
        object.__setattr__(self, "pump_path_phase_rad", dict(self.pump_path_phase_rad))
        object.__setattr__(self, "coef_tables", dict(self.coef_tables))
        object.__setattr__(self, "provenance", dict(self.provenance))

    # ----------------------------------------------------------------------------------------------

    @property
    def theta_sq0_rad(self) -> float:
        """Nominal squeezing angle of ``cov`` at φ_p = 0: θ_sq,0 = φ_OPA,0/2 + π/2 (conv §4)."""
        return 0.5 * self.phi_opa0_rad + 0.5 * np.pi

    @property
    def n_modes(self) -> int:
        """Number of modes in the basis (1 for single_mode and per sideband for cw_sideband)."""
        return 1 if self.modes is None else self.modes.shape[0]

    def lo_mode_cov(self, lo_mode: Any) -> np.ndarray:
        """2×2 covariance (q, p) of the LO mode (cw_sideband: (n_f, 2, 2)); see :meth:`project`."""
        if self.basis in ("single_mode", "cw_sideband"):
            if lo_mode is not None:
                raise ValueError(f"{self.basis}: lo_mode must be None (the basis is already one mode per entry)")
            return np.array(self.cov)
        if lo_mode is None:
            raise ValueError("temporal: lo_mode u(T) on the grid is required")
        u = np.asarray(lo_mode)
        if u.shape != self.grid.shape:
            raise ValueError(f"temporal: lo_mode must have shape {self.grid.shape}, got {u.shape}")
        w = np.gradient(self.grid)
        norm = float(np.sum(np.abs(u) ** 2 * w))
        if not norm > 0.0:
            raise ValueError("temporal: lo_mode has zero norm")
        u = u / np.sqrt(norm)  # strong-LO homodyne measures the shape only
        c = (self.modes * (u.conj() * w)).sum(axis=1)  # c_n = ∫ u* f_n dT
        cap = float(np.sum(np.abs(c) ** 2))
        if cap > 1.0 + 1e-6:
            raise ValueError(f"temporal: Σ|c_n|² = {cap} > 1 (modes not orthonormal on this grid?)")
        cap = min(cap, 1.0)
        # â_LO = Σ c_n â_n  ⇒  q_LO = Σ (Re c q_n − Im c p_n),  p_LO = Σ (Im c q_n + Re c p_n)
        a = np.block([[c.real[None, :], -c.imag[None, :]], [c.imag[None, :], c.real[None, :]]])
        return a @ self.cov @ a.T + (1.0 - cap) * 0.5 * np.eye(2)

    def project(self, lo_mode: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Project onto the LO mode → (V_−, V_+, θ_off) [variances in conv §2 units, θ_off rad].

        For ``cw_sideband`` each output is an array over ``grid``; otherwise 0-d. ``lo_mode`` is a
        temporal mode function on ``grid`` (temporal) or None (cw_sideband / single_mode).
        θ_off = θ_min − θ_sq,0 ∈ (−π/2, π/2] (module docstring).
        """
        vm, vp, th = gaussian.noise_levels_from_cov(self.lo_mode_cov(lo_mode))
        return np.asarray(vm), np.asarray(vp), np.asarray(wrap_half_pi(np.asarray(th) - self.theta_sq0_rad))

    # ----------------------------------------------------------------------------------------------

    @classmethod
    def single_mode_reference(
        cls,
        *,
        loss: float,
        alpha_per_w: float,
        pump_w: float,
        phi_opa0_rad: float,
        provenance: Mapping[str, Any] | None = None,
    ) -> SourceOutput:
        """Prompt §7-1 single-mode reference (Ref A), N = 1 (protected list "Ref A blocks").

        R± = L + (1 − L)·exp(±2√(αP)) (shot-noise units), i.e. a squeezed vacuum with r = √(αP)
        followed by a pure loss η = 1 − L (T12 §1.2). Returned as ``single_mode`` with
        V± = R±·V_SN at the nominal angle θ_sq,0 (θ_off = 0). Phase jitter is *not* included: apply
        :func:`cranberry.offline.gaussian.jitter_average` to the projected (V_−, V_+).

        Args:
            loss: total effective loss L ∈ [0, 1] (lumped, D48).
            alpha_per_w: SHG efficiency α [1/W] (Ref A: 906 %/W = 9.06 /W).
            pump_w: pump power P [W].
            phi_opa0_rad: φ_OPA,0 (conv §4; required, D34).
        """
        if not (0.0 <= loss <= 1.0):
            raise ValueError(f"loss must be in [0, 1], got {loss}")
        if alpha_per_w < 0 or pump_w < 0:
            raise ValueError("alpha_per_w and pump_w must be >= 0")
        r = float(np.sqrt(alpha_per_w * pump_w))
        prov = {
            "model": "single_mode_reference",
            "formula": "R± = L + (1 − L) exp(±2 sqrt(alpha P)) [prompt §7-1, Ref A]",
            "loss": float(loss),
            "alpha_per_w": float(alpha_per_w),
            "pump_w": float(pump_w),
            "r": r,
        }
        prov.update(provenance or {})
        theta0 = 0.5 * float(phi_opa0_rad) + 0.5 * np.pi
        cov = gaussian.pure_loss(gaussian.squeezed_vacuum(r, theta0), 1.0 - float(loss))
        return cls(basis="single_mode", grid=None, cov=cov, phi_opa0_rad=phi_opa0_rad, provenance=prov)


def _ascending(grid: Any, what: str) -> np.ndarray:
    if grid is None:
        raise ValueError(f"{what} is required")
    g = np.array(grid, dtype=float)
    if g.ndim != 1 or g.size < 1 or (g.size > 1 and not np.all(np.diff(g) > 0)):
        raise ValueError(f"{what} must be a non-empty strictly ascending 1-D array")
    g.flags.writeable = False
    return g
