"""``mode_overlap_table``: temporal-mode overlap of the squeezed pulse and its LO under a timing error
and a dispersion mismatch (T03 §2.4–§2.7, §8; T17 §5.2). Owner: WP-B11.

Offline (O) block. Product: Table ``<name>.overlap`` on axes (``dt_s``, ``dL_m``) with values
(V_−, V_+, θ_off, η_t) — what the time loop interpolates from Δt_k(t) (T03 §8.2; ``Table.interp``).

D51-6 (δL_k vs δτ_k): the two axes are independent. ``dt_s`` is the residual timing δτ_k = τ_s − τ_LO
(squeezed later → positive; µm-scale path errors after the delay trim); ``dL_m`` is the length
mismatch that sets the dispersion difference φ2 = β2·δL, φ3 = β3·δL (cm–100 m scale). The squeezed side
gets ũ_s = ũ e^{i(ΩΔt + φ2Ω²/2 + φ3Ω³/6)} (conv §8); the LO is the source's optimal mode (matched at
Δt = δL = 0), so each entry is ``source.project(U†u_LO)`` and η_t = |∫u_LO* U u_LO dT|². For a single
mode θ_off = arg κ (T03 §2.5: the GDD turns timing into a squeezing-angle offset the lock cannot see).

η_t bookkeeping (M5-09): V_± already contain η_t (and θ_off the mode-shape rotation), i.e. the
``projection`` convention of ``lo_temporal_shape``; the ``eta_t`` column is informational. A preset that
interpolates V_± must not register an η_t loss for the same LO as well (``lo_temporal_shape`` with
``eta_t_carrier: projection`` registers its element disabled).

Fast timing jitter (T03 §2.6, §8.2): with ``jitter_fast`` on, the LO-mode covariance is averaged over
Δt ~ N(dt, σ_t,fast²) (Gauss–Hermite, ``n_jitter`` nodes) before (V_−, V_+, θ_off) are formed — the
measurement is linear in that covariance. Slow jitter belongs to the time loop (no double counting).

Closed forms used by the validation tests (Gaussian mode, intensity ∝ e^{−T²/T0²}):

* :func:`gaussian_overlap`: κ = (1 − iχ)^{−1/2} exp[−Δt²(1 + iχ)/(4T0²(1 + χ²))], χ = φ2/(2T0²)
* :func:`jitter_averaged_eta`: ⟨η_t⟩ = (1 + σ²/T0²)^{−1/2} exp[−μ²/(2(T0² + σ²))]
* :func:`net_dispersion`: (φ2, φ3) of the squeezed path minus the LO path (identical paths → 0, V9)
* :func:`timing_after_phase_lock`: Δτ_res = Δτ_th(1 − R_τ,st/R_τ,th) (TEC correction → 0, V10)

Params (D34): ``source`` (block spec of a temporal source, or {gaussian: {T0_s, r, n_points,
window_s, phi_opa0_rad}}), ``dt_grid_s``, ``dL_grid_m`` (ascending), ``beta2_s2_per_m``,
``beta3_s3_per_m``, ``pad`` (zero padding factor of the FFT), ``toggles`` {dispersion (off → φ2 = φ3 =
0), jitter_fast}; jitter_fast on: ``jitter_fast_s``, ``n_jitter``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from functools import cached_property
from typing import Any

import numpy as np

from cranberry.blocks.source import temporal_modes as tm
from cranberry.core.block import BaseBlock, Block, BlockSpec, require
from cranberry.core.registry import make_block, register_block
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput, wrap_half_pi
from cranberry.offline.tables import Table

__all__ = ["ModeOverlapTable", "make_mode_overlap_table", "gaussian_overlap", "gaussian_mode", "gaussian_source",
           "jitter_averaged_eta", "net_dispersion", "timing_after_phase_lock", "source_from_param"]

KIND = "mode_overlap_table"


def gaussian_overlap(dt_s, phi2_s2, t0_s) -> np.ndarray:
    """κ for equal Gaussian modes (module docstring; T03 §2.5)."""
    chi = np.asarray(phi2_s2, dtype=float) / (2 * t0_s**2)
    dt = np.asarray(dt_s, dtype=float)
    return (1 - 1j * chi) ** -0.5 * np.exp(-(dt**2) * (1 + 1j * chi) / (4 * t0_s**2 * (1 + chi**2)))


def gaussian_mode(grid: np.ndarray, t0_s: float) -> np.ndarray:
    """u(T) = (πT0²)^{−1/4} e^{−T²/(2T0²)} (intensity ∝ e^{−T²/T0²}, ∫u² = 1)."""
    return (np.pi * t0_s**2) ** -0.25 * np.exp(-0.5 * (np.asarray(grid) / t0_s) ** 2)


def gaussian_source(t0_s: float, r: float, n_points: int, window_s: float, phi_opa0_rad: float) -> SourceOutput:
    """Single-mode temporal ``SourceOutput``: squeezed vacuum (r) in a Gaussian mode (T03 §2.4, N = 1)."""
    grid = (np.arange(int(n_points)) - int(n_points) // 2) * (window_s / int(n_points))
    u = gaussian_mode(grid, t0_s)
    u = u / np.sqrt(np.sum(u**2) * np.gradient(grid))
    cov = gaussian.squeezed_vacuum(r, 0.5 * phi_opa0_rad + 0.5 * np.pi)
    return SourceOutput(basis="temporal", grid=grid, cov=cov, phi_opa0_rad=phi_opa0_rad, modes=u[None, :],
                        provenance={"model": "gaussian_single_mode", "T0_s": t0_s, "r": r})


def source_from_param(obj: Mapping[str, Any]) -> SourceOutput:
    """``source`` param of the B11 offline blocks: a block spec (``{block, name?, params}``) of a source
    with ``source_output``, or ``{gaussian: {...}}`` (:func:`gaussian_source`)."""
    if "gaussian" in obj:
        g = obj["gaussian"]
        require(g, "T0_s", "r", "n_points", "window_s", "phi_opa0_rad", block="source.gaussian")
        return gaussian_source(float(g["T0_s"]), float(g["r"]), int(g["n_points"]), float(g["window_s"]),
                               float(g["phi_opa0_rad"]))
    return make_block(BlockSpec.from_obj(obj)).source_output


def jitter_averaged_eta(mu_s, sigma_s, t0_s) -> np.ndarray:
    """⟨η_t⟩ over Δt ~ N(μ, σ²) for equal Gaussian modes, φ2 = 0 (T03 §2.6)."""
    mu, sg = np.asarray(mu_s, dtype=float), np.asarray(sigma_s, dtype=float)
    return (1 + sg**2 / t0_s**2) ** -0.5 * np.exp(-(mu**2) / (2 * (t0_s**2 + sg**2)))


def net_dispersion(path_s: Iterable[Mapping[str, float]], path_lo: Iterable[Mapping[str, float]]) -> tuple[float, float]:
    """(φ2, φ3) = Σ_squeezed (β2 L, β3 L) − Σ_LO (β2 L, β3 L); segments {L_m, beta2_s2_per_m, beta3_s3_per_m}."""
    def tot(path):
        p = list(path)
        return (sum(s["L_m"] * s["beta2_s2_per_m"] for s in p), sum(s["L_m"] * s["beta3_s3_per_m"] for s in p))
    (a2, a3), (b2, b3) = tot(path_s), tot(path_lo)
    return a2 - b2, a3 - b3


def timing_after_phase_lock(dtau_th_s: float, r_tau_st: float, r_tau_th: float, tec: bool) -> float:
    """Residual timing after an ideal phase lock cancels the thermal phase with a stretcher (T03 §2.7):
    Δτ_th (1 − R_τ,st/R_τ,th); a TEC acts through the same mechanism and leaves 0."""
    return 0.0 if tec else dtau_th_s * (1.0 - r_tau_st / r_tau_th)


class ModeOverlapTable(BaseBlock):
    """The ``mode_overlap_table`` block (module docstring)."""

    toggle_names = ("dispersion", "jitter_fast")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "source", "dt_grid_s", "dL_grid_m", "beta2_s2_per_m", "beta3_s3_per_m", "pad", block=name)
        super().__init__(KIND, name, params)
        if self._toggles["jitter_fast"]:
            require(params, "jitter_fast_s", "n_jitter", block=name)
        for k in ("dt_grid_s", "dL_grid_m"):
            a = np.asarray(params[k], dtype=float)
            if a.ndim != 1 or a.size < 1 or np.any(np.diff(a) <= 0):
                raise ValueError(f"block {name!r}: {k} must be a non-empty ascending list")

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    @cached_property
    def source(self) -> SourceOutput:
        so = source_from_param(self.params["source"])
        if so.basis != "temporal":
            raise ValueError(f"block {self.name!r}: source must be temporal (CW has no timing state, T03 §6)")
        return so

    @cached_property
    def lo_mode(self) -> np.ndarray:
        return tm.optimal_lo_mode(self.source)[0]

    def _phis(self, dl_m: float) -> tuple[float, float]:
        if not self._toggles["dispersion"]:
            return 0.0, 0.0
        return float(self.params["beta2_s2_per_m"]) * dl_m, float(self.params["beta3_s3_per_m"]) * dl_m

    def _point(self, dt_s: float, dl_m: float) -> tuple[np.ndarray, float]:
        """(2×2 LO-mode covariance, η_t) at one timing / length point (no jitter average)."""
        so, u, pad = self.source, self.lo_mode, int(self.params["pad"])
        p2, p3 = self._phis(dl_m)
        lo = tm.transform_mode(u, so.grid, dt_s=-dt_s, phi2_s2=-p2, phi3_s3=-p3, pad=pad)  # U† u_LO
        us = tm.transform_mode(u, so.grid, dt_s=dt_s, phi2_s2=p2, phi3_s3=p3, pad=pad)
        return so.lo_mode_cov(lo), abs(tm.overlap(u, us, so.grid)) ** 2

    def evaluate(self, dt_s: float, dl_m: float) -> tuple[float, float, float, float]:
        """(V_−, V_+, θ_off, η_t) at (Δt, δL), with the fast-jitter average if enabled."""
        if self._toggles["jitter_fast"]:
            x, w = np.polynomial.hermite_e.hermegauss(int(self.params["n_jitter"]))
            w = w / w.sum()
            pts = [self._point(dt_s + float(self.params["jitter_fast_s"]) * xi, dl_m) for xi in x]
            cov = sum(wi * c for wi, (c, _) in zip(w, pts, strict=True))
            eta = float(sum(wi * e for wi, (_, e) in zip(w, pts, strict=True)))
        else:
            cov, eta = self._point(dt_s, dl_m)
        vm, vp, th = gaussian.noise_levels_from_cov(cov)
        return vm, vp, float(wrap_half_pi(th - self.source.theta_sq0_rad)), eta

    @cached_property
    def table(self) -> Table:
        dts = np.asarray(self.params["dt_grid_s"], dtype=float)
        dls = np.asarray(self.params["dL_grid_m"], dtype=float)
        vals = np.array([[self.evaluate(a, b) for b in dls] for a in dts])
        return Table(f"{self.name}.overlap", {"dt_s": dts, "dL_m": dls}, vals, "V_minus,V_plus [vac=1/2]; theta_off [rad]; eta_t [1]",
                     {"block": self.name, "columns": ["V_minus", "V_plus", "theta_off", "eta_t"],
                      "toggles": dict(self._toggles), "sources": ["T03 §2.4–§2.6", "D51-6"]})


@register_block(KIND)
def make_mode_overlap_table(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``mode_overlap_table`` block (module docstring)."""
    return ModeOverlapTable(params, name)
