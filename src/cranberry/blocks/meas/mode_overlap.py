"""``mode_overlap``: spatial / temporal / LO-shape mode-matching efficiency as one loss element (G06). Owner: WP-B5.

Offline. Mode mismatch with a single-mode source is a pure loss η = |⟨u_s|u_LO⟩|² (T12 §1.3, T03
§2.4). One block instance registers one ``LossElement`` of kind ``mode_overlap_spatial`` /
``mode_overlap_temporal`` / ``lo_shape`` at its ``segment`` or ``observer``; the ``LossRegistry``
rejects a second element of the same kind on the same path pair (VT-F3), so each overlap is counted
once. Models (``model``):

* spatial ``gaussian_beam`` (T12 §1.3, waists at the same position):
  η = [2w₁w₂/(w₁² + w₂²)]² · exp[−2d²/(w₁² + w₂²)] · exp[−(πnθ)²·2w₁²w₂²/(λ²(w₁² + w₂²))]
* spatial ``visibility``: η = V²(P₁ + P₂)²/(4P₁P₂) (T12 §1.3)
* temporal ``gaussian_pulse`` (T03 §2.5, squeezed delayed by Δt with GDD mismatch φ₂, T₀ = 1/e
  intensity half-width): κ = (1 − iχ)^{−1/2} exp[−Δt²(1 + iχ)/(4T₀²(1 + χ²))], χ = φ₂/(2T₀²),
  η = |κ|²; arg κ is the angle offset θ_off (``theta_off_rad``; not a loss — for the overlap tables)
* any kind ``eta``: η given (e.g. LO-shape overlaps of T17: Gaussian 0.951, rectangle 0.864, with
  the pulsed-OPO optimal mode of WP-B11)

Numerical helpers for sampled modes (conv §8 grids): :func:`temporal_overlap` κ = ∫u₁*u₂ dT and
:func:`best_shape_overlap` (LO shape ``gauss`` / ``rect`` / ``one_sided_exp`` optimized against a
sampled mode, T17 §6.3).

Toggle ``enabled`` (MODE_MISMATCH_ENABLE, T12 §8.1): off ⇒ the element is recorded disabled (η = 1).
Params (all required, D34): overlap (spatial | temporal | lo_shape), model, model parameters,
segment, observer (exactly one non-null), applies_to, source, toggles.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
import scipy.optimize

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["ModeOverlap", "make_mode_overlap", "OVERLAP_KIND", "MODELS", "gaussian_beam_overlap",
           "overlap_from_visibility", "visibility_from_overlap", "gaussian_pulse_kappa", "temporal_overlap",
           "lo_shape", "best_shape_overlap", "LO_SHAPES"]

#: block ``overlap`` → LossElement kind
OVERLAP_KIND: dict[str, str] = {"spatial": "mode_overlap_spatial", "temporal": "mode_overlap_temporal",
                                "lo_shape": "lo_shape"}
MODELS: dict[str, tuple[str, ...]] = {"spatial": ("gaussian_beam", "visibility", "eta"),
                                      "temporal": ("gaussian_pulse", "eta"), "lo_shape": ("eta",)}
LO_SHAPES: tuple[str, ...] = ("gauss", "rect", "one_sided_exp")


def gaussian_beam_overlap(w1_m: float, w2_m: float, d_m: float, theta_rad: float, wavelength_m: float,
                          n: float) -> float:
    """η_mm of two Gaussian beams (offset d, tilt θ, medium index n; T12 §1.3)."""
    s = w1_m**2 + w2_m**2
    return ((2.0 * w1_m * w2_m / s) ** 2 * math.exp(-2.0 * d_m**2 / s)
            * math.exp(-((math.pi * n * theta_rad) ** 2) * 2.0 * w1_m**2 * w2_m**2 / (wavelength_m**2 * s)))


def overlap_from_visibility(v_vis: float, p1_w: float, p2_w: float) -> float:
    """η_mm = V²(P₁ + P₂)²/(4P₁P₂) (T12 §1.3)."""
    return v_vis**2 * (p1_w + p2_w) ** 2 / (4.0 * p1_w * p2_w)


def visibility_from_overlap(eta_mm: float, p1_w: float, p2_w: float) -> float:
    """Fringe visibility V = 2√(P₁P₂)/(P₁ + P₂)·√η_mm."""
    return 2.0 * math.sqrt(p1_w * p2_w) / (p1_w + p2_w) * math.sqrt(eta_mm)


def gaussian_pulse_kappa(dt_s: float, phi2_s2: float, t0_s: float) -> complex:
    """κ of identical Gaussian pulses with delay Δt and GDD mismatch φ₂ (T03 §2.5)."""
    chi = phi2_s2 / (2.0 * t0_s**2)
    return (1.0 - 1j * chi) ** -0.5 * cmath.exp(-(dt_s**2) * (1.0 + 1j * chi) / (4.0 * t0_s**2 * (1.0 + chi**2)))


def temporal_overlap(u1: np.ndarray, u2: np.ndarray, dt_s: float) -> complex:
    """κ = ∫u₁*(T) u₂(T) dT on a uniform grid (modes normalized here, ∫|u|² dT = 1)."""
    u1, u2 = np.asarray(u1, dtype=complex), np.asarray(u2, dtype=complex)
    n1, n2 = np.sum(np.abs(u1) ** 2) * dt_s, np.sum(np.abs(u2) ** 2) * dt_s
    return complex(np.sum(np.conj(u1) * u2) * dt_s / math.sqrt(n1 * n2))


def lo_shape(shape: str, t_s: np.ndarray, a: float, b: float) -> np.ndarray:
    """LO amplitude u(t) (unnormalized): gauss (centre a, intensity FWHM b), rect (start a, width b),
    one_sided_exp (start a, field decay rate b [1/s])."""
    t = np.asarray(t_s, dtype=float)
    if shape == "gauss":
        return np.exp(-2.0 * math.log(2.0) * ((t - a) / b) ** 2)
    if shape == "rect":
        return ((t >= a) & (t < a + b)).astype(float)
    if shape == "one_sided_exp":
        return np.where(t >= a, np.exp(-b * np.clip(t - a, 0.0, None)), 0.0)
    raise ValueError(f"shape must be one of {LO_SHAPES}")


def best_shape_overlap(u_mode: np.ndarray, t_s: np.ndarray, shape: str, x0: tuple[float, float]) -> dict[str, Any]:
    """Maximize |κ|² of a real LO ``shape`` against the sampled mode ``u_mode`` (T17 §6.3); x0 = (a, b)."""
    dt = float(t_s[1] - t_s[0])
    scale = np.array([abs(x0[0]) or dt, abs(x0[1]) or dt])

    def neg(x: np.ndarray) -> float:
        u = lo_shape(shape, t_s, *(x * scale))
        return 1.0 if not np.any(u) else -abs(temporal_overlap(u_mode, u, dt)) ** 2

    res = scipy.optimize.minimize(neg, np.asarray(x0, dtype=float) / scale, method="Nelder-Mead",
                                  options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 4000})
    a, b = res.x * scale
    return {"eta": -float(res.fun), "a": float(a), "b": float(b)}


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class ModeOverlap(BaseBlock):
    """``mode_overlap`` (module docstring)."""

    KIND = "mode_overlap"
    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "overlap", "model", "segment", "observer", "applies_to", "source", "toggles", block=name)
        super().__init__(self.KIND, name, params)
        p = params
        if p["overlap"] not in OVERLAP_KIND:
            raise ValueError(f"block {name!r}: overlap must be one of {tuple(OVERLAP_KIND)}")
        if p["model"] not in MODELS[p["overlap"]]:
            raise ValueError(f"block {name!r}: model for {p['overlap']} must be one of {MODELS[p['overlap']]}")
        if (p["segment"] is None) == (p["observer"] is None):
            raise ValueError(f"block {name!r}: give exactly one of segment / observer")
        self.theta_off_rad = 0.0
        m = p["model"]
        if m == "eta":
            require(p, "eta", block=name)
            self.eta = float(p["eta"])
        elif m == "gaussian_beam":
            require(p, "w1_m", "w2_m", "d_m", "theta_rad", "wavelength_m", "n", block=name)
            self.eta = gaussian_beam_overlap(*(float(p[k]) for k in ("w1_m", "w2_m", "d_m", "theta_rad",
                                                                    "wavelength_m", "n")))
        elif m == "visibility":
            require(p, "v_vis", "p1_w", "p2_w", block=name)
            self.eta = overlap_from_visibility(float(p["v_vis"]), float(p["p1_w"]), float(p["p2_w"]))
        else:
            require(p, "dt_s", "phi2_s2", "t0_s", block=name)
            k = gaussian_pulse_kappa(float(p["dt_s"]), float(p["phi2_s2"]), float(p["t0_s"]))
            self.eta, self.theta_off_rad = abs(k) ** 2, cmath.phase(k)
        if not 0.0 <= self.eta <= 1.0 + 1e-12:
            raise ValueError(f"block {name!r}: overlap {self.eta} not in [0, 1]")
        self.eta = min(self.eta, 1.0)

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def losses(self) -> list[LossElement]:
        p = self.params
        return [LossElement(id=f"{self.name}.overlap", eta=self.eta, kind=OVERLAP_KIND[p["overlap"]],
                            applies_to=set(p["applies_to"]), source=str(p["source"]), enabled=self._toggles["enabled"],
                            segment=None if p["segment"] is None else str(p["segment"]),
                            observer=None if p["observer"] is None else str(p["observer"]))]


@register_block(ModeOverlap.KIND)
def make_mode_overlap(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``mode_overlap`` (module docstring)."""
    return ModeOverlap(params, name)
