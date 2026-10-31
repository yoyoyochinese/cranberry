"""``lo_temporal_shape``: temporal waveform of the LO and its overlap with the squeezed mode (T17 §6.3,
T03 §2.4, D45, D50). Owner: WP-B11.

Offline (O) block. The squeezed state comes from ``source`` (a nested temporal source spec or a
Gaussian single mode, ``fib.mode_overlap.source_from_param``); its optimal LO mode u_opt is the
reference. Shapes (``shape``, no default, D34):

* ``optimal``       u_LO = u_opt (pulse shaped by an EOM, η_t = 1)
* ``gauss``         e^{−(t − c)²/(2T0²)} (intensity ∝ e^{−(t−c)²/T0²}); params ``center_s``, ``t0_s``
* ``rect``          window [start, start + width]; params ``start_s``, ``width_s``
* ``one_sided_exp`` e^{−γ(t − d)} for t ≥ d; params ``delay_s``, ``decay_rate_rad_s`` γ (κ_c: the
                    free decay of the cavity)
* ``cw_digital``    CW LO + BHD with a single-pole bandwidth B + digital weight w = u_opt: the measured
                    mode is u_eff(s) = Σ_{j≥0} k_j w(s + jδt), k_j = (1 − e^{−2πBδt}) e^{−2πBjδt} (exact
                    for the bin-averaged first-order response; B → ∞ gives u_opt). Param ``bhd_bandwidth_hz``

Shapes are bin-averaged (32 sub-samples per bin) so that edges move continuously; η_t uses the
exact ∫|LO|² (the out-of-bin-space part of a sharp edge sees vacuum, also added to ``projection``). ``fit: true`` maximises
η_t over the shape parameters (gauss: centre and width; rect: start and width; one_sided_exp: delay
with γ fixed), starting from the moments of u_opt — the setting of the T17 §6.3 numbers, as corrected
by D60-5 / D70-2 (0.955 Gaussian, 0.866 rectangular, 0.84 one-sided exponential; T17 printed 0.951 /
0.864 from a coarser search); ``fit: false`` uses the given values.

Products: ``lo_mode``, ``eta_t`` = |∫u_LO* u_opt dT|², ``projection`` → (V_−, V_+, θ_off), ``losses()``.

η_t bookkeeping is single-source (M5-09): ``eta_t_carrier`` (D34, no default) says which product
carries η_t, and the other one is built without it, so a consumer cannot count it twice:

* ``projection``: ``projection`` = source.project(u_LO) (η_t, the mode-shape θ_off and the out-of-span
  vacuum included — exact also for a multimode source); the ``lo_shape`` element is registered
  **disabled** (η = 1).
* ``registry``: the ``lo_shape`` element (η_t, BHD only, at ``observer``) is enabled and ``projection``
  is source.project(u_opt) (η_t-free, θ_off = 0); the LossRegistry applies η_t once. Single-mode reading:
  for a single-mode source both carriers give the same V (V7), for a multimode source the ``projection``
  carrier is exact (K ≈ 1.1 for the OPO presets, T17 §5.2) and the mode-shape θ_off is lost here.

``mode_overlap_table`` follows the ``projection`` convention (its V_± include η_t; its η_t column is
informational and must not be registered as a loss as well).

Params (D34): ``source``, ``shape``, ``fit``, shape params as above (``decay_rate_rad_s`` always for
one_sided_exp), ``observer``, ``source_tag``, ``eta_t_carrier`` (projection | registry), ``toggles``
{shape_mismatch}: off ⇒ u_LO = u_opt and the element is registered disabled (η = 1).
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cranberry.blocks.fib.mode_overlap import source_from_param
from cranberry.blocks.source import temporal_modes as tm
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement
from cranberry.offline.source_output import SourceOutput

__all__ = ["LoTemporalShape", "make_lo_temporal_shape", "SHAPES", "ETA_T_CARRIERS"]

KIND = "lo_temporal_shape"
SHAPES = ("optimal", "gauss", "rect", "one_sided_exp", "cw_digital")
ETA_T_CARRIERS = ("projection", "registry")
_SHAPE_PARAMS = {"optimal": (), "gauss": ("center_s", "t0_s"), "rect": ("start_s", "width_s"),
                 "one_sided_exp": ("delay_s",), "cw_digital": ("bhd_bandwidth_hz",)}
_SUB = 32


class LoTemporalShape(BaseBlock):
    """The ``lo_temporal_shape`` block (module docstring)."""

    toggle_names = ("shape_mismatch",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "source", "shape", "fit", "observer", "source_tag", "eta_t_carrier", block=name)
        if params["shape"] not in SHAPES:
            raise ValueError(f"block {name!r}: shape must be one of {SHAPES}, got {params['shape']!r}")
        if params["eta_t_carrier"] not in ETA_T_CARRIERS:
            raise ValueError(f"block {name!r}: eta_t_carrier must be one of {ETA_T_CARRIERS}, "
                             f"got {params['eta_t_carrier']!r}")
        super().__init__(KIND, name, params)
        shape = params["shape"]
        if shape == "one_sided_exp":
            require(params, "decay_rate_rad_s", block=name)
        if shape == "cw_digital" or not params["fit"]:
            require(params, *_SHAPE_PARAMS[shape], block=name)

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    @cached_property
    def source(self) -> SourceOutput:
        so = source_from_param(self.params["source"])
        if so.basis != "temporal":
            raise ValueError(f"block {self.name!r}: source must be temporal")
        return so

    @cached_property
    def u_opt(self) -> np.ndarray:
        return tm.optimal_lo_mode(self.source)[0]

    def _binned(self, f) -> tuple[np.ndarray, float]:
        t = self.source.grid
        dt = t[1] - t[0]
        v = f(t[:, None] + ((np.arange(_SUB) + 0.5) / _SUB - 0.5)[None, :] * dt)
        return v.mean(axis=1), float(np.sum(v**2) / _SUB * dt)

    def waveform(self, shape: str, x: tuple[float, ...]) -> tuple[np.ndarray, float]:
        """(bin averages of the LO waveform on the source grid, its exact ∫|LO|² dT) for parameters ``x``."""
        dt = self.source.grid[1] - self.source.grid[0]
        if shape == "optimal":
            return self.u_opt, float(np.sum(np.abs(self.u_opt) ** 2) * dt)
        if shape == "gauss":
            c, t0 = x
            return self._binned(lambda t: np.exp(-0.5 * ((t - c) / t0) ** 2))
        if shape == "rect":
            a, w = x
            return self._binned(lambda t: ((t >= a) & (t < a + abs(w))).astype(float))
        if shape == "one_sided_exp":
            (d,), g = x, float(self.params["decay_rate_rad_s"])
            return self._binned(lambda t: np.where(t >= d, np.exp(-g * np.maximum(t - d, 0.0)), 0.0))
        a = np.exp(-2 * np.pi * float(x[0]) * dt)
        u, out, acc = self.u_opt, np.empty_like(self.u_opt), 0j
        for i in range(u.size - 1, -1, -1):  # u_eff[n] = (1 − a) Σ_j a^j u[n + j]
            acc = u[i] + a * acc
            out[i] = (1 - a) * acc
        return out, float(np.sum(np.abs(out) ** 2) * dt)

    def eta_of(self, wf: tuple[np.ndarray, float]) -> float:
        """|∫LO* u_opt dT|² / ∫|LO|² dT with the exact norm (the part of a sharp-edged LO outside the
        bin space sees vacuum; using the norm of the bin averages would overestimate η by ~δt/τ)."""
        w, nrm = wf
        return abs(tm.overlap(w, self.u_opt, self.source.grid)) ** 2 / nrm if nrm > 0 else 0.0

    @cached_property
    def shape_params(self) -> tuple[float, ...]:
        p, shape = self.params, self.params["shape"]
        if shape == "optimal":
            return ()
        if shape == "cw_digital" or not p["fit"]:
            return tuple(float(p[k]) for k in _SHAPE_PARAMS[shape])
        t, u = self.source.grid, self.u_opt
        dt = t[1] - t[0]
        w2 = np.abs(u) ** 2 * dt
        c = float(np.sum(t * w2))
        s = float(np.sqrt(np.sum((t - c) ** 2 * w2)))
        cost = lambda x: -self.eta_of(self.waveform(shape, tuple(x)))  # noqa: E731
        if shape == "one_sided_exp":
            r = minimize_scalar(lambda d: cost((d,)), bounds=(t[0], c + 3 * s), method="bounded", options={"xatol": 1e-6 * s})
            return (float(r.x),)
        x0 = (c, s) if shape == "gauss" else (c - 1.5 * s, 3 * s)
        r = minimize(cost, x0, method="Nelder-Mead", options={"xatol": 1e-6 * s, "fatol": 1e-12, "maxiter": 4000})
        return tuple(float(v) for v in r.x)

    @cached_property
    def _wf(self) -> tuple[np.ndarray, float]:
        if not self._toggles["shape_mismatch"]:
            return self.waveform("optimal", ())
        return self.waveform(self.params["shape"], self.shape_params)

    @property
    def lo_mode(self) -> np.ndarray:
        """In-span LO mode on the source grid (unit norm)."""
        w = self._wf[0].astype(complex)
        g = self.source.grid
        return w / np.sqrt(np.sum(np.abs(w) ** 2) * (g[1] - g[0]))

    @property
    def span_fraction(self) -> float:
        """∫|bin-averaged LO|² / ∫|LO|² ≤ 1: the rest of the LO measures vacuum."""
        g = self.source.grid
        return float(np.sum(np.abs(self._wf[0]) ** 2) * (g[1] - g[0]) / self._wf[1])

    @property
    def eta_t(self) -> float:
        return self.eta_of(self._wf)

    @property
    def projection(self) -> tuple[float, float, float]:
        """(V_−, V_+, θ_off). ``eta_t_carrier: projection``: seen by the LO, including the out-of-span
        vacuum (pure loss ``span_fraction``); ``registry``: seen by u_opt (η_t-free, module docstring)."""
        if self.params["eta_t_carrier"] == "registry":
            vm, vp, th = self.source.project(self.u_opt)
            return float(vm), float(vp), float(th)
        vm, vp, th = self.source.project(self.lo_mode)
        f = self.span_fraction
        return float(f * vm + 0.5 * (1 - f)), float(f * vp + 0.5 * (1 - f)), float(th)

    def losses(self) -> list[LossElement]:
        on = self._toggles["shape_mismatch"] and self.params["eta_t_carrier"] == "registry"
        return [LossElement(id=f"{self.name}.lo_shape", eta=min(self.eta_t, 1.0), kind="lo_shape", applies_to={"bhd"},
                            source=str(self.params["source_tag"]), enabled=on, observer=str(self.params["observer"]))]


@register_block(KIND)
def make_lo_temporal_shape(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``lo_temporal_shape`` block (module docstring)."""
    return LoTemporalShape(params, name)
