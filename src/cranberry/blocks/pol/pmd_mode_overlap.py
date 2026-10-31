"""``pmd_mode_overlap``: PMD temporal-mode penalty of an SMF segment for pulsed squeezed light / LO
(T09 §1.7, §3.6, §4.2; E40). Owner: WP-B7. Offline (O) block: nothing is done in the time loop.

    τ_DGD = D_PMD·√L                     (mean DGD, long fiber, random coupling; SMF-28 Ultra
                                          D_PMD ≤ 0.1 ps/√km, link design 0.04 ps/√km, UNVERIFIED)
    σ_t   = τ_p / (2√(2 ln 2))           (rms intensity width of a Gaussian pulse, τ_p = FWHM)

Worst case: the pulse is split equally between the two principal states (PSPs, group delays ±τ/2
about the mean); the analyzer is the output state at the carrier. The projected mode is
f(T) = ½[u(T − τ/2) + u(T + τ/2)] and η = |⟨u_ref|f⟩|² (polarization loss included), with
⟨u(T)|u(T − s)⟩ = exp(−s²/(8σ_t²)):

    timing_reference: psp       η = [(1 + exp(−τ²/(8σ_t²)))/2]²   (T09 §1.7 [AE]: reference mode aligned
                                                                    with one PSP; 100 fs: 4.2 % at 63 m,
                                                                    32 % at 630 m; 1 ps, 630 m: 0.44 %)
    timing_reference: centered  η = exp(−τ²/(16σ_t²))             (reference re-timed to the mean delay,
                                                                    as a timing alignment does; 2.2 % / 20 %)

``dphi_cw_pulse_rad`` is T09's [AE] estimate (τ/(2σ_t))²/8 of the CW-probe vs pulse-averaged phase
difference (0 for the symmetric worst-case split above; it needs an asymmetric U(Ω)).
``pulse_projection`` evaluates any frequency-dependent Jones matrix U(Ω) on a pulse (conv §8 Fourier
convention) and is the numerical check of both η forms (and of T09 test 7, PM two-point crosstalk).
The loss is registered on ``segment`` with kind ``fiber`` (not ``mode_overlap_temporal``, which the
timing / dispersion overlap of the same path pair owns, G06).

ponytail: only the worst equal split and the √L (long-fiber) DGD; add the PSP split fraction or the
short-fiber limit B_g·L/c when a measured U(Ω) or a short SMF needs it.

Toggle ``pol_pmd``: off ⇒ η = 1 (element kept in the record, disabled).

Params (all required, D34): ``D_PMD_s_per_sqrt_m``, ``L_m``, ``tau_p_fwhm_s``, ``timing_reference``
(psp | centered), ``segment``, ``applies_to``, ``source``, ``toggles``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import numpy as np

from cranberry.blocks.loss.per_budget import _noop
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["PmdModeOverlap", "make_pmd_mode_overlap", "pulse_projection", "TIMING_REFERENCES"]

KIND = "pmd_mode_overlap"
TIMING_REFERENCES: tuple[str, ...] = ("psp", "centered")
_FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def pulse_projection(u: np.ndarray, dt: float, U_of_omega: Callable[[np.ndarray], np.ndarray],  # noqa: N803
                     e_in, g) -> np.ndarray:
    """Projected output mode f(T) = g^† · F⁻¹[U(Ω) e_in Ã(Ω)] of the input mode u(T) (grid step dt).

    conv §8: Ã(Ω_k) = N·dt·ifft(u)[k], Ω_k = 2π·fftfreq(N, dt), A = fft(Ã)/(N·dt); a delay τ is
    Ã·e^{+iΩτ}. ``U_of_omega(Ω)`` returns shape (N, 2, 2). ⟨u|f⟩ = Σ conj(u)·f·dt.
    """
    n = len(u)
    om = 2.0 * np.pi * np.fft.fftfreq(n, dt)
    a = n * dt * np.fft.ifft(u)
    v = np.einsum("kij,j->ki", U_of_omega(om), np.asarray(e_in, dtype=complex))
    ft = (np.conj(np.asarray(g, dtype=complex)) @ v.T) * a
    return np.fft.fft(ft) / (n * dt)


class PmdModeOverlap(BaseBlock):
    """The ``pmd_mode_overlap`` block (module docstring)."""

    toggle_names = ("pol_pmd",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "D_PMD_s_per_sqrt_m", "L_m", "tau_p_fwhm_s", "timing_reference", "segment", "applies_to",
                "source", block=name)
        super().__init__(KIND, name, params)
        if params["timing_reference"] not in TIMING_REFERENCES:
            raise ValueError(f"block {name!r}: timing_reference must be one of {TIMING_REFERENCES}")
        if min(float(params["L_m"]), float(params["tau_p_fwhm_s"])) <= 0 or float(params["D_PMD_s_per_sqrt_m"]) < 0:
            raise ValueError(f"block {name!r}: need L_m > 0, tau_p_fwhm_s > 0, D_PMD_s_per_sqrt_m ≥ 0")
        self.tau_dgd_s = float(params["D_PMD_s_per_sqrt_m"]) * math.sqrt(float(params["L_m"]))
        self.sigma_t_s = float(params["tau_p_fwhm_s"]) * _FWHM_TO_SIGMA
        x = self.tau_dgd_s ** 2 / (8.0 * self.sigma_t_s ** 2)
        self._eta = ((1.0 + math.exp(-x)) / 2.0) ** 2 if params["timing_reference"] == "psp" else math.exp(-x / 2.0)

    @property
    def eta(self) -> float:
        """Worst-case PMD mode-overlap efficiency (1 with pol_pmd off)."""
        return self._eta if self._toggles["pol_pmd"] else 1.0

    @property
    def dphi_cw_pulse_rad(self) -> float:
        """T09 [AE] estimate (τ_DGD/(2σ_t))²/8 of the CW vs pulse phase difference (0 when off)."""
        return (self.tau_dgd_s / (2.0 * self.sigma_t_s)) ** 2 / 8.0 if self._toggles["pol_pmd"] else 0.0

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def losses(self) -> list[LossElement]:
        p = self.params
        return [LossElement(id=f"{self.name}.pmd", eta=self._eta, kind="fiber", applies_to=set(p["applies_to"]),
                            source=str(p["source"]), enabled=self._toggles["pol_pmd"], segment=str(p["segment"]))]


@register_block(KIND)
def make_pmd_mode_overlap(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``pmd_mode_overlap`` block (module docstring)."""
    return PmdModeOverlap(params, name)
