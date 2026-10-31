"""``spectral_filter``: passive spectral filter on the squeezed beam, e.g. the single-mode filter before
the GBS (D44, T03 §5.4, T17 §2.7). Owner: WP-B11.

Offline (O) block. Amplitude transmission t(Ω) (zero phase for ``gauss`` / ``rect``, causal single
pole for ``lorentz``; B = FWHM of |t|² [Hz], f = Ω/2π − ``center_hz``)::

    gauss    t = exp(−2 ln2 (f/B)²)
    rect     t = 1 for |f| ≤ B/2, else 0
    lorentz  t = γ/(γ − iΩ'),  γ = πB          (h(T) ∝ e^{−γT}θ(T) with the conv §8 sign)

:meth:`SpectralFilter.apply` maps a ``SourceOutput`` through it: â′ = A â + vacuum with
A_mn = ∫ f_m* (h ∗ f_n) dT (temporal basis, FFT on a zero-padded grid; components leaving
span{f_n} are lost, so A is a contraction) and real form R = [[Re A, −Im A], [Im A, Re A]]:
cov′ = R cov Rᵀ + (𝟙 − R Rᵀ)/2 (conv §2). ``cw_sideband``: the ±f_s pair sees the pure loss |t(f_s)|²
(centred filters only). ``single_mode`` is a CW mode at f = 0: unchanged (|t(0)|² = 1).

The flat insertion loss ``insertion_eta`` is one registered ``component`` element on ``segment`` (G06);
the shape-dependent part is what :meth:`apply` returns (it is not registered, so it is not counted
twice).

Params (D34): ``shape`` (gauss | rect | lorentz), ``bandwidth_hz``, ``center_hz``, ``insertion_eta``,
``segment``, ``applies_to``, ``source_tag``, ``pad``, ``toggles`` {enabled}: off ⇒ ``apply`` is the
identity and the element is registered disabled.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.core import fft
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement
from cranberry.offline.source_output import SourceOutput

__all__ = ["SpectralFilter", "make_spectral_filter", "SHAPES"]

KIND = "spectral_filter"
SHAPES = ("gauss", "rect", "lorentz")


class SpectralFilter(BaseBlock):
    """The ``spectral_filter`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "shape", "bandwidth_hz", "center_hz", "insertion_eta", "segment", "applies_to", "source_tag", "pad",
                block=name)
        if params["shape"] not in SHAPES:
            raise ValueError(f"block {name!r}: shape must be one of {SHAPES}, got {params['shape']!r}")
        if not float(params["bandwidth_hz"]) > 0 or not 0.0 <= float(params["insertion_eta"]) <= 1.0:
            raise ValueError(f"block {name!r}: need bandwidth_hz > 0 and 0 <= insertion_eta <= 1")
        super().__init__(KIND, name, params)

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def t_of(self, omega) -> np.ndarray:
        """Amplitude transmission t(Ω) [Ω in rad/s, sideband frequency, conv §1]."""
        p = self.params
        b = float(p["bandwidth_hz"])
        w = np.asarray(omega, dtype=float) - 2 * np.pi * float(p["center_hz"])
        f = w / (2 * np.pi)
        if p["shape"] == "gauss":
            return np.exp(-2 * np.log(2) * (f / b) ** 2).astype(complex)
        if p["shape"] == "rect":
            return (np.abs(f) <= 0.5 * b).astype(complex)
        g = np.pi * b
        return g / (g - 1j * w)

    def apply(self, so: SourceOutput) -> SourceOutput:
        """Filtered ``SourceOutput`` (module docstring)."""
        if not self._toggles["enabled"] or so.basis == "single_mode":
            return so
        if so.basis == "cw_sideband":
            if float(self.params["center_hz"]) != 0.0:
                raise ValueError(f"block {self.name!r}: cw_sideband needs a centred filter (center_hz = 0)")
            eta = (np.abs(self.t_of(2 * np.pi * so.grid)) ** 2)[:, None, None]
            return _replace(so, eta * so.cov + (1.0 - eta) * 0.5 * np.eye(2))  # pure loss per sideband (conv §2)
        g, modes = so.grid, so.modes
        dt = g[1] - g[0]
        n, m = g.size, g.size * int(self.params["pad"])
        a = np.zeros((modes.shape[0], n + 2 * m), dtype=complex)
        a[:, m : m + n] = modes
        t0 = float(g[0]) - m * dt
        om = fft.omega_grid(a.shape[1], dt)
        filt = fft.ifft_agrawal(fft.fft_agrawal(a, dt, t0=t0) * self.t_of(om), dt, t0=t0)[:, m : m + n]
        amat = (modes.conj() * np.gradient(g)) @ filt.T  # A_mn = ∫ f_m* (h ∗ f_n)
        r = np.block([[amat.real, -amat.imag], [amat.imag, amat.real]])
        cov = r @ so.cov @ r.T + 0.5 * (np.eye(r.shape[0]) - r @ r.T)
        return _replace(so, 0.5 * (cov + cov.T))

    def losses(self) -> list[LossElement]:
        return [LossElement(id=f"{self.name}.il", eta=float(self.params["insertion_eta"]), kind="component",
                            applies_to=set(self.params["applies_to"]), source=str(self.params["source_tag"]),
                            enabled=self._toggles["enabled"], segment=str(self.params["segment"]))]


def _replace(so: SourceOutput, cov: np.ndarray) -> SourceOutput:
    return SourceOutput(basis=so.basis, grid=so.grid, cov=cov, phi_opa0_rad=so.phi_opa0_rad, modes=so.modes,
                        pump_path_phase_rad=so.pump_path_phase_rad, gain_table=so.gain_table, coef_tables=so.coef_tables,
                        tail_into_next_slot=so.tail_into_next_slot, provenance={**so.provenance, "spectral_filter": True})


@register_block(KIND)
def make_spectral_filter(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``spectral_filter`` block (module docstring)."""
    return SpectralFilter(params, name)
