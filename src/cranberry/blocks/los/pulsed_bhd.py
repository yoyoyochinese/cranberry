"""``pulsed_bhd_kernel``: pulse-resolved / gated / integrating balanced homodyne readout (T13 §1.8,
T11 §3.7, N10, D45, D46). Owner: WP-B11 (stand-in: WP-B5's DET family had no pulsed kernel on dev).

Offline (O) block. Per-pulse estimate X_j = ∫w(t − jT_rep)(h ∗ i_−)(t) dt (T13 §1.8): the charge of
pulse m enters estimate j with

    c_m = ∫ w(t) h(t + m T_rep) dt / ∫ w(t) h(t) dt          (m ≥ 0, c_0 = 1; h causal)

w = rectangular gate [0, ``gate_s``] from the pulse arrival; h = first-order low pass (bandwidth B,
c_m = e^{−2πB m/f_rep} for any gate) or second order (B = f_0, ``q``; ringing → sign changes, by
quadrature). Adjacent-slot correlation of an i.i.d. pulse train ρ_adj = Σ c_m c_{m+1}/Σ c_m²; the
normalised per-pulse variance is unchanged when every slot holds the same state (T13 §1.8).

Electronic noise (T13 §1.3): clearance C = η_q N_LO/σ_el² (σ_el in electrons rms per gate), registered
as one ``electronic`` loss element η_el = C/(1 + C) (BHD only, at ``observer``).

Readout modes (``mode``, D45):

* ``pulse_resolved`` / ``gated``: R_meas = R (with the crosstalk numbers above); ``gated`` also checks
  D46/D50 — the gate must be at least the mode length (``mode_length_s``), else a CHK warning
  ``gate_shorter_than_mode`` (``chk_log``; the computation continues, D34).
* ``integrating``: the LO covers K = ``n_slots_integrated`` slots and only one is squeezed (N10):
  R_meas = (R + K − 1)/K (K = 32, −10 dB → −0.12 dB, T11 §7.7).

Params (D34): ``response`` (first_order | second_order), ``bandwidth_hz``, ``q`` (second_order),
``gate_s``, ``rep_rate_hz``, ``m_max``, ``n_lo_photons``, ``sigma_el_electrons``,
``quantum_efficiency`` (enters C only; the QE loss itself is ``detector_efficiency``'s), ``mode``,
``n_slots_integrated`` (integrating), ``mode_length_s`` (gated), ``observer``, ``source_tag``,
``toggles`` {crosstalk (off → c_m = δ_m0), electronic_noise (off → element disabled)}.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np
from scipy.integrate import quad

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["PulsedBhdKernel", "make_pulsed_bhd_kernel", "MODES"]

KIND = "pulsed_bhd_kernel"
MODES = ("pulse_resolved", "gated", "integrating")
RESPONSES = ("first_order", "second_order")


class PulsedBhdKernel(BaseBlock):
    """The ``pulsed_bhd_kernel`` block (module docstring)."""

    toggle_names = ("crosstalk", "electronic_noise")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "response", "bandwidth_hz", "gate_s", "rep_rate_hz", "m_max", "n_lo_photons", "sigma_el_electrons",
                "quantum_efficiency", "mode", "observer", "source_tag", block=name)
        if params["response"] not in RESPONSES or params["mode"] not in MODES:
            raise ValueError(f"block {name!r}: response in {RESPONSES}, mode in {MODES}")
        super().__init__(KIND, name, params)
        if params["response"] == "second_order":
            require(params, "q", block=name)
        require(params, *{"integrating": ("n_slots_integrated",), "gated": ("mode_length_s",)}.get(params["mode"], ()),
                block=name)
        if not (float(params["bandwidth_hz"]) > 0 and float(params["gate_s"]) > 0 and float(params["rep_rate_hz"]) > 0):
            raise ValueError(f"block {name!r}: bandwidth_hz, gate_s, rep_rate_hz must be > 0")
        self.chk_log = ChkLog()
        if params["mode"] == "gated" and float(params["gate_s"]) < float(params["mode_length_s"]):
            self.chk_log.warn("gate_shorter_than_mode", name,
                              f"gate {params['gate_s']} s < mode length {params['mode_length_s']} s (D46/D50)",
                              gate_s=float(params["gate_s"]), mode_length_s=float(params["mode_length_s"]))

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _h(self, t: float) -> float:
        p = self.params
        w0 = 2 * np.pi * float(p["bandwidth_hz"])
        if t < 0:
            return 0.0
        if p["response"] == "first_order":
            return w0 * np.exp(-w0 * t)
        z = 0.5 / float(p["q"])
        if z < 1.0:
            wd = w0 * np.sqrt(1 - z * z)
            return w0 / np.sqrt(1 - z * z) * np.exp(-z * w0 * t) * np.sin(wd * t)
        if z == 1.0:
            return w0 * w0 * t * np.exp(-w0 * t)
        s = w0 * np.sqrt(z * z - 1)
        return w0 / (2 * np.sqrt(z * z - 1)) * (np.exp((-z * w0 + s) * t) - np.exp((-z * w0 - s) * t))

    @cached_property
    def c_m(self) -> np.ndarray:
        """Crosstalk coefficients c_0..c_{m_max} (module docstring)."""
        p = self.params
        m = np.arange(int(p["m_max"]) + 1)
        if not self._toggles["crosstalk"]:
            return (m == 0).astype(float)
        trep, g = 1.0 / float(p["rep_rate_hz"]), float(p["gate_s"])
        if p["response"] == "first_order":
            return np.exp(-2 * np.pi * float(p["bandwidth_hz"]) * m * trep)
        w = np.array([quad(lambda t, k=k: self._h(t + k * trep), 0.0, g, limit=400)[0] for k in m])
        return w / w[0]

    @property
    def adjacent_correlation(self) -> float:
        c = self.c_m
        return float(np.sum(c[:-1] * c[1:]) / np.sum(c**2))

    @property
    def clearance(self) -> float:
        p = self.params
        return float(p["quantum_efficiency"]) * float(p["n_lo_photons"]) / float(p["sigma_el_electrons"]) ** 2

    def measured_r(self, r):
        """Measured normalised variance R_meas from the per-pulse R (module docstring; no electronic noise,
        which is the registered element)."""
        if self.params["mode"] == "integrating":
            k = int(self.params["n_slots_integrated"])
            return (np.asarray(r, dtype=float) + k - 1) / k
        return np.asarray(r, dtype=float)

    def losses(self) -> list[LossElement]:
        c = self.clearance
        return [LossElement(id=f"{self.name}.electronic", eta=c / (1 + c), kind="electronic", applies_to={"bhd"},
                            source=str(self.params["source_tag"]), enabled=self._toggles["electronic_noise"],
                            observer=str(self.params["observer"]))]


@register_block(KIND)
def make_pulsed_bhd_kernel(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``pulsed_bhd_kernel`` block (module docstring)."""
    return PulsedBhdKernel(params, name)
