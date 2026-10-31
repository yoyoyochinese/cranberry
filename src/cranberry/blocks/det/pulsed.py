"""``pulsed_bhd_kernel``: per-pulse BHD estimator coefficients (T13 §1.8, §4.3; T11 §3.7). Owner: WP-B5.

Offline. Pulse j's estimate X_j = ∫ w(t − jT)(h ∗ i_−)(t) dt; the charge of pulse j − m enters with

    c_m = ∫_0^W h(t + mT) dt / ∫_0^W h(t) dt,   m = 0..n_taps,  T = 1/f_rep            (c_0 = 1)

for a boxcar window w of length W ≤ T starting at the pulse (causal h ⇒ c_m = 0 for m < 0).
Responses: ``first_order`` (bandwidth B, h ∝ e^{−2πBt}: c_m = e^{−2πBm/f_rep} exactly, any W) or ``second_order`` (f_0, Q; ringing gives sign-changing c_m), evaluated with
``scipy.signal.step``. Adjacent-pulse correlation of a uniform independent train
ρ_adj = Σ c_m c_{m+1} / Σ c_m² (first order: c_1). Per-pulse clearance C_pulse = η_q N_LO / σ_el²
(σ_el in electrons rms); its equivalent loss C/(1 + C) is registered by ``detector_efficiency``
(``clearance_lin`` = C_pulse), not here, so it is counted once.

Pulse shape (M5-11): the formula treats each pulse's charge as a δ function at t = jT (the window
starts at the pulse). That is exact for ps pulses but **underestimates c_m for ns modes**: a pulse-pumped
OPO mode (Borealis type, τ_eff 14.5 ns, centroid 6.9 ns after the pump, tail ~20 ns) at f_rep ≈ 45 MHz
(T = 22 ns) has a duration comparable to T, so charge arriving late in the slot leaks more into the next
estimate (at 6 MHz the effect is negligible). Optional ``pulse_charge_profile`` {t_s, weight} (the
charge arrival profile of one pulse relative to its slot start, e.g. t_s = the source grid and weight =
|u_LO(t)|² of the measured mode) replaces the δ by that profile:

    c_m ∝ Σ_k q_k [τ(mT − s_k) − τ(mT − s_k + W)],   τ(x) = 1 − s(x) (= 1 for x ≤ 0),

normalised to c_0 = 1 (absent or null ⇒ δ pulses, the formula above; same "check-style" optional key as
``mode_length_s``).

Toggles: ``pulsed_crosstalk`` (off ⇒ c = [1]), ``circuit_noise`` (off ⇒ C_pulse = inf).

Consistency check (D46/D50, CHK — never stops the run): if ``mode_length_s`` is given (not null) and
``window_s`` < ``mode_length_s``, a ``gate_shorter_than_mode`` warning goes to ``chk_log``. It is a
check, not a modelling choice, so leaving the key out simply skips it (D34 concerns modelling choices).

Params (all required, D34): f_rep_hz, response {model: first_order, bandwidth_hz} |
{model: second_order, f0_hz, q}, window_s, n_taps, sigma_el_e, n_lo_photons, qe, toggles. Optional:
mode_length_s (check), pulse_charge_profile (ns modes, above).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
import scipy.signal

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block

__all__ = ["PulsedBhdKernel", "make_pulsed_bhd_kernel", "RESPONSE_MODELS", "c1_first_order", "rho_adjacent",
           "crosstalk_coefficients"]

RESPONSE_MODELS: tuple[str, ...] = ("first_order", "second_order")


def c1_first_order(bandwidth_hz: float, f_rep_hz: float) -> float:
    """c_1 = exp(−2πB/f_rep) for a first-order response and a short window (T13 §1.8)."""
    return math.exp(-2.0 * math.pi * bandwidth_hz / f_rep_hz)


def rho_adjacent(c: Any) -> float:
    """Apparent adjacent-pulse correlation of an independent uniform train: Σ c_m c_{m+1} / Σ c_m²."""
    c = np.asarray(c, dtype=float)
    return float(np.dot(c[:-1], c[1:]) / np.dot(c, c))


def _tail(response: Mapping[str, Any], t: np.ndarray) -> np.ndarray:
    """1 − s(t), s = unit-DC-gain step response (integral of h), at increasing times t ≥ 0."""
    if response["model"] == "first_order":
        return np.exp(-2.0 * math.pi * float(response["bandwidth_hz"]) * t)  # exact, no cancellation
    w0 = 2.0 * math.pi * float(response["f0_hz"])
    sys = scipy.signal.lti([w0 * w0], [1.0, w0 / float(response["q"]), w0 * w0])
    return 1.0 - scipy.signal.step(sys, T=t)[1]


def crosstalk_coefficients(response: Mapping[str, Any], f_rep_hz: float, window_s: float, n_taps: int,
                           profile: Mapping[str, Any] | None = None) -> np.ndarray:
    """c_m, m = 0..n_taps (module docstring; ``profile`` = ``pulse_charge_profile`` or None for δ pulses)."""
    if response.get("model") not in RESPONSE_MODELS:
        raise ValueError(f"response.model must be one of {RESPONSE_MODELS}")
    period = 1.0 / f_rep_hz
    if not 0.0 < window_s <= period:
        raise ValueError(f"window_s must be in (0, 1/f_rep], got {window_s}")
    m = np.arange(n_taps + 1)
    sk = np.zeros(1) if profile is None else np.asarray(profile["t_s"], dtype=float)
    qk = np.ones(1) if profile is None else np.asarray(profile["weight"], dtype=float)
    if sk.shape != qk.shape or sk.ndim != 1 or np.any(qk < 0) or not qk.sum() > 0:
        raise ValueError("pulse_charge_profile needs 1-D t_s and non-negative weight of the same length")
    lag = m[:, None] * period - sk[None, :]
    ends = np.concatenate((lag, lag + window_s)).ravel()
    # one increasing grid (scipy.signal.step needs increasing times ≥ 0); τ(x ≤ 0) = 1 (causal h)
    t, inv = np.unique(np.maximum(ends, 0.0), return_inverse=True)
    tail = np.where(ends <= 0.0, 1.0, _tail(response, t)[inv]).reshape(2, m.size, sk.size)
    num = (tail[0] - tail[1]) @ qk
    return num / num[0]


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class PulsedBhdKernel(BaseBlock):
    """``pulsed_bhd_kernel`` (module docstring)."""

    KIND = "pulsed_bhd_kernel"
    toggle_names = ("pulsed_crosstalk", "circuit_noise")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "f_rep_hz", "response", "window_s", "n_taps", "sigma_el_e", "n_lo_photons", "qe",
                "toggles", block=name)
        super().__init__(self.KIND, name, params)
        p = params
        self._c = crosstalk_coefficients(p["response"], float(p["f_rep_hz"]), float(p["window_s"]), int(p["n_taps"]),
                                         p.get("pulse_charge_profile"))
        self.chk_log = ChkLog()
        mode_len = p.get("mode_length_s")
        if mode_len is not None and float(p["window_s"]) < float(mode_len):
            self.chk_log.warn("gate_shorter_than_mode", name,
                              f"gate {p['window_s']} s < mode length {mode_len} s (D46/D50)",
                              window_s=float(p["window_s"]), mode_length_s=float(mode_len))

    @property
    def c(self) -> np.ndarray:
        """c_m, m = 0..n_taps ([1] with ``pulsed_crosstalk`` off)."""
        return self._c.copy() if self._toggles["pulsed_crosstalk"] else np.ones(1)

    @property
    def clearance_pulse(self) -> float:
        """C_pulse = η_q N_LO / σ_el² (inf with ``circuit_noise`` off or σ_el = 0)."""
        s = float(self.params["sigma_el_e"])
        if not self._toggles["circuit_noise"] or s == 0.0:
            return math.inf
        return float(self.params["qe"]) * float(self.params["n_lo_photons"]) / (s * s)

    @property
    def rho_adj(self) -> float:
        c = self.c
        return rho_adjacent(c) if c.size > 1 else 0.0

    def table(self) -> dict[str, Any]:
        """Offline table handed to ``pulse_sample`` / the record."""
        return {"c": self.c.tolist(), "clearance_pulse": self.clearance_pulse, "rho_adj": self.rho_adj}

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state


@register_block(PulsedBhdKernel.KIND)
def make_pulsed_bhd_kernel(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``pulsed_bhd_kernel`` (module docstring)."""
    return PulsedBhdKernel(params, name)
