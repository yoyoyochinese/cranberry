"""Balanced (homodyne) detector models: ``detector_efficiency`` (O) and ``bhd_frontend`` (M). Owner: WP-B5.

Notation and equations: T13 (notes/bhd-detection.md) §1, conv §2, §5. Normalized variances R are in
shot-noise units (vacuum R = 1 = V/V_SN, conv §2).

Frequency response of the difference channel (T13 §1.4), normalized to 1 at DC::

    H(f)     = H_TIA(f) · H_dip(f)
    H_TIA(f) = 1/(1 + i f/f_3dB)                                (tia.model = first_order)
             = 1/(1 − (f/f_0)² + i f/(Q f_0))                   (tia.model = second_order)
    H_dip(f) = (1 + r e^{−i2πfτ})/(1 + r)                       min at f = 1/(2τ), depth (1 − r)/(1 + r)

Circuit (electronic) noise S_el is input-referred and added after H (T13 §2), so the clearance is
C(f) = |H(f)|² C_0 with C_0 = S_SN/S_el = 2eℜP_LO/i_n² (or C_0 = C_ref · P_LO/P_ref, C ∝ P_LO). A
measurement normalized to its own LO-only trace is a pure loss η_e = C/(1 + C) (conv §5, T13 §1.3).
Imbalance (T13 §1.5): x(f) = |δ|²(1 + N_ex), N_ex = η_q P_LO S_RIN/(2hν), |δ|² = 1/(4·CMRR_lin) from
the Thorlabs CMRR definition (|H_1 − H_2|/|H_1| = 10^(−CMRR/20)); a pure loss η_bal = 1/(1 + x).
QE is a pure loss η_q = ℜhν/e (T13 §1.2). η_det = η_q η_e η_bal multiplies the optical η.

``detector_efficiency`` (offline): registers these three as ``LossElement`` records at the BHD
observer (G06; kinds ``detector`` / ``electronic`` / ``detector``), evaluated at the measurement
frequency f_s. It is the only place the detector enters the metric chain (``LossRegistry`` →
``quadrature_for_metric``). Electronic noise has no ``bands()`` (D56-5).

``bhd_frontend`` (measurement stage): the raw electrical outputs for the ESA emulation (T04 §5) —
signal / LO-only / dark PSDs of the RF output, the DC voltage with offset and clipping, the MONITOR
voltage. **Feed it R before the detector's electronic and imbalance elements** (R = optical η incl.
QE); otherwise they are counted twice. Coupler sign s (conv §5): v_DC = s·G ℜ (P_1 − P_2) + V_off.
Load: ``50ohm`` halves G and the swing (PDB manual p.13, p.14 note 1).

Toggles (protected list; off ⇒ identity, i.e. the element is recorded disabled / the term is 0):
``qe``, ``circuit_noise``, ``imbalance_leak``, ``freq_response``, ``dip``, ``subtract_dark``
(on ⇒ dark-noise-corrected: η_e treated as 1 / dark PSD subtracted), ``saturation`` (frontend only).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
import scipy.constants as sc

from cranberry.core.block import BaseBlock, Block, MissingParameterError, require
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = [
    "DetectorEfficiency",
    "BhdFrontend",
    "make_detector_efficiency",
    "make_bhd_frontend",
    "TIA_MODELS",
    "LOADS",
    "device_params",
    "qe_from_responsivity",
    "eta_electronic",
    "r_meas_electronic",
    "imbalance_from_cmrr",
    "n_ex_from_rin",
    "eta_imbalance",
    "r_meas_imbalance",
    "h_tia",
    "h_dip",
    "monitor_voltage",
]

TIA_MODELS: tuple[str, ...] = ("first_order", "second_order")
LOADS: tuple[str, ...] = ("high_z", "50ohm")


# ------------------------------------------------------------------------------------------------
# formulas (T13 §1)


def qe_from_responsivity(responsivity_a_per_w: float, wavelength_m: float) -> float:
    """η_q = ℜ h ν / e (T13 §1.2)."""
    return float(responsivity_a_per_w) * sc.h * sc.c / (float(wavelength_m) * sc.e)


def eta_electronic(clearance_lin: Any) -> Any:
    """η_e = C/(1 + C) (conv §5); C = inf → 1."""
    c = np.asarray(clearance_lin, dtype=float)
    out = np.where(np.isinf(c), 1.0, c / (1.0 + np.where(np.isinf(c), 0.0, c)))
    return float(out) if out.ndim == 0 else out


def r_meas_electronic(r: Any, clearance_lin: Any) -> Any:
    """R_meas = (C R + 1)/(C + 1): measurement normalized to its own LO-only trace (T13 §1.3)."""
    c = np.asarray(clearance_lin, dtype=float)
    return (c * np.asarray(r, dtype=float) + 1.0) / (c + 1.0)


def imbalance_from_cmrr(cmrr_lin: float) -> float:
    """|δ|² = 1/(4·CMRR_lin), CMRR_lin = 10^(CMRR_dB/10) (T13 §1.5, Thorlabs definition)."""
    return 0.25 / float(cmrr_lin)


def n_ex_from_rin(p_lo_w: float, s_rin_per_hz: Any, qe: float, wavelength_m: float) -> Any:
    """LO excess amplitude noise in shot-noise units N_ex = η_q P_LO S_RIN/(2hν) (one-sided RIN, T13 §1.5)."""
    return qe * p_lo_w * np.asarray(s_rin_per_hz, dtype=float) * wavelength_m / (2.0 * sc.h * sc.c)


def eta_imbalance(delta_abs2: float, n_ex: Any) -> Any:
    """η_bal = 1/(1 + x), x = |δ|²(1 + N_ex) (T13 §1.5)."""
    return 1.0 / (1.0 + delta_abs2 * (1.0 + np.asarray(n_ex, dtype=float)))


def r_meas_imbalance(r: Any, delta_abs2: float, n_ex: Any) -> Any:
    """R_meas = (R + x)/(1 + x) (T13 §1.5)."""
    x = delta_abs2 * (1.0 + np.asarray(n_ex, dtype=float))
    return (np.asarray(r, dtype=float) + x) / (1.0 + x)


def h_tia(f_hz: Any, tia: Mapping[str, Any]) -> np.ndarray:
    """H_TIA(f) of ``tia`` = {model: first_order, f3db_hz} | {model: second_order, f0_hz, q}."""
    f = np.asarray(f_hz, dtype=float)
    if tia["model"] == "first_order":
        return 1.0 / (1.0 + 1j * f / float(tia["f3db_hz"]))
    if tia["model"] == "second_order":
        u = f / float(tia["f0_hz"])
        return 1.0 / (1.0 - u * u + 1j * u / float(tia["q"]))
    raise ValueError(f"tia.model must be one of {TIA_MODELS}, got {tia['model']!r}")


def h_dip(f_hz: Any, r: float, tau_s: float) -> np.ndarray:
    """H_dip(f) = (1 + r e^{−i2πfτ})/(1 + r) (T13 §1.4; Ref C dips: τ = 1/(2 f_dip))."""
    return (1.0 + r * np.exp(-2j * np.pi * np.asarray(f_hz, dtype=float) * tau_s)) / (1.0 + r)


def monitor_voltage(p_w: Any, conversion_v_per_w: float, v_max_v: float) -> Any:
    """MONITOR output V = conversion·P clipped to [0, v_max] (PDB manual p.8, p.14: 10 V/mW, 10 V → ~1 mW)."""
    return np.clip(conversion_v_per_w * np.asarray(p_w, dtype=float), 0.0, v_max_v)


def device_params(params: Mapping[str, Any], sections: tuple[str, ...], block: str) -> dict[str, Any]:
    """Spec params over the values of ``devices/<params['device']>.yaml`` sections (``device: null`` → params).

    Device values (null leaves dropped) fill keys the spec does not give; the spec always wins (D34).
    """
    require(params, "device", block=block)
    merged: dict[str, Any] = {}
    if params["device"] is not None:
        from cranberry.hardware.devices import device_values, load_device

        dev = device_values(load_device(str(params["device"])))
        for s in sections:
            merged.update({k: v for k, v in dev.get(s, {}).items() if v is not None})
    merged.update(params)
    return merged


def _one_of(p: Mapping[str, Any], keys: tuple[str, ...], block: str, spec: Mapping[str, Any] | None = None) -> str:
    """The one key of ``keys`` given; keys of the spec itself take priority over device-filled ones."""
    given = [k for k in keys if spec is not None and spec.get(k) is not None] or [k for k in keys if p.get(k) is not None]
    if len(given) != 1:
        raise MissingParameterError(f"block {block!r}: give exactly one of {keys} (D34), got {given}")
    return given[0]


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class _DetBase(BaseBlock):
    """Shared parameter handling (QE, clearance, imbalance, response) of the two BHD blocks."""

    DEVICE_SECTIONS: tuple[str, ...] = ("photodiode", "rf", "bhd")

    def _setup(self, kind: str, name: str | None, params: Mapping[str, Any]) -> None:
        name = name or kind
        self.p = device_params(params, self.DEVICE_SECTIONS, name)
        require(self.p, "p_lo_w", "toggles", block=name)
        super().__init__(kind, name, params)
        p, t = self.p, self._toggles
        self.p_lo_w = float(p["p_lo_w"])
        # responsivity (frontend currents) / QE (loss)
        self.responsivity = None if p.get("responsivity_a_per_w") is None else float(p["responsivity_a_per_w"])
        # clearance C_0 at p_lo_w, low frequency
        if t.get("circuit_noise", False):
            k = _one_of(p, ("clearance_lin", "i_n_a_rthz"), name, params)
            if k == "clearance_lin":
                require(p, "p_lo_ref_w", block=name)
                self.c0 = float(p["clearance_lin"]) * self.p_lo_w / float(p["p_lo_ref_w"])
            else:
                if self.responsivity is None:
                    raise MissingParameterError(f"block {name!r}: i_n_a_rthz needs responsivity_a_per_w")
                self.c0 = 2.0 * sc.e * self.responsivity * self.p_lo_w / float(p["i_n_a_rthz"]) ** 2
        else:
            self.c0 = math.inf
        # imbalance x = |δ|²(1 + N_ex)
        self.delta_abs2, self.n_ex = 0.0, 0.0
        if t.get("imbalance_leak", False):
            k = _one_of(p, ("imbalance_abs2", "cmrr_lin"), name, params)
            self.delta_abs2 = float(p["imbalance_abs2"]) if k == "imbalance_abs2" else imbalance_from_cmrr(p["cmrr_lin"])
            k = _one_of(p, ("n_ex", "s_rin_per_hz"), name, params)
            if k == "n_ex":
                self.n_ex = float(p["n_ex"])
            else:
                require(p, "wavelength_m", block=name)
                self.n_ex = float(n_ex_from_rin(self.p_lo_w, p["s_rin_per_hz"], self.eta_q_for_rin, p["wavelength_m"]))
        if t.get("freq_response", False):
            require(p, "tia", block=name)
            if p["tia"].get("model") not in TIA_MODELS:
                raise ValueError(f"block {name!r}: tia.model must be one of {TIA_MODELS}")
        if t.get("dip", False):
            require(p, "dip_r", "dip_tau_s", block=name)

    @property
    def eta_q_for_rin(self) -> float:
        """QE used in N_ex (1 if not known to this block)."""
        return 1.0

    def h(self, f_hz: Any) -> np.ndarray:
        """Complex response H(f) of the difference channel (toggles ``freq_response``, ``dip``)."""
        f = np.asarray(f_hz, dtype=float)
        out = np.ones(f.shape, dtype=complex)
        if self._toggles.get("freq_response", False):
            out = out * h_tia(f, self.p["tia"])
        if self._toggles.get("dip", False):
            out = out * h_dip(f, float(self.p["dip_r"]), float(self.p["dip_tau_s"]))
        return out

    def clearance(self, f_hz: Any) -> Any:
        """C(f) = |H(f)|² C_0 (inf with ``circuit_noise`` off)."""
        out = np.abs(self.h(f_hz)) ** 2 * self.c0 if math.isfinite(self.c0) else np.full(np.shape(f_hz), math.inf)
        return float(out) if np.ndim(out) == 0 else out

    @property
    def x_leak(self) -> float:
        """x = |δ|²(1 + N_ex) (0 with ``imbalance_leak`` off)."""
        return self.delta_abs2 * (1.0 + self.n_ex)

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state


class DetectorEfficiency(_DetBase):
    """``detector_efficiency`` (module docstring). Offline: registers QE, electronic and imbalance losses.

    Params (D34; ``device: <name>`` fills unspecified keys from sections photodiode / rf / bhd)::

        device, observer, applies_to, source, f_s_hz, p_lo_w, toggles
        qe | (responsivity_a_per_w, wavelength_m)                 # toggle qe
        clearance_lin + p_lo_ref_w | i_n_a_rthz (+ responsivity)  # toggle circuit_noise
        imbalance_abs2 | cmrr_lin;  n_ex | s_rin_per_hz (+ wavelength_m)   # toggle imbalance_leak
        tia: {model, f3db_hz | (f0_hz, q)}                        # toggle freq_response
        dip_r, dip_tau_s                                          # toggle dip
    """

    KIND = "detector_efficiency"
    toggle_names = ("qe", "circuit_noise", "imbalance_leak", "freq_response", "dip", "subtract_dark")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        p0 = device_params(params, self.DEVICE_SECTIONS, name or self.KIND)
        require(p0, "observer", "applies_to", "source", "f_s_hz", block=name or self.KIND)
        k = _one_of(p0, ("qe", "responsivity_a_per_w"), name or self.KIND, params)
        if k == "qe":
            self._eta_q = float(p0["qe"])
        else:
            require(p0, "wavelength_m", block=name or self.KIND)
            self._eta_q = qe_from_responsivity(p0["responsivity_a_per_w"], p0["wavelength_m"])
        if not 0.0 <= self._eta_q <= 1.0:
            raise ValueError(f"block {name!r}: quantum efficiency {self._eta_q} not in [0, 1]")
        self._setup(self.KIND, name, params)
        self.f_s_hz = float(self.p["f_s_hz"])

    @property
    def eta_q_for_rin(self) -> float:
        return self._eta_q

    @property
    def eta_q(self) -> float:
        return self._eta_q

    def breakdown(self, f_hz: Any | None = None) -> dict[str, Any]:
        """η_q, η_e, η_bal, η_det at ``f_hz`` (default f_s), with the toggles applied (off ⇒ 1)."""
        f = self.f_s_hz if f_hz is None else f_hz
        t = self._toggles
        eq = self._eta_q if t["qe"] else 1.0
        ee = 1.0 if (t["subtract_dark"] or not t["circuit_noise"]) else eta_electronic(self.clearance(f))
        eb = 1.0 / (1.0 + self.x_leak)
        return {"eta_q": eq, "eta_e": ee, "eta_bal": eb, "eta_det": eq * np.asarray(ee) * eb,
                "clearance_lin": self.clearance(f), "x_leak": self.x_leak}

    def losses(self) -> list[LossElement]:
        p, t = self.p, self._toggles
        common = dict(applies_to=set(p["applies_to"]), source=str(p["source"]), observer=str(p["observer"]))
        c = self.clearance(self.f_s_hz)
        return [
            LossElement(id=f"{self.name}.qe", eta=self._eta_q, kind="detector", enabled=t["qe"], **common),
            LossElement(id=f"{self.name}.electronic", eta=float(eta_electronic(c)), kind="electronic",
                        enabled=t["circuit_noise"] and not t["subtract_dark"], **common),
            LossElement(id=f"{self.name}.imbalance", eta=1.0 / (1.0 + self.x_leak), kind="detector",
                        enabled=t["imbalance_leak"], **common),
        ]


class BhdFrontend(_DetBase):
    """``bhd_frontend`` (module docstring). Measurement stage: RF-output PSDs, DC / MONITOR voltages.

    Params (D34; ``device: <name>`` fills keys from sections photodiode / rf / monitor / bhd)::

        device, p_lo_w, responsivity_a_per_w, gain_v_per_a, load (high_z | 50ohm),
        load_divisor_50ohm, coupler_sign (±1), v_off_v, v_sat_v, toggles
        clearance_lin + p_lo_ref_w | i_n_a_rthz                    # toggle circuit_noise
        imbalance_abs2 | cmrr_lin;  n_ex | s_rin_per_hz + wavelength_m   # toggle imbalance_leak
        tia, dip_r, dip_tau_s                                      # toggles freq_response, dip
        monitor_v_per_w, monitor_v_max_v                           # MONITOR (optional, high-Z)
    """

    KIND = "bhd_frontend"
    DEVICE_SECTIONS = ("photodiode", "rf", "monitor", "bhd")
    toggle_names = ("circuit_noise", "imbalance_leak", "freq_response", "dip", "subtract_dark", "saturation")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        self._setup(self.KIND, name, params)
        p = self.p
        require(p, "responsivity_a_per_w", "gain_v_per_a", "load", "coupler_sign", "v_off_v", block=self.name)
        if p["load"] not in LOADS:
            raise ValueError(f"block {self.name!r}: load must be one of {LOADS}")
        if p["coupler_sign"] not in (1, -1):
            raise ValueError(f"block {self.name!r}: coupler_sign must be +1 or -1")
        div = 1.0
        if p["load"] == "50ohm":
            require(p, "load_divisor_50ohm", block=self.name)
            div = float(p["load_divisor_50ohm"])
        self.gain = float(p["gain_v_per_a"]) / div
        if self._toggles["saturation"]:
            require(p, "v_sat_v", block=self.name)
        self.v_sat = None if p.get("v_sat_v") is None else float(p["v_sat_v"]) / div

    @property
    def s_shot_a2_per_hz(self) -> float:
        """One-sided shot-noise PSD of the difference current S_SN = 2eℜP_LO [A²/Hz] (T13 §1.1)."""
        return 2.0 * sc.e * self.responsivity * self.p_lo_w

    def psd(self, f_hz: Any, r: Any) -> dict[str, Any]:
        """RF-output PSDs [V²/Hz] for signal R(f) (shot-noise units, before electronic / imbalance loss).

        Returns ``signal``, ``shot`` (LO only, signal port vacuum), ``dark`` (LO off), and
        ``r_meas`` = (signal − d)/(shot − d) with d = dark if ``subtract_dark`` else 0 (T13 §1.3, §1.5).
        """
        h2 = np.abs(self.h(f_hz)) ** 2
        g2s = self.gain**2 * self.s_shot_a2_per_hz
        dark = np.zeros_like(h2) if not math.isfinite(self.c0) else g2s / self.c0 * np.ones_like(h2)
        x = self.x_leak
        signal = g2s * h2 * (np.asarray(r, dtype=float) + x) + dark
        shot = g2s * h2 * (1.0 + x) + dark
        d = dark if self._toggles["subtract_dark"] else 0.0
        return {"signal": signal, "shot": shot, "dark": dark, "r_meas": (signal - d) / (shot - d)}

    def samples(self, r: float, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        """White Monte-Carlo records (per-sample variance = low-frequency PSD) of signal, LO-only and dark."""
        s = self.psd(0.0, r)
        return {k: math.sqrt(float(s[k])) * rng.standard_normal(n) if s[k] > 0 else np.zeros(n)
                for k in ("signal", "shot", "dark")}

    def dc_voltage(self, p1_w: Any, p2_w: Any) -> Any:
        """RF OUTPUT DC v = s·G ℜ (P_1 − P_2) + V_off, clipped to ±v_sat (toggle ``saturation``)."""
        v = (self.p["coupler_sign"] * self.gain * self.responsivity * (np.asarray(p1_w) - np.asarray(p2_w))
             + float(self.p["v_off_v"]))
        return np.clip(v, -self.v_sat, self.v_sat) if self._toggles["saturation"] else v

    def monitor(self, p_w: Any) -> Any:
        """MONITOR± voltage of one photodiode (high-Z; saturates at v_max, ~1 mW for the PDB450C)."""
        require(self.p, "monitor_v_per_w", "monitor_v_max_v", block=self.name)
        return monitor_voltage(p_w, float(self.p["monitor_v_per_w"]), float(self.p["monitor_v_max_v"]))


@register_block(DetectorEfficiency.KIND)
def make_detector_efficiency(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``detector_efficiency`` (module docstring)."""
    return DetectorEfficiency(params, name)


@register_block(BhdFrontend.KIND)
def make_bhd_frontend(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``bhd_frontend`` (module docstring)."""
    return BhdFrontend(params, name)
