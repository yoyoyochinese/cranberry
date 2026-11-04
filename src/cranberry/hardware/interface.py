"""Abstract controller interface (protected list; prompt §5.3, PLAN §2.10, T16 §4.6, G16-h). Owner: WP-F5.

The same optimisation loop (BoTorch, Phase 6) drives either the simulation (``SimController``) or real
hardware (``MokuController`` / ``DSC1Controller``, written when the hardware is used, G16-h) through
this interface only. A driver loop looks the same for every backend::

    ctl.configure(config)                    # once: structure, device presets, wiring
    for x in candidates:                     # BoTorch
        applied = ctl.set_params(x)          # canonical names → clipped/quantised values actually used
        ctl.reset_integrators(); ctl.engage()
        data = ctl.acquire(duration_s, ["ctrl.lock1.u", "v.pdb450c"])
        objective = f(data); ctl.disengage()
    record = ctl.log_record()                # JSON for the run record (prompt §6)

In-lock identification (Phase 6 P7, D85-7; backward-compatible additions to the D55 API)::

    ctl.engage()                             # every loop (unchanged call)
    ctl.disengage("lock0.fast")              # P7 literal: only the fast stage opens, the slow stage holds
    res = ctl.measure_loop_gain("lock0", "stepped_sine", freqs_hz, amplitude=0.05, settle_s=5e-3, dwell_s=25e-3)
    res["l_hat"], res["coherence"], res["phase_margin_deg"]
    data = ctl.acquire(0.1, ["v.adc0", "ctrl.exc0.x"], injections_on=True)   # in-loop record, excitation on

Names
    Runtime parameters use the canonical, device-agnostic names ``<block name>.<key>`` of the
    simulator's runtime params (``core.layout``), e.g. ``"lock1.pid.kp"``, ``"lock1.pid.ki_per_s"``,
    ``"lock1.pid.kd_s"`` (keys of a PID: ``hardware.mapping.CANONICAL_KEYS``). Signals use State slot
    names (``ctrl.lock1.u``, ``v.pdb450c``, ``lock.lockD.state``; conv §9). A hardware wrapper keeps a
    table from these names to its physical channels / probe points, so callers never see device names.

Units and signs
    SI with unit suffixes in names (G16-i): times [s], frequencies [Hz], voltages [V]. PID gains are
    canonical (K_p [1], K_i [1/s], K_d [s], error e = v − S); converting to the device's own
    representation (DSC1 e = S − v and entry units, Moku crossovers) happens inside the backend through
    ``hardware.mapping`` (conv §9). Transfer functions use the Laplace convention s = +i·2πf.

Ranges
    Values outside what the backend can realise are clipped/rounded to the nearest realisable value
    with a ``hardware.mapping.DeviceRangeWarning``; the applied values are returned (never silently
    changed).

Backend differences (documented, not hidden)
    Hardware runs continuously: ``acquire`` records the live system. ``SimController`` has no
    persistent physical state: each ``acquire`` runs the spec from its initial state with the current
    parameters and engage/ramp settings (a reproducible experiment per call). See ``hardware.sim``.
    ``measure_loop_gain`` is one locked acquisition with the excitation on in both: the simulator records
    it and runs ``meas.loop_identification.identify``; Moku:Delta maps it to its FRA (``hardware.moku``).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = ["ControllerInterface"]


class ControllerInterface(ABC):
    """Device-agnostic controller API (PLAN §2.10, T16 §4.6). See the module docstring for conventions."""

    @abstractmethod
    def configure(self, config: Mapping[str, Any]) -> None:
        """Apply a full controller configuration: loops (which controller drives which output and its
        device preset), ramp generators, injection points, and their wiring. Device-agnostic values only
        (T16 §4.6 item 1); the accepted keys are documented by each backend. Leaves every loop
        disengaged and every ramp disabled."""

    @abstractmethod
    def set_params(self, params: Mapping[str, float]) -> dict[str, float]:
        """Set runtime parameters (canonical names, e.g. ``"lock1.pid.kp"``). Returns the values actually
        applied (after device range clipping / quantisation; out-of-range values warn and are clipped).
        Unknown names raise ``KeyError``."""

    @abstractmethod
    def get_params(self) -> dict[str, float]:
        """Current runtime parameters (canonical names → applied values)."""

    @abstractmethod
    def engage(self, loop_id: str | None = None) -> None:
        """Close loop ``loop_id``; None (the original call) closes every configured loop. A stage of a
        stacked loop that must be switched on its own (P7: fast stage open, slow stage holding) is
        configured as a loop of its own (D85-7). Unknown ids raise ``KeyError``."""

    @abstractmethod
    def disengage(self, loop_id: str | None = None) -> None:
        """Open loop ``loop_id`` (None: every loop; outputs hold or return to offset, device-dependent,
        documented per device). Unknown ids raise ``KeyError``."""

    @abstractmethod
    def reset_integrators(self) -> None:
        """Zero all integrator states."""

    @abstractmethod
    def ramp(self, channel: str, *, amplitude: float, period_s: float, offset: float, enable: bool) -> None:
        """Configure the ramp generator on output ``channel`` (lock acquisition / phase scan). Waveform is
        device-dependent (DSC1 sawtooth, DSC1 manual p.7; Moku ramps/triangle, Moku p.45). ``amplitude``
        and ``offset`` are in the output's unit (V); out-of-range values warn and are clipped."""

    @abstractmethod
    def acquire(self, duration_s: float, signals: Sequence[str],
                injections_on: bool | Sequence[str] = False) -> dict[str, np.ndarray]:
        """Record ``signals`` for ``duration_s``; returns {name: float64 samples} plus ``"t_s"`` (common
        time axis; the sample rate is ``1 / diff(t_s)``). Unknown signal names raise ``KeyError``.
        ``injections_on``: False (the original behaviour) records with every injection / excitation off;
        True switches every configured one on, a sequence only the named points, at their current
        runtime parameters (``set_params``), so the in-loop signals can be recorded under excitation."""

    @abstractmethod
    def measure_transfer_function(self, injection_point: str, freqs_hz: np.ndarray) -> np.ndarray:
        """Swept-sine transfer function at ``injection_point`` (complex, one value per frequency,
        Laplace convention H(i·2πf)): ratio of the response to the reference signal configured for that
        point, both demodulated at the injected frequency. Used by VT-13 (prompt §7-7)."""

    def measure_loop_gain(self, loop_id: str, method: str, freqs_hz: Sequence[float] | None = None, *,
                          band_hz: Sequence[float] | None = None, amplitude: float,
                          **settings: Any) -> dict[str, Any]:
        """In-lock open-loop gain of loop ``loop_id`` from one locked acquisition with an excitation
        (Phase 6 P7, D85-7). The loop stays as engaged by the caller; the excitation point configured
        for the loop decides what is measured (summing point: L̂ = −resp/ref; actuator input with the
        stage open: the plant, L̂ = −P̂·C).

        Args:
            method: excitation waveform ``stepped_sine`` | ``multisine`` | ``chirp`` | ``prbs`` | ``white``
                (``blocks.ctrl.excitation``).
            freqs_hz: frequencies of a stepped sine / multisine, or None with ``band_hz`` (+ ``n_freqs``).
            band_hz: [f_lo, f_hi] (chirp sweep, analysis band of prbs / white, or the range of
                ``n_freqs`` log-spaced frequencies).
            amplitude: √2 × RMS of the excitation in the unit of the injection point.
            settings: timing and analysis keys of the backend (``settle_s``, ``dwell_s``, ``period_s``,
                ``n_periods``, ``t_start_s``, ...; documented per backend).

        Returns:
            dict with at least ``f_hz``, ``l_hat`` (complex L̂(i2πf), Laplace convention), ``coherence``,
            ``s_hat``, ``t_hat``, ``f_u_hz``, ``phase_margin_deg``, ``f_180_hz``, ``gain_margin_db``
            (``meas.loop_identification`` conventions). Backend extras are allowed.

        Not abstract, so backends written before D85-7 stay valid; they raise ``NotImplementedError``.
        """
        raise NotImplementedError(f"{type(self).__name__} does not implement measure_loop_gain (D85-7)")

    @abstractmethod
    def status(self) -> dict[str, Any]:
        """Lock state, saturation flags, errors (JSON-serialisable)."""

    @abstractmethod
    def log_record(self) -> dict[str, Any]:
        """JSON-serialisable snapshot of configuration + parameters for the run record (prompt §6;
        hardware adds model and firmware)."""
