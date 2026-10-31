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
    def engage(self) -> None:
        """Close the loop(s)."""

    @abstractmethod
    def disengage(self) -> None:
        """Open the loop(s) (outputs hold or return to offset, device-dependent, documented per device)."""

    @abstractmethod
    def reset_integrators(self) -> None:
        """Zero all integrator states."""

    @abstractmethod
    def ramp(self, channel: str, *, amplitude: float, period_s: float, offset: float, enable: bool) -> None:
        """Configure the ramp generator on output ``channel`` (lock acquisition / phase scan). Waveform is
        device-dependent (DSC1 sawtooth, DSC1 manual p.7; Moku ramps/triangle, Moku p.45). ``amplitude``
        and ``offset`` are in the output's unit (V); out-of-range values warn and are clipped."""

    @abstractmethod
    def acquire(self, duration_s: float, signals: Sequence[str]) -> dict[str, np.ndarray]:
        """Record ``signals`` for ``duration_s``; returns {name: float64 samples} plus ``"t_s"`` (common
        time axis; the sample rate is ``1 / diff(t_s)``). Unknown signal names raise ``KeyError``."""

    @abstractmethod
    def measure_transfer_function(self, injection_point: str, freqs_hz: np.ndarray) -> np.ndarray:
        """Swept-sine transfer function at ``injection_point`` (complex, one value per frequency,
        Laplace convention H(i·2πf)): ratio of the response to the reference signal configured for that
        point, both demodulated at the injected frequency. Used by VT-13 (prompt §7-7)."""

    @abstractmethod
    def status(self) -> dict[str, Any]:
        """Lock state, saturation flags, errors (JSON-serialisable)."""

    @abstractmethod
    def log_record(self) -> dict[str, Any]:
        """JSON-serialisable snapshot of configuration + parameters for the run record (prompt §6;
        hardware adds model and firmware)."""
