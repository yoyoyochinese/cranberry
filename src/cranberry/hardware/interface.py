"""Abstract controller interface (protected list; prompt §5.3, PLAN §2.10, G16-h). Owner: WP-F5.

The same optimisation loop (BoTorch, Phase 6) drives either the simulation (``SimController``) or real
hardware (``MokuController`` / ``DSC1Controller``, written when the hardware is used). Parameters are
the device-agnostic canonical ones (``hardware.mapping.PidCoefficients``); mapping to a device's own
representation happens in ``hardware.mapping``.

Units: times [s], frequencies [Hz], voltages [V]; ``acquire`` returns sampled signals as float64
arrays with their sample rate.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = ["ControllerInterface"]


class ControllerInterface(ABC):
    """Device-agnostic controller API (PLAN §2.10)."""

    @abstractmethod
    def configure(self, config: Mapping[str, Any]) -> None:
        """Apply a full controller configuration (structure, device preset, wiring)."""

    @abstractmethod
    def set_params(self, params: Mapping[str, float]) -> dict[str, float]:
        """Set runtime parameters (canonical names, e.g. ``"lock1.pid.kp"``). Returns the values actually
        applied (after device range clipping / quantisation; out-of-range values warn and are clipped)."""

    @abstractmethod
    def get_params(self) -> dict[str, float]:
        """Current runtime parameters (canonical names)."""

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
        """Configure the triangular ramp generator on output ``channel`` (lock acquisition / phase scan)."""

    @abstractmethod
    def acquire(self, duration_s: float, signals: Sequence[str]) -> dict[str, np.ndarray]:
        """Record ``signals`` for ``duration_s``; returns {name: samples} plus ``"t_s"``."""

    @abstractmethod
    def measure_transfer_function(self, injection_point: str, freqs_hz: np.ndarray) -> np.ndarray:
        """Sine-injection transfer function at ``injection_point`` (complex, one value per frequency);
        used by VT-13."""

    @abstractmethod
    def status(self) -> dict[str, Any]:
        """Lock state, saturation flags, errors."""

    @abstractmethod
    def log_record(self) -> dict[str, Any]:
        """JSON-serialisable snapshot of configuration + parameters for the run record."""
