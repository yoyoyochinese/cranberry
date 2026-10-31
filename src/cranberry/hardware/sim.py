"""``SimController``: the ControllerInterface implemented on the time-loop simulation (PLAN §2.10).
Owner: WP-F5 (connects to ``runner.run`` once WP-F1 lands).

``acquire`` runs the simulation and returns the recorded series; ``measure_transfer_function`` injects
sines (VT-13).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from cranberry.hardware.interface import ControllerInterface

__all__ = ["SimController"]

_MSG = "hardware.sim.SimController is implemented by WP-F5"


class SimController(ControllerInterface):
    """Simulation-backed controller."""

    def __init__(self, spec: Any):
        raise NotImplementedError(_MSG)

    def configure(self, config: Mapping[str, Any]) -> None:
        raise NotImplementedError(_MSG)

    def set_params(self, params: Mapping[str, float]) -> dict[str, float]:
        raise NotImplementedError(_MSG)

    def get_params(self) -> dict[str, float]:
        raise NotImplementedError(_MSG)

    def engage(self) -> None:
        raise NotImplementedError(_MSG)

    def disengage(self) -> None:
        raise NotImplementedError(_MSG)

    def reset_integrators(self) -> None:
        raise NotImplementedError(_MSG)

    def ramp(self, channel: str, *, amplitude: float, period_s: float, offset: float, enable: bool) -> None:
        raise NotImplementedError(_MSG)

    def acquire(self, duration_s: float, signals: Sequence[str]) -> dict[str, np.ndarray]:
        raise NotImplementedError(_MSG)

    def measure_transfer_function(self, injection_point: str, freqs_hz: np.ndarray) -> np.ndarray:
        raise NotImplementedError(_MSG)

    def status(self) -> dict[str, Any]:
        raise NotImplementedError(_MSG)

    def log_record(self) -> dict[str, Any]:
        raise NotImplementedError(_MSG)
