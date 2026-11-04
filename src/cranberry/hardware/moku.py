"""Moku:Delta mapping of the in-lock identification (D85-7): a documented adapter stub, no network calls.

Facts used (``notes/fact-p6-thermo-optic-moku.md`` §2, ``refs/manuals/MokuDelta-Specifications.pdf`` v26-0817):

* The Frequency Response Analyzer (FRA) runs in Multi-instrument Mode (MiM, p79–80) next to the PID Controller;
  its In/In1 modes are on p30. The FRA sweep is a stepped sine (settling + averaging per point), i.e. the
  ``stepped_sine`` of ``blocks.ctrl.excitation``.
* MiM wiring cannot add two sources into one destination. The sum e + x is formed **inside the PID Controller**
  by its 4-input control matrix (coefficients −20…+20, p64–65): row of the loop's output = 1·In_err + g·In_exc,
  with In_exc the FRA output routed in MiM. PID Controller: 4 in / 4 out, "400 kHz bandwidth, 30° lag" (p64–65).
* The Laser Lock Box documents no injection point (p43–44): the loop must run in the PID Controller for this
  procedure.
* **UNVERIFIED: the unit of the control-matrix coefficient in the API** (the Liquid Instruments API documents a
  2-channel form in dB, the PDF a linear coefficient). g enters the reduction below directly (1 dB = 12 %), so g
  must be confirmed on the bench (e.g. set 1, apply a DC offset on In_exc with the loop open, read the PID input
  monitor) before L̂ is trusted.

Signal flow and reduction (summing point, the ``minus_l`` relation of ``meas.loop_identification``)::

    error e ──► Moku In1 ──┬──────────────────────────► PID matrix: 1·In1 + g·FRA ──► PID ──► Out1 ──► actuator
                           └─► FRA input (MiM route)        ▲
                                FRA output x ───────────────┘   (MiM route into the PID's second matrix input)

The sum e + g·x is internal to the PID, so it cannot be recorded; the FRA measures H = e/x in its In/Out
normalisation (In1 mode: the input over the FRA output). With e = −L(e + g x):

    L̂ = −H / (g + H),   S = (g + H)/g = 1/(1 + L̂),   T = −H/g = L̂/(1 + L̂)          (:func:`reduce_summing_point`)

The other waveforms (multisine, chirp, PRBS, white) are not FRA modes; on Moku:Delta they would need the
Waveform Generator + Data Logger in MiM with the same control-matrix sum, then ``identify`` offline (UNVERIFIED:
whether the Data Logger can record the internal FRA / generator output at the needed rate).

``MokuDeltaController`` keeps the canonical-name table and produces the instrument settings
(:meth:`MokuDeltaController.loop_gain_plan`); every method that would talk to the device raises
``NotImplementedError`` (G16-h: the driver is written when the hardware is used).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from cranberry.hardware.interface import ControllerInterface

__all__ = ["MokuDeltaController", "reduce_summing_point", "MATRIX_RANGE", "GAIN_UNIT_STATUS"]

#: Control-matrix coefficient range of the PID Controller (p64–65).
MATRIX_RANGE: tuple[float, float] = (-20.0, 20.0)
#: The API unit of the coefficient is not confirmed (module docstring).
GAIN_UNIT_STATUS = "UNVERIFIED: dB (API docs) vs linear (PDF p64-65); confirm on the bench"


def reduce_summing_point(h_e_over_x: Any, g: float) -> dict[str, np.ndarray]:
    """L̂, S, T from the FRA ratio H = e/x with control-matrix gain ``g`` (module docstring)."""
    h = np.asarray(h_e_over_x, dtype=complex)
    if g == 0:
        raise ValueError("control-matrix gain g must be non-zero")
    return {"l_hat": -h / (g + h), "s_hat": (g + h) / g, "t_hat": -h / g}


def _stub(what: str) -> NotImplementedError:
    return NotImplementedError(f"MokuDeltaController.{what}: adapter stub, no device I/O (D85-7, G16-h)")


class MokuDeltaController(ControllerInterface):
    """Adapter stub (module docstring).

    Args:
        channels: canonical name → Moku resource, e.g. ``{"v.adc0": "In1", "ctrl.pid0f.u": "Out1"}``.
        matrix_gain: the intended linear coefficient g of the FRA input in the PID control matrix
            (its API encoding is UNVERIFIED, :data:`GAIN_UNIT_STATUS`).
    """

    def __init__(self, channels: Mapping[str, str], *, matrix_gain: float = 1.0):
        lo, hi = MATRIX_RANGE
        if not lo <= matrix_gain <= hi or matrix_gain == 0:
            raise ValueError(f"matrix_gain must be non-zero in [{lo}, {hi}], got {matrix_gain!r}")
        self.channels = dict(channels)
        self.matrix_gain = float(matrix_gain)

    def loop_gain_plan(self, loop_id: str, method: str, freqs_hz: Sequence[float], *, amplitude: float,
                       settle_s: float, dwell_s: float, error_in: str, output: str) -> dict[str, Any]:
        """Instrument settings for ``measure_loop_gain`` at the summing point (no device I/O)."""
        if method != "stepped_sine":
            raise NotImplementedError(f"method {method!r} is not an FRA mode (module docstring)")
        f = np.sort(np.asarray(freqs_hz, dtype=float))
        if f.size < 2 or f[0] <= 0:
            raise ValueError("freqs_hz needs >= 2 positive frequencies")
        return {
            "loop_id": loop_id,
            "mim": {"instruments": ["PID Controller", "Frequency Response Analyzer"],
                    "routes": [(self.channels[error_in], "PID.In1"), (self.channels[error_in], "FRA.In1"),
                               ("FRA.Out1", "PID.In2"), ("PID.Out1", self.channels[output])]},
            "pid_control_matrix": {"row": "PID.Out1", "coefficients": {"In1": 1.0, "In2": self.matrix_gain},
                                   "gain_unit": GAIN_UNIT_STATUS},
            "fra": {"mode": "In1 (input / output normalisation)", "start_hz": float(f[0]), "stop_hz": float(f[-1]),
                    "points": int(f.size), "spacing": "log", "amplitude_vpp": 2.0 * float(amplitude),
                    "settling_s": float(settle_s), "averaging_s": float(dwell_s - settle_s)},
            "reduction": "L = -H/(g + H), S = (g + H)/g, T = -H/g with H = e/x (reduce_summing_point)",
        }

    # ------------------------------------------------------------ ControllerInterface (no device I/O)

    def configure(self, config):
        raise _stub("configure")

    def set_params(self, params):
        raise _stub("set_params")

    def get_params(self):
        raise _stub("get_params")

    def engage(self, loop_id=None):
        raise _stub("engage")

    def disengage(self, loop_id=None):
        raise _stub("disengage")

    def reset_integrators(self):
        raise _stub("reset_integrators")

    def ramp(self, channel, *, amplitude, period_s, offset, enable):
        raise _stub("ramp")

    def acquire(self, duration_s, signals, injections_on=False):
        raise _stub("acquire")

    def measure_transfer_function(self, injection_point, freqs_hz):
        raise _stub("measure_transfer_function")

    def measure_loop_gain(self, loop_id, method, freqs_hz=None, *, band_hz=None, amplitude, **settings):
        raise _stub("measure_loop_gain (use loop_gain_plan + reduce_summing_point)")

    def status(self):
        return {"backend": "moku_delta_stub", "connected": False, "matrix_gain": self.matrix_gain,
                "gain_unit": GAIN_UNIT_STATUS}

    def log_record(self):
        return {"controller": "MokuDeltaController (stub)", "channels": self.channels,
                "matrix_gain": self.matrix_gain, "gain_unit": GAIN_UNIT_STATUS}
