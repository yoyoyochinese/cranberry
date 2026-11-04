"""Moku:Delta mapping of the in-lock identification (D85-7): a documented adapter stub, no network calls.

Facts used (``notes/fact-p6-thermo-optic-moku.md`` §2, ``refs/manuals/MokuDelta-Specifications.pdf`` v26-0817):

* The Frequency Response Analyzer (FRA) runs in Multi-instrument Mode (MiM, p79–80) next to the PID Controller;
  its measurement modes are on p30: In/Out (dB), In/In1 (dB) and In (dBm). H = e/x (the FRA input over its own
  output) is the **In/Out** mode. The same page documents a Normalization (reference sweep) that divides out the
  cable / front-end response of the FRA path; it does not see the control-matrix path (g below). The FRA sweep is a stepped sine (settling + averaging per point), i.e. the
  ``stepped_sine`` of ``blocks.ctrl.excitation``.
* MiM wiring cannot add two sources into one destination. The sum e + x is formed **inside the PID Controller**
  by its 4-input control matrix (coefficients −20…+20, p64–65): row of the loop's output = 1·In_err + g·In_exc,
  with In_exc the FRA output routed in MiM. PID Controller: 4 in / 4 out, "400 kHz bandwidth, 30° lag" (p64–65).
* The Laser Lock Box documents no injection point (p43–44): the loop must run in the PID Controller for this
  procedure.
* **UNVERIFIED: the unit of the control-matrix coefficient in the API** (the Liquid Instruments API documents a
  2-channel form in dB, the PDF a linear coefficient; a dB / linear mix-up gives errors of 10–20× in g).

Calibration of g (Phase 5 review I-1). g is the complex response from the FRA output to the PID input through the
MiM route and the control matrix, not only the matrix coefficient. It must be calibrated **as a complex function of
frequency across the whole identification band**: g(f) = (PID input monitor)/(FRA output) with the loop OPEN (PID
output disabled or the actuator disconnected), with the same FRA settings, and passed to
:func:`reduce_summing_point` as an array on the identification frequencies (a DC check alone is not enough). The
error of an assumed g′ against the true g follows from (1 + L̂)/L̂ = (g′/g)(1 + L)/L, so it is not a fixed 12 % per
dB: for |L| ≲ 1 a 1 dB error in |g| gives ≈ 9–13 % in |L̂|, it grows with |L| (refD model, +1 dB: −16 % at |L| = 7,
f_u +9 %, PM −7°, GM −0.8 dB), and for |L| → ∞ L̂ tends to 1/(g′/g − 1) (8.2 for +1 dB) instead of growing with L,
so the low-frequency high-gain points are the most affected. A pure route delay of 10–100 ns is negligible
(≤ 0.1 %) below 100 kHz but belongs to g(f) at higher frequencies (reviews/scripts/p6_wp6i_x_checks.py part 1).

Signal flow and reduction (summing point, the ``minus_l`` relation of ``meas.loop_identification``)::

    error e ──► Moku In1 ──┬──────────────────────────► PID matrix: 1·In1 + g·FRA ──► PID ──► Out1 ──► actuator
                           └─► FRA input (MiM route)        ▲
                                FRA output x ───────────────┘   (MiM route into the PID's second matrix input)

The sum e + g·x is internal to the PID, so it cannot be recorded; the FRA measures H = e/x in its In/Out
mode (the input over the FRA output, spec p30). With e = −L(e + g x):

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


def reduce_summing_point(h_e_over_x: Any, g: Any) -> dict[str, np.ndarray]:
    """L̂, S, T from the FRA ratio H = e/x with the injection gain ``g``: a scalar or, as calibrated (module
    docstring), a complex array g(f) on the same frequencies as ``h_e_over_x``."""
    h = np.asarray(h_e_over_x, dtype=complex)
    g = np.asarray(g, dtype=complex)
    if g.ndim and g.shape != h.shape:
        raise ValueError(f"g(f) must have the shape of H {h.shape}, got {g.shape}")
    if np.any(g == 0):
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
            "fra": {"mode": "In/Out", "normalization": "optional reference sweep of the FRA path (spec p30)",
"start_hz": float(f[0]), "stop_hz": float(f[-1]),
                    "points": int(f.size), "spacing": "log", "amplitude_vpp": 2.0 * float(amplitude),
                    "settling_s": float(settle_s), "averaging_s": float(dwell_s - settle_s)},
            "g_calibration": "complex g(f) = PID input monitor / FRA output, loop open, same sweep (module docstring)",
            "reduction": "L = -H/(g + H), S = (g + H)/g, T = -H/g with H = e/x (reduce_summing_point, g(f) array)",
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
