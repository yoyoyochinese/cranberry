"""IIR filters as second-order sections: ``iir``, ``lpf``, ``hpf``, ``band_split`` (PLAN §4.3 CTRL; T16 §4.3–4.4).

All four run the same cascade of biquads (direct form II transposed, a0 = 1), coefficients fixed at
build time (structural, h = rate_div / f_sim):

* ``iir``: ``sos`` given ([[b0, b1, b2, a0, a1, a2], ...], scipy layout, normalised so a0 = 1).
* ``lpf`` / ``hpf``: Butterworth of ``order`` at ``f_c_hz`` designed by ``method``:
  ``bilinear`` (scipy ``butter(..., fs)``: bilinear transform prewarped at f_c, exact −3 dB at f_c) or
  ``zoh`` (analogue Butterworth, zero-order-hold equivalent, ``scipy.signal.cont2discrete``).
* ``band_split``: the ``lpf`` design on ``in`` → ``out_low``, and ``out_high = in − out_low``
  (exactly complementary: low + high = in at every sample; parallel band-split loops, T16 §4.4).

Params: ``in``, ``out`` (``band_split``: ``out_low``, ``out_high``), ``unit`` (unit tag of the outputs),
design keys above. Toggle ``enabled``: False makes the block an identity (``band_split``: low = in,
high = 0). Slots: outputs, ``ctrl.<name>.z`` (section states, shape (n_sections, 2)).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
from scipy import signal

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block

__all__ = ["make_iir", "make_lpf", "make_hpf", "make_band_split", "butter_sos"]


def _step(S, I, P, i_in, i_out, i_high, z, sos, enabled):  # noqa: N803
    x0 = S[i_in]
    x = x0
    if enabled:
        for s in range(sos.shape[0]):
            z1 = S[z + 2 * s]
            y = sos[s, 0] * x + z1
            S[z + 2 * s] = sos[s, 1] * x - sos[s, 4] * y + S[z + 2 * s + 1]
            S[z + 2 * s + 1] = sos[s, 2] * x - sos[s, 5] * y
            x = y
    S[i_out] = x
    if i_high >= 0:
        S[i_high] = x0 - x


def butter_sos(btype: str, order: int, f_c_hz: float, fs_hz: float, method: str) -> np.ndarray:
    """Butterworth SOS (a0 = 1) for ``method`` ``bilinear`` (prewarped at f_c) or ``zoh``."""
    if not (0 < f_c_hz < fs_hz / 2):
        raise ValueError(f"f_c_hz must be in (0, fs/2) = (0, {fs_hz / 2:g}), got {f_c_hz!r}")
    if method == "bilinear":
        return signal.butter(int(order), f_c_hz, btype=btype, fs=fs_hz, output="sos")
    if method == "zoh":
        z, p, k = signal.butter(int(order), 2 * np.pi * f_c_hz, btype=btype, analog=True, output="zpk")
        ad, bd, cd, dd, _ = signal.cont2discrete(signal.zpk2ss(z, p, k), 1.0 / fs_hz, method="zoh")
        num, den = signal.ss2tf(ad, bd, cd, dd)
        num = num.ravel()
        num[np.abs(num) < 1e-12 * np.max(np.abs(num))] = 0.0  # ZOH: leading coefficient is 0 up to rounding
        d = len(num) - len(np.trim_zeros(num, "f"))
        # scipy's tf2sos drops leading numerator zeros (the ZOH delay): keep them as pure z⁻¹ sections
        delay = np.tile([0.0, 1.0, 0.0, 1.0, 0.0, 0.0], (d, 1))
        return np.vstack([signal.tf2sos(num[d:], den), delay])
    raise ValueError(f"method must be 'bilinear' or 'zoh', got {method!r}")


class _Sos(KitBlock):
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("enabled",)
    OUTS = ("out",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "unit", *self.OUTS, *self.DESIGN, block=name)
        self.z_name = f"ctrl.{name}.z"
        super().__init__(params, name, reads=(params["in"],), writes=(*(params[k] for k in self.OUTS), self.z_name))
        self.slot_units = {params[k]: str(params["unit"]) for k in self.OUTS}
        self.sos: np.ndarray | None = None

    def design(self, fs_hz: float) -> np.ndarray:  # pragma: no cover - overridden
        raise NotImplementedError

    def init(self, layout, rng):
        self.sos = np.ascontiguousarray(self.design(layout.f_sim_hz / self.rate_div), dtype=np.float64)
        if self.sos.ndim != 2 or self.sos.shape[1] != 6 or not np.allclose(self.sos[:, 3], 1.0):
            raise ValueError(f"block {self.name!r}: sos must be (n, 6) with a0 = 1")
        super().init(layout, rng)

    def slots(self):
        u = str(self.params["unit"])
        return [*((self.params[k], (), np.float64, u, 0, 0.0) for k in self.OUTS),
                (self.z_name, (self.sos.shape[0], 2), np.float64, u, 0, 0.0)]

    def make_args(self, layout):
        outs = [self.o(layout, self.params[k]) for k in self.OUTS]
        return (self.o(layout, self.params["in"]), outs[0], outs[1] if len(outs) > 1 else -1,
                self.o(layout, self.z_name), self.sos, self.toggles["enabled"])

    def freq_response(self, f_hz: Any) -> np.ndarray:
        """Discrete response of the (low / single) output at f_hz, z = e^{+i2πfh} (needs ``init``)."""
        _, h = signal.sosfreqz(self.sos, worN=np.atleast_1d(np.asarray(f_hz, float)), fs=1.0 / self.h)
        return h if self.toggles["enabled"] else np.ones_like(h)


class Iir(_Sos):
    KIND, DESIGN = "iir", ("sos",)

    def design(self, fs_hz):
        sos = np.asarray(self.params["sos"], dtype=np.float64).reshape(-1, 6)
        return sos / sos[:, 3:4]


class Lpf(_Sos):
    KIND, DESIGN, BTYPE = "lpf", ("order", "f_c_hz", "method"), "lowpass"

    def design(self, fs_hz):
        return butter_sos(self.BTYPE, self.params["order"], float(self.params["f_c_hz"]), fs_hz, self.params["method"])


class Hpf(Lpf):
    KIND, BTYPE = "hpf", "highpass"


class BandSplit(Lpf):
    KIND, OUTS = "band_split", ("out_low", "out_high")


@register_block("iir")
def make_iir(params: Mapping[str, Any], *, name: str | None = None) -> Iir:
    """Factory of the ``iir`` block (module docstring)."""
    return Iir(params, name)


@register_block("lpf")
def make_lpf(params: Mapping[str, Any], *, name: str | None = None) -> Lpf:
    """Factory of the ``lpf`` block (module docstring)."""
    return Lpf(params, name)


@register_block("hpf")
def make_hpf(params: Mapping[str, Any], *, name: str | None = None) -> Hpf:
    """Factory of the ``hpf`` block (module docstring)."""
    return Hpf(params, name)


@register_block("band_split")
def make_band_split(params: Mapping[str, Any], *, name: str | None = None) -> BandSplit:
    """Factory of the ``band_split`` block (module docstring)."""
    return BandSplit(params, name)
