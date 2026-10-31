"""``adc`` / ``dac``: quantisation, range, averaging, white noise, output divider (PLAN §4.3 CTRL; T16 §4.1).

Step (h = rate_div / f_sim)::

    x   = mean of the last ``avg_samples`` inputs                       (adc only; DSC1 2-point average)
    x  += noise_vrms · w,  w ~ N(0, 1)                                   (toggle noise; DrawBuffer, G11)
    q   = (hi − lo) / 2^bits,  c = (lo + hi)/2                           (range_v = [lo, hi])
    y   = c + q · clip(floor((x − c)/q + ½), −2^(bits−1), 2^(bits−1) − 1)   (toggle quantization;
          two's-complement codes, so the top value is hi − q)
          else clip(x, lo, hi)                                           (toggle saturation clips either way)
    y  *= R_load / (R_load + R_out)                                      (dac only; R_load = inf → 1)

``noise_vrms`` is the per-sample standard deviation (white up to the block's Nyquist frequency); a
band-limited spec value (DSC1 0.3 mV rms over 100 Hz–250 kHz) must be converted by the preset.
The input / output analogue bandwidths are separate ``lpf`` blocks (composition).

Params (structural; from the spec or ``device.adc`` / ``device.dac``, D34): ``in``, ``out``, ``device``,
``bits``, ``range_v``, ``noise_vrms``; adc ``avg_samples``; dac ``r_out_ohm`` (device
``output_impedance_ohm``), ``r_load_ohm``. Toggles: ``quantization``, ``noise``, ``saturation``.
Slots: ``out`` (V), ``ctrl.<name>.sat`` (int, ±1 when clipped), adc ``ctrl.<name>.hist``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock, device_section
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["Converter", "make_adc", "make_dac"]


def _step(S, I, P, i_in, i_out, i_sat, hist, n_avg, lo, hi, levels, sigma, gain, buf, idx,  # noqa: N803
          quant_on, noise_on, sat_on):
    x = S[i_in]
    if n_avg > 1:
        acc = x
        for j in range(n_avg - 1, 0, -1):
            S[hist + j] = S[hist + j - 1]
            acc += S[hist + j]
        S[hist] = x
        x = acc / n_avg
    if noise_on:
        x += sigma * S[buf + I[idx]]
        I[idx] += 1
    sat = 0
    if quant_on:
        q = (hi - lo) / levels
        c = 0.5 * (lo + hi)
        code = float(math.floor((x - c) / q + 0.5))
        if sat_on:
            if code > 0.5 * levels - 1.0:
                code, sat = 0.5 * levels - 1.0, 1
            elif code < -0.5 * levels:
                code, sat = -0.5 * levels, -1
        y = c + code * q
    else:
        y = x
        if sat_on:
            if y > hi:
                y, sat = hi, 1
            elif y < lo:
                y, sat = lo, -1
    S[i_out] = y * gain
    I[i_sat] = sat


class Converter(KitBlock):
    STEP = staticmethod(_step)
    STEP_NB = numba.njit(_step)
    toggle_names = ("quantization", "noise", "saturation")
    DEVICE_KEYS: dict[str, str] = {}
    EXTRA: tuple[str, ...] = ()

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "out", block=name)
        dev = device_section(params, self.KIND, name)
        dev = {mine: dev[theirs] for mine, theirs in self.DEVICE_KEYS.items() if theirs in dev}
        self.sat_name, self.hist_name = f"ctrl.{name}.sat", f"ctrl.{name}.hist"
        super().__init__(params, name, reads=(params["in"],), writes=(params["out"], self.sat_name), device=dev)
        self.need("bits", "range_v", "noise_vrms", *self.EXTRA)
        lo, hi = (float(v) for v in self.values["range_v"])
        if not hi > lo or int(self.values["bits"]) < 1:
            raise ValueError(f"block {self.name!r}: need range_v[1] > range_v[0] and bits >= 1")
        self.draws = DrawBuffer(self.name, 1, "standard_normal")
        if self.toggles["noise"]:
            self.writes = (*self.writes, *self.draws.slot_names)
        if self.n_avg() > 1:
            self.writes = (*self.writes, self.hist_name)
        self.slot_units = {params["out"]: "V"}

    def n_avg(self) -> int:
        return 1

    def gain(self) -> float:
        return 1.0

    def slots(self):
        out = [(self.params["out"], (), np.float64, "V", 0, 0.0), (self.sat_name, (), np.int64, "1", 0, 0)]
        if self.n_avg() > 1:
            out.append((self.hist_name, (self.n_avg(),), np.float64, "V", 0, 0.0))
        return out

    def init(self, layout, rng):
        super().init(layout, rng)
        if self.toggles["noise"]:
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self.toggles["noise"]:
            self.draws.refill(state, n0, n_samples)

    def make_args(self, layout):
        if self.toggles["noise"]:
            self.draws.bind(layout)
        lo, hi = (float(v) for v in self.values["range_v"])
        t = self.toggles
        return (self.o(layout, self.params["in"]), self.o(layout, self.params["out"]), self.o(layout, self.sat_name),
                self.o(layout, self.hist_name if self.n_avg() > 1 else None), self.n_avg(), lo, hi,
                float(2 ** int(self.values["bits"])), float(self.values["noise_vrms"]), self.gain(),
                self.draws.buf_off, self.draws.idx_off, t["quantization"], t["noise"], t["saturation"])

    @property
    def lsb_v(self) -> float:
        """Quantisation step q = (hi − lo) / 2^bits [V] (before the output divider)."""
        lo, hi = (float(v) for v in self.values["range_v"])
        return (hi - lo) / 2 ** int(self.values["bits"])


class Adc(Converter):
    KIND = "adc"
    DEVICE_KEYS = {"bits": "bits", "range_v": "range_v", "avg_samples": "avg_samples"}
    EXTRA = ("avg_samples",)

    def n_avg(self) -> int:
        return int(self.values["avg_samples"])


class Dac(Converter):
    KIND = "dac"
    DEVICE_KEYS = {"bits": "bits", "range_v": "range_v", "r_out_ohm": "output_impedance_ohm"}
    EXTRA = ("r_out_ohm", "r_load_ohm")

    def gain(self) -> float:
        rl, ro = float(self.values["r_load_ohm"]), float(self.values["r_out_ohm"])
        return 1.0 if math.isinf(rl) else rl / (rl + ro)


@register_block("adc")
def make_adc(params: Mapping[str, Any], *, name: str | None = None) -> Adc:
    """Factory of the ``adc`` block (module docstring)."""
    return Adc(params, name)


@register_block("dac")
def make_dac(params: Mapping[str, Any], *, name: str | None = None) -> Dac:
    """Factory of the ``dac`` block (module docstring)."""
    return Dac(params, name)
