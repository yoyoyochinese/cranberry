"""Signal helpers: ``adder``, ``ramp``, ``sine_injection`` (PLAN §4.3 CTRL; SimController contract, WP-F5).

* ``adder``: out = Σ_k gains[k] · ins[k]. Runtime ``gains`` (array). Params ``ins``, ``out``, ``unit``.
  No toggle (a sum has no effect to switch off; set a gain to 0).
* ``ramp``: out = offset_v + amplitude_v · w(t/period_s), w ∈ [−1, 1]; ``waveform`` ``sawtooth``
  (DSC1 p.7, p.20: from offset − amplitude to offset + amplitude, then jump back) or ``triangle``
  (Moku scan). Runtime ``amplitude_v``, ``period_s``, ``offset_v``; toggle ``enabled`` (False: out = 0).
  Params ``out``, ``waveform``.
* ``sine_injection``: out = in + amplitude · sin(2π freq_hz t) (``in`` may be null: out is the sine).
  Runtime ``freq_hz``, ``amplitude`` (unit of the injected signal); toggle ``enabled`` (False: out = in,
  identity). Params ``in``, ``out``, ``unit``. ``SimController.measure_transfer_function`` uses it.

t = n / f_sim is the time of the current sample (runner slot ``t``).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block

__all__ = ["make_adder", "make_ramp", "make_sine_injection", "WAVEFORMS"]

WAVEFORMS: dict[str, int] = {"sawtooth": 0, "triangle": 1}


def _adder_step(S, I, P, i_out, ins, pg):  # noqa: N803
    acc = 0.0
    for j in range(ins.shape[0]):
        acc += P[pg + j] * S[ins[j]]
    S[i_out] = acc


def _ramp_step(S, I, P, i_t, i_out, pa, pp, po, wave, enabled):  # noqa: N803
    if not enabled:
        S[i_out] = 0.0
        return
    x = S[i_t] / P[pp]
    ph = x - math.floor(x)
    w = 2.0 * ph - 1.0 if wave == 0 else 1.0 - 4.0 * abs(ph - 0.5)
    S[i_out] = P[po] + P[pa] * w


def _sine_step(S, I, P, i_t, i_in, i_out, pf, pa, enabled):  # noqa: N803
    y = S[i_in] if i_in >= 0 else 0.0
    if enabled:
        y += P[pa] * math.sin(2.0 * math.pi * P[pf] * S[i_t])
    S[i_out] = y


class Adder(KitBlock):
    KIND = "adder"
    STEP = staticmethod(_adder_step)
    STEP_NB = numba.njit(_adder_step)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "ins", "gains", "out", "unit", block=name)
        if not params["ins"] or len(params["ins"]) != len(params["gains"]):
            raise ValueError(f"block {name!r}: ins and gains need equal non-zero length")
        super().__init__(params, name, reads=tuple(params["ins"]), writes=(params["out"],), runtime=("gains",))
        self.slot_units = {params["out"]: str(params["unit"])}

    def slots(self):
        return [(self.params["out"], (), np.float64, str(self.params["unit"]), 0, 0.0)]

    def make_args(self, layout):
        ins = np.array([self.o(layout, s) for s in self.params["ins"]], dtype=np.int64)
        return (self.o(layout, self.params["out"]), ins, self.p(layout, "gains"))


class Ramp(KitBlock):
    KIND = "ramp"
    STEP = staticmethod(_ramp_step)
    STEP_NB = numba.njit(_ramp_step)
    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "out", "waveform", block=name)
        if params["waveform"] not in WAVEFORMS:
            raise ValueError(f"block {name!r}: waveform must be one of {sorted(WAVEFORMS)}")
        super().__init__(params, name, reads=("t",), writes=(params["out"],),
                         runtime=("amplitude_v", "period_s", "offset_v"))
        self.slot_units = {params["out"]: "V"}

    def slots(self):
        return [(self.params["out"], (), np.float64, "V", 0, 0.0)]

    def make_args(self, layout):
        return (self.o(layout, "t"), self.o(layout, self.params["out"]), self.p(layout, "amplitude_v"),
                self.p(layout, "period_s"), self.p(layout, "offset_v"), WAVEFORMS[self.params["waveform"]],
                self.toggles["enabled"])


class SineInjection(KitBlock):
    KIND = "sine_injection"
    STEP = staticmethod(_sine_step)
    STEP_NB = numba.njit(_sine_step)
    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "out", "unit", block=name)
        reads = ("t",) if params["in"] is None else ("t", params["in"])
        super().__init__(params, name, reads=reads, writes=(params["out"],), runtime=("freq_hz", "amplitude"))
        self.slot_units = {params["out"]: str(params["unit"])}

    def slots(self):
        return [(self.params["out"], (), np.float64, str(self.params["unit"]), 0, 0.0)]

    def make_args(self, layout):
        return (self.o(layout, "t"), self.o(layout, self.params["in"]), self.o(layout, self.params["out"]),
                self.p(layout, "freq_hz"), self.p(layout, "amplitude"), self.toggles["enabled"])


@register_block("adder")
def make_adder(params: Mapping[str, Any], *, name: str | None = None) -> Adder:
    """Factory of the ``adder`` block (module docstring)."""
    return Adder(params, name)


@register_block("ramp")
def make_ramp(params: Mapping[str, Any], *, name: str | None = None) -> Ramp:
    """Factory of the ``ramp`` block (module docstring)."""
    return Ramp(params, name)


@register_block("sine_injection")
def make_sine_injection(params: Mapping[str, Any], *, name: str | None = None) -> SineInjection:
    """Factory of the ``sine_injection`` block (module docstring)."""
    return SineInjection(params, name)
