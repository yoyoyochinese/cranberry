"""Modulators: ``eom_phase``, ``dither_modulator``, ``aom_shift`` (PLAN §4.3 MOD; T16 §3.2, §3.5; T14; conv §4).

* ``eom_phase``: fibre EOM as a phase actuator, δψ = polarity · π u / V_π (``gain_mode: vpi``, runtime
  ``v_pi_v`` = V_π,LF at servo frequencies) or ``direct``; optional first-order bandwidth
  (``poles_hz``: [] for a static EOM, GHz devices; toggle ``bandwidth``); input clip
  ``drive_min_v``/``drive_max_v`` (toggle ``saturation``). Insertion loss is registered through
  ``losses()`` (G06, prompt §5.5) on ``segment``; on the squeezed path it caps squeezing (T16 §3.2).
  Same step as :class:`cranberry.blocks.act.pzt.PhaseActuator`. RAM is an error-signal effect (WP-B3).
* ``dither_modulator``: δψ = m sin(2π f_m t + φ) on a segment (``m_rad``, ``f_m_hz``, ``phase_rad``
  runtime). ``dither_mode``: ``timeseries`` writes it; ``analytic`` writes 0 and the effect enters the
  metric through χ = J0(2m) (conv §6; :meth:`DitherModulator.chi`, used by WP-F4). Toggle ``enabled``
  (False: 0). An electrical dither (Ref C: added to the EOM drive) is a ``sine_injection``.
* ``aom_shift`` (conv §4, T16 §3.5): order m = ±1 (structural ``order``). ``input: frequency``:
  dψ/dt = −m 2π δf_RF with δf_RF = ``k_hz_per_v`` · u (VCO), accumulated without wrapping (no range);
  ``input: phase``: δψ = −m φ_RF with φ_RF = ``k_rad_per_v`` · u. ``polarity`` keeps the physical sign
  as a calibration knob. Diffraction efficiency → ``insertion_loss_db`` (never on the squeezed path).
  Toggle ``enabled`` (False: 0) and ``insertion_loss``.

All three write ``act.<name>.dpsi`` [rad] and ``drv.<segment>.<driver>`` (split over ``sub``; ``segment``
may be null for no driver output).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
from scipy import special

from cranberry.blocks.act.pzt import PhaseActuator
from cranberry.blocks.ctrl._kit import DrvTarget, KitBlock, drv_write, drv_write_nb, loss_elements
from cranberry.core.block import require
from cranberry.core.registry import register_block

__all__ = ["EomPhase", "DitherModulator", "AomShift", "make_eom_phase", "make_dither_modulator", "make_aom_shift"]

_TWO_PI = 2.0 * math.pi
_NO_DRV = (0, 1, 0, 0, 0, -1, 1.0)


class EomPhase(PhaseActuator):
    KIND = "eom_phase"
    ALLOWED = ("vpi", "direct")
    DYN = "poles_hz"
    DYN_TOGGLE = "bandwidth"
    HYSTERESIS = False
    toggle_names = ("bandwidth", "saturation", "insertion_loss")


class _DrvSource(KitBlock):
    """A source writing a phase into act.<name>.dpsi and optionally a segment driver."""

    def __init__(self, params, name, reads, runtime):
        name = name or self.KIND
        require(params, "segment", block=name)
        self.drv = None if params["segment"] is None else DrvTarget(params, name, "rad", split=True)
        self.out = f"act.{name}.dpsi"
        writes = (self.out,) + (() if self.drv is None else (self.drv.slot,))
        super().__init__(params, name, reads=reads, writes=writes, runtime=runtime)

    def init(self, layout, rng):
        if self.drv is not None:
            self.drv.check_rate(self.rate_div)
        super().init(layout, rng)

    def slots(self):
        return [(self.out, (), np.float64, "rad", 0, 0.0)] + ([] if self.drv is None else [self.drv.slot_def()])

    def drv_args(self, layout):
        return (self.drv is not None, *(_NO_DRV if self.drv is None else self.drv.args(layout)))


def _make_dither_step(write):
    def step(S, I, P, i_t, i_out, pm, pf, pph, on, has_drv, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale):  # noqa: N803
        d = P[pm] * math.sin(_TWO_PI * P[pf] * S[i_t] + P[pph]) if on else 0.0
        S[i_out] = d
        if has_drv:
            write(S, I, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale, d)

    return step


class DitherModulator(_DrvSource):
    KIND = "dither_modulator"
    STEP = staticmethod(_make_dither_step(drv_write))
    STEP_NB = numba.njit(_make_dither_step(drv_write_nb))
    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        require(params, "dither_mode", block=name or self.KIND)
        if params["dither_mode"] not in ("timeseries", "analytic"):
            raise ValueError("dither_mode must be 'timeseries' or 'analytic'")
        super().__init__(params, name, ("t",), ("m_rad", "f_m_hz", "phase_rad"))

    def make_args(self, layout):
        on = self.toggles["enabled"] and self.params["dither_mode"] == "timeseries"
        return (self.o(layout, "t"), self.o(layout, self.out), self.p(layout, "m_rad"), self.p(layout, "f_m_hz"),
                self.p(layout, "phase_rad"), on, *self.drv_args(layout))

    def chi(self, params=None) -> float:
        """χ = J0(2m) of the dither (conv §6), 1 when disabled; m from ``params`` (runtime) or the spec."""
        if not self.toggles["enabled"]:
            return 1.0
        m = float(self.params["m_rad"] if params is None else params[self.param_name("m_rad")])
        return float(special.j0(2.0 * m))


def _make_aom_step(write):
    def step(S, I, P, i_in, i_out, pk, ppol, m, freq, h, on, has_drv,  # noqa: N803
             d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale):
        if not on:
            S[i_out] = 0.0
        elif freq:
            S[i_out] = S[i_out] - m * _TWO_PI * P[ppol] * P[pk] * S[i_in] * h
        else:
            S[i_out] = -m * P[ppol] * P[pk] * S[i_in]
        if has_drv:
            write(S, I, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale, S[i_out])

    return step


class AomShift(_DrvSource):
    KIND = "aom_shift"
    STEP = staticmethod(_make_aom_step(drv_write))
    STEP_NB = numba.njit(_make_aom_step(drv_write_nb))
    toggle_names = ("enabled", "insertion_loss")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        nm = name or self.KIND
        require(params, "in", "order", "input", block=nm)
        if params["order"] not in (1, -1) or params["input"] not in ("frequency", "phase"):
            raise ValueError(f"block {nm!r}: order must be +1/-1 and input 'frequency' | 'phase'")
        gain = "k_hz_per_v" if params["input"] == "frequency" else "k_rad_per_v"
        self.gain_key = gain
        super().__init__(params, name, (params["in"],), (gain, "polarity"))
        self.slot_units = {params["in"]: "V"}

    def make_args(self, layout):
        return (self.o(layout, self.params["in"]), self.o(layout, self.out), self.p(layout, self.gain_key),
                self.p(layout, "polarity"), float(self.params["order"]), self.params["input"] == "frequency",
                self.h, self.toggles["enabled"], *self.drv_args(layout))

    def losses(self):
        return loss_elements(self, self.toggles["insertion_loss"])


@register_block("eom_phase")
def make_eom_phase(params: Mapping[str, Any], *, name: str | None = None) -> EomPhase:
    """Factory of the ``eom_phase`` block (module docstring)."""
    return EomPhase(params, name)


@register_block("dither_modulator")
def make_dither_modulator(params: Mapping[str, Any], *, name: str | None = None) -> DitherModulator:
    """Factory of the ``dither_modulator`` block (module docstring)."""
    return DitherModulator(params, name)


@register_block("aom_shift")
def make_aom_shift(params: Mapping[str, Any], *, name: str | None = None) -> AomShift:
    """Factory of the ``aom_shift`` block (module docstring)."""
    return AomShift(params, name)
