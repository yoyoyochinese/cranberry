"""``demodulator`` and ``rf_reference`` (T14 §2.2–§2.3, §3.6, §12.1). Owner: WP-B3.

``demodulator``
    Reference 2cos(2π f_ref t + φ), φ = ``phi_d_rad`` + φ_RF(t) (``in_phase``, a slot written by
    ``rf_reference`` — normally through a ``delay_buffer`` of τ_RF; toggle
    ``rf_reference_phase_noise``). Input from a generator (``_common`` signal representation):

    * baseband (``rf_rate_mode`` off, input shape (2,) = (a, b)): I = a cos φ − b sin φ,
      Q = a sin φ + b cos φ (the exact DC terms of the RF products);
    * RF rate (``rf_rate_mode`` on, scalar input x): I = 2x cos(2π f_ref t + φ), Q = 2x sin(...).

    Then ``lpf_order`` (structural, ≥ 0) cascaded one-pole low-pass sections y += α(x − y),
    α = 1 − exp(−2π f_lpf h), on I and Q; ``out_i`` = I + V_off (toggle ``demod_offset``), ``out_q`` = Q.
    White noise of PSD S at the input gives PSD 2S on I (T14 §2.2).

``rf_reference``
    RF phase φ_RF(t) of an AOM drive / demodulation reference: white phase noise (one-sided
    ``s_phi_white_rad2_hz``) plus white frequency noise (``s_nu_white_hz2_hz``, a random walk,
    S_φ = S_ν/f², conv §7; per-step increment variance 2π² S_ν h). Writes ``out`` [rad] (toggle
    ``phase_noise``, off = zero source). A delay τ_RF between the AOM and the demodulator is a
    ``delay_buffer`` block; their difference then has S = 4 sin²(π f τ_RF) S_φRF (T14 §3.6). Declares
    S_φRF above its Nyquist frequency to the ``BandLedger`` on ``band_targets`` (G05).
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.err._common import ErrBlock, helpers, register
from cranberry.core.block import require
from cranberry.offline.band_ledger import BandDecl

__all__ = ["Demodulator", "RfReference", "make_demodulator", "make_rf_reference"]

TWO_PI = 2.0 * math.pi
_DP = ("f_ref_hz", "phi_d_rad", "v_off_v", "lpf_hz")
C_DEM = namedtuple("C_DEM", ("s_in", "s_phase", "s_i", "s_q", "s_lpf", "t", "buf", "idx", "fb", "mode", "law",
                             "f_rf", "f_demod_offset", "f_rf_reference_phase_noise", *(f"p_{k}" for k in _DP)))
C_RFR = namedtuple("C_RFR", ("s_out", "s_walk", "t", "buf", "idx", "fb", "mode", "law", "f_phase_noise",
                             "p_s_phi_white_rad2_hz", "p_s_nu_white_hz2_hz"))


def _build(jit: bool):
    draw = helpers(jit).draw
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def dem_step(S, Z, I, P, c):  # noqa: N803
        phi = P[c.p_phi_d_rad]
        if c.f_rf_reference_phase_noise:
            phi += S[c.s_phase]
        if c.f_rf:
            arg = TWO_PI * P[c.p_f_ref_hz] * S[c.t] + phi
            x = S[c.s_in]
            vi, vq = 2.0 * x * math.cos(arg), 2.0 * x * math.sin(arg)
        else:
            a, b = S[c.s_in], S[c.s_in + 1]
            cp, sp = math.cos(phi), math.sin(phi)
            vi, vq = a * cp - b * sp, a * sp + b * cp
        n = c.mode  # LPF order
        if n > 0:
            alpha = 1.0 - math.exp(-TWO_PI * P[c.p_lpf_hz] / c.fb)
            for j in range(n):
                S[c.s_lpf + j] += alpha * (vi - S[c.s_lpf + j])
                vi = S[c.s_lpf + j]
                S[c.s_lpf + n + j] += alpha * (vq - S[c.s_lpf + n + j])
                vq = S[c.s_lpf + n + j]
        S[c.s_i] = vi + (P[c.p_v_off_v] if c.f_demod_offset else 0.0)
        S[c.s_q] = vq

    @dec
    def rfr_step(S, Z, I, P, c):  # noqa: N803
        if not c.f_phase_noise:
            S[c.s_out] = 0.0
            return
        h = 1.0 / c.fb
        S[c.s_walk] += math.sqrt(2.0 * math.pi**2 * P[c.p_s_nu_white_hz2_hz] * h) * draw(S, I, c.buf, c.idx)
        S[c.s_out] = S[c.s_walk] + math.sqrt(0.5 * P[c.p_s_phi_white_rad2_hz] * c.fb) * draw(S, I, c.buf, c.idx)

    return dem_step, rfr_step


_PY = _build(False)
_NB = _build(True)


class Demodulator(ErrBlock):
    KIND, CONST, STEP = "demodulator", C_DEM, (_PY[0], _NB[0])

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "in_phase", "out_i", "out_q", "lpf_order", "toggles", block=name)
        order = params["lpf_order"]
        if isinstance(order, bool) or int(order) != order or order < 0:
            raise ValueError(f"block {name!r}: lpf_order must be an int >= 0, got {order!r}")
        tg = params["toggles"]
        if tg.get("rf_reference_phase_noise") and params["in_phase"] is None:
            raise ValueError(f"block {name!r}: rf_reference_phase_noise is on but in_phase is null")
        self._rf = bool(tg.get("rf_rate_mode"))
        self.setup(
            self.KIND, name, params,
            inputs={"in": (params["in"], "V"),
                    "phase": (params["in_phase"] if tg.get("rf_reference_phase_noise") else None, "rad")},
            outputs={"i": (params["out_i"], (), "V"), "q": (params["out_q"], (), "V"),
                     "lpf": (f"dem.{name}.lpf", (max(2 * int(order), 1),), "V")},
            runtime=_DP, toggles=("rf_rate_mode", "demod_offset", "rf_reference_phase_noise"),
            n_draws=0, mode=int(order),
        )

    @property
    def flags(self) -> dict[str, bool]:
        f = {k: v for k, v in self.toggles.items() if k != "rf_rate_mode"}
        f["rf"] = self._rf
        return f

    def bind(self, layout):
        shape = layout.slot(self.params["in"]).shape
        want = () if self._rf else (2,)
        if tuple(shape) != want:
            raise ValueError(f"block {self.name!r}: input {self.params['in']!r} has shape {shape}, expected {want} "
                             f"(rf_rate_mode={self._rf}; generator and demodulator must agree)")
        super().bind(layout)


class RfReference(ErrBlock):
    KIND, CONST, STEP = "rf_reference", C_RFR, (_PY[1], _NB[1])

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "out", "band_targets", "toggles", block=name)
        on = bool(params["toggles"].get("phase_noise"))
        self.setup(
            self.KIND, name, params, inputs={},
            outputs={"out": (params["out"], (), "rad"), "walk": (f"rf.{name}.walk", (), "rad")},
            runtime=("s_phi_white_rad2_hz", "s_nu_white_hz2_hz"), toggles=("phase_noise",),
            n_draws=2 if on else 0,
        )

    def bands(self) -> list[BandDecl]:
        if not self.toggles["phase_noise"]:
            return []
        if self._fb is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        s_phi, s_nu = float(self.params["s_phi_white_rad2_hz"]), float(self.params["s_nu_white_hz2_hz"])
        return [BandDecl(source=self.name, quantity="phase", targets=tuple(self.params["band_targets"]),
                         f_lo_hz=0.0, f_hi_hz=0.5 * self._fb, unit="rad",
                         psd_outside=lambda f: s_phi + s_nu / np.asarray(f, dtype=float) ** 2)]


make_demodulator = register("demodulator", Demodulator)
make_rf_reference = register("rf_reference", RfReference)
