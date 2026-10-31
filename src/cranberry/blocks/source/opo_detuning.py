"""``opo_detuning_noise`` and ``opo_length_lock``: OPO cavity detuning Δ(t) in the time loop (T02 §2.4,
§6). Owner: WP-B6.

Both are time-loop (T) blocks writing a detuning [rad/s] (conv §4 / T02 §1: Δ = ω_cav − ω0). The usual
chain is ``opo_detuning_noise`` (free-running Δ_free) → ``opo_length_lock`` (residual Δ_s) →
``drv.<opo segment>.detuning`` → ``segment_phase``, whose OPO segment × beam coefficients
(``opo_cw_spectrum.coef_entries``: squeezed −1/(κ(1 + x²)), probe frame −1/(2κ), probe transit −1/κ)
turn Δ_s into the squeezing-angle and probe-phase rotations; the lock combination then carries the
residual (Δ/κ)(1 − x²)/(2(1 + x²)) (VT-23). Without a length lock, wire ``opo_detuning_noise.out`` to
the drv slot directly.

Δ = ω_cav − ω_L: the laser frequency is a detuning too (Phase 5 M3B-02). A preset adds −2π·``laser.<s>.dnu``
to Δ_free (an ``adder`` before ``opo_length_lock``), so a lock whose light follows the laser (r = 1)
cancels it below f_ugf only; the part above the block Nyquist frequency belongs to the laser's ledger
pair (θ_sq lags φ_L by 1/(κ(1 + x²)), refB.yaml header), not to this block's ``detuning`` row.

``opo_detuning_noise`` (T02 §2.4, §5)
    Δ_free = −(ω0/l)·δl + (∂Δ/∂T)·ΔT,  ω0 = 2πν0, l = round-trip optical path length.

    * δl [m]: round-trip optical path length noise, a two-pole h/f² section (``length_noise``:
      ``h_m2_hz``, ``f_lo_hz``, ``f_c_hz``; S_δl ≈ h/f² for f_lo ≪ f ≪ f_c, acoustic / mount noise) plus
      band-pass peaks (``peaks``: [{f_hz, q, area_m2}]), exact discretisation of ``noise._lgs``, one RNG
      sub-stream per toggle (G11). Sign: δl > 0 lowers the resonance, Δ < 0 (T02 §2.4).
    * ΔT [K]: a temperature slot written by another block (``thermal``: ``temperature_slot``,
      ``ddelta_dT_rad_s_per_k`` = −2πν0(α_L + (1/n)dn/dT) for a monolithic cavity, Ref C; T02 §2.4).

    Toggles {length_noise, peaks, thermal}; all off → Δ_free = 0. Ledger (G05): the δl PSD above the
    block Nyquist frequency, as quantity ``detuning`` [rad/s] on ``segment`` (the thermal part is
    declared by the block that writes the temperature).

``opo_length_lock`` (effective loop; Ref B does not document its OPO length lock, T02 §6.1)
    The lock light sees Δ_lock = r·Δ_free + a (r = ``lock_light_ratio``: 1 for a lock beam in the same
    mode, ≠ 1 when it sits at another frequency / polarization, e.g. a Borealis-type lock beam with a
    different dn/dT, T17). Per step h (pure integrator, unity-gain frequency f_ugf)::

        e   = r·Δ_free + a + w − Δ_off            (w: sensor noise, one-sided PSD S_w [(rad/s)²/Hz])
        Δ_s = Δ_free + a                          → ``out``
        a  ← a − 2π f_ugf h · e

    so |Δ_s/Δ_free| ≈ f/f_ugf below f_ugf and Δ_s → (1 − r)Δ_free + Δ_off in steady state.
    # ponytail: one integrator stands in for sensor + servo + PZT; compose the ERR / CTRL / ACT blocks
    # (B3 / B4) when a documented OPO lock (e.g. Liu PDH, Borealis lock beam) needs its real dynamics.
    Runtime params ``f_ugf_hz``, ``lock_light_ratio``, ``offset_rad_s``, ``sensor_psd``. Slots:
    ``out`` (Δ_s), ``opo.<name>.act`` (a), ``opo.<name>.err`` (e). Toggles {lock (off → Δ_s = Δ_free,
    identity), sensor_noise}. No ledger entry (in-loop sensor noise, cf. D56-5).
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.noise._lgs import EMPTY_BANK, Bank, bank_step, bank_step_nb, resonance, two_pole
from cranberry.blocks.noise.laser import resonance_psd
from cranberry.blocks.noise.mech import two_pole_psd
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer
from cranberry.offline.band_ledger import BandDecl

__all__ = ["OpoDetuningNoise", "OpoLengthLock", "make_opo_detuning_noise", "make_opo_length_lock"]

NC = namedtuple("NC", "out th temp kth b1 t1 b2 t2")
LC = namedtuple("LC", "inp out act err lock sens buf idx h pf pr po ps")


def _build(jit: bool):
    bs = bank_step_nb if jit else bank_step

    def noise_step(S, I, c):  # noqa: N803
        S[c.out] = c.kth * S[c.temp] if c.th else 0.0
        if c.b1.on:
            bs(S, I, c.b1, c.t1)
        if c.b2.on:
            bs(S, I, c.b2, c.t2)

    def lock_step(S, I, P, c):  # noqa: N803
        d = S[c.inp]
        if not c.lock:
            S[c.out] = d
            return
        a = S[c.act]
        w = 0.0
        if c.sens:
            w = math.sqrt(0.5 * P[c.ps] / c.h) * S[c.buf + I[c.idx]]
            I[c.idx] += 1
        e = P[c.pr] * d + a + w - P[c.po]
        S[c.out] = d + a
        S[c.err] = e
        S[c.act] = a - 2.0 * math.pi * P[c.pf] * c.h * e

    if jit:
        return numba.njit(noise_step), numba.njit(lock_step)
    return noise_step, lock_step


(_NOISE_PY, _LOCK_PY), (_NOISE_NB, _LOCK_NB) = _build(False), _build(True)


class OpoDetuningNoise(BaseBlock):
    """``opo_detuning_noise`` (module docstring)."""

    toggle_names = ("length_noise", "peaks", "thermal")
    KIND = "opo_detuning_noise"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "out", "segment", "nu0_hz", "round_trip_m", "toggles", block=name)
        tg = params["toggles"]
        if tg.get("length_noise"):
            require(params, "h_m2_hz", "f_lo_hz", "f_c_hz", block=name)
        if tg.get("peaks"):
            require(params, "peaks", block=name)
        reads: tuple[str, ...] = ()
        if tg.get("thermal"):
            require(params, "temperature_slot", "ddelta_dT_rad_s_per_k", block=name)
            reads = (str(params["temperature_slot"]),)
        if not (float(params["nu0_hz"]) > 0 and float(params["round_trip_m"]) > 0):
            raise ValueError(f"block {name!r}: need nu0_hz > 0 and round_trip_m > 0")
        super().__init__(self.KIND, name, params, reads=reads, writes=(str(params["out"]),))
        self.slot_units = {str(params["out"]): "rad/s"}
        if reads:
            self.slot_units[reads[0]] = "K"
        self.gain = -2.0 * math.pi * float(params["nu0_hz"]) / float(params["round_trip_m"])  # Δ per δl
        self.banks: dict[str, Bank] = {}
        self._h: float | None = None

    def _sections(self) -> dict[str, list]:
        p, tg, out = self.params, self._toggles, {}
        if tg["length_noise"]:
            out["length_noise"] = [two_pole(float(p["h_m2_hz"]), float(p["f_lo_hz"]), float(p["f_c_hz"]))]
        if tg["peaks"] and p["peaks"]:
            out["peaks"] = [resonance(float(k["f_hz"]), float(k["q"]), float(k["area_m2"])) for k in p["peaks"]]
        return out

    def init(self, layout, rng):
        super().init(layout, rng)
        self._h = h = block_dt_s(layout.f_sim_hz, self.rate_div)
        out = str(self.params["out"])
        layout.add(out, (), np.float64, "rad/s")
        extra: list[str] = []
        for label, secs in self._sections().items():
            b = self.banks.setdefault(label, Bank(self.name, label))
            for sec in secs:
                b.add(sec, [(out, 0, self.gain)])
            extra += b.init(layout, rng, self.rate_div, h)
        self.writes = (*self.writes, *extra)

    def bind(self, layout):
        th = self._toggles["thermal"]
        bk = [(self.banks[k].bind(layout), self.banks[k].scratch_off(layout)) if k in self.banks else (EMPTY_BANK, -1)
              for k in ("length_noise", "peaks")]
        self._c = NC(layout.slot(str(self.params["out"])).offset, th,
                     layout.slot(str(self.params["temperature_slot"])).offset if th else 0,
                     float(self.params["ddelta_dT_rad_s_per_k"]) if th else 0.0, bk[0][0], bk[0][1], bk[1][0], bk[1][1])
        c, nb = self._c, _NOISE_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        for b in self.banks.values():
            b.refill(state, n0, n_samples)

    def __call__(self, state, params):
        _NOISE_PY(state.S, state.I, self._c)
        return state

    def psd_delta(self, f) -> np.ndarray:
        """One-sided PSD of the length-noise part of Δ_free [(rad/s)²/Hz] (analytic)."""
        p, tg = self.params, self._toggles
        f = np.asarray(f, dtype=float)
        s = np.zeros_like(f)
        if tg["length_noise"]:
            s = s + two_pole_psd(f, float(p["h_m2_hz"]), float(p["f_lo_hz"]), float(p["f_c_hz"]))
        if tg["peaks"]:
            for k in p["peaks"]:
                s = s + resonance_psd(f, float(k["f_hz"]), float(k["q"]), float(k["area_m2"]))
        return self.gain**2 * s

    def bands(self) -> list[BandDecl]:
        if self._h is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        on = self._toggles["length_noise"] or self._toggles["peaks"]
        return [BandDecl(self.name, "detuning", (str(self.params["segment"]),), 0.0, 0.5 / self._h, "rad/s",
                         psd_outside=self.psd_delta if on else None, variance_outside=None if on else 0.0)]


class OpoLengthLock(BaseBlock):
    """``opo_length_lock`` (module docstring)."""

    toggle_names = ("lock", "sensor_noise")
    KIND = "opo_length_lock"
    RUNTIME = ("f_ugf_hz", "lock_light_ratio", "offset_rad_s", "sensor_psd")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "out", "toggles", block=name)
        self.act, self.err = f"opo.{name}.act", f"opo.{name}.err"
        super().__init__(self.KIND, name, params, reads=(str(params["in"]),),
                         writes=(str(params["out"]), self.act, self.err), runtime_params=self.RUNTIME)
        if float(params["sensor_psd"]) < 0 or float(params["f_ugf_hz"]) < 0:
            raise ValueError(f"block {name!r}: need sensor_psd >= 0 and f_ugf_hz >= 0")
        self.slot_units = {str(params["in"]): "rad/s", str(params["out"]): "rad/s"}
        self.draws: DrawBuffer | None = None

    def init(self, layout, rng):
        super().init(layout, rng)
        self._h = block_dt_s(layout.f_sim_hz, self.rate_div)
        if 2.0 * math.pi * float(self.params["f_ugf_hz"]) * self._h >= 1.0:
            raise ValueError(f"block {self.name!r}: f_ugf_hz too high for the block rate (2π f_ugf h >= 1)")
        for s in (str(self.params["out"]), self.act, self.err):
            layout.add(s, (), np.float64, "rad/s")
        if self._toggles["lock"] and self._toggles["sensor_noise"]:
            self.draws = DrawBuffer(self.name, 1, "standard_normal")
            self.draws.init(layout, rng, self.rate_div)
            self.writes = (*self.writes, *self.draws.slot_names)

    def bind(self, layout):
        if self.draws is not None:
            self.draws.bind(layout)
        o = lambda s: layout.slot(s).offset  # noqa: E731
        q = lambda k: layout.param(self.param_name(k)).offset  # noqa: E731
        d = self.draws
        self._c = LC(o(str(self.params["in"])), o(str(self.params["out"])), o(self.act), o(self.err),
                     self._toggles["lock"], d is not None, d.buf_off if d else -1, d.idx_off if d else -1, self._h,
                     q("f_ugf_hz"), q("lock_light_ratio"), q("offset_rad_s"), q("sensor_psd"))
        c, nb = self._c, _LOCK_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, P, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        if self.draws is not None:
            self.draws.refill(state, n0, n_samples)

    def __call__(self, state, params):
        _LOCK_PY(state.S, state.I, params.P, self._c)
        return state


@register_block(OpoDetuningNoise.KIND)
def make_opo_detuning_noise(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opo_detuning_noise`` block (module docstring)."""
    return OpoDetuningNoise(params, name)


@register_block(OpoLengthLock.KIND)
def make_opo_length_lock(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opo_length_lock`` block (module docstring)."""
    return OpoLengthLock(params, name)
