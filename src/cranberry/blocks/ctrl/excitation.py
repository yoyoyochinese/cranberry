"""``excitation``: broadband / swept excitation added to a loop signal while the loop stays locked (Phase 6
P7, D85-7; notes/lit-kashiwazaki2023.md §5.3 N7-1). CTRL block of the time loop.

    out = in + x(t)          (``in`` may be null: out = x)

One run covers every frequency (the loop is never reopened), unlike ``sine_injection`` +
``SimController.measure_transfer_function`` (one run from the initial state per frequency).
Extra slots: ``ctrl.<name>.x`` (x itself, the instrument of the estimate), ``ctrl.<name>.f_hz`` (current
frequency: the step frequency of a stepped sine, the instantaneous frequency of a chirp, 0 otherwise) and
``ctrl.<name>.step`` (int64: the step index of a stepped sine, the period index of a periodic waveform, 0
for white noise; −1 before ``t_start_s`` or with the toggle off; n_steps after the last step).

Time origin: the excitation starts at base sample n₀ = ⌈t_start_s·f_sim⌉ (τ = (n − n₀)/f_sim). Before it
x = 0 (lock acquisition).

**Amplitude convention.** ``amplitude`` (runtime) is √2 × RMS of x, i.e. the amplitude of a sine of the
same power, for every waveform: a sine of amplitude A; a multisine of N tones, each of amplitude A/√N; a
PRBS of levels ±A/√2; white Gaussian noise of standard deviation A/√2. Equal ``amplitude`` = equal
excitation power, whatever the waveform. (The peak of a multisine is crest × A/√2, crest ≈ 1.5–2 with
Schroeder phases; check it against the fringe range.)

Waveforms (``waveform``, structural):

* ``stepped_sine``: step k lasts ``dwell_s`` (settle included); x = A sin(2π f_k τ_k), τ_k from the step
  start. Each f_k is snapped to the grid 1/T_meas, T_meas = dwell_s − settle_s, so the measurement window
  [settle_s, dwell_s) of every step holds an integer number of periods (leakage-free lock-in when T_meas
  is a multiple of the record step). Duplicates after snapping are dropped. After the last step x = 0.
* ``multisine``: tones on the grid 1/``period_s`` (snapped, duplicates dropped), Schroeder phases
  φ_k = −π k(k − 1)/N (k = 1..N; Schroeder 1970, low crest factor for a flat amplitude spectrum).
  Periodic; one period is tabulated at the block rate (period_s·f_block must be an integer).
* ``chirp``: periodic sweep f_start → f_stop over ``period_s`` (``spacing`` lin: f = f₀ + (f₁ − f₀)τ/T;
  log: f = f₀ (f₁/f₀)^{τ/T}), x = A sin(phase), phase = 2π∫f dτ from the period start.
* ``prbs``: maximal-length LFSR sequence of order ``prbs_order`` (period 2ⁿ − 1 bits), each bit held
  ``bit_s`` (bit_s·f_block must be an integer); levels ±A/√2. Tabulated.
* ``white``: Gaussian, standard deviation A/√2, one draw per step of the block (``DrawBuffer``, G11).

Frequency list (``stepped_sine`` / ``multisine``): either ``freqs_hz`` (list) or ``f_start_hz``,
``f_stop_hz``, ``n_freqs``, ``spacing`` (``log`` | ``lin``), and ``amp_f_corner_hz``: step / tone k gets
the amplitude × g_k, g_k = max(1, f_k/amp_f_corner_hz) (``.inf``: flat). Above the loop's unity-gain
frequency |T| ≈ |L| falls about as 1/f, so a corner near f_u keeps the response (−T·x) above the loop
noise while the phase excursion at the fringe (≈ |T|·|x|/k_d) stays about constant. With shaping the
amplitude convention above holds for the flat part (g_k = 1). Every other key is required only by the
waveform that uses it (D34: no defaults). ``settle_s`` is required for every waveform: the analysis
(``loop_identification``) discards [t_start, t_start + settle_s) and, for a stepped sine, the first
``settle_s`` of every step.

Params: ``in``, ``out``, ``unit``, ``waveform``, ``settle_s`` and the waveform keys above; runtime
``amplitude``, ``t_start_s``; toggle ``enabled`` (False: out = in, x = 0, step = −1, no draws).
:func:`excitation_schedule` returns the frequencies / timing the analysis needs (same snapping).
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
from cranberry.core.rng import DrawBuffer

__all__ = ["Excitation", "make_excitation", "EXC_WAVEFORMS", "excitation_schedule", "prbs_bits", "schroeder_phases"]

EXC_WAVEFORMS: dict[str, int] = {"stepped_sine": 0, "multisine": 1, "chirp": 2, "prbs": 3, "white": 4}
_SPACINGS = ("log", "lin")
_WAVE_KEYS: dict[str, tuple[str, ...]] = {
    "stepped_sine": ("dwell_s", "amp_f_corner_hz"),
    "multisine": ("period_s", "amp_f_corner_hz"),
    "chirp": ("f_start_hz", "f_stop_hz", "spacing", "period_s"),
    "prbs": ("prbs_order", "bit_s"),
    "white": (),
}
#: Fibonacci LFSR feedback taps (1-based bit positions) of maximal-length sequences, order 2–16
#: (Xilinx XAPP052 table).
_PRBS_TAPS: dict[int, tuple[int, ...]] = {
    2: (2, 1), 3: (3, 2), 4: (4, 3), 5: (5, 3), 6: (6, 5), 7: (7, 6), 8: (8, 6, 5, 4), 9: (9, 5), 10: (10, 7),
    11: (11, 9), 12: (12, 6, 4, 1), 13: (13, 4, 3, 1), 14: (14, 5, 3, 1), 15: (15, 14), 16: (16, 15, 13, 4),
}


def prbs_bits(order: int) -> np.ndarray:
    """One period (2ⁿ − 1 values) of the maximal-length sequence of ``order`` as ±1."""
    if order not in _PRBS_TAPS:
        raise ValueError(f"prbs_order must be in {sorted(_PRBS_TAPS)}, got {order!r}")
    n = (1 << order) - 1
    reg = [1] * order
    out = np.empty(n)
    for i in range(n):
        out[i] = 1.0 if reg[-1] else -1.0
        fb = 0
        for tp in _PRBS_TAPS[order]:
            fb ^= reg[tp - 1]
        reg = [fb] + reg[:-1]
    return out


def schroeder_phases(n: int) -> np.ndarray:
    """φ_k = −π k(k − 1)/N, k = 1..N."""
    k = np.arange(1, n + 1, dtype=float)
    return -np.pi * k * (k - 1) / n


def _freq_list(p: Mapping[str, Any], where: str) -> np.ndarray:
    if "freqs_hz" in p and p["freqs_hz"] is not None:
        f = np.asarray(p["freqs_hz"], dtype=float).ravel()
    else:
        require(p, "f_start_hz", "f_stop_hz", "n_freqs", "spacing", block=where)
        if p["spacing"] not in _SPACINGS:
            raise ValueError(f"block {where!r}: spacing must be one of {_SPACINGS}, got {p['spacing']!r}")
        lo, hi, n = float(p["f_start_hz"]), float(p["f_stop_hz"]), int(p["n_freqs"])
        f = np.geomspace(lo, hi, n) if p["spacing"] == "log" else np.linspace(lo, hi, n)
    if f.size == 0 or not np.all(f > 0):
        raise ValueError(f"block {where!r}: frequencies must be > 0, got {f!r}")
    return f


def _snap(f: np.ndarray, grid_s: float) -> np.ndarray:
    m = np.maximum(np.round(f * grid_s), 1.0)
    _, first = np.unique(m, return_index=True)
    return m[np.sort(first)] / grid_s


def _gains(f: np.ndarray, p: Mapping[str, Any], where: str) -> np.ndarray:
    fc = float(p["amp_f_corner_hz"])
    if not fc > 0:
        raise ValueError(f"block {where!r}: amp_f_corner_hz must be > 0 (.inf: flat)")
    return np.maximum(1.0, f / fc)


def excitation_schedule(params: Mapping[str, Any], where: str = "excitation") -> dict[str, Any]:
    """Waveform, frequencies (snapped as the block does) and timing of an ``excitation`` (module
    docstring). Keys: ``waveform``, ``settle_s``, ``t_start_s``, ``amplitude`` and per waveform
    ``freqs_hz`` + ``gains`` + ``dwell_s`` (stepped_sine), ``freqs_hz`` + ``gains`` + ``period_s`` (multisine),
    ``f_start_hz`` /
    ``f_stop_hz`` / ``spacing`` / ``period_s`` (chirp), ``prbs_order`` / ``bit_s`` / ``period_s`` (prbs)."""
    require(params, "waveform", "settle_s", block=where)
    w = params["waveform"]
    if w not in EXC_WAVEFORMS:
        raise ValueError(f"block {where!r}: waveform must be one of {sorted(EXC_WAVEFORMS)}, got {w!r}")
    require(params, *_WAVE_KEYS[w], block=where)
    settle = float(params["settle_s"])
    if settle < 0:
        raise ValueError(f"block {where!r}: settle_s must be >= 0")
    out: dict[str, Any] = {"waveform": w, "settle_s": settle, "t_start_s": float(params.get("t_start_s", 0.0)),
                           "amplitude": float(params.get("amplitude", 0.0))}
    if w == "stepped_sine":
        dwell = float(params["dwell_s"])
        if not dwell > settle:
            raise ValueError(f"block {where!r}: dwell_s must exceed settle_s (measurement window > 0)")
        out.update(dwell_s=dwell, freqs_hz=_snap(_freq_list(params, where), dwell - settle))
        out["gains"] = _gains(out["freqs_hz"], params, where)
    elif w == "multisine":
        per = float(params["period_s"])
        if not per > 0:
            raise ValueError(f"block {where!r}: period_s must be > 0")
        out.update(period_s=per, freqs_hz=_snap(_freq_list(params, where), per))
        out["gains"] = _gains(out["freqs_hz"], params, where)
    elif w == "chirp":
        f0, f1, per = float(params["f_start_hz"]), float(params["f_stop_hz"]), float(params["period_s"])
        if not (f0 > 0 and f1 > 0 and per > 0) or params["spacing"] not in _SPACINGS:
            raise ValueError(f"block {where!r}: chirp needs f_start_hz, f_stop_hz, period_s > 0 and spacing in {_SPACINGS}")
        if params["spacing"] == "log" and f0 == f1:
            raise ValueError(f"block {where!r}: a log chirp needs f_start_hz != f_stop_hz")
        out.update(f_start_hz=f0, f_stop_hz=f1, spacing=params["spacing"], period_s=per)
    elif w == "prbs":
        order, bit = int(params["prbs_order"]), float(params["bit_s"])
        if order not in _PRBS_TAPS or not bit > 0:
            raise ValueError(f"block {where!r}: prbs_order in {sorted(_PRBS_TAPS)} and bit_s > 0 required")
        out.update(prbs_order=order, bit_s=bit, period_s=((1 << order) - 1) * bit)
    return out


def _samples(x_s: float, rate_hz: float, what: str, where: str) -> int:
    n = x_s * rate_hz
    k = int(round(n))
    if k < 1 or abs(n - k) > 1e-6 * max(1.0, n):
        raise ValueError(f"block {where!r}: {what} × block rate = {n:.9g} must be a positive integer")
    return k


def _exc_step(S, I, P, i_n, i_in, i_out, i_x, i_f, i_step, pa, pts, wave, fs, rdiv,  # noqa: N803
              freqs, gains, dwell_n, table, f0, f1, per_n, logsp, buf, idx, enabled):
    y = S[i_in] if i_in >= 0 else 0.0
    w = 0.0
    if wave == 4 and enabled:
        w = S[buf + I[idx]]
        I[idx] += 1
    x = 0.0
    f = 0.0
    step = -1
    if enabled:
        n0 = int(math.ceil(P[pts] * fs - 1e-6))
        rel = I[i_n] - n0
        if rel >= 0:
            a = P[pa]
            if wave == 0:
                k = rel // dwell_n
                if k < freqs.shape[0]:
                    f = freqs[k]
                    x = a * gains[k] * math.sin(2.0 * math.pi * f * ((rel - k * dwell_n) / fs))
                step = k if k < freqs.shape[0] else freqs.shape[0]
            elif wave == 1 or wave == 3:
                j = rel // rdiv
                x = a * table[j % table.shape[0]]
                step = j // table.shape[0]
            elif wave == 2:
                step = rel // per_n
                tau = (rel - step * per_n) / fs
                tp = per_n / fs
                if logsp:
                    r = f1 / f0
                    q = r ** (tau / tp)
                    f = f0 * q
                    ph = 2.0 * math.pi * f0 * tp * (q - 1.0) / math.log(r)
                else:
                    f = f0 + (f1 - f0) * tau / tp
                    ph = 2.0 * math.pi * (f0 * tau + 0.5 * (f1 - f0) * tau * tau / tp)
                x = a * math.sin(ph)
            else:
                x = a * w * 0.7071067811865476
                step = 0
    S[i_out] = y + x
    S[i_x] = x
    S[i_f] = f
    I[i_step] = step


class Excitation(KitBlock):
    """The ``excitation`` block (module docstring)."""

    KIND = "excitation"
    STEP = staticmethod(_exc_step)
    STEP_NB = numba.njit(_exc_step)
    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "out", "unit", "waveform", "settle_s", block=name)
        self.schedule = excitation_schedule(params, name)
        self.x_name, self.f_name, self.step_name = (f"ctrl.{name}.{k}" for k in ("x", "f_hz", "step"))
        reads = ("n",) if params["in"] is None else ("n", params["in"])
        super().__init__(params, name, reads=reads, writes=(params["out"], self.x_name, self.f_name, self.step_name),
                         runtime=("amplitude", "t_start_s"))
        self.draws = DrawBuffer(self.name, 1, "standard_normal")
        if self._white():
            self.writes = (*self.writes, *self.draws.slot_names)
        u = str(params["unit"])
        self.slot_units = {params["out"]: u, self.x_name: u}

    def _white(self) -> bool:
        return self.schedule["waveform"] == "white" and self.toggles["enabled"]

    def slots(self):
        u = str(self.params["unit"])
        return [(self.params["out"], (), np.float64, u, 0, 0.0), (self.x_name, (), np.float64, u, 0, 0.0),
                (self.f_name, (), np.float64, "Hz", 0, 0.0), (self.step_name, (), np.int64, "1", 0, -1)]

    def init(self, layout, rng):
        super().init(layout, rng)
        if self._white():
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self._white():
            self.draws.refill(state, n0, n_samples)

    def table(self, f_block_hz: float) -> np.ndarray:
        """One period of the tabulated waveform (multisine, prbs) at the block rate, unit amplitude."""
        sc = self.schedule
        if sc["waveform"] == "multisine":
            n_tab = _samples(sc["period_s"], f_block_hz, "period_s", self.name)
            m = np.round(sc["freqs_hz"] * sc["period_s"])
            if np.any(m >= n_tab / 2):
                raise ValueError(f"block {self.name!r}: multisine tone at or above the block Nyquist frequency")
            j = np.arange(n_tab)[:, None]
            ph = schroeder_phases(m.size)[None, :]
            g = np.asarray(sc["gains"])[None, :]
            return (g * np.cos(2 * np.pi * m[None, :] * j / n_tab + ph)).sum(axis=1) / math.sqrt(m.size)
        if sc["waveform"] == "prbs":
            nb = _samples(sc["bit_s"], f_block_hz, "bit_s", self.name)
            return np.repeat(prbs_bits(sc["prbs_order"]), nb) * (1.0 / math.sqrt(2.0))
        return np.zeros(1)

    def make_args(self, layout):
        if self._white():
            self.draws.bind(layout)
        sc, fs = self.schedule, float(layout.f_sim_hz)
        f_blk = fs / self.rate_div
        freqs, gains, dwell_n, per_n = np.zeros(1), np.ones(1), 1, 1
        if sc["waveform"] == "stepped_sine":
            freqs, gains = np.asarray(sc["freqs_hz"], dtype=float), np.asarray(sc["gains"], dtype=float)
            dwell_n = _samples(sc["dwell_s"], fs, "dwell_s", self.name)
            if np.any(freqs >= f_blk / 2):
                raise ValueError(f"block {self.name!r}: stepped-sine frequency at or above the block Nyquist frequency")
        if sc["waveform"] == "chirp":
            per_n = _samples(sc["period_s"], fs, "period_s", self.name)
        return (self.o(layout, "n"), self.o(layout, self.params["in"]), self.o(layout, self.params["out"]),
                self.o(layout, self.x_name), self.o(layout, self.f_name), self.o(layout, self.step_name),
                self.p(layout, "amplitude"), self.p(layout, "t_start_s"), EXC_WAVEFORMS[sc["waveform"]], fs,
                int(self.rate_div), freqs, gains, dwell_n, self.table(f_blk), float(sc.get("f_start_hz", 1.0)),
                float(sc.get("f_stop_hz", 1.0)), per_n, sc.get("spacing") == "log", self.draws.buf_off,
                self.draws.idx_off, self.toggles["enabled"])


@register_block("excitation")
def make_excitation(params: Mapping[str, Any], *, name: str | None = None) -> Excitation:
    """Factory of the ``excitation`` block (module docstring)."""
    return Excitation(params, name)
