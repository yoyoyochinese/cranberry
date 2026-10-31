"""ESA emulation: ``esa_zero_span`` and ``esa_swept_spectrum`` (T04 §5, G16-f). Owner: WP-B5.

Measurement stage (after the time loop); the time-loop kernel is a no-op. Random numbers come from
the ``rng`` argument of ``evaluate`` — pass ``core.rng.stream(run seed, block name)`` (G11).

Zero span (T04 §5.1), complex baseband at f_esa = ``oversample``·ENBW, ENBW = ``enbw_factor``·RBW:

1. z(t): unit complex white Gaussian noise → RBW filter (``rbw_shape``: ``rect`` of width ENBW, or
   ``gaussian`` |H|² = exp(−π f²/ENBW²)), normalized so E|z|² = 1;
2. detected quantity d(t): ``power`` P = (V′(t)/V_SN)|z|², ``log`` 10 log10 P, ``voltage`` √P;
3. one-pole video filter (−3 dB at VBW) on d(t); ``detector`` ``sample`` | ``average`` per display
   point; ``n_traces`` averaged in the detection domain; converted to dB;
4. normalized by the mean of a shot-noise trace (V′ = V_SN) made with the same settings and noise
   model (independent noise).

``radiometer_noise`` off replaces |z|² by its expectation: power E|z|² = 1, log E[10 log10|z|²] =
−10γ/ln 10 = −2.507 dB, voltage E|z| = √π/2 (−1.049 dB) (AN 1303).

Radiometer scatter (analytic, :func:`radiometer_std_db`, VBW ≪ RBW): the video filter's two-sided
noise bandwidth is πVBW, so Var(P̄)/P̄² = πVBW/B_stat with B_stat = (∫|H|²)²/∫|H|⁴ = ENBW (rect) or
√2·ENBW (gaussian). For log averaging the variance is (10/ln 10)²·πVBW·∫Li₂(|ρ(τ)|²)dτ, i.e. the
power-domain value times :data:`LOG_VAR_RATIO` (1.325 rect, 1.341 gaussian), not π²/6 = 1.645: the
π²/6 factor of T04 §5.3 (AN 1303) holds only for independent samples. See WP-B5 report.

Swept spectrum (T04 §5.5): each bin S_k = ∫ |H_RBW(f − f_k)|² S(f) df / ENBW for the signal, LO-only
(shot) and dark PSDs given as callables (e.g. ``bhd_frontend.psd``); radiometer noise per bin with
N_eff = B_stat/(πVBW) (video filter settled in each bin; ponytail: no inter-bin smearing for fast
sweeps — add a swept-time model if sweep_time < n_points/VBW matters). Outputs (a) S/S_shot per bin
and, with ``subtract_circuit_noise``, (b) N = (S − S_dark)/(S_shot − S_dark) (Ref A Fig. 4(b)).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import numba
import numpy as np
import scipy.signal

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.fft import fft_agrawal, ifft_agrawal, omega_grid
from cranberry.core.registry import register_block
from cranberry.core.units import V_SN

__all__ = ["EsaZeroSpan", "EsaSweptSpectrum", "make_esa_zero_span", "make_esa_swept_spectrum", "AVERAGE_TYPES",
           "RBW_SHAPES", "DETECTORS", "LOG_BIAS_DB", "VOLTAGE_BIAS_DB", "LOG_VAR_RATIO", "radiometer_std_db",
           "one_pole"]

AVERAGE_TYPES: tuple[str, ...] = ("power", "log", "voltage")
RBW_SHAPES: tuple[str, ...] = ("rect", "gaussian")
DETECTORS: tuple[str, ...] = ("sample", "average")
_DB = 10.0 / math.log(10.0)
#: E[10 log10 X], X ~ Exp(1): −10γ/ln 10 (AN 1303: log average 2.51 dB low)
LOG_BIAS_DB: float = -_DB * 0.5772156649015329
#: 20 log10 E[√X] = 20 log10(√π/2) (AN 1303: voltage average 1.05 dB low)
VOLTAGE_BIAS_DB: float = 20.0 * math.log10(math.sqrt(math.pi) / 2.0)
#: ∫Li₂(|ρ|²)dτ / ∫|ρ|²dτ for the RBW shapes (computed numerically; checked in the unit tests)
LOG_VAR_RATIO: dict[str, float] = {"rect": 1.3248108, "gaussian": 1.3414873}
#: B_stat / ENBW
_BSTAT: dict[str, float] = {"rect": 1.0, "gaussian": math.sqrt(2.0)}


def radiometer_std_db(rbw_hz: float, vbw_hz: float, enbw_factor: float, rbw_shape: str, average: str) -> float:
    """Analytic standard deviation [dB] of a settled trace of stationary Gaussian noise (module docstring)."""
    rel2 = math.pi * vbw_hz / (_BSTAT[rbw_shape] * enbw_factor * rbw_hz)
    if average == "power":
        return _DB * math.sqrt(rel2)
    if average == "log":
        return _DB * math.sqrt(rel2 * LOG_VAR_RATIO[rbw_shape])
    raise ValueError("radiometer_std_db: average must be 'power' or 'log'")


def one_pole(x: np.ndarray, vbw_hz: float, fs_hz: float, x0: float) -> np.ndarray:
    """One-pole low-pass (−3 dB at VBW) y[n] = y[n−1] + a(x[n] − y[n−1]), a = 1 − e^{−2πVBW/fs}, y[−1] = x0."""
    a = -math.expm1(-2.0 * math.pi * vbw_hz / fs_hz)
    y, _ = scipy.signal.lfilter([a], [1.0, a - 1.0], x, zi=[(1.0 - a) * x0])
    return y


def _check(p: Mapping[str, Any], name: str) -> None:
    for key, allowed in (("average", AVERAGE_TYPES), ("rbw_shape", RBW_SHAPES)):
        if p[key] not in allowed:
            raise ValueError(f"block {name!r}: {key} must be one of {allowed}, got {p[key]!r}")
    for key in ("rbw_hz", "vbw_hz", "enbw_factor"):
        if not float(p[key]) > 0.0:
            raise ValueError(f"block {name!r}: {key} must be > 0")


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class _EsaBase(BaseBlock):
    toggle_names = ("radiometer_noise",)

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    @property
    def enbw_hz(self) -> float:
        return float(self.params["enbw_factor"]) * float(self.params["rbw_hz"])

    def _rbw_power(self, f_hz: np.ndarray) -> np.ndarray:
        """|H_RBW(f)|² with ∫|H|² df = ENBW and |H(0)|² = 1."""
        b = self.enbw_hz
        if self.params["rbw_shape"] == "rect":
            return (np.abs(f_hz) <= 0.5 * b).astype(float)
        return np.exp(-math.pi * (f_hz / b) ** 2)


class EsaZeroSpan(_EsaBase):
    """``esa_zero_span`` (module docstring).

    Params (all required, D34): f_s_hz (centre, recorded), rbw_hz, vbw_hz, enbw_factor, rbw_shape,
    average, detector, oversample, sweep_time_s, n_points, n_traces, toggles {radiometer_noise}.
    """

    KIND = "esa_zero_span"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "f_s_hz", "rbw_hz", "vbw_hz", "enbw_factor", "rbw_shape", "average", "detector",
                "oversample", "sweep_time_s", "n_points", "n_traces", "toggles", block=name)
        super().__init__(self.KIND, name, params)
        _check(params, name)
        if params["detector"] not in DETECTORS:
            raise ValueError(f"block {name!r}: detector must be one of {DETECTORS}")
        self.fs_hz = float(params["oversample"]) * self.enbw_hz
        self.n = int(round(float(params["sweep_time_s"]) * self.fs_hz))
        if not 1 <= int(params["n_points"]) <= self.n:
            raise ValueError(f"block {name!r}: need 1 <= n_points <= samples per sweep ({self.n})")

    def _z2(self, rng: np.random.Generator) -> np.ndarray:
        """|z|² of RBW-filtered unit complex Gaussian noise, E|z|² = 1 (circular filtering)."""
        dt = 1.0 / self.fs_hz
        w = (rng.standard_normal(self.n) + 1j * rng.standard_normal(self.n)) / math.sqrt(2.0)
        h2 = self._rbw_power(omega_grid(self.n, dt) / (2.0 * math.pi))
        h = np.sqrt(h2 / h2.mean())
        return np.abs(ifft_agrawal(fft_agrawal(w, dt) * h, dt)) ** 2

    def _detected(self, p_rel: np.ndarray, rng: np.random.Generator | None) -> np.ndarray:
        avg = self.params["average"]
        if self._toggles["radiometer_noise"]:
            if rng is None:
                raise ValueError(f"block {self.name!r}: radiometer_noise needs an rng (core.rng.stream)")
            z2 = self._z2(rng)
            return {"power": p_rel * z2, "log": _DB * np.log(p_rel * z2), "voltage": np.sqrt(p_rel * z2)}[avg]
        return {"power": p_rel, "log": _DB * np.log(p_rel) + LOG_BIAS_DB,
                "voltage": np.sqrt(p_rel) * math.sqrt(math.pi) / 2.0}[avg]

    def _trace(self, p_rel: np.ndarray, rng: np.random.Generator | None) -> np.ndarray:
        """Display trace (detection domain) of one sweep, averaged over n_traces."""
        npts = int(self.params["n_points"])
        n_tau = min(self.n, max(1, int(self.fs_hz / (2.0 * math.pi * float(self.params["vbw_hz"])))))
        acc = np.zeros(npts)
        for _ in range(int(self.params["n_traces"])):
            d = self._detected(p_rel, rng)
            y = one_pole(d, float(self.params["vbw_hz"]), self.fs_hz, float(d[:n_tau].mean()))
            if self.params["detector"] == "sample":
                acc += y[np.linspace(0, self.n - 1, npts).round().astype(int)]
            else:
                acc += np.array([b.mean() for b in np.array_split(y, npts)])
        return acc / int(self.params["n_traces"])

    def _to_db(self, d: np.ndarray) -> np.ndarray:
        avg = self.params["average"]
        return d if avg == "log" else (_DB * np.log(d) if avg == "power" else 2.0 * _DB * np.log(d))

    def evaluate(self, v_prime: Any, *, dt_s: float | None = None, rng: np.random.Generator | None = None) -> dict:
        """Zero-span trace of V′ (conv units, vacuum 1/2): scalar, or V′(t) sampled every ``dt_s`` from t = 0.

        Returns ``t_s`` (display points), ``trace_db`` (dB re the mean shot-noise trace), ``raw_db``,
        ``shot_db``, ``mean_db`` / ``std_db`` of ``trace_db``.
        """
        t = np.arange(self.n) / self.fs_hz
        v = np.asarray(v_prime, dtype=float)
        if v.ndim == 0:
            p_rel = np.full(self.n, float(v) / V_SN)
        else:
            if dt_s is None:
                raise ValueError("evaluate: V′(t) needs dt_s")
            p_rel = np.interp(t, np.arange(v.size) * dt_s, v) / V_SN
        if not np.all(p_rel > 0.0):
            raise ValueError("evaluate: V′ must be > 0")
        raw = self._to_db(self._trace(p_rel, rng))
        shot = self._to_db(self._trace(np.ones(self.n), rng))
        trace = raw - shot.mean()
        npts = int(self.params["n_points"])
        return {"t_s": np.linspace(0.0, t[-1], npts), "trace_db": trace, "raw_db": raw, "shot_db": shot,
                "mean_db": float(trace.mean()), "std_db": float(trace.std()), "f_s_hz": float(self.params["f_s_hz"])}


class EsaSweptSpectrum(_EsaBase):
    """``esa_swept_spectrum`` (module docstring).

    Params (all required, D34): f_start_hz, f_stop_hz, n_points, rbw_hz, vbw_hz, enbw_factor,
    rbw_shape, average, n_sub (integration points per bin over ±3 ENBW), toggles {radiometer_noise,
    subtract_circuit_noise}.
    """

    KIND = "esa_swept_spectrum"
    toggle_names = ("radiometer_noise", "subtract_circuit_noise")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "f_start_hz", "f_stop_hz", "n_points", "rbw_hz", "vbw_hz", "enbw_factor", "rbw_shape",
                "average", "n_sub", "toggles", block=name)
        super().__init__(self.KIND, name, params)
        _check(params, name)
        self.f_hz = np.linspace(float(params["f_start_hz"]), float(params["f_stop_hz"]), int(params["n_points"]))

    def _bins(self, psd: Callable[[np.ndarray], Any]) -> np.ndarray:
        b = self.enbw_hz
        n_sub = int(self.params["n_sub"])
        df = np.linspace(-3.0 * b, 3.0 * b, n_sub) if n_sub > 1 else np.zeros(1)  # n_sub = 1: no RBW smoothing
        w = self._rbw_power(df)
        w = w / w.sum()
        f = np.abs(self.f_hz[:, None] + df[None, :])  # one-sided PSD: fold negative frequencies
        return (np.broadcast_to(np.asarray(psd(f), dtype=float), f.shape) * w).sum(axis=1)

    def _measure(self, s: np.ndarray, rng: np.random.Generator | None) -> np.ndarray:
        """Displayed value per bin, back in linear power units (bias of the average type removed)."""
        if not self._toggles["radiometer_noise"]:
            return s
        if rng is None:
            raise ValueError(f"block {self.name!r}: radiometer_noise needs an rng (core.rng.stream)")
        p = self.params
        if p["average"] == "power":
            k = _BSTAT[p["rbw_shape"]] * self.enbw_hz / (math.pi * float(p["vbw_hz"]))
            return s * rng.gamma(k, 1.0 / k, s.shape)
        sd = radiometer_std_db(float(p["rbw_hz"]), float(p["vbw_hz"]), float(p["enbw_factor"]), p["rbw_shape"],
                               "log" if p["average"] == "log" else "power")
        # ponytail: voltage average uses the power-domain scatter (Gaussian in dB)
        return s * 10.0 ** (sd * rng.standard_normal(s.shape) / 10.0)

    def evaluate(self, signal_psd: Callable, shot_psd: Callable, dark_psd: Callable,
                 *, rng: np.random.Generator | None = None) -> dict:
        """Swept spectra of the signal / LO-only / dark PSDs (callables of f [Hz], any common unit).

        Returns ``f_hz``, ``signal_db`` / ``shot_db`` / ``dark_db`` (10 log10 of the bin powers),
        ``normalized_db`` = signal − shot (a) and ``corrected_db`` (b, None unless
        ``subtract_circuit_noise``; NaN where radiometer noise makes S − S_dark ≤ 0, as on a real ESA).
        """
        s, sh, dk = (self._measure(self._bins(fn), rng) for fn in (signal_psd, shot_psd, dark_psd))

        def db(x: np.ndarray) -> np.ndarray:
            with np.errstate(invalid="ignore", divide="ignore"):
                return _DB * np.log(np.where(x > 0, x, np.nan))

        out = {"f_hz": self.f_hz, "signal_db": db(s), "shot_db": db(sh), "dark_db": db(dk) if np.all(dk > 0) else None,
               "normalized_db": db(s / sh), "corrected_db": None}
        if self._toggles["subtract_circuit_noise"]:
            out["corrected_db"] = db((s - dk) / (sh - dk))
        return out


@register_block(EsaZeroSpan.KIND)
def make_esa_zero_span(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``esa_zero_span`` (module docstring)."""
    return EsaZeroSpan(params, name)


@register_block(EsaSweptSpectrum.KIND)
def make_esa_swept_spectrum(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``esa_swept_spectrum`` (module docstring)."""
    return EsaSweptSpectrum(params, name)
