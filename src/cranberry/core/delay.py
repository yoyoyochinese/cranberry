"""Delay buffers: integer-sample ring + fractional part (protected list; PLAN §2.6, G09). Owner: WP-F1.

The delay buffer is the only legal feedback path in a pipeline: a controller writes ``ctrl.<x>.u``,
a ``delay_buffer`` block delays it by τ, and the actuator reads the delayed slot. Loop delay is the
sum of its parts (controller pure delay from the device preset + ZOH h/2 when not discretised exactly
+ optical propagation n_g·L/c from actuator to detector + optional extra; PLAN §2.6).

**Loop-delay accounting.** A feedback loop closes where a block reads a slot written later in the
execution order (the actuator output read by the plant / segment drivers): that read sees the
previous sample, i.e. the loop closure itself contributes exactly one base sample 1/f_sim (PLAN
§2.5: "actuators act on the next sample"). The total modelled loop delay is therefore
``tau_s + 1/f_sim`` (+ the delays of the blocks' own discretisations); a preset wanting a total τ_L
sets ``tau_s = τ_L − 1/f_sim``. ``core.pipeline.check_pipeline`` enforces the "exactly once": a loop
that can only close through two or more back edges (a ``delay_buffer`` ordered before the controller
it delays) would add k/f_sim with k ≥ 2 and is rejected at build time (Phase 5 F1-1). Blocks held at
``rate_div = r`` add (r − 1)/(2 f_sim) more (``core.rates``, "Hold delay").

Delay arithmetic. With the block's own step h = rate_div / f_sim, D = τ/h samples:

* D integer (within 1e−9): pure ring read y[n] = x[n − D] for either method;
* ``"linear"``  : D = k + μ, μ ∈ (0, 1): y[n] = (1 − μ) x[n − k] + μ x[n − k − 1];
* ``"thiran1"`` : first-order Thiran all-pass after an integer delay k:
  D = k + d with d ∈ [0.5, 1.5) (d = D for D < 0.5), a = (1 − d)/(1 + d),
  H(z) = z^−k (a + z^−1)/(1 + a z^−1), i.e. v[n] = a u[n] + u[n−1] − a v[n−1] with u[n] = x[n−k].
  Unit magnitude; maximally flat group delay d at DC (Thiran 1971; Välimäki & Laakso 2000).
  **Use thiran1 only for D ≥ 0.5** and when keeping the magnitude at 1 is the point: for D < 0.5 the
  pole −a → −1 (D = 0.01: a = 0.980, the impulse response rings at Nyquist for ~35 samples) and its
  phase error in the loop band exceeds linear's. ``"linear"`` is the default choice (both reference
  presets use it; D57-7 / D65-4). A thiran1 design with 0 < D < 0.5 emits a :class:`ChkWarning`
  (code ``thiran1_short_delay``; Phase 5 F1-6).

Ring semantics (also used for ``laser.*`` / ``drv.*`` emission-time histories, WP-F6): the head is
the row of the newest sample; a push advances the head (mod ring_len) and writes that row; the value
d samples older is row ``(head − d) % ring_len``. The ring starts filled with zeros (zero history).
:func:`ring_push` / :func:`ring_at` are pure-Python functions; :data:`ring_push_nb` /
:data:`ring_at_nb` are the same source compiled with numba for kernels.

Block kind ``delay_buffer`` (factory :func:`make_delay_buffer`), params (all required, D34):

    in      : str    slot name to delay (a scalar float64 slot)
    out     : str    slot name written (a scalar float64 slot registered by this block; its unit
                     tag is copied from ``in`` unless ``unit`` is given)
    tau_s   : float  delay [s], >= 0 (structural: it sets the ring length)
    method  : "linear" | "thiran1"
    unit    : str    optional unit tag of ``out`` (needed only if ``in`` is registered later)

Internal slots: ``delay.<name>.ring`` (ring of length k + 2) and ``delay.<name>.ap`` (Thiran state).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.pipeline import DELAY_BUFFER_KIND
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block

__all__ = [
    "DELAY_METHODS",
    "DelayDesign",
    "design_delay",
    "interp_variance_deficit",
    "ring_push",
    "ring_at",
    "ring_push_nb",
    "ring_at_nb",
    "delay_step",
    "delay_step_nb",
    "DelayLine",
    "DelayBuffer",
    "make_delay_buffer",
]

DELAY_METHODS: tuple[str, ...] = ("linear", "thiran1")

_MODE_INT, _MODE_LINEAR, _MODE_THIRAN = 0, 1, 2
_INT_TOL = 1e-9


class DelayDesign:
    """Integer/fractional decomposition of a delay of ``delay_samples`` (see module docstring).

    Attributes:
        delay_samples: D = τ / h.
        method: "linear" | "thiran1".
        mode: 0 integer, 1 linear, 2 Thiran.
        k: integer part read from the ring.
        mu: linear-interpolation weight (mode 1), else 0.
        a: Thiran coefficient (mode 2), else 0.
        ring_len: k + 2.
    """

    def __init__(self, delay_samples: float, method: str, *, source: str = "delay"):
        if method not in DELAY_METHODS:
            raise ValueError(f"delay method must be one of {DELAY_METHODS}, got {method!r}")
        if not (np.isfinite(delay_samples) and delay_samples >= 0):
            raise ValueError(f"delay must be finite and >= 0, got {delay_samples!r} samples")
        D = float(delay_samples)
        self.delay_samples = D
        self.method = method
        self.mu = 0.0
        self.a = 0.0
        if abs(D - round(D)) <= _INT_TOL:
            self.mode, self.k = _MODE_INT, int(round(D))
        elif method == "linear":
            self.mode, self.k = _MODE_LINEAR, int(math.floor(D))
            self.mu = D - self.k
        else:
            self.mode = _MODE_THIRAN
            if D < 0.5:
                ChkLog().warn("thiran1_short_delay", source,
                              f"thiran1 with D = {D:.4g} samples < 0.5: pole at {-(1.0 - D) / (1.0 + D):.4f} rings at "
                              "Nyquist; use method 'linear' (core.delay docstring)", delay_samples=D)
            self.k = 0 if D < 0.5 else int(math.floor(D - 0.5))
            d = D - self.k
            self.a = (1.0 - d) / (1.0 + d)
        self.ring_len = self.k + 2

    def freq_response(self, f_hz: Any, fs_hz: float) -> np.ndarray:
        """Exact frequency response H(e^{jωh}) of this discrete delay at ``f_hz`` for step rate ``fs_hz``."""
        z1 = np.exp(-2j * np.pi * np.asarray(f_hz, dtype=float) / fs_hz)  # z^-1
        h = z1**self.k
        if self.mode == _MODE_LINEAR:
            h = h * ((1.0 - self.mu) + self.mu * z1)
        elif self.mode == _MODE_THIRAN:
            h = h * (self.a + z1) / (1.0 + self.a * z1)
        return h


def design_delay(tau_s: float, step_hz: float, method: str, *, source: str = "delay") -> DelayDesign:
    """:class:`DelayDesign` for a delay ``tau_s`` at step rate ``step_hz`` (= f_sim / rate_div);
    ``source`` names the owner in CHK warnings."""
    return DelayDesign(tau_s * step_hz, method, source=source)


def interp_variance_deficit(
    channels: Any,
    step_hz: float,
    psd: Any,
    *,
    f_lo_hz: float = 0.0,
    n_f: int = 4096,
) -> float:
    """In-band variance the linear fractional delay misses (D65-4, Phase 5 F6-1):

        ΔVar = ∫_{f_lo}^{f_N} (|H_exact(f)|² − |H_lin(f)|²) S(f) df,   f_N = step_hz / 2,

    summed over independent ``channels``. A channel is a list of taps ``(delay_s, weight)`` that see the
    *same* driver sequence (coherent: H = Σ w·e^{−i2πf·d} exactly, Σ w·H_lin(d) with the ``"linear"``
    design of :class:`DelayDesign`); different channels are independent drivers with the same one-sided
    PSD ``psd(f)`` [unit²/Hz] (callable, vectorised). Midpoint rule on ``n_f`` points (f = 0 is never
    evaluated, so 1/f² spectra work with ``f_lo_hz = 0``).

    It is an **offline correction, never applied automatically**: when a preset needs the total in-band
    variance of a path-phase combination whose driver is flat up to f_N, it declares ΔVar in the
    ``BandLedger`` (e.g. ``BandDecl(source=..., quantity="phase", targets=(...), f_lo_hz=...,
    f_hi_hz=f_N, unit="rad", variance_outside=ΔVar)`` on a dedicated target), alongside the
    out-of-band part. Spectral conclusions should instead be drawn at f ≤ f_N/4 (≤ 5 % error) or
    f ≤ f_N/10 (1 %), where no correction is needed (D65-4). ``path.segments.subsegment_channels``
    builds the channels of a beam combination on a segment.
    """
    if not (step_hz > 0 and np.isfinite(step_hz)):
        raise ValueError(f"step_hz must be finite and > 0, got {step_hz!r}")
    f_n = 0.5 * float(step_hz)
    if not (0.0 <= f_lo_hz < f_n):
        raise ValueError(f"need 0 <= f_lo_hz < f_N = {f_n:g}, got {f_lo_hz!r}")
    df = (f_n - f_lo_hz) / int(n_f)
    f = f_lo_hz + (np.arange(int(n_f)) + 0.5) * df
    s = np.asarray(psd(f), dtype=float) * np.ones_like(f)
    total = 0.0
    for taps in channels:
        h_ex = np.zeros(f.shape, dtype=complex)
        h_li = np.zeros(f.shape, dtype=complex)
        for d_s, w in taps:
            h_ex += float(w) * np.exp(-2j * np.pi * f * float(d_s))
            h_li += float(w) * design_delay(float(d_s), step_hz, "linear").freq_response(f, step_hz)
        total += float(np.sum((np.abs(h_ex) ** 2 - np.abs(h_li) ** 2) * s) * df)
    return total


def ring_push(S, I, off, ring_len, head, x):  # noqa: N803 (flat-array names, G02)
    """Push ``x`` as the newest sample of the ring at ``S[off : off + ring_len]``; returns the new head."""
    h = I[head] + 1
    if h == ring_len:
        h = 0
    I[head] = h
    S[off + h] = x
    return h


def ring_at(S, off, ring_len, h, d):  # noqa: N803
    """Value ``d`` (0 <= d < ring_len) samples older than the newest, given head ``h``."""
    j = h - d
    if j < 0:
        j += ring_len
    return S[off + j]


ring_push_nb = numba.njit(inline="always")(ring_push)
ring_at_nb = numba.njit(inline="always")(ring_at)


def _make_delay_step(push, at):
    def delay_step(S, I, x, off, ring_len, head, ap, mode, k, mu, a):  # noqa: N803
        """Push ``x`` and return the delayed output (module docstring); ``ap`` = Thiran state index in S."""
        h = push(S, I, off, ring_len, head, x)
        if mode == 0:
            return at(S, off, ring_len, h, k)
        if mode == 1:
            return (1.0 - mu) * at(S, off, ring_len, h, k) + mu * at(S, off, ring_len, h, k + 1)
        v = a * at(S, off, ring_len, h, k) + at(S, off, ring_len, h, k + 1) - a * S[ap]
        S[ap] = v
        return v

    return delay_step


#: Pure-Python delay step ``delay_step(S, I, x, off, ring_len, head, ap, mode, k, mu, a) -> y``.
delay_step = _make_delay_step(ring_push, ring_at)
#: The same step compiled for kernels.
delay_step_nb = numba.njit(inline="always")(_make_delay_step(ring_push_nb, ring_at_nb))


class DelayLine:
    """Stand-alone delay line (Python), same arithmetic as the ``delay_buffer`` block.

    Args:
        tau_s: delay [s] (>= 0).
        f_sim_hz: sample rate of the writer [Hz].
        method: one of :data:`DELAY_METHODS`.
    """

    def __init__(self, tau_s: float, f_sim_hz: float, method: str):
        self.design = design_delay(tau_s, f_sim_hz, method)
        self.tau_s = float(tau_s)
        self.f_sim_hz = float(f_sim_hz)
        L = self.design.ring_len
        self._S = np.zeros(L + 1)  # ring rows 0..L-1, Thiran state at L
        self._I = np.zeros(1, dtype=np.int64)
        self._y = 0.0

    def push(self, x: float) -> None:
        """Write the newest sample (and advance the fractional-delay state)."""
        d = self.design
        self._y = float(delay_step(self._S, self._I, float(x), 0, d.ring_len, 0, d.ring_len, d.mode, d.k, d.mu, d.a))

    def read(self) -> float:
        """Value delayed by tau_s relative to the newest sample."""
        return self._y

    def process(self, x: Any) -> np.ndarray:
        """Push every sample of ``x`` and return the delayed sequence."""
        out = np.empty(len(x))
        for i, xi in enumerate(x):
            self.push(xi)
            out[i] = self._y
        return out

    def freq_response(self, f_hz: Any) -> np.ndarray:
        return self.design.freq_response(f_hz, self.f_sim_hz)


class DelayBuffer(BaseBlock):
    """The ``delay_buffer`` block (module docstring)."""

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or DELAY_BUFFER_KIND
        require(params, "in", "out", "tau_s", "method", block=name)
        if params["method"] not in DELAY_METHODS:
            raise ValueError(f"block {name!r}: method must be one of {DELAY_METHODS}, got {params['method']!r}")
        if not (np.isfinite(float(params["tau_s"])) and float(params["tau_s"]) >= 0):
            raise ValueError(f"block {name!r}: tau_s must be finite and >= 0, got {params['tau_s']!r}")
        self.ring_name = f"delay.{name}.ring"
        self.ap_name = f"delay.{name}.ap"
        super().__init__(
            DELAY_BUFFER_KIND,
            name,
            params,
            reads=(params["in"],),
            writes=(params["out"], self.ring_name, self.ap_name),
        )
        self.design: DelayDesign | None = None

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        src, out = self.params["in"], self.params["out"]
        if "unit" in self.params:
            unit = str(self.params["unit"])
        elif src in layout:
            unit = layout.slot(src).unit
        else:
            raise ValueError(
                f"block {self.name!r}: input slot {src!r} is not registered yet; give params['unit'] "
                "or order its writer before this delay_buffer"
            )
        self.slot_units = {src: unit, out: unit}
        self.design = design_delay(float(self.params["tau_s"]), 1.0 / block_dt_s(layout.f_sim_hz, self.rate_div),
                                   self.params["method"], source=self.name)
        layout.add(out, (), np.float64, unit)
        layout.add(self.ring_name, (), np.float64, unit, ring_len=self.design.ring_len)
        layout.add(self.ap_name, (), np.float64, unit)

    def _offsets(self, layout):
        ring = layout.slot(self.ring_name)
        for nm in (self.params["in"], self.params["out"]):
            sl = layout.slot(nm)
            if sl.array != "S" or sl.shape != () or sl.ring_len:
                raise TypeError(f"block {self.name!r}: slot {nm!r} must be a scalar float64 slot")
        return (layout.slot(self.params["in"]).offset, layout.slot(self.params["out"]).offset,
                ring.offset, ring.ring_len, ring.head, layout.slot(self.ap_name).offset)

    def bind(self, layout) -> None:
        self._off = self._offsets(layout)
        i_in, i_out, off, L, head, ap = self._off
        d = self.design
        mode, k, mu, a = d.mode, d.k, d.mu, d.a

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            S[i_out] = delay_step_nb(S, I, S[i_in], off, L, head, ap, mode, k, mu, a)

        self.kernel = kernel

    def __call__(self, state, params):
        i_in, i_out, off, L, head, ap = self._off
        d = self.design
        S = state.S
        S[i_out] = delay_step(S, state.I, float(S[i_in]), off, L, head, ap, d.mode, d.k, d.mu, d.a)
        return state


@register_block(DELAY_BUFFER_KIND)
def make_delay_buffer(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``delay_buffer`` block (see module docstring)."""
    return DelayBuffer(params, name)
