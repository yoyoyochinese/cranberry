"""Shared machinery of the ERR blocks (WP-B3): wiring, sensor noise, dither harmonics, A± helpers.

Signal representation (T14 §2.1, G16-b)
    A generator writes loop-rate voltages. A carrier at signed frequency f_c is written in one of
    two forms, chosen by the toggle ``rf_rate_mode`` (generator and ``demodulator`` must agree):

    * baseband (``rf_rate_mode`` off): ``out`` has shape (2,) = (a, b), the quadrature amplitudes of
      v(t) = DC + a·cos(2π f_c t) + b·sin(2π f_c t) + (other harmonics, removed by the LPF). The
      dither of a dither-type lock is applied analytically (``dither_mode: analytic`` on the
      modulator), the carrier phase reference is the dither sin(2π f_m t).
    * RF rate (``rf_rate_mode`` on): ``out`` is the scalar instantaneous voltage v(t), evaluated
      from the instantaneous beam phases (which then carry a time-series dither,
      ``dither_mode: timeseries``). Needs f_sim > 2·|f_c|; used for cross-validation (T14 V1, V2, V6).

    ``demodulator`` (2cos(2π f_ref t + φ_d), LPF) turns either form into the same I/Q:
    I = a cos φ_d − b sin φ_d, Q = a sin φ_d + b cos φ_d.

Sensor noise (T14 §2.2)
    White noise of one-sided PSD S_v [V²/Hz] near the carrier becomes, after demodulation with
    2cos(...), white noise of PSD 2·S_v. Baseband: a and b each get PSD 2·S_v (per-sample std
    √(S_v f_b)); RF rate and DC outputs: PSD S_v (std √(S_v f_b / 2)); f_b = f_sim / rate_div.
    In a dither lock S_v depends on the instantaneous phase, i.e. the noise is cyclostationary at
    f_m; demodulation then weights it: a gets 2⟨2cos²θ S_v(θ)⟩, b gets 2⟨2sin²θ S_v(θ)⟩ (θ = 2π f_m t).
    The baseband mode uses these weights (computed in the same quadrature, with S_v = law′(p)²·2Q for
    the detector law in use), so it equals the RF-rate mode for every law. The size of the weight
    depends on the law (D65-1, M1-A-04): for Ref C at the squeezing point it is 1.2257 × S_v(Δ0) for
    the SQUARE law only (sensor-limited residual 3.24 → 3.58 mrad, continuous first-order loop with
    b_loop = 10 kHz); for the LOG law S_v ∝ 2Q/P² is constant over a flat band, the weight is exactly 1
    and the change enters through the slope instead (0.905 × the stationary T14 slope, 1.222 in PSD):
    4.11 mrad with κ_det = 1.32 (T14 §6.3's 4.2 multiplied the rms by κ_det; stationary value 3.72).
    Above the block's Nyquist frequency sensor noise is not represented; it never reaches the lock
    (filtered by the loop), so each generator declares a zero-outside sensor band to the
    ``BandLedger`` (G05) and adds nothing out of band itself.

Dither harmonics (T14 §4–§6)
    For v(θ) = law(k, (d + c1 cos(x1 + β sinθ) + c2 cos(x2 + β sinθ))·(1 + a_RAM sin(θ + φ_a)) + d_un),
    θ = 2π f_m t, the DC and first-harmonic coefficients are computed by the trapezoid rule over
    one dither period with :data:`N_QUAD` points. The integrand is periodic and analytic, so this is
    exact to rounding for |β| ≲ 3 (Fourier content beyond harmonic N ~ J_N(β)); it reproduces the
    Bessel forms (e.g. −2 J1(β) sin x) without a Bessel routine in the kernel and works for any
    detector law (square / linear envelope / log).

Reference vs compiled
    Every step function is written once inside a ``_build(jit)`` factory and instantiated twice:
    pure Python (``jit=False``, the reference ``block(state, params)``) and ``numba.njit``
    (``jit=True``, the kernel). Both run the same equations (VT-44: agree to 1e−12).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from types import SimpleNamespace
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, require
from cranberry.core.rates import block_dt_s
from cranberry.core.rng import DrawBuffer
from cranberry.offline.band_ledger import BandDecl
from cranberry.offline.tables import Table

__all__ = [
    "E_CHARGE_C",
    "N_QUAD",
    "DETECTOR_LAWS",
    "ErrBlock",
    "helpers",
    "amplitude_gains_opa",
    "amplitude_gains_opo",
    "amplitude_gains_from_table",
    "register",
]

#: elementary charge [C] (shot noise S_I = 2 e I, one-sided)
E_CHARGE_C = 1.602176634e-19
#: quadrature points per dither period (see module docstring)
N_QUAD = 32
#: detector laws v = law(k, P): square (k·P), linear envelope (k·√P), log (k·10·log10 P) [T14 §6]
DETECTOR_LAWS: tuple[str, ...] = ("square", "linear", "log")


def helpers(jit: bool) -> SimpleNamespace:
    """Shared step helpers, pure Python (``jit=False``) or numba-compiled (``jit=True``)."""
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def law(code, k, p):
        if code == 0:
            return k * p
        if code == 1:
            return k * math.sqrt(p)
        return k * 10.0 * math.log10(p)

    @dec
    def noise(code, k, p, cx, n0, n1, nq0, nq1, nq2):
        """S_v(p) = n0 + n1·p + law′(p)²·(nq0 + nq1·cx + nq2·cx²) [V²/Hz] (PD shot/electronic; radiometer)."""
        s = n0 + n1 * p
        if nq0 != 0.0 or nq1 != 0.0 or nq2 != 0.0:
            if code == 0:
                dv = k
            elif code == 1:
                dv = 0.5 * k / math.sqrt(p)
            else:
                dv = 10.0 / math.log(10.0) * k / p
            s += dv * dv * (nq0 + nq1 * cx + nq2 * cx * cx)
        return s

    @dec
    def harmonic(code, k, d, c1, x1, c2, x2, beta, ram_a, ram_phi, n0, n1, nq0, nq1, nq2, d_un):
        """Over one dither period θ of p(θ) = (d + c1 cos(x1 + β sinθ) + c2 cos(x2 + β sinθ))(1 + a sin(θ + φ_a))
        + d_un (d_un: part not modulated by RAM, e.g. the detector's dark noise), v = law(k, p): returns
        (dc, a, b, s_a, s_b) with v ≈ dc + a cosθ + b sinθ and the cyclostationary noise weights
        s_a = ⟨2cos²θ S_v(θ)⟩, s_b = ⟨2sin²θ S_v(θ)⟩ (S_v from :func:`noise`, cx = cos(x1 + β sinθ))."""
        dc = 0.0
        a = 0.0
        b = 0.0
        sa = 0.0
        sb = 0.0
        for j in range(N_QUAD):
            th = 2.0 * math.pi * j / N_QUAD
            s = math.sin(th)
            co = math.cos(th)
            cx = math.cos(x1 + beta * s)
            p = (d + c1 * cx + c2 * math.cos(x2 + beta * s)) * (1.0 + ram_a * math.sin(th + ram_phi)) + d_un
            v = law(code, k, p)
            dc += v
            a += v * co
            b += v * s
            sv = noise(code, k, p, cx, n0, n1, nq0, nq1, nq2)
            sa += 2.0 * co * co * sv
            sb += 2.0 * s * s * sv
        return dc / N_QUAD, 2.0 * a / N_QUAD, 2.0 * b / N_QUAD, sa / N_QUAD, sb / N_QUAD

    @dec
    def draw(S, I, buf, idx):
        w = S[buf + I[idx]]
        I[idx] += 1
        return w

    @dec
    def emit(S, I, c, out, a, b, s_a, s_b):
        """Write a carrier output: baseband (a, b) + noise (PSD 2 s_a, 2 s_b) or RF/DC scalar a + noise (PSD s_a)."""
        if c.f_rf:
            if c.buf >= 0:
                a += math.sqrt(s_a * c.fb * 0.5) * draw(S, I, c.buf, c.idx)
            S[out] = a
        else:
            if c.buf >= 0:
                a += math.sqrt(s_a * c.fb) * draw(S, I, c.buf, c.idx)
                b += math.sqrt(s_b * c.fb) * draw(S, I, c.buf, c.idx)
            S[out] = a
            S[out + 1] = b

    return SimpleNamespace(law=law, noise=noise, harmonic=harmonic, draw=draw, emit=emit)


# --------------------------------------------------------------------------------------------------
# A± (amplitude gains of the two quadratures, T14 §1) — inputs to parametric_beat_detector / bhd_probe_beat
# --------------------------------------------------------------------------------------------------


def amplitude_gains_opa(r: float) -> tuple[float, float]:
    """Single-pass OPA (Ref A): A± = e^{±r}, r = √(αP) (T14 §1)."""
    return math.exp(r), math.exp(-r)


def amplitude_gains_opo(x: float) -> tuple[float, float]:
    """OPO below threshold, quasi-static (Ref B): A+ = 1/(1−x), A− = 1/(1+x), x = √(P/P_th) (T02, D51-10)."""
    if not 0.0 <= x < 1.0:
        raise ValueError(f"OPO pump parameter x must be in [0, 1), got {x}")
    return 1.0 / (1.0 - x), 1.0 / (1.0 + x)


def amplitude_gains_from_table(table: Table, **coords: Any) -> tuple[float, float]:
    """(A+, A−) from an F3 ``Table`` whose entries have value_shape (2,) = (A+, A−) (SourceOutput.gain_table)."""
    if table.value_shape != (2,):
        raise ValueError(f"gain table {table.name!r}: value_shape must be (2,) = (A+, A−), got {table.value_shape}")
    ap, am = (float(v) for v in table.interp(**coords))
    return ap, am


# --------------------------------------------------------------------------------------------------
# block base
# --------------------------------------------------------------------------------------------------


class ErrBlock(BaseBlock):
    """Wiring shared by ERR blocks.

    Subclasses set ``CONST`` (namedtuple type: fields ``s_<role>`` slot offsets, ``p_<key>`` P
    offsets, ``f_<flag>`` bools, plus ``t``, ``buf``, ``idx``, ``fb``, ``mode``, ``law``) and
    ``STEP = (py_step, nb_step)``, and call :meth:`setup` from ``__init__``.
    """

    CONST: Any = None
    STEP: tuple[Any, Any] = (None, None)

    def setup(
        self,
        kind: str,
        name: str | None,
        params: Mapping[str, Any],
        *,
        inputs: Mapping[str, tuple[str | None, str]],
        outputs: Mapping[str, tuple[str | None, tuple[int, ...], str]],
        runtime: Iterable[str],
        toggles: Iterable[str],
        n_draws: int,
        mode: int = 0,
        law: int = 0,
    ) -> None:
        """inputs: role → (slot or None, unit); outputs: role → (slot or None, shape, unit)."""
        self._inputs = {r: s for r, (s, _) in inputs.items()}
        self._outputs = {r: (s, shp, u) for r, (s, shp, u) in outputs.items()}
        reads = tuple(s for s in self._inputs.values() if s is not None)
        writes = tuple(s for s, _, _ in self._outputs.values() if s is not None)
        super().__init__(kind, name, params, reads=reads, writes=writes, runtime_params=tuple(runtime),
                         toggle_names=tuple(toggles))
        self.slot_units = {s: u for s, u in inputs.values() if s is not None}
        self.slot_units.update({s: u for s, _, u in outputs.values() if s is not None})
        self._n_draws = int(n_draws)
        self._mode, self._law = int(mode), int(law)
        self.draws: DrawBuffer | None = None
        self._fb: float | None = None

    @property
    def flags(self) -> dict[str, bool]:
        """Structural flags passed to the step as ``f_<name>`` (default: the toggles)."""
        return dict(self.toggles)

    def init(self, layout, rng):
        super().init(layout, rng)
        for slot, shape, unit in self._outputs.values():
            if slot is not None:
                layout.add(slot, shape, np.float64, unit)
        self._fb = 1.0 / block_dt_s(layout.f_sim_hz, self.rate_div)
        if self._n_draws > 0:
            self.draws = DrawBuffer(self.name, self._n_draws, "standard_normal")
            self.draws.init(layout, rng, self.rate_div)
            self.writes = (*self.writes, *self.draws.slot_names)

    def bind(self, layout):
        vals: dict[str, Any] = {}
        for field in self.CONST._fields:
            if field.startswith(("s_", "p_")):
                vals[field] = -1
        for role, slot in self._inputs.items():
            vals[f"s_{role}"] = layout.slot(slot).offset if slot is not None else -1
        for role, (slot, _, _) in self._outputs.items():
            vals[f"s_{role}"] = layout.slot(slot).offset if slot is not None else -1
        for key in self.runtime_params:
            vals[f"p_{key}"] = layout.param(self.param_name(key)).offset
        for flag, on in self.flags.items():
            vals[f"f_{flag}"] = bool(on)
        if self.draws is not None:
            self.draws.bind(layout)
            vals["buf"], vals["idx"] = self.draws.buf_off, self.draws.idx_off
        else:
            vals["buf"], vals["idx"] = -1, -1
        vals.update(t=layout.slot("t").offset, fb=float(self._fb), mode=self._mode, law=self._law)
        unknown = set(vals) - set(self.CONST._fields)
        if unknown:
            raise RuntimeError(f"{self.kind}: constants {sorted(unknown)} not in {self.CONST.__name__}")
        c = self.CONST(**{f: vals.get(f, False) for f in self.CONST._fields})
        self._c = c
        py_step, nb_step = self.STEP
        self._py = py_step

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb_step(S, Z, I, P, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        if self.draws is not None:
            self.draws.refill(state, n0, n_samples)

    def __call__(self, state, params):
        self._py(state.S, state.Z, state.I, params.P, self._c)
        return state

    def bands(self) -> list[BandDecl]:
        """Sensor-noise band: represented up to the block Nyquist frequency, nothing added outside (G05)."""
        if self.draws is None:
            return []
        if self._fb is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        targets = tuple(s for s, _, _ in self._outputs.values() if s is not None)
        return [BandDecl(source=self.name, quantity="sensor", targets=targets, f_lo_hz=0.0,
                         f_hi_hz=0.5 * self._fb, unit="V", variance_outside=0.0)]


def choice(params: Mapping[str, Any], key: str, options: tuple[str, ...], block: str) -> int:
    """Index of the required structural option ``params[key]`` in ``options``."""
    require(params, key, block=block)
    if params[key] not in options:
        raise ValueError(f"block {block!r}: {key} must be one of {options}, got {params[key]!r}")
    return options.index(params[key])


def register(kind: str, cls: type) -> Any:
    """Register ``cls(params, name)`` as factory ``make_<kind>``."""
    from cranberry.core.registry import register_block

    @register_block(kind)
    def factory(params, *, name=None):
        return cls(params, name)

    factory.__name__ = factory.__qualname__ = f"make_{kind}"
    return factory
