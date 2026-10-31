"""Slow drift sources: ``segment_thermal``, ``crystal_thermal``, ``pd_offset_drift`` (T08 §2.2–§2.4, T01 §2.5,
T15 §2.4). Owner: WP-B2.

All three share one engine: a common input u(t) and first-order lags per output,

    u(t)  = c0 [offset] + r·t [drift] + A·sin(2π t/T_per) [periodic] + w(t) [random walk, S_w = b/f²]
    τ_i dx_i/dt = u − x_i          (τ_i = 0: x_i = u)

``segment_thermal`` (T08 eq. 2.5–2.6): u = ambient temperature change of one thermal zone (the old
``ambient_temperature`` is this u; one block per zone), x_i = ΔT of segment i with its own τ_th,i,
written to ``drv.<seg>.dT`` [K] (shape (N_sub,) with every entry equal when the topology places dT
``per_sub``, else scalar). Segments of one block are fully correlated at low frequency; unequal τ_th
make them differ (T08 §2.4). Toggles {drift, random, periodic}.

``crystal_thermal`` (T01 §2.5, D42 (a)): the same for one crystal segment (scalar ``drv.<segment>.dT``,
pump-path phase through the topology coefficient) plus the phase mismatch for the OPA's rotation
δθ_mm (B1): Δkℓ/2 = 2.784 (δT_0 + ΔT)/ΔT_FWHM, i.e. ``opa.<name>.dk`` [1/m] = 5.568 (δT_0 + ΔT)/(ℓ ΔT_FWHM),
δT_0 = T_set − T_pm (``detune_k``). Toggles {drift, random, periodic}.

``pd_offset_drift`` (T15 §2.4 item 1): electrical DC offset V_off = c0 + r·t + w(t) [V] on ``out``
(no lag). Toggles {offset, drift, random}.

Discretisation: exact for an input linear between samples (first-order hold): x_{n+1} = a x_n +
(1 − a) u_n + (u_{n+1} − u_n)(1 − τ(1 − a)/h), a = e^{−h/τ}; ramp responses are exact (T08 V7, V10).
Random walk: per-step increment std √(2π² b h) (one-sided S_w = b/f²). States start at x = u(0).

Ledger (G05): the represented band is [0, f_N]; above it the random-walk tail
b/f²·|1/(1 + i2πfτ)|² is declared (``variance_outside = 0`` when ``random`` is off).
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer
from cranberry.offline.band_ledger import BandDecl
from cranberry.path.segments import select_ids
from cranberry.path.topology import as_topology

__all__ = ["SegmentThermal", "CrystalThermal", "PdOffsetDrift", "make_segment_thermal", "make_crystal_thermal",
           "make_pd_offset_drift", "k_T_rad_per_m_K", "DK_HALF_PER_FWHM"]

#: Δkℓ/2 at the half maximum of sinc²(Δkℓ/2) (T01 §2.5)
DK_HALF_PER_FWHM = 2.784


def k_T_rad_per_m_K(wavelength_m: float, dn_dT_per_k: float, n: float, alpha_l_per_k: float) -> float:
    """Thermal phase coefficient k_T = (2π/λ)(dn/dT + n α_L) [rad/(m·K)] (T08 eq. 2.2; > 0, conv §4)."""
    return 2.0 * math.pi / wavelength_m * (dn_dT_per_k + n * alpha_l_per_k)


DC = namedtuple("DC", "t h c0 r amp per f_rand sig wbuf widx w nout xs oo size a b dk kdk det")


def _build(jit: bool):
    def u_of(c, t, w):
        return c.c0 + c.r * t + c.amp * math.sin(2.0 * math.pi * t / c.per) + w

    if jit:
        u_of = numba.njit(u_of)

    def step(S, I, c):  # noqa: N803
        t = S[c.t]
        u0 = u_of(c, t, S[c.w])
        if c.f_rand:
            S[c.w] += c.sig * S[c.wbuf + I[c.widx]]
            I[c.widx] += 1
        u1 = u_of(c, t + c.h, S[c.w])
        for i in range(c.nout):
            x = S[c.xs + i]
            for j in range(c.size[i]):
                S[c.oo[i] + j] = x
            if c.dk >= 0 and i == 0:
                S[c.dk] = c.kdk * (c.det + x)
            S[c.xs + i] = c.a[i] * x + (1.0 - c.a[i]) * u0 + (u1 - u0) * c.b[i]

    return numba.njit(step) if jit else step


_STEP_PY, _STEP_NB = _build(False), _build(True)


class _DriftLag(BaseBlock):
    """Engine (module docstring). Subclasses set ``outputs`` [(slot, shape, unit, tau_s)] before init."""

    KIND = ""
    UNIT = "1"
    #: toggle → spec keys it needs
    NEED: dict[str, tuple[str, ...]] = {}

    def __init__(self, params, name, outputs, extra_writes=()):
        name = name or self.KIND
        require(params, "toggles", block=name)
        for t, keys in self.NEED.items():
            if params["toggles"].get(t):
                require(params, *keys, block=name)
        self.outputs = outputs
        super().__init__(self.KIND, name, params, reads=("t",), writes=(*(o[0] for o in outputs), *extra_writes),
                         toggle_names=tuple(self.NEED))
        self.slot_units = {o[0]: o[2] for o in outputs}
        self._h: float | None = None

    def _g(self, key: str, toggle: str) -> float:
        return float(self.params[key]) if self.toggles.get(toggle) else 0.0

    def init(self, layout, rng):
        super().init(layout, rng)
        self._h = block_dt_s(layout.f_sim_hz, self.rate_div)
        c0 = self._g("offset_v", "offset")
        for slot, shape, unit, _ in self.outputs:
            layout.add(slot, shape, np.float64, unit, init=c0)
        self.state = (f"rng.{self.name}.lag", f"rng.{self.name}.walk")
        layout.add(self.state[0], (len(self.outputs),), np.float64, "1", init=c0)
        layout.add(self.state[1], (), np.float64, "1")
        self.writes = (*self.writes, *self.state)
        self.draws = None
        if self.toggles.get("random"):
            self.draws = DrawBuffer(self.name, 1, "standard_normal")
            self.draws.init(layout, rng, self.rate_div)
            self.writes = (*self.writes, *self.draws.slot_names)

    def _crystal(self, layout) -> tuple[int, float, float]:
        return -1, 0.0, 0.0

    def bind(self, layout):
        h, off = self._h, lambda s: layout.slot(s).offset  # noqa: E731
        if self.draws is not None:
            self.draws.bind(layout)
        taus = [o[3] for o in self.outputs]
        a = np.array([math.exp(-h / t) if t > 0 else 0.0 for t in taus])
        b = np.array([1.0 - t * (1.0 - ai) / h if t > 0 else 1.0 for t, ai in zip(taus, a)])
        periodic = self.toggles.get("periodic", False)
        self._c = DC(
            off("t"), h, self._g("offset_v", "offset"), self._g(self.RATE, "drift"),
            self._g(self.AMP, "periodic"), float(self.params[self.PER]) if periodic else 1.0,
            bool(self.toggles.get("random")), math.sqrt(2.0 * math.pi**2 * self._g(self.WALK, "random") * h),
            self.draws.buf_off if self.draws else -1, self.draws.idx_off if self.draws else -1, off(self.state[1]),
            len(self.outputs), off(self.state[0]), np.array([off(o[0]) for o in self.outputs], dtype=np.int64),
            np.array([int(np.prod(o[1], dtype=np.int64)) for o in self.outputs], dtype=np.int64), a, b,
            *self._crystal(layout),
        )
        c, nb = self._c, _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        if self.draws is not None:
            self.draws.refill(state, n0, n_samples)

    def __call__(self, state, params):
        _STEP_PY(state.S, state.I, self._c)
        return state

    def bands(self) -> list[BandDecl]:
        if self._h is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        f_n, b = 0.5 / self._h, self._g(self.WALK, "random")
        out = []
        for slot, _, unit, tau in self.outputs:
            target = self._target(slot)
            if b > 0.0:
                psd = lambda f, tau=tau: b / np.asarray(f) ** 2 / (1.0 + (2.0 * np.pi * np.asarray(f) * tau) ** 2)
                out.append(BandDecl(self.name, self.QUANTITY, (target,), 0.0, f_n, unit, psd_outside=psd))
            else:
                out.append(BandDecl(self.name, self.QUANTITY, (target,), 0.0, f_n, unit, variance_outside=0.0))
        return out

    def _target(self, slot: str) -> str:
        return slot.split(".")[1]  # drv.<seg>.dT → segment id


class SegmentThermal(_DriftLag):
    KIND, QUANTITY = "segment_thermal", "dT"
    RATE, AMP, PER, WALK = "r_T_k_per_s", "periodic_amp_k", "periodic_period_s", "b_T_k2_hz"
    NEED = {"drift": ("r_T_k_per_s",), "random": ("b_T_k2_hz",), "periodic": ("periodic_amp_k", "periodic_period_s")}

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "topology", "segments", "tau_th_s", block=name)
        topo = as_topology(params["topology"])
        outs = []
        for s in select_ids(params["segments"], topo.segments, "segments", name):
            seg = topo.segments[s]
            if "dT" not in seg.drivers:
                raise ValueError(f"block {name!r}: segment {s!r} declares no dT driver")
            shape = (seg.n_sub,) if seg.drivers["dT"][0] == "per_sub" else ()
            outs.append((f"drv.{s}.dT", shape, "K", float(params["tau_th_s"][s])))
        super().__init__(params, name, outs)


class CrystalThermal(_DriftLag):
    KIND, QUANTITY = "crystal_thermal", "dT"
    RATE, AMP, PER, WALK = SegmentThermal.RATE, SegmentThermal.AMP, SegmentThermal.PER, SegmentThermal.WALK
    NEED = SegmentThermal.NEED

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "segment", "tau_th_s", "crystal_length_m", "dT_fwhm_k", "detune_k", block=name)
        self.dk_slot = f"opa.{name}.dk"
        super().__init__(params, name, [(f"drv.{params['segment']}.dT", (), "K", float(params["tau_th_s"]))],
                         (self.dk_slot,))
        self.slot_units[self.dk_slot] = "1/m"

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.dk_slot, (), np.float64, "1/m")

    def _crystal(self, layout):
        k = 2.0 * DK_HALF_PER_FWHM / (float(self.params["crystal_length_m"]) * float(self.params["dT_fwhm_k"]))
        return layout.slot(self.dk_slot).offset, k, float(self.params["detune_k"])


class PdOffsetDrift(_DriftLag):
    KIND, QUANTITY = "pd_offset_drift", "offset"
    RATE, AMP, PER, WALK = "r_v_per_s", "", "", "b_v2_hz"
    NEED = {"offset": ("offset_v",), "drift": ("r_v_per_s",), "random": ("b_v2_hz",)}

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "out", block=name)
        super().__init__(params, name, [(params["out"], (), "V", 0.0)])

    def _target(self, slot: str) -> str:
        return slot


def _register(cls):
    @register_block(cls.KIND)
    def factory(params: Mapping[str, Any], *, name: str | None = None) -> Block:
        return cls(params, name)

    factory.__name__ = factory.__qualname__ = f"make_{cls.KIND}"
    return factory


make_segment_thermal = _register(SegmentThermal)
make_crystal_thermal = _register(CrystalThermal)
make_pd_offset_drift = _register(PdOffsetDrift)
