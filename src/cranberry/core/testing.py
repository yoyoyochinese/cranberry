"""Toy blocks for framework tests and benchmarks (WP-F1). Not physics blocks.

Importing this module registers the kinds ``toy_noise_ar1``, ``toy_sine``, ``toy_sum``,
``toy_gain``, ``toy_pid``, ``toy_sat`` and ``toy_integrator`` (``core.registry.ensure_loaded`` does
not import it, so production registries stay clean). Each block has a compiled kernel and a
pure-Python reference step written independently from the same equations, follows the block
contract (all parameters required, D34; runtime parameters in P), and is small enough to read in
one glance. The real CTRL / ACT / NOISE blocks are WP-B2 / WP-B4.

=================  =================================================================  ===========
kind               step (h = rate_div / f_sim)                                        runtime P
=================  =================================================================  ===========
toy_noise_ar1      out ← a·out + b·w,  w ~ N(0, 1) from a DrawBuffer (G11)             a, b
toy_sine           out ← amp·sin(2π·freq_hz·t)                                         amp, freq_hz
toy_sum            out ← Σ_k gains[k]·ins[k]                                           gains
toy_gain           out ← gain·in                                                       gain
toy_pid            e = setpoint − in; I ← I + ki·h·e; out ← kp·e + I + kd·(e − e₋)/h   kp, ki, kd,
                                                                                        setpoint
toy_sat            out ← clip(gain·in, −limit, +limit)                                 gain, limit
toy_integrator     out ← out + h·k·in                                                  k
=================  =================================================================  ===========

Every toy block has the toggle ``enabled``; off makes it a zero source (noise, sine) or an identity
(gain, sum of one input, sat) or a hold (pid, integrator).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["ToyBlock", "toy"]


class ToyBlock(BaseBlock):
    """Common wiring: ``in``/``ins`` read, ``out`` written (registered here with ``unit``)."""

    toggle_names = ("enabled",)
    RUNTIME: tuple[str, ...] = ()
    KIND = ""

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "out", "unit", block=name)
        ins = self._inputs(params, name)
        self.state_slots = self._state_slots(name)
        super().__init__(self.KIND, name, params, reads=tuple(ins), writes=(params["out"], *self.state_slots),
                         runtime_params=self.RUNTIME)
        self.slot_units = {params["out"]: str(params["unit"])}

    def _inputs(self, params, name) -> list[str]:
        require(params, "in", block=name)
        return [params["in"]]

    def _state_slots(self, name) -> tuple[str, ...]:
        return ()

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"], (), np.float64, str(self.params["unit"]))
        for s in self.state_slots:
            layout.add(s, (), np.float64, "1")
        self.h = block_dt_s(layout.f_sim_hz, self.rate_div)

    def bind(self, layout):
        self.o = {nm: layout.slot(nm).offset for nm in (*self.reads, *self.writes)}
        self.p = {k: layout.param(self.param_name(k)).offset for k in self.runtime_params}
        self.kernel = self._kernel(layout)

    def _kernel(self, layout):  # pragma: no cover - overridden
        raise NotImplementedError


class NoiseAR1(ToyBlock):
    KIND, RUNTIME = "toy_noise_ar1", ("a", "b")

    def _inputs(self, params, name):
        return []

    def init(self, layout, rng):
        super().init(layout, rng)
        self.draws = DrawBuffer(self.name, 1, "standard_normal")
        self.draws.init(layout, rng, self.rate_div)
        self.writes = (*self.writes, *self.draws.slot_names)

    def bind(self, layout):
        self.draws.bind(layout)
        super().bind(layout)

    def prepare(self, state, params, n0, n_samples):
        if self.toggles["enabled"]:
            self.draws.refill(state, n0, n_samples)

    def _kernel(self, layout):
        o, ia, ib = self.o[self.params["out"]], self.p["a"], self.p["b"]
        buf, idx = self.draws.buf_off, self.draws.idx_off
        if not self.toggles["enabled"]:
            @numba.njit
            def k_off(S, Z, I, P):  # noqa: N803
                S[o] = 0.0
            return k_off

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            w = S[buf + I[idx]]
            I[idx] += 1
            S[o] = P[ia] * S[o] + P[ib] * w
        return k

    def __call__(self, state, params):
        out = self.params["out"]
        if not self.toggles["enabled"]:
            state[out] = 0.0
            return state
        w = self.draws.take(state)
        state[out] = float(params[self.param_name("a")]) * float(state[out]) + float(params[self.param_name("b")]) * w
        return state


class Sine(ToyBlock):
    KIND, RUNTIME = "toy_sine", ("amp", "freq_hz")

    def _inputs(self, params, name):
        return ["t"]

    def _kernel(self, layout):
        o, it, ia, iff = self.o[self.params["out"]], layout.slot("t").offset, self.p["amp"], self.p["freq_hz"]
        on = self.toggles["enabled"]

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            S[o] = P[ia] * math.sin(2.0 * math.pi * P[iff] * S[it]) if on else 0.0
        return k

    def __call__(self, state, params):
        amp, f = float(params[self.param_name("amp")]), float(params[self.param_name("freq_hz")])
        state[self.params["out"]] = amp * math.sin(2.0 * math.pi * f * state.t) if self.toggles["enabled"] else 0.0
        return state


class Sum(ToyBlock):
    KIND, RUNTIME = "toy_sum", ("gains",)

    def _inputs(self, params, name):
        require(params, "ins", "gains", block=name)
        if len(params["ins"]) != len(params["gains"]) or not params["ins"]:
            raise ValueError(f"block {name!r}: ins and gains need equal non-zero length")
        return list(params["ins"])

    def init(self, layout, rng):
        BaseBlock.init(self, layout, rng)  # gains has shape (k,)
        layout.add(self.params["out"], (), np.float64, str(self.params["unit"]))
        self.h = block_dt_s(layout.f_sim_hz, self.rate_div)

    def _kernel(self, layout):
        o = self.o[self.params["out"]]
        src = np.array([layout.slot(s).offset for s in self.params["ins"]], dtype=np.int64)
        g0, n_in, on = self.p["gains"], len(src), self.toggles["enabled"]

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            acc = 0.0
            for j in range(n_in):
                acc += (P[g0 + j] if on else 1.0) * S[src[j]]
            S[o] = acc
        return k

    def __call__(self, state, params):
        g = params[self.param_name("gains")]
        acc = 0.0
        for j, s in enumerate(self.params["ins"]):
            acc += (float(g[j]) if self.toggles["enabled"] else 1.0) * float(state[s])
        state[self.params["out"]] = acc
        return state


class Gain(ToyBlock):
    KIND, RUNTIME = "toy_gain", ("gain",)

    def _kernel(self, layout):
        i, o, ig, on = self.o[self.params["in"]], self.o[self.params["out"]], self.p["gain"], self.toggles["enabled"]

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            S[o] = (P[ig] if on else 1.0) * S[i]
        return k

    def __call__(self, state, params):
        g = float(params[self.param_name("gain")]) if self.toggles["enabled"] else 1.0
        state[self.params["out"]] = g * float(state[self.params["in"]])
        return state


class PID(ToyBlock):
    KIND, RUNTIME = "toy_pid", ("kp", "ki", "kd", "setpoint")

    def _state_slots(self, name):
        return (f"ctrl.{name}.integ", f"ctrl.{name}.eprev")

    def _kernel(self, layout):
        i, o = self.o[self.params["in"]], self.o[self.params["out"]]
        ii, ie = self.o[self.state_slots[0]], self.o[self.state_slots[1]]
        kp, ki, kd, sp = self.p["kp"], self.p["ki"], self.p["kd"], self.p["setpoint"]
        h, on = self.h, self.toggles["enabled"]

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            if on:
                e = P[sp] - S[i]
                S[ii] += P[ki] * h * e
                S[o] = P[kp] * e + S[ii] + P[kd] * (e - S[ie]) / h
                S[ie] = e
        return k

    def __call__(self, state, params):
        if not self.toggles["enabled"]:
            return state
        g = lambda key: float(params[self.param_name(key)])  # noqa: E731
        integ, eprev = self.state_slots
        e = g("setpoint") - float(state[self.params["in"]])
        state[integ] = float(state[integ]) + g("ki") * self.h * e
        state[self.params["out"]] = g("kp") * e + float(state[integ]) + g("kd") * (e - float(state[eprev])) / self.h
        state[eprev] = e
        return state


class Sat(ToyBlock):
    KIND, RUNTIME = "toy_sat", ("gain", "limit")

    def _kernel(self, layout):
        i, o, ig, il, on = (self.o[self.params["in"]], self.o[self.params["out"]], self.p["gain"], self.p["limit"],
                            self.toggles["enabled"])

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            y = P[ig] * S[i]
            if on:
                y = min(max(y, -P[il]), P[il])
            S[o] = y
        return k

    def __call__(self, state, params):
        y = float(params[self.param_name("gain")]) * float(state[self.params["in"]])
        if self.toggles["enabled"]:
            lim = float(params[self.param_name("limit")])
            y = min(max(y, -lim), lim)
        state[self.params["out"]] = y
        return state


class Integrator(ToyBlock):
    KIND, RUNTIME = "toy_integrator", ("k",)

    def _kernel(self, layout):
        i, o, ik, h, on = self.o[self.params["in"]], self.o[self.params["out"]], self.p["k"], self.h, self.toggles["enabled"]

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            if on:
                S[o] += h * P[ik] * S[i]
        return k

    def __call__(self, state, params):
        if self.toggles["enabled"]:
            out = self.params["out"]
            state[out] = float(state[out]) + self.h * float(params[self.param_name("k")]) * float(state[self.params["in"]])
        return state


def _register(cls):
    @register_block(cls.KIND)
    def factory(params, *, name=None):
        return cls(params, name)

    factory.__name__ = factory.__qualname__ = f"make_{cls.KIND}"
    return factory


make_toy_noise_ar1 = _register(NoiseAR1)
make_toy_sine = _register(Sine)
make_toy_sum = _register(Sum)
make_toy_gain = _register(Gain)
make_toy_pid = _register(PID)
make_toy_sat = _register(Sat)
make_toy_integrator = _register(Integrator)


def toy(kind: str, name: str, rate_div: int = 1, *, enabled: bool = True, **params: Any) -> dict[str, Any]:
    """Pipeline dict element for a toy block: ``toy("toy_gain", "g", in_="v.x", out="v.y", unit="V", gain=2.0)``.

    A trailing underscore in a keyword is dropped (``in_`` → ``in``).
    """
    p = {k.rstrip("_"): v for k, v in params.items()}
    p["toggles"] = {"enabled": enabled}
    obj: dict[str, Any] = {"block": kind, "name": name, "params": p}
    if rate_div != 1:
        obj["rate_div"] = rate_div
    return obj
