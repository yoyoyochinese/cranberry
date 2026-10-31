"""Test-support driver blocks for the path model (WP-F6). Not physics blocks.

Importing this module registers ``pathtest_table``, ``pathtest_white`` and ``pathtest_record`` (``core.registry.
ensure_loaded`` does not import it). Both have a compiled kernel and a pure-Python step.

=================  ==================================================================================
kind               step
=================  ==================================================================================
pathtest_table     out ← table[n] for n < n_rows, else 0 (table = slot ``<out>_table`` of shape
                   (n_rows, *shape), filled by the test before running; impulse / step schedules)
pathtest_white     out ← sigma·w, w ~ N(0, 1) of the out shape, from a DrawBuffer (G11)
pathtest_record    row n of rec.<name>.buf ← the ``ins`` scalar slots (n < n_rows)
=================  ==================================================================================

Helpers: ``pathtest`` / ``record`` (pipeline elements), ``run_and_record`` (build + fill tables + run),
``seg_spec`` / ``step_spec`` / ``beam_spec`` / ``coef_spec`` / ``topo_spec`` / ``path_pipeline``
(topology specs for tests; presets must still give every value explicitly, D34).

Params: ``out`` (slot), ``shape`` (list, [] for scalar), ``unit``; ``n_rows`` (table) / ``sigma``
(white, runtime).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["TableDriver", "WhiteDriver", "Recorder", "pathtest", "record", "run_and_record", "seg_spec", "step_spec",
           "beam_spec", "coef_spec", "topo_spec", "path_pipeline"]


class _Driver(BaseBlock):
    KIND = ""
    RUNTIME: tuple[str, ...] = ()

    def __init__(self, params: Mapping[str, Any], name: str | None, extra_writes=()):
        name = name or self.KIND
        require(params, "out", "shape", "unit", *self.RUNTIME, block=name)
        self.shape = tuple(int(x) for x in params["shape"])
        self.size = int(np.prod(self.shape, dtype=np.int64)) if self.shape else 1
        super().__init__(self.KIND, name, params, reads=("n",), writes=(params["out"], *extra_writes),
                         runtime_params=self.RUNTIME)

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"], self.shape, np.float64, str(self.params["unit"]))


class TableDriver(_Driver):
    KIND = "pathtest_table"

    def __init__(self, params, name):
        require(params, "out", "n_rows", block=name or self.KIND)
        self.table = f"{params['out']}_table"
        super().__init__(params, name, (self.table,))
        self.n_rows = int(params["n_rows"])

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.table, (self.n_rows, *self.shape), np.float64, str(self.params["unit"]))

    def bind(self, layout):
        o, tab, n_off = layout.slot(self.params["out"]).offset, layout.slot(self.table).offset, layout.slot("n").offset
        size, rows = self.size, self.n_rows

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            r = I[n_off]
            for j in range(size):
                S[o + j] = S[tab + r * size + j] if r < rows else 0.0

        self.kernel = kernel

    def __call__(self, state, params):
        r = state.n
        state[self.params["out"]] = state[self.table][r] if r < self.n_rows else 0.0
        return state


class WhiteDriver(_Driver):
    KIND, RUNTIME = "pathtest_white", ("sigma",)

    def init(self, layout, rng):
        super().init(layout, rng)
        self.draws = DrawBuffer(self.name, self.size, "standard_normal")
        self.draws.init(layout, rng, self.rate_div)
        self.writes = (*self.writes, *self.draws.slot_names)

    def prepare(self, state, params, n0, n_samples):
        self.draws.refill(state, n0, n_samples)

    def bind(self, layout):
        self.draws.bind(layout)
        o, size = layout.slot(self.params["out"]).offset, self.size
        buf, idx, isg = self.draws.buf_off, self.draws.idx_off, layout.param(self.param_name("sigma")).offset

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            for j in range(size):
                S[o + j] = P[isg] * S[buf + I[idx]]
                I[idx] += 1

        self.kernel = kernel

    def __call__(self, state, params):
        sigma = float(params[self.param_name("sigma")])
        vals = [sigma * self.draws.take(state) for _ in range(self.size)]
        state[self.params["out"]] = np.array(vals).reshape(self.shape) if self.shape else vals[0]
        return state


class Recorder(BaseBlock):
    """``pathtest_record``: row n of ``rec.<name>.buf`` (shape (n_rows, len(ins))) ← the ``ins`` scalar
    slots at sample n (n < n_rows). Place it after the writers of ``ins``."""

    KIND = "pathtest_record"

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "ins", "n_rows", block=name)
        self.buf = f"rec.{name}.buf"
        self.n_rows = int(params["n_rows"])
        super().__init__(self.KIND, name, params, reads=(*params["ins"], "n"), writes=(self.buf,))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.buf, (self.n_rows, len(self.params["ins"])), np.float64, "1")

    def bind(self, layout):
        src = np.array([layout.slot(s).offset for s in self.params["ins"]], dtype=np.int64)
        buf, n_off, k, rows = layout.slot(self.buf).offset, layout.slot("n").offset, len(src), self.n_rows

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            r = I[n_off]
            if r < rows:
                for j in range(k):
                    S[buf + r * k + j] = S[src[j]]

        self.kernel = kernel

    def __call__(self, state, params):
        r = state.n
        if r < self.n_rows:
            state[self.buf][r] = [np.asarray(state[s]).item() for s in self.params["ins"]]  # scalar or 1-element slot
        return state

    def series(self, state) -> dict[str, np.ndarray]:
        """{slot: recorded series} (copies)."""
        buf = np.array(state[self.buf])
        return {s: buf[:, j] for j, s in enumerate(self.params["ins"])}


def _register(cls):
    @register_block(cls.KIND)
    def factory(params, *, name=None):
        return cls(params, name)

    factory.__name__ = factory.__qualname__ = f"make_{cls.KIND}"
    return factory


make_pathtest_table = _register(TableDriver)
make_pathtest_white = _register(WhiteDriver)
make_pathtest_record = _register(Recorder)


def record(ins, n_rows: int, name: str = "rec") -> dict[str, Any]:
    """Pipeline dict element for a ``pathtest_record`` block."""
    return {"block": "pathtest_record", "name": name, "params": {"ins": list(ins), "n_rows": int(n_rows)}}


def run_and_record(pipe, fs_hz: float, n_steps: int, tables: Mapping[str, Any] | None = None, *,
                   compiled: bool = False, seed: int = 0, external=(), setup=None, values=None):
    """Build ``pipe`` (its last element a ``pathtest_record`` block), fill ``tables`` ({slot: array}) and
    ``values`` ({slot: value}), run ``n_steps`` samples and return ({slot: series}, state, params)."""
    from cranberry.core.compile import run_compiled, run_reference
    from cranberry.core.layout import SlotLayout
    from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
    from cranberry.core.state import Params, State

    lay = SlotLayout(fs_hz)
    if setup is not None:
        setup(lay)
    blk = build(spec_from_obj(pipe), lay, seed, external=external)
    st, pa = State(lay), Params(lay)
    for k, v in (tables or {}).items():
        st[k] = v
    for k, v in (values or {}).items():
        st[k] = v
    (run_compiled if compiled else run_reference)(blk, st, pa, n_steps)
    rec = [leaf for leaf in iter_leaves(blk) if isinstance(leaf, Recorder)][-1]
    return rec.series(st), st, pa


def pathtest(kind: str, name: str, out: str, shape=(), unit: str = "1", **params: Any) -> dict[str, Any]:
    """Pipeline dict element for a test driver block."""
    return {"block": kind, "name": name, "params": {"out": out, "shape": list(shape), "unit": unit, **params}}


# ------------------------------------------------------------------------------------------------
# topology-spec builders for tests (test conveniences only: presets give every value explicitly)


def seg_spec(id: str, L_m: float, *, n_sub: int = 1, drivers: Mapping[str, Any] | None = None, pol_model: str = "off",
             pol_params: Mapping[str, Any] | None = None, spatial: str = "lumped", n_g: float = 1.468,
             actuator_sub: int | None = None, kind: str = "pm") -> dict[str, Any]:
    """A segment entry of a topology spec."""
    out = {"id": id, "kind": kind, "L_m": L_m, "n_g": n_g, "beta2_s2_per_m": 0.0, "beta3_s3_per_m": 0.0,
           "zone": "z0", "n_sub": n_sub, "spatial": spatial, "pol_model": pol_model,
           "pol_params": dict(pol_params or {}), "drivers": dict(drivers or {})}
    if actuator_sub is not None:
        out["actuator_sub"] = actuator_sub
    return out


def step_spec(segment: str, direction: int = 1, *, axis: str | None = None, sub_in: int | None = None,
              sub_out: int | None = None) -> dict[str, Any]:
    """A path step of a beam spec."""
    out: dict[str, Any] = {"segment": segment, "direction": direction}
    if axis is not None:
        out["axis"] = axis
    if sub_in is not None:
        out["sub_in"] = sub_in
    if sub_out is not None:
        out["sub_out"] = sub_out
    return out


def beam_spec(id: str, path: list, observer: str, *, source: str = "las", jones_in=(1, 0), jones_out=None,
              freq_offset_hz: float = 0.0, power_w: float = 1e-3) -> dict[str, Any]:
    """A beam entry of a topology spec."""
    out = {"id": id, "source": source, "path": list(path), "freq_offset_hz": freq_offset_hz,
           "jones_in": list(jones_in), "power_w": power_w, "observer": observer}
    if jones_out is not None:
        out["jones_out"] = list(jones_out)
    return out


def coef_spec(segment: str, driver: str, beam: str, phase: Any, group_delay: Any = 0.0) -> list[dict[str, Any]]:
    """Phase and group-delay coefficient entries for one (segment, driver, beam)."""
    return [{"segment": segment, "driver": driver, "beam": beam, "effect": "phase", "value": phase},
            {"segment": segment, "driver": driver, "beam": beam, "effect": "group_delay", "value": group_delay}]


def topo_spec(segments, beams, observers, coef, locks=(), ports=()) -> dict[str, Any]:
    """A topology spec; ``observers`` may be ids (kind detector) or full entries."""
    obs = [o if isinstance(o, Mapping) else {"id": o, "kind": "detector"} for o in observers]
    return {"segments": list(segments), "beams": list(beams), "observers": obs, "locks": list(locks),
            "ports": list(ports), "coef": list(coef)}


def path_pipeline(topo: Mapping[str, Any], sources: Mapping[str, str], *, phase: bool = True,
                  group_delay: bool = True, pm_axis_drift: bool = True, emission_time: bool = True,
                  polarization: bool = True, locks: bool = False) -> list[dict[str, Any]]:
    """``segment_phase`` → ``beam_assembler`` (→ ``lock_error``) pipeline elements for ``topo``."""
    out = [{"block": "segment_phase", "name": "sp", "params": {
                "topology": topo, "segments": "all",
                "toggles": {"phase": phase, "group_delay": group_delay, "pm_axis_drift": pm_axis_drift}}},
           {"block": "beam_assembler", "name": "ba", "params": {
                "topology": topo, "beams": "all", "sources": dict(sources),
                "toggles": {"emission_time": emission_time, "polarization": polarization}}}]
    if locks:
        out.append({"block": "lock_error", "name": "le", "params": {"topology": topo, "locks": "all"}})
    return out
