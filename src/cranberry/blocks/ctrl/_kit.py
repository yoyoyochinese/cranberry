"""Shared plumbing of the WP-B4 blocks (CTRL / ACT / MOD). Not a block module.

One step function per block kind, written once in numba-compatible Python with the signature
``step(S, I, P, *args)``. The pure-Python reference ``__call__`` runs it interpreted on the State
arrays; the kernel is the same source compiled (``numba.njit``) and closed over the same ``args``
(the core.delay pattern). ``args`` holds only slot / parameter offsets, structural constants and
structural coefficient arrays (D55 rule 3: kernels capture no runtime values or random numbers).

Common conventions of the B4 blocks
    * Wiring by slot name: ``in`` (read) and ``out`` (written, registered by the block).
    * Runtime parameters are float entries of P; every one of them must be in the spec params
      (``SimController.get_params`` reads them from there).
    * ``device: <name> | null``: structural device facts (bits, ranges, discretisation, ...) are read
      with ``device_values(load_device(name))`` and may be overridden by the spec (D34, prompt §5.3).
    * Segment drivers (G08, PLAN §2.3.2, F6 slot names): an actuator / modulator writes
      ``drv.<segment>.<driver>`` of shape ``(n_sub,)``, non-zero only on subsegments ``sub = [lo, hi)``.
      Driver values are *intensive* (per-subsegment densities), matching F6's ``per_sub`` placement,
      which reads drv[m]/N_sub (D61): phase-like drivers (stretcher, eom, ...) write
      δψ·N_sub/(hi − lo) on ``[lo, hi)``, so that F6's Σ drv[m]/N_sub is the actuator's δψ (with a
      1-element slot, N_sub = 1, this is just δψ); temperature (TEC) is copied to every subsegment of
      the range. ``ring_len > 0`` makes the slot a history ring (head = newest row, ``core.delay``)
      for emission-time reads by ``beam_assembler`` (then the block must run at ``rate_div = 1``).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, MissingParameterError, require
from cranberry.core.rates import block_dt_s

__all__ = ["KitBlock", "make_kernel", "device_section", "DrvTarget", "drv_write", "drv_write_nb", "loss_elements",
           "LOSS_KEYS"]


def make_kernel(step_nb, args: tuple):
    """Compiled kernel ``kernel(S, Z, I, P)`` calling ``step_nb(S, I, P, *args)``."""

    @numba.njit
    def kernel(S, Z, I, P):  # noqa: N803
        step_nb(S, I, P, *args)

    return kernel


def device_section(params: Mapping[str, Any], section: str, block: str) -> dict[str, Any]:
    """Plain values of ``devices/<params['device']>.yaml`` section ``section`` ({} for ``device: null``)."""
    require(params, "device", block=block)
    if params["device"] is None:
        return {}
    from cranberry.hardware.devices import device_values, load_device

    dev = device_values(load_device(str(params["device"])))
    if section not in dev:
        raise KeyError(f"block {block!r}: device {params['device']!r} has no section {section!r}")
    return {k: v for k, v in dev[section].items() if v is not None}


class KitBlock(BaseBlock):
    """Leaf block driven by one step function (module docstring).

    Subclasses set ``KIND``, ``STEP`` (pure Python) and ``STEP_NB`` (compiled), list their slots in
    :meth:`slots` and build ``args`` in :meth:`make_args`.

    Attributes:
        values: spec params merged over device values (structural reads use this).
        h: own time step rate_div / f_sim [s] (set in ``init``).
    """

    KIND = ""
    STEP = None
    STEP_NB = None

    def __init__(self, params: Mapping[str, Any], name: str | None, *, reads: Iterable[str], writes: Iterable[str],
                 runtime: Iterable[str] = (), device: Mapping[str, Any] | None = None):
        name = name or self.KIND
        self.values: dict[str, Any] = {**(device or {}), **dict(params)}
        super().__init__(self.KIND, name, params, reads=tuple(reads), writes=tuple(writes), runtime_params=runtime)

    def need(self, *keys: str) -> None:
        """Raise MissingParameterError unless every key is in spec params or device values (D34)."""
        missing = [k for k in keys if k not in self.values]
        if missing:
            raise MissingParameterError(f"missing required parameter(s) for block {self.name!r}: {', '.join(missing)} "
                                        "(D34: no code defaults; give it in the spec or the device preset)")

    # slots: (name, shape, dtype, unit, ring_len, init)
    def slots(self) -> list[tuple]:
        return []

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        self.h = block_dt_s(layout.f_sim_hz, self.rate_div)
        for nm, shape, dtype, unit, ring, init in self.slots():
            layout.add(nm, shape, dtype, unit, ring, init=init)

    def make_args(self, layout) -> tuple:  # pragma: no cover - overridden
        raise NotImplementedError

    def p(self, layout, key: str) -> int:
        """Offset in P of runtime parameter ``key``."""
        return layout.param(self.param_name(key)).offset

    @staticmethod
    def o(layout, slot: str | None) -> int:
        """Offset of slot (−1 for None)."""
        return -1 if slot is None else layout.slot(slot).offset

    def bind(self, layout) -> None:
        self.args = self.make_args(layout)
        self.kernel = make_kernel(self.STEP_NB, self.args)

    def __call__(self, state, params):
        type(self).STEP(state.S, state.I, params.P, *self.args)
        return state


class DrvTarget:
    """Segment-driver output ``drv.<segment>.<driver>`` (module docstring). Params: ``segment``,
    ``driver``, ``n_sub``, ``sub`` = [lo, hi), ``ring_len``."""

    KEYS = ("segment", "driver", "n_sub", "sub", "ring_len")

    def __init__(self, params: Mapping[str, Any], block: str, unit: str, split: bool):
        require(params, *self.KEYS, block=block)
        self.slot = f"drv.{params['segment']}.{params['driver']}"
        self.n_sub = int(params["n_sub"])
        lo, hi = (int(v) for v in params["sub"])
        if not (self.n_sub >= 1 and 0 <= lo < hi <= self.n_sub):
            raise ValueError(f"block {block!r}: need n_sub >= 1 and 0 <= sub[0] < sub[1] <= n_sub, got "
                             f"n_sub={params['n_sub']!r}, sub={params['sub']!r}")
        self.lo, self.hi = lo, hi
        self.ring_len = int(params["ring_len"])
        if self.ring_len < 0:
            raise ValueError(f"block {block!r}: ring_len must be >= 0")
        self.unit = unit
        # intensive density: F6 per_sub reads drv[m]/N_sub, so Σ over [lo, hi) gives δψ (D61)
        self.scale = self.n_sub / (hi - lo) if split else 1.0
        self.block = block

    def slot_def(self) -> tuple:
        return (self.slot, (self.n_sub,), np.float64, self.unit, self.ring_len, 0.0)

    def check_rate(self, rate_div: int) -> None:
        if self.ring_len and rate_div != 1:
            raise ValueError(f"block {self.block!r}: a driver history ring (ring_len > 0) needs rate_div = 1")

    def args(self, layout) -> tuple:
        """(off, n_sub, lo, hi, ring_len, head, scale) for :func:`drv_write`."""
        s = layout.slot(self.slot)
        return (s.offset, self.n_sub, self.lo, self.hi, self.ring_len, s.head, self.scale)


def drv_write(S, I, off, n_sub, lo, hi, ring_len, head, scale, value):  # noqa: N803
    """Write ``value * scale`` into subsegments [lo, hi) (pushing a ring row first if ring_len > 0)."""
    base = off
    if ring_len > 0:
        r = I[head] + 1
        if r == ring_len:
            r = 0
        I[head] = r
        base = off + r * n_sub
    v = value * scale
    for j in range(lo, hi):
        S[base + j] = v


drv_write_nb = numba.njit(inline="always")(drv_write)

#: Params of the optional insertion-loss registration (G06).
LOSS_KEYS = ("insertion_loss_db", "loss_applies_to", "loss_source")


def loss_elements(block: BaseBlock, toggle_on: bool) -> list:
    """``[LossElement]`` from ``insertion_loss_db`` (None → no element) on ``params['segment']`` (G06).

    ``insertion_loss_lin`` (loss factor 10^(IL/10) ≥ 1, or None; D56-1) may replace
    ``insertion_loss_db``: YAML presets cannot give a null under a dB name, and a numeric
    ``insertion_loss_db`` is renamed to ``insertion_loss_lin`` on load (``config.db``; WP-I2 addition).

    Registered whatever the beam: the LossRegistry decides whether it lies on a squeezed / LO / probe
    path (prompt §5.5). ``toggle_on`` False keeps the record with η treated as 1.
    """
    import math

    from cranberry.offline.loss_registry import LossElement

    prm = block.params
    if "insertion_loss_lin" in prm and "insertion_loss_db" not in prm:
        il = prm["insertion_loss_lin"]
        il_db = None if il is None else 10.0 * math.log10(float(il))
    else:
        require(prm, "insertion_loss_db", block=block.name)
        il_db = None if prm["insertion_loss_db"] is None else float(prm["insertion_loss_db"])
    if il_db is None:
        return []
    require(prm, "loss_applies_to", "loss_source", "segment", block=block.name)
    return [LossElement.from_il_db(f"{block.name}.il", il_db, kind="component",
                                   applies_to=frozenset(prm["loss_applies_to"]), source=str(prm["loss_source"]),
                                   enabled=bool(toggle_on), segment=str(prm["segment"]))]
