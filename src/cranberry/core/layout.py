"""Slot registry: flat typed arrays + named slots (G02, PLAN §2.3.2, §4.2).

The dynamic ``State`` of a time-loop run is three flat arrays plus one for runtime parameters:

    S : float64[n_float]      real-valued slots
    Z : complex128[n_complex] complex slots (Jones matrices/vectors, I/Q pairs, ...)
    I : int64[n_int]          integer slots (sample counter, flags, ring heads, lock state)
    P : float64[n_param]      runtime parameters (``set_params``-changeable, e.g. PID gains)

A ``SlotLayout`` assigns every named slot a contiguous region of one of these arrays. Compiled
kernels (``@njit`` closures, WP-F1) capture the integer offsets as constants; the pure-Python
reference implementation and the recorder read the same memory through named views
(``state["ctrl.lock1.u"]``, ``state.ctrl.lock1.u``; see ``core.state``).

Slot naming scheme (contract):

* A slot name is a dotted path ``<group>.<owner>[.<field>...]`` matching ``[A-Za-z][A-Za-z0-9_]*``
  per component. ``<group>`` must be one of :data:`SLOT_GROUPS` (the State field table, PLAN
  §2.3.2); ``t`` and ``n`` are the only single-component names and are always present (written by
  the runner). Examples:

  =======================  ===========================================  ==========  ======
  name                     meaning                                       array       unit
  =======================  ===========================================  ==========  ======
  ``t``                    time of the current sample                    S           s
  ``n``                    sample index                                  I           1
  ``laser.dx1.phi``        laser phase φ_L (+ ring: emission-time hist)  S           rad
  ``drv.arm1.dT``          per-subsegment driver (shape (N_sub,), ring)  S           K
  ``psi.arm1.pr``          segment × beam phase ψ                        S           rad
  ``tau_g.arm1.pr``           segment group delay τ_g(t)                    S           s
  ``J.arm1``               segment Jones matrix, shape (2, 2)            Z           1
  ``beam.lo.phi``          beam phase at its observer                    S           rad
  ``v.pdb450c``            a voltage                                     S           V
  ``dem.l1.i``             demodulated I                                 S           V
  ``ctrl.lock1.u``         controller output                             S           V
  ``delay.lock1``          delay-buffer ring                             S           V
  ``lock.lockD.state``     lock state / event counter                    I           1
  ``metric.port1.delta``   metric (i) sample                             S           rad
  ``rng.n1.buf``           pre-drawn random numbers of block ``n1``      S           1
  ``rec.recorder.buf``     runner's decimated recording buffer           S           1
  =======================  ===========================================  ==========  ======

  ``<owner>`` is the block/segment/beam/lock id (block names may themselves contain dots, e.g.
  ``ctrl.lock1.pid.u``). Adding a group to :data:`SLOT_GROUPS` is a one-line additive change;
  propose it in the WP report.

* Ring buffers (history for emission-time evaluation, delay lines; PLAN §2.6): ``ring_len > 0``
  allocates storage of shape ``(ring_len, *shape)`` and one int64 head index (not a named slot;
  read via ``Slot.head`` / ``State.ring_head``). Head semantics (``core.delay``): the head is the
  row holding the **newest** sample; a push advances the head by one (mod ``ring_len``) and then
  writes that row; the value ``d`` samples older than the newest is row ``(head - d) % ring_len``.
  The head starts at 0 and every row starts at ``init``.

* Runtime parameters are named ``<block name>.<param key>`` (e.g. ``lock1.pid.kp``) and live in P.
  Only float64 runtime parameters exist; anything that changes array sizes or kernel structure is a
  structural parameter (G02) and is fixed when the block is made.

Lifecycle: blocks call :meth:`SlotLayout.add` / :meth:`SlotLayout.add_param` in ``Block.init``;
the builder then calls :meth:`SlotLayout.freeze`; after that ``State(layout)`` / ``Params(layout)``
allocate the arrays and ``Block.bind(layout)`` resolves offsets.

``f_sim_hz`` (the base sample rate, G10) is carried by the layout because it is the one build-time
global every block needs to discretise itself; a block with ``rate_div = r`` steps at f_sim/r.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

__all__ = ["SLOT_GROUPS", "RESERVED_SLOTS", "Slot", "SlotLayout", "LayoutFrozenError"]

#: Allowed first components of slot names (PLAN §2.3.2 State fields, plus ``rf`` from §4.3 NOISE,
#: ``rng`` for pre-drawn random-number buffers (``core.rng.DrawBuffer``) and ``rec`` for the
#: runner's recorder buffers; the last two added by WP-F1).
SLOT_GROUPS: frozenset[str] = frozenset(
    {
        "t", "n", "laser", "drv", "psi", "tau_g", "J", "beam", "pump", "opa", "opo", "rf",
        "v", "dem", "ctrl", "delay", "act", "timing", "leak", "lock", "metric",
        "rng", "rec",
    }
)
#: Slots that every layout has, written by the runner.
RESERVED_SLOTS: tuple[str, ...] = ("t", "n")

_COMPONENT = r"[A-Za-z][A-Za-z0-9_]*"
_NAME_RE = re.compile(rf"^{_COMPONENT}(\.{_COMPONENT})*$")

_DTYPE_ARRAY: dict[np.dtype, Literal["S", "Z", "I"]] = {
    np.dtype(np.float64): "S",
    np.dtype(np.complex128): "Z",
    np.dtype(np.int64): "I",
}


class LayoutFrozenError(RuntimeError):
    """Raised when a frozen layout is modified, or an unfrozen one is used to allocate."""


@dataclass(frozen=True)
class Slot:
    """One named region of a flat array.

    Attributes:
        name: dotted slot name (see module docstring).
        array: which flat array holds it: "S" (float64), "Z" (complex128), "I" (int64), "P" (params).
        shape: logical shape of one sample (``()`` for a scalar).
        unit: unit tag (``core.units.UNIT_TAGS``; "1" = dimensionless).
        ring_len: 0 for a plain slot, else number of history rows.
        offset: start index in the flat array.
        size: number of elements (``prod(storage_shape)``).
        head: index into I of the ring head, or -1 if ``ring_len == 0``.
        init: initial fill value (plain slots) / default value (params).
    """

    name: str
    array: Literal["S", "Z", "I", "P"]
    shape: tuple[int, ...]
    unit: str
    ring_len: int
    offset: int
    size: int
    head: int
    init: Any

    @property
    def storage_shape(self) -> tuple[int, ...]:
        """``(ring_len, *shape)`` for a ring slot, else ``shape``."""
        return (self.ring_len, *self.shape) if self.ring_len else self.shape

    @property
    def stop(self) -> int:
        return self.offset + self.size


def _check_name(name: str, *, slot: bool) -> None:
    if not isinstance(name, str) or not _NAME_RE.match(name):
        raise ValueError(f"invalid name {name!r}: dotted path of [A-Za-z][A-Za-z0-9_]* components")
    if slot:
        parts = name.split(".")
        if parts[0] not in SLOT_GROUPS:
            raise ValueError(f"slot {name!r}: unknown group {parts[0]!r}; allowed {sorted(SLOT_GROUPS)}")
        if len(parts) == 1 and name not in RESERVED_SLOTS:
            raise ValueError(f"slot {name!r}: needs <group>.<owner>[...]")


def _norm_shape(shape: int | tuple[int, ...]) -> tuple[int, ...]:
    shp = (shape,) if isinstance(shape, (int, np.integer)) else tuple(shape)
    if any(int(d) != d or d < 1 for d in shp):
        raise ValueError(f"shape entries must be positive ints, got {shape!r}")
    return tuple(int(d) for d in shp)


class SlotLayout:
    """Build-time registry of State slots and runtime parameters (G02).

    Args:
        f_sim_hz: base sample rate of the time loop [Hz] (``RunSpec.run.f_sim_hz``).
    """

    def __init__(self, f_sim_hz: float):
        if not (np.isfinite(f_sim_hz) and f_sim_hz > 0):
            raise ValueError(f"f_sim_hz must be finite and > 0, got {f_sim_hz!r}")
        self.f_sim_hz = float(f_sim_hz)
        self._slots: dict[str, Slot] = {}
        self._params: dict[str, Slot] = {}
        self._sizes = {"S": 0, "Z": 0, "I": 0, "P": 0}
        self._frozen = False
        self.add("t", (), np.float64, "s")
        self.add("n", (), np.int64, "1")

    @property
    def dt_s(self) -> float:
        """Base time step 1/f_sim [s]."""
        return 1.0 / self.f_sim_hz

    @property
    def frozen(self) -> bool:
        return self._frozen

    def _alloc(self, array: str, size: int) -> int:
        off = self._sizes[array]
        self._sizes[array] += size
        return off

    def add(
        self,
        name: str,
        shape: int | tuple[int, ...] = (),
        dtype: Any = np.float64,
        unit: str = "1",
        ring_len: int = 0,
        *,
        init: Any = 0,
    ) -> Slot:
        """Register slot ``name`` and return its :class:`Slot`.

        Args:
            name: dotted slot name (module docstring). Must be unique.
            shape: per-sample shape; ``()`` for a scalar.
            dtype: float64, complex128 or int64 (selects S, Z or I).
            unit: unit tag (``core.units.UNIT_TAGS``).
            ring_len: history length; > 0 makes a ring slot with storage ``(ring_len, *shape)``.
            init: initial value filled at allocation.
        """
        if self._frozen:
            raise LayoutFrozenError(f"cannot add slot {name!r}: layout is frozen")
        _check_name(name, slot=True)
        if name in self._slots:
            raise ValueError(f"slot {name!r} already registered")
        dt = np.dtype(dtype)
        if dt not in _DTYPE_ARRAY:
            raise TypeError(f"slot {name!r}: dtype must be float64, complex128 or int64, got {dt}")
        if int(ring_len) != ring_len or ring_len < 0:
            raise ValueError(f"slot {name!r}: ring_len must be a non-negative int")
        shp = _norm_shape(shape)
        array = _DTYPE_ARRAY[dt]
        size = int(np.prod((ring_len, *shp) if ring_len else shp, dtype=np.int64))
        offset = self._alloc(array, size)
        head = self._alloc("I", 1) if ring_len else -1
        slot = Slot(name, array, shp, str(unit), int(ring_len), offset, size, head, init)
        self._slots[name] = slot
        return slot

    def add_param(self, name: str, default: Any, shape: int | tuple[int, ...] = (), unit: str = "1") -> Slot:
        """Register runtime parameter ``name`` (``<block name>.<key>``) in P with value ``default``.

        ``default`` is the value from the block spec (preset), not a code default (D34).
        """
        if self._frozen:
            raise LayoutFrozenError(f"cannot add param {name!r}: layout is frozen")
        _check_name(name, slot=False)
        if name in self._params:
            raise ValueError(f"param {name!r} already registered")
        shp = _norm_shape(shape)
        size = int(np.prod(shp, dtype=np.int64))
        offset = self._alloc("P", size)
        slot = Slot(name, "P", shp, str(unit), 0, offset, size, -1, default)
        self._params[name] = slot
        return slot

    def freeze(self) -> None:
        """Fix the layout; afterwards offsets are final and arrays may be allocated."""
        self._frozen = True

    def slot(self, name: str) -> Slot:
        try:
            return self._slots[name]
        except KeyError:
            raise KeyError(f"no slot {name!r}") from None

    def param(self, name: str) -> Slot:
        try:
            return self._params[name]
        except KeyError:
            raise KeyError(f"no runtime param {name!r}") from None

    def __contains__(self, name: object) -> bool:
        return name in self._slots

    def has_param(self, name: str) -> bool:
        return name in self._params

    def slot_names(self) -> tuple[str, ...]:
        return tuple(self._slots)

    def param_names(self) -> tuple[str, ...]:
        return tuple(self._params)

    @property
    def sizes(self) -> dict[str, int]:
        """Element counts of the flat arrays: {"S", "Z", "I", "P"}."""
        return dict(self._sizes)

    def _require_frozen(self) -> None:
        if not self._frozen:
            raise LayoutFrozenError("layout must be frozen before allocation")

    def allocate(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Allocate (S, Z, I) filled with each slot's ``init`` (ring heads start at 0)."""
        self._require_frozen()
        arrays = {
            "S": np.zeros(self._sizes["S"], dtype=np.float64),
            "Z": np.zeros(self._sizes["Z"], dtype=np.complex128),
            "I": np.zeros(self._sizes["I"], dtype=np.int64),
        }
        for s in self._slots.values():
            arrays[s.array][s.offset : s.stop] = s.init
        return arrays["S"], arrays["Z"], arrays["I"]

    def allocate_params(self) -> np.ndarray:
        """Allocate P filled with each parameter's spec value."""
        self._require_frozen()
        p = np.zeros(self._sizes["P"], dtype=np.float64)
        for s in self._params.values():
            p[s.offset : s.stop] = np.asarray(s.init, dtype=np.float64).reshape(-1)
        return p
