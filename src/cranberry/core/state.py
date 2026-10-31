"""Dynamic State and runtime Params: named views over flat arrays (G02, PLAN §2.3.2).

``State`` owns the flat arrays S (float64), Z (complex128), I (int64) allocated from a frozen
``SlotLayout``. ``Params`` owns P (float64 runtime parameters). Both expose writable named views:

    state["ctrl.lock1.u"]            -> numpy view of shape slot.storage_shape (0-d for scalars)
    state["ctrl.lock1.u"] = 0.3      -> in-place write
    state.ctrl.lock1.u               -> same view via attribute path
    params["lock1.pid.kp"]           -> view into P

Step contract (prompt §5.2, G02): every block is ``step(state, params) -> state``; the step updates
the arrays **in place** and returns the **same** State object (no copies in the loop). Compiled
kernels receive the raw arrays ``(state.S, state.Z, state.I, params.P)`` instead.

Time (written by the runner only): ``state.n`` is the sample index, ``state.t`` = n / f_sim [s].
"""

from __future__ import annotations

from typing import Any

import numpy as np

from cranberry.core.layout import LayoutFrozenError, Slot, SlotLayout

__all__ = ["State", "Params"]


def _view(arr: np.ndarray, slot: Slot) -> np.ndarray:
    return arr[slot.offset : slot.stop].reshape(slot.storage_shape)


class _SlotPath:
    """Attribute-path proxy: ``state.ctrl.lock1.u`` resolves to the view of ``ctrl.lock1.u``."""

    __slots__ = ("_owner", "_prefix")

    def __init__(self, owner: State, prefix: str):
        object.__setattr__(self, "_owner", owner)
        object.__setattr__(self, "_prefix", prefix)

    def __getattr__(self, attr: str) -> Any:
        name = f"{self._prefix}.{attr}"
        if name in self._owner.layout:
            return self._owner.view(name)
        return _SlotPath(self._owner, name)

    def __setattr__(self, attr: str, value: Any) -> None:
        self._owner[f"{self._prefix}.{attr}"] = value

    def __repr__(self) -> str:
        return f"<slot path {self._prefix!r}>"


class State:
    """Dynamic simulation state: flat arrays S, Z, I plus named views (G02).

    Args:
        layout: a frozen ``SlotLayout``.
    """

    def __init__(self, layout: SlotLayout):
        if not layout.frozen:
            raise LayoutFrozenError("State requires a frozen SlotLayout")
        self.layout = layout
        self.S, self.Z, self.I = layout.allocate()

    def _arr(self, slot: Slot) -> np.ndarray:
        return {"S": self.S, "Z": self.Z, "I": self.I}[slot.array]

    def view(self, name: str) -> np.ndarray:
        """Writable view of slot ``name`` with shape ``slot.storage_shape``."""
        slot = self.layout.slot(name)
        return _view(self._arr(slot), slot)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.view(name)

    def __setitem__(self, name: str, value: Any) -> None:
        self.view(name)[...] = value

    def __getattr__(self, attr: str) -> Any:
        # only called for names not found normally (S, Z, I, layout, methods are found first)
        if attr.startswith("_"):
            raise AttributeError(attr)
        layout = self.__dict__.get("layout")
        if layout is None:
            raise AttributeError(attr)
        if attr in layout:
            return self.view(attr)
        return _SlotPath(self, attr)

    @property
    def n(self) -> int:
        """Current sample index."""
        return int(self.I[self.layout.slot("n").offset])

    @property
    def t(self) -> float:
        """Current time [s]."""
        return float(self.S[self.layout.slot("t").offset])

    def ring_head(self, name: str) -> int:
        """Head index of ring slot ``name`` (semantics defined by WP-F1, core.delay)."""
        slot = self.layout.slot(name)
        if slot.head < 0:
            raise ValueError(f"slot {name!r} is not a ring slot")
        return int(self.I[slot.head])

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(S, Z, I) for compiled kernels."""
        return self.S, self.Z, self.I

    def copy(self) -> State:
        """Deep copy (same layout object, copied arrays). For snapshots/tests, not for the loop."""
        new = object.__new__(State)
        new.layout = self.layout
        new.S, new.Z, new.I = self.S.copy(), self.Z.copy(), self.I.copy()
        return new


class Params:
    """Runtime parameters P (float64) with named views; changed by ``set_params`` (G02, BoTorch).

    Args:
        layout: a frozen ``SlotLayout``.
    """

    def __init__(self, layout: SlotLayout):
        if not layout.frozen:
            raise LayoutFrozenError("Params requires a frozen SlotLayout")
        self.layout = layout
        self.P = layout.allocate_params()

    def view(self, name: str) -> np.ndarray:
        slot = self.layout.param(name)
        return _view(self.P, slot)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.view(name)

    def __setitem__(self, name: str, value: Any) -> None:
        self.view(name)[...] = value

    def set(self, values: dict[str, Any]) -> None:
        """Set several parameters; all names are validated before any value is written."""
        slots = [self.layout.param(k) for k in values]
        for slot, v in zip(slots, values.values()):
            _view(self.P, slot)[...] = v

    def as_dict(self) -> dict[str, float | list]:
        """Plain-Python snapshot {name: float | nested list} (for logs and ``get_params``)."""
        return {k: self.view(k).tolist() for k in self.layout.param_names()}

    def copy(self) -> Params:
        new = object.__new__(Params)
        new.layout = self.layout
        new.P = self.P.copy()
        return new
