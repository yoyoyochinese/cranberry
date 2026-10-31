"""Block contract: the factory interface of the time-loop pipeline (protected list; PLAN §2.4, G02).

Factory
    Every block kind ``<kind>`` has a factory ``make_<kind>(params, *, name=None) -> Block``
    registered with ``@register_block("<kind>")`` (``core.registry``). YAML ``{block: <kind>, ...}``
    and the Python call produce the same object (prompt §5.2). ``name`` defaults to the kind.

Block
    A block is callable as the pure-Python reference step ``block(state, params) -> state`` (updates
    ``state`` in place and returns the same object) and carries its metadata:

    =================  ========================================================================
    ``kind``           registry key, e.g. ``"pid"``
    ``name``           unique instance name, e.g. ``"lock1.pid"`` (also the RNG stream key, G11)
    ``spec``           ``BlockSpec`` (kind, name, params, rate_div): exact YAML round trip
    ``reads``          State slot names read (wired by name from params ``in``/``in_*``)
    ``writes``         State slot names written (``out``/``out_*``); declared in ``init``
    ``rate_div``       step only when ``n % rate_div == 0``; values are held in between (G10)
    ``runtime_params`` keys of ``params`` that are runtime parameters (P array, ``set_params``)
    ``toggles``        {toggle name: bool}, structural on/off flags (protected list)
    ``kernel``         compiled step ``kernel(S, Z, I, P) -> None`` or None (Python only); set
                       by ``bind``; must agree with ``__call__`` to 1e−12 (VT-44)
    ``init``           declare written slots / runtime params into the (unfrozen) layout and keep
                       the RNG stream ``core.rng.stream(root, name)`` passed by the builder
    ``bind``           after ``layout.freeze()``: resolve slot offsets, build ``kernel``
    ``losses()``       ``list[LossElement]`` registered into the ``LossRegistry`` (G06)
    ``bands()``        ``list[BandDecl]`` declared to the ``BandLedger`` (G05)
    =================  ========================================================================

    Optional members (WP-F1; looked up with ``getattr``, so blocks without them stay valid):

    =================  ========================================================================
    ``slot_units``     {slot name: unit tag} the block expects for slots it reads or writes;
                       compared with the layout's tag by ``check_pipeline`` (check 4)
    ``prepare``        ``prepare(state, params, n0, n_samples)``: Python hook run before each
                       piece of <= ``core.rates.PIECE_SAMPLES`` base samples starting at n0
                       (open-loop work, e.g. pre-drawing random numbers, G11). Must consume
                       exactly what the piece needs so results are chunk-length invariant.
    =================  ========================================================================

Kernel rules (G02, ``core.compile``)
    ``kernel(S, Z, I, P)`` may capture only slot / parameter offsets, structural parameters,
    ``f_sim`` and ``rate_div`` as constants. Anything seed-dependent or runtime-changeable lives in
    S / Z / I / P, because a compiled loop is reused across builds with the same structure
    (sweeps, seed ensembles). Kernels must not use ``np.random``; random numbers come from a
    ``core.rng.DrawBuffer`` filled in ``prepare``.

Lifecycle (builder = ``core.pipeline.build``, WP-F1)::

    blk = make_<kind>(params, name=...)          # structural params fixed, validated
    blk.init(layout, core.rng.stream(root, blk.name))
    layout.freeze(); state = State(layout); params = Params(layout)
    blk.bind(layout)
    loop: if state.n % blk.rate_div == 0: blk(state, params)   # or blk.kernel(S, Z, I, P)

Parameters (G02, D34)
    * **Structural** parameters (array sizes, filter orders, number of subsegments, discretisation
      method, toggles, slot wiring) are fixed by the factory and change the kernel. Changing them
      means rebuilding.
    * **Runtime** parameters (gains, setpoints, ...) are float64 entries of P named
      ``<block name>.<key>``; listed in ``runtime_params``; changeable while running.
    * **No code defaults for modelling choices (D34).** Every option and toggle must come from the
      spec (preset / overlay). A factory raises :class:`MissingParameterError` (via
      :func:`require`) instead of silently choosing a value. ``rate_div = 1`` (no decimation) is the
      only spec-level fallback.

Toggles
    A block lists its toggle names in ``toggle_names``; the spec gives every one of them in
    ``params["toggles"]`` (no defaults, D34). A toggle switched off must make the block an identity
    or a zero source, with a unit test (PLAN §7.3). ``RunSpec.toggles`` entries
    ``"<block name>.<toggle>"`` override spec values at build time (used by ``ablate``).

Sign / unit conventions inside blocks follow notes/conventions.md; device-specific conventions
(e.g. DSC1 e = S − v) are converted inside the block (conv §9). Device values come from
``device: <name>`` + ``devices/*.yaml``, never from code (prompt §5.3).
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from cranberry.core.layout import SlotLayout
    from cranberry.core.state import Params, State
    from cranberry.offline.band_ledger import BandDecl
    from cranberry.offline.loss_registry import LossElement

__all__ = [
    "KernelFn",
    "BlockSpec",
    "Block",
    "BaseBlock",
    "MissingParameterError",
    "require",
    "KIND_RE",
]

#: Compiled step on flat arrays: ``kernel(S, Z, I, P) -> None``, slot offsets captured as constants.
KernelFn = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None]

#: Block kinds are lowercase identifiers (they are also YAML keys and ``make_<kind>`` suffixes).
KIND_RE = re.compile(r"^[a-z][a-z0-9_]*$")


class MissingParameterError(KeyError):
    """A required spec parameter (modelling choice, D34) is missing."""


def require(params: Mapping[str, Any], *keys: str, block: str = "") -> None:
    """Raise :class:`MissingParameterError` listing every key of ``keys`` missing from ``params``."""
    missing = [k for k in keys if k not in params]
    if missing:
        where = f" for block {block!r}" if block else ""
        raise MissingParameterError(f"missing required parameter(s){where}: {', '.join(missing)} (D34: no code defaults)")


@dataclass(frozen=True)
class BlockSpec:
    """Pure-data description of one leaf block; the YAML form is ``to_obj()``.

    YAML / dict form (exact round trip, VT-F1)::

        {block: <kind>, name: <name>, params: {...}}          # + rate_div: <int> when != 1

    Attributes:
        kind: registry key.
        name: unique instance name.
        params: structural + runtime parameters as given in the config (after overlay/callable
            resolution by WP-F2). Treated as immutable.
        rate_div: decimation factor (>= 1).
    """

    kind: str
    name: str
    params: Mapping[str, Any] = field(default_factory=dict)
    rate_div: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.kind, str) or not KIND_RE.match(self.kind):
            raise ValueError(f"invalid block kind {self.kind!r}")
        if not isinstance(self.name, str) or not self.name:
            raise ValueError(f"block name must be a non-empty str, got {self.name!r}")
        if isinstance(self.rate_div, bool) or int(self.rate_div) != self.rate_div or self.rate_div < 1:
            raise ValueError(f"rate_div must be an int >= 1, got {self.rate_div!r}")

    def to_obj(self) -> dict[str, Any]:
        obj: dict[str, Any] = {"block": self.kind, "name": self.name, "params": dict(self.params)}
        if self.rate_div != 1:
            obj["rate_div"] = int(self.rate_div)
        return obj

    @classmethod
    def from_obj(cls, obj: Mapping[str, Any]) -> BlockSpec:
        unknown = set(obj) - {"block", "name", "params", "rate_div"}
        if unknown:
            raise ValueError(f"unknown block keys {sorted(unknown)} in {dict(obj)!r}")
        if "block" not in obj:
            raise ValueError(f"block spec needs a 'block' key: {dict(obj)!r}")
        kind = obj["block"]
        return cls(kind, obj.get("name", kind), dict(obj.get("params") or {}), obj.get("rate_div", 1))


@runtime_checkable
class Block(Protocol):
    """Structural type of every pipeline element, leaf or composite (see module docstring)."""

    kind: str
    name: str
    reads: tuple[str, ...]
    writes: tuple[str, ...]
    rate_div: int
    runtime_params: tuple[str, ...]
    kernel: KernelFn | None

    @property
    def spec(self) -> Any: ...  # BlockSpec for leaves; SeqSpec / ParallelSpec / StackSpec for composites

    @property
    def toggles(self) -> Mapping[str, bool]: ...

    def __call__(self, state: State, params: Params) -> State: ...

    def init(self, layout: SlotLayout, rng: np.random.Generator) -> None: ...

    def bind(self, layout: SlotLayout) -> None: ...

    def losses(self) -> list[LossElement]: ...

    def bands(self) -> list[BandDecl]: ...


class BaseBlock(ABC):
    """Convenience base for leaf blocks implementing the bookkeeping of the :class:`Block` contract.

    Subclasses implement ``__call__`` (reference step) and usually extend ``init`` (declare written
    slots with ``layout.add``; call ``super().init`` to register runtime params) and ``bind``
    (resolve offsets, set ``self.kernel``), and override ``losses`` / ``bands`` when applicable.

    Args:
        kind: registry key.
        name: instance name (None → ``kind``).
        params: the spec params (kept verbatim for ``spec``).
        reads / writes: State slot names (wiring).
        rate_div: decimation factor.
        runtime_params: keys of ``params`` that are runtime (P) parameters; each must be present.
        toggle_names: names of this block's toggles; ``params["toggles"]`` must give all of them.
    """

    #: subclasses may override; also accepted as constructor args
    toggle_names: tuple[str, ...] = ()

    def __init__(
        self,
        kind: str,
        name: str | None,
        params: Mapping[str, Any],
        *,
        reads: Iterable[str] = (),
        writes: Iterable[str] = (),
        rate_div: int = 1,
        runtime_params: Iterable[str] = (),
        toggle_names: Iterable[str] | None = None,
    ):
        self._spec = BlockSpec(kind, name or kind, dict(params), rate_div)
        self.kind = self._spec.kind
        self.name = self._spec.name
        self.params: Mapping[str, Any] = self._spec.params
        self.reads = tuple(reads)
        self.writes = tuple(writes)
        self.rate_div = int(rate_div)
        self.runtime_params = tuple(runtime_params)
        require(self.params, *self.runtime_params, block=self.name)
        if toggle_names is not None:
            self.toggle_names = tuple(toggle_names)
        self._toggles = self._resolve_toggles(self.params.get("toggles", {}))
        self.kernel: KernelFn | None = None
        self.rng: np.random.Generator | None = None
        #: expected unit tags of read / written slots (check 4 of ``check_pipeline``); optional
        self.slot_units: dict[str, str] = {}

    def _resolve_toggles(self, given: Mapping[str, Any]) -> dict[str, bool]:
        given = dict(given or {})
        unknown = set(given) - set(self.toggle_names)
        if unknown:
            raise ValueError(f"block {self.name!r}: unknown toggle(s) {sorted(unknown)}; known {list(self.toggle_names)}")
        missing = [t for t in self.toggle_names if t not in given]
        if missing:
            raise MissingParameterError(f"block {self.name!r}: toggle(s) {missing} not set in params['toggles'] (D34)")
        for k, v in given.items():
            if not isinstance(v, (bool, np.bool_)):
                raise TypeError(f"block {self.name!r}: toggle {k!r} must be bool, got {v!r}")
        return {k: bool(given[k]) for k in self.toggle_names}

    @property
    def spec(self) -> BlockSpec:
        return self._spec

    @property
    def toggles(self) -> Mapping[str, bool]:
        return dict(self._toggles)

    def set_rate_div(self, rate_div: int) -> None:
        """Set the decimation factor (used by ``make_block`` from ``BlockSpec.rate_div``); keeps ``spec`` in sync."""
        self._spec = BlockSpec(self.kind, self.name, self.params, rate_div)
        self.rate_div = self._spec.rate_div

    def param_name(self, key: str) -> str:
        """Runtime-parameter name in P: ``<block name>.<key>``."""
        return f"{self.name}.{key}"

    @abstractmethod
    def __call__(self, state: State, params: Params) -> State:
        """Reference step: update ``state`` in place and return the same object."""

    def init(self, layout: SlotLayout, rng: np.random.Generator) -> None:
        """Register runtime params (spec values) and keep ``rng``. Subclasses add their slots."""
        self.rng = rng
        for key in self.runtime_params:
            value = np.asarray(self.params[key], dtype=np.float64)
            layout.add_param(self.param_name(key), value, value.shape)

    def bind(self, layout: SlotLayout) -> None:
        """Resolve offsets and set ``self.kernel``. Default: Python-only block (kernel None)."""
        self.kernel = None

    def prepare(self, state: State, params: Params, n0: int, n_samples: int) -> None:
        """Python hook before each piece of ``n_samples`` base samples starting at ``n0``. Default: no-op."""

    def losses(self) -> list[LossElement]:
        return []

    def bands(self) -> list[BandDecl]:
        return []

    def __rshift__(self, other: Block) -> Block:
        from cranberry.core.pipeline import Seq

        return Seq.chain(self, other)

    def __or__(self, other: Block) -> Block:
        from cranberry.core.pipeline import Parallel

        return Parallel.join(self, other)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} kind={self.kind!r} name={self.name!r}>"
