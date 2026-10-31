"""Nested pipeline composition: ``Seq``, ``Parallel``, ``Stack`` (protected list; PLAN §2.5, G09).

Semantics (G09):

* ``Seq([a, b, c])`` (also ``a >> b >> c``): run in order within one sample; ``b`` sees what ``a``
  wrote in the same sample.
* ``Parallel([p, q])`` (also ``p | q``): branches reading the same inputs. Their write sets must be
  disjoint (checked at build time by :func:`check_pipeline`); they run sequentially in declaration
  order, so the result does not depend on the order.
* ``Stack(stages=[c1, c2, c3], actuators=[a1, a2, a3])``: stacked cascade (Ref C, T16 §4.4): stage
  c_k drives actuator a_k and also feeds stage c_{k+1}. Pure sugar, expanded to
  ``Seq([c1, Parallel([a1, Seq([c2, Parallel([a2, Seq([c3, a3])])])])])``
  = ``c1 >> (a1 | (c2 >> (a2 | (c3 >> a3))))``.
* Wiring is by slot name (block params ``in`` / ``out``), not by position.
* Feedback is allowed only through a ``delay_buffer`` block (``core.delay``, kind
  :data:`DELAY_BUFFER_KIND`); a zero-delay cycle within one sample is an algebraic loop and is
  rejected by :func:`check_pipeline`.
* Composites are Blocks too: ``reads`` / ``writes`` are unions; ``losses()`` / ``bands()`` are
  concatenations; nesting depth is unlimited; a composite may live in its own YAML file.

Build-time checks (:func:`check_pipeline`, WP-F1): (1) every read slot has a writer (or is a
reserved slot ``t`` / ``n``), (2) no same-sample cycle except through ``delay_buffer``,
(3) parallel branches write disjoint slots, (4) unit tags of wired slots agree.

Pipeline YAML / dict form (round trip via :func:`spec_to_obj` / :func:`spec_from_obj`, VT-F1)::

    [ ... ]                                     # SeqSpec without a name (the top-level pipeline)
    {seq: [ ... ], name: <str>}                 # named SeqSpec
    {parallel: [ ... ], name: <str>?}           # ParallelSpec
    {stack: {stages: [ ... ], actuators: [ ... ]}, name: <str>?}   # StackSpec
    {block: <kind>, name: <str>, params: {...}, rate_div: <int>?}  # BlockSpec (core.block)

Composites hold no slots, runtime params or RNG of their own: the builder initialises leaves
individually (``iter_leaves``) with ``core.rng.stream(root, leaf.name)`` (G11).
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Union

from cranberry.core.block import Block, BlockSpec
from cranberry.core.rates import should_step

if TYPE_CHECKING:
    import numpy as np

    from cranberry.core.layout import SlotLayout
    from cranberry.core.state import Params, State
    from cranberry.offline.band_ledger import BandDecl
    from cranberry.offline.loss_registry import LossElement

__all__ = [
    "DELAY_BUFFER_KIND",
    "SeqSpec",
    "ParallelSpec",
    "StackSpec",
    "PipelineSpec",
    "Composite",
    "Seq",
    "Parallel",
    "Stack",
    "spec_from_obj",
    "spec_to_obj",
    "from_spec",
    "iter_leaves",
    "build",
    "check_pipeline",
]

#: Block kind of the delay buffer, the only legal feedback path (G09).
DELAY_BUFFER_KIND = "delay_buffer"


@dataclass(frozen=True)
class SeqSpec:
    children: tuple[PipelineSpec, ...]
    name: str | None = None


@dataclass(frozen=True)
class ParallelSpec:
    children: tuple[PipelineSpec, ...]
    name: str | None = None


@dataclass(frozen=True)
class StackSpec:
    stages: tuple[PipelineSpec, ...]
    actuators: tuple[PipelineSpec, ...]
    name: str | None = None


PipelineSpec = Union[BlockSpec, SeqSpec, ParallelSpec, StackSpec]


def spec_from_obj(obj: Any) -> PipelineSpec:
    """Parse the plain dict/list pipeline form (module docstring) into specs. Tags are WP-F2's job."""
    if isinstance(obj, (list, tuple)):
        return SeqSpec(tuple(spec_from_obj(o) for o in obj))
    if not isinstance(obj, Mapping):
        raise TypeError(f"pipeline element must be a list or mapping, got {type(obj).__name__}: {obj!r}")
    keys = set(obj)
    if "block" in keys:
        return BlockSpec.from_obj(obj)
    for key, cls in (("seq", SeqSpec), ("parallel", ParallelSpec)):
        if key in keys:
            if keys - {key, "name"}:
                raise ValueError(f"unknown keys {sorted(keys - {key, 'name'})} in {key} element")
            return cls(tuple(spec_from_obj(o) for o in obj[key]), obj.get("name"))
    if "stack" in keys:
        if keys - {"stack", "name"}:
            raise ValueError(f"unknown keys {sorted(keys - {'stack', 'name'})} in stack element")
        body = obj["stack"]
        if not isinstance(body, Mapping) or set(body) != {"stages", "actuators"}:
            raise ValueError("stack element must be {stack: {stages: [...], actuators: [...]}}")
        return StackSpec(
            tuple(spec_from_obj(o) for o in body["stages"]),
            tuple(spec_from_obj(o) for o in body["actuators"]),
            obj.get("name"),
        )
    raise ValueError(f"cannot parse pipeline element {dict(obj)!r}: expected block/seq/parallel/stack")


def spec_to_obj(spec: PipelineSpec) -> Any:
    """Inverse of :func:`spec_from_obj` (unnamed SeqSpec → plain list)."""
    if isinstance(spec, BlockSpec):
        return spec.to_obj()
    if isinstance(spec, SeqSpec):
        items = [spec_to_obj(c) for c in spec.children]
        return items if spec.name is None else {"seq": items, "name": spec.name}
    if isinstance(spec, ParallelSpec):
        obj: dict[str, Any] = {"parallel": [spec_to_obj(c) for c in spec.children]}
    elif isinstance(spec, StackSpec):
        obj = {
            "stack": {
                "stages": [spec_to_obj(c) for c in spec.stages],
                "actuators": [spec_to_obj(c) for c in spec.actuators],
            }
        }
    else:
        raise TypeError(f"not a pipeline spec: {spec!r}")
    if spec.name is not None:
        obj["name"] = spec.name
    return obj


def _union(seqs: Sequence[Sequence[str]]) -> tuple[str, ...]:
    out: dict[str, None] = {}
    for s in seqs:
        out.update(dict.fromkeys(s))
    return tuple(out)


class Composite:
    """Common base of Seq / Parallel / Stack (a Block whose step runs its children)."""

    kind: str = ""

    def __init__(self, children: Sequence[Block], name: str | None = None):
        if not children:
            raise ValueError(f"{type(self).__name__} needs at least one child")
        for c in children:
            if not isinstance(c, Block):
                raise TypeError(f"{type(self).__name__} child is not a Block: {c!r}")
        self.children: tuple[Block, ...] = tuple(children)
        self._name = name
        self.rate_div = 1
        self.runtime_params: tuple[str, ...] = ()
        self.kernel = None  # composed by core.compile (WP-F1)

    @property
    def name(self) -> str:
        return self._name if self._name is not None else self.kind

    @property
    def reads(self) -> tuple[str, ...]:
        return _union([c.reads for c in self.children])

    @property
    def writes(self) -> tuple[str, ...]:
        return _union([c.writes for c in self.children])

    @property
    def toggles(self) -> Mapping[str, bool]:
        out: dict[str, bool] = {}
        for leaf in iter_leaves(self):
            out.update({f"{leaf.name}.{k}": v for k, v in leaf.toggles.items()})
        return out

    def __call__(self, state: State, params: Params) -> State:
        n = state.n
        for c in self.children:
            if should_step(n, c.rate_div):
                c(state, params)
        return state

    def init(self, layout: SlotLayout, rng: np.random.Generator) -> None:
        """No-op: composites own nothing; the builder initialises leaves (see module docstring)."""

    def bind(self, layout: SlotLayout) -> None:
        for c in self.children:
            c.bind(layout)

    def losses(self) -> list[LossElement]:
        return [e for c in self.children for e in c.losses()]

    def bands(self) -> list[BandDecl]:
        return [b for c in self.children for b in c.bands()]

    def __rshift__(self, other: Block) -> Block:
        return Seq.chain(self, other)

    def __or__(self, other: Block) -> Block:
        return Parallel.join(self, other)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} name={self._name!r} children={[c.name for c in self.children]}>"


class Seq(Composite):
    """Sequential composition (``a >> b``)."""

    kind = "seq"

    @classmethod
    def chain(cls, a: Block, b: Block) -> Seq:
        """``a >> b``; an unnamed Seq on either side is flattened."""
        left = list(a.children) if isinstance(a, Seq) and a._name is None else [a]
        right = list(b.children) if isinstance(b, Seq) and b._name is None else [b]
        return cls(left + right)

    @property
    def spec(self) -> SeqSpec:
        return SeqSpec(tuple(c.spec for c in self.children), self._name)


class Parallel(Composite):
    """Parallel branches with disjoint writes (``p | q``), executed in declaration order."""

    kind = "parallel"

    @classmethod
    def join(cls, a: Block, b: Block) -> Parallel:
        """``a | b``; an unnamed Parallel on either side is flattened."""
        left = list(a.children) if isinstance(a, Parallel) and a._name is None else [a]
        right = list(b.children) if isinstance(b, Parallel) and b._name is None else [b]
        return cls(left + right)

    @property
    def spec(self) -> ParallelSpec:
        return ParallelSpec(tuple(c.spec for c in self.children), self._name)


class Stack(Composite):
    """Stacked cascade (Ref C): sugar for nested Seq/Parallel (module docstring).

    Args:
        stages: controller stages c_1..c_n (each a Block, typically a Seq like ``adc >> pid``).
        actuators: actuator chains a_1..a_n (typically ``delay_buffer >> actuator``).
        name: optional name.
    """

    kind = "stack"

    def __init__(self, stages: Sequence[Block], actuators: Sequence[Block], name: str | None = None):
        if len(stages) != len(actuators) or not stages:
            raise ValueError("Stack needs equally many (>= 1) stages and actuators")
        self.stages: tuple[Block, ...] = tuple(stages)
        self.actuators: tuple[Block, ...] = tuple(actuators)
        expanded = self.expand()
        super().__init__([expanded], name)

    def expand(self) -> Seq:
        """``c1 >> (a1 | (c2 >> (a2 | ... (c_n >> a_n))))`` as explicit Seq / Parallel objects."""
        inner: Block = Seq([self.stages[-1], self.actuators[-1]])
        for c, a in zip(reversed(self.stages[:-1]), reversed(self.actuators[:-1])):
            inner = Seq([c, Parallel([a, inner])])
        return inner  # type: ignore[return-value]

    @property
    def spec(self) -> StackSpec:
        return StackSpec(
            tuple(s.spec for s in self.stages), tuple(a.spec for a in self.actuators), self._name
        )


def from_spec(spec: PipelineSpec) -> Block:
    """Construct the block tree for ``spec`` via the registry (no layout, no RNG yet)."""
    from cranberry.core.registry import make_block

    if isinstance(spec, BlockSpec):
        return make_block(spec)
    if isinstance(spec, SeqSpec):
        return Seq([from_spec(c) for c in spec.children], spec.name)
    if isinstance(spec, ParallelSpec):
        return Parallel([from_spec(c) for c in spec.children], spec.name)
    if isinstance(spec, StackSpec):
        return Stack([from_spec(c) for c in spec.stages], [from_spec(c) for c in spec.actuators], spec.name)
    raise TypeError(f"not a pipeline spec: {spec!r}")


def iter_leaves(block: Block) -> Iterator[Block]:
    """Depth-first leaves (non-composites) in execution order."""
    if isinstance(block, Composite):
        for c in block.children:
            yield from iter_leaves(c)
    else:
        yield block


def build(pipeline_spec: PipelineSpec, layout: SlotLayout, root_seed: int) -> Block:
    """Build a runnable pipeline: ``from_spec`` → ``leaf.init(layout, stream(root_seed, leaf.name))``
    for every leaf → :func:`check_pipeline` → ``layout.freeze()`` → ``bind``. Implemented by WP-F1."""
    raise NotImplementedError("core.pipeline.build is implemented by WP-F1")


def check_pipeline(block: Block, layout: SlotLayout) -> None:
    """The four build-time checks of the module docstring; raises on violation. Implemented by WP-F1."""
    raise NotImplementedError("core.pipeline.check_pipeline is implemented by WP-F1")
