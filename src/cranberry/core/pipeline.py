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

Build-time checks (:func:`check_pipeline`, WP-F1; all violations are collected and raised together
as :class:`PipelineCheckError`):

0. wiring: leaf names are unique; every slot a leaf reads or writes is registered in the layout.
1. every read slot has a writer in the pipeline (or is a reserved slot ``t`` / ``n``, or is listed
   in ``external``: written from outside the pipeline, e.g. by a test harness).
2. feedback only through ``delay_buffer``. Execution order is the depth-first leaf order. A read
   is *forward* when a writer of the slot runs earlier in the same sample, *own state* when the
   reader is the only writer, and a *back edge* otherwise (it sees the previous sample's value;
   the loop closure contributes one sample of delay, see ``core.delay``). With edges
   writer → reader, every cycle must pass through a ``delay_buffer`` leaf (no zero-delay /
   algebraic loop), and every back edge must close a cycle; one that does not is a
   *read-before-write* ordering error (an unmodelled one-sample delay). Every back edge must also
   close a cycle on which it is the **only** back edge (fewest back edges on any closing cycle,
   0-1 BFS, = 1): a loop whose every closing path needs two or more back edges (e.g. a
   ``delay_buffer`` placed before the controller it delays: ``plant → sensor → d → pid → act``)
   would add k ≥ 2 samples at the closure, total loop delay ``tau_s + k/f_sim`` instead of the
   ``tau_s + 1/f_sim`` of D55-1, and is rejected naming k and the extra delay (Phase 5 F1-1).
   Cycles that merely traverse two loops which each close with one back edge (coupled locks) are
   allowed: each of their back edges also closes a one-back-edge loop of its own.
3. ``Parallel`` branches are independent: their write sets are disjoint and no branch reads a slot
   written by a sibling branch (so the result does not depend on the declaration order, G09).
4. unit tags: for every slot in a leaf's optional ``slot_units`` mapping, the expected tag equals
   the tag registered in the layout.

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

from collections import deque
from collections.abc import Iterable, Iterator, Mapping, Sequence
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
    "PipelineCheckError",
    "apply_toggle_overrides",
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
        """Bind every child and keep ``layout`` (``core.compile.compile_loop`` reads it from the root)."""
        self.layout = layout
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


class PipelineCheckError(ValueError):
    """One or more build-time checks failed; ``violations`` lists every message."""

    def __init__(self, violations: Sequence[str]):
        self.violations = list(violations)
        super().__init__("pipeline check failed:\n  - " + "\n  - ".join(self.violations))


def _walk(block: Block) -> Iterator[Block]:
    yield block
    if isinstance(block, Composite):
        for c in block.children:
            yield from _walk(c)


def _leaf_writes(block: Block) -> set[str]:
    return {w for leaf in iter_leaves(block) for w in leaf.writes}


def _leaf_reads(block: Block) -> set[str]:
    return {r for leaf in iter_leaves(block) for r in leaf.reads}


def check_pipeline(block: Block, layout: SlotLayout, *, external: Iterable[str] = ()) -> None:
    """Run the build-time checks of the module docstring on a pipeline whose leaves are initialised
    (all slots registered in ``layout``); raise :class:`PipelineCheckError` listing every violation.

    Args:
        block: root of the pipeline.
        layout: the layout the leaves were initialised into (frozen or not).
        external: slot names written from outside the pipeline (exempt from check 1).
    """
    from cranberry.core.layout import RESERVED_SLOTS

    leaves = list(iter_leaves(block))
    exempt = set(external) | set(RESERVED_SLOTS)
    errs: list[str] = []

    # (0) wiring: unique names, registered slots
    counts: dict[str, int] = {}
    for leaf in leaves:
        counts[leaf.name] = counts.get(leaf.name, 0) + 1
    errs += [f"duplicate leaf name {n!r} ({c}x): names key RNG streams and runtime params"
             for n, c in counts.items() if c > 1]
    for leaf in leaves:
        for role, names in (("reads", leaf.reads), ("writes", leaf.writes)):
            errs += [f"{leaf.name!r} {role} unregistered slot {s!r}" for s in names if s not in layout]

    # (3) parallel branches are independent
    for node in _walk(block):
        if not isinstance(node, Parallel):
            continue
        wsets = [_leaf_writes(c) for c in node.children]
        rsets = [_leaf_reads(c) for c in node.children]
        for i, ci in enumerate(node.children):
            for j, cj in enumerate(node.children):
                if i == j:
                    continue
                if j > i and wsets[i] & wsets[j]:
                    errs.append(f"parallel {node.name!r}: branches {ci.name!r} and {cj.name!r} both write "
                                f"{sorted(wsets[i] & wsets[j])}")
                if rsets[i] & wsets[j]:
                    errs.append(f"parallel {node.name!r}: branch {ci.name!r} reads {sorted(rsets[i] & wsets[j])} "
                                f"written by sibling branch {cj.name!r} (result would depend on branch order)")

    # (4) unit tags
    for leaf in leaves:
        for slot, unit in (getattr(leaf, "slot_units", None) or {}).items():
            if slot in layout and layout.slot(slot).unit != unit:
                errs.append(f"{leaf.name!r} expects unit {unit!r} for slot {slot!r}, "
                            f"layout has {layout.slot(slot).unit!r}")

    # (1) every read has a writer; (2) dataflow graph, edges writer -> reader
    writers: dict[str, list[int]] = {}
    for i, leaf in enumerate(leaves):
        for w in dict.fromkeys(leaf.writes):
            writers.setdefault(w, []).append(i)
    edges: dict[int, set[int]] = {i: set() for i in range(len(leaves))}
    back: list[tuple[int, int, str]] = []
    back_pairs: set[tuple[int, int]] = set()
    for i, leaf in enumerate(leaves):
        for r in dict.fromkeys(leaf.reads):
            ws = writers.get(r, [])
            others = [j for j in ws if j != i]
            if not others:
                if not ws and r not in exempt:
                    errs.append(f"{leaf.name!r} reads {r!r}, which no block writes (check 1)")
                continue  # own state, external or reserved
            prior = [j for j in others if j < i]
            if prior:
                edges[max(prior)].add(i)
            else:
                j = max(ws)
                edges[j].add(i)
                back.append((j, i, r))
                back_pairs.add((j, i))

    def reachable(src: int, dst: int, skip: set[int]) -> bool:
        stack, seen = [src], {src}
        while stack:
            u = stack.pop()
            if u == dst:
                return True
            for v in edges[u]:
                if v not in seen and v not in skip:
                    seen.add(v)
                    stack.append(v)
        return False

    def fewest_back_edges(src: int, dst: int) -> int | None:
        """Fewest back edges on a path src → dst (0-1 BFS over the dataflow graph), None if unreachable."""
        dist, dq = {src: 0}, deque([src])
        while dq:
            u = dq.popleft()
            for v in edges[u]:
                w = int((u, v) in back_pairs)
                if v not in dist or dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    (dq.append if w else dq.appendleft)(v)
        return dist.get(dst)

    delay_nodes = {i for i, leaf in enumerate(leaves) if leaf.kind == DELAY_BUFFER_KIND}
    f_sim = getattr(layout, "f_sim_hz", None)
    for j, i, r in back:
        if not reachable(i, j, set()):
            errs.append(f"read-before-write: {leaves[i].name!r} reads {r!r} before its writer {leaves[j].name!r} "
                        "runs (it would see the previous sample); reorder the pipeline or route feedback "
                        "through a delay_buffer")
        elif j not in delay_nodes and i not in delay_nodes and reachable(i, j, delay_nodes):
            errs.append(f"zero-delay feedback loop: {leaves[i].name!r} reads {r!r} from {leaves[j].name!r} on a "
                        "cycle without a delay_buffer (algebraic loop, G09)")
        else:
            k = (fewest_back_edges(i, j) or 0) + 1  # + the closing edge j → i itself
            if k > 1:
                extra = f" = {(k - 1) / f_sim:.4g} s" if f_sim else ""
                errs.append(f"feedback loop through {leaves[j].name!r} → {leaves[i].name!r} (via {r!r}) closes with {k} "
                            f"back edges: total loop delay tau_s + {k}/f_sim, {k - 1} sample(s){extra} more than the "
                            "tau_s + 1/f_sim of D55-1; place the delay_buffer after the controller it delays (just "
                            "before the actuator) so that the loop closes exactly once (Phase 5 F1-1)")
    if errs:
        raise PipelineCheckError(errs)


def build(
    pipeline_spec: PipelineSpec | Block,
    layout: SlotLayout,
    root_seed: Any,
    *,
    external: Iterable[str] = (),
) -> Block:
    """Build a runnable pipeline.

    Steps: ``from_spec`` (skipped when a never-built Block tree is passed) →
    ``leaf.init(layout, core.rng.stream(root_seed, leaf.name))`` for every leaf in execution order →
    :func:`check_pipeline` → ``layout.freeze()`` → ``bind``; the root keeps ``layout`` as
    ``root.layout``.

    Args:
        pipeline_spec: pipeline spec tree, or a Block tree made by factories that was never built.
        layout: an unfrozen ``SlotLayout`` (may already hold externally written slots).
        root_seed: run seed (int >= 0 or SeedSequence, G11).
        external: slots written from outside the pipeline (see :func:`check_pipeline`).

    Returns:
        The root Block; allocate ``State(layout)`` / ``Params(layout)`` to run it.
    """
    from cranberry.core.rng import stream

    block = pipeline_spec if isinstance(pipeline_spec, Block) else from_spec(pipeline_spec)
    names = [leaf.name for leaf in iter_leaves(block)]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        raise PipelineCheckError([f"duplicate leaf name {n!r}: names key RNG streams and runtime params"
                                  for n in dups])
    for leaf in iter_leaves(block):
        leaf.init(layout, stream(root_seed, leaf.name))
    check_pipeline(block, layout, external=external)
    layout.freeze()
    block.bind(layout)
    if not isinstance(block, Composite):
        block.layout = layout  # type: ignore[attr-defined]
    return block


def apply_toggle_overrides(spec: PipelineSpec, toggles: Mapping[str, bool]) -> PipelineSpec:
    """Return ``spec`` with ``RunSpec.toggles`` entries ``"<leaf name>.<toggle>": bool`` written into
    the matching leaf's ``params["toggles"]`` (used by ``run``; ``ablate`` builds on it). Unknown
    leaves raise KeyError; unknown toggle names are rejected by the block factory."""
    import dataclasses

    leaf_names: set[str] = set()

    def collect(s: PipelineSpec) -> None:
        if isinstance(s, BlockSpec):
            leaf_names.add(s.name)
        elif isinstance(s, (SeqSpec, ParallelSpec)):
            for c in s.children:
                collect(c)
        elif isinstance(s, StackSpec):
            for c in (*s.stages, *s.actuators):
                collect(c)

    collect(spec)
    by_leaf: dict[str, dict[str, bool]] = {}
    for key, val in toggles.items():
        leaf, _, tog = key.rpartition(".")
        if leaf not in leaf_names:
            raise KeyError(f"toggle override {key!r}: no leaf named {leaf!r}")
        by_leaf.setdefault(leaf, {})[tog] = bool(val)

    def walk(s: PipelineSpec) -> PipelineSpec:
        if isinstance(s, BlockSpec):
            if s.name not in by_leaf:
                return s
            params = dict(s.params)
            params["toggles"] = {**dict(params.get("toggles") or {}), **by_leaf[s.name]}
            return dataclasses.replace(s, params=params)
        if isinstance(s, (SeqSpec, ParallelSpec)):
            return dataclasses.replace(s, children=tuple(walk(c) for c in s.children))
        if isinstance(s, StackSpec):
            return dataclasses.replace(s, stages=tuple(walk(c) for c in s.stages),
                                       actuators=tuple(walk(c) for c in s.actuators))
        raise TypeError(f"not a pipeline spec: {s!r}")

    return walk(spec) if by_leaf else spec
