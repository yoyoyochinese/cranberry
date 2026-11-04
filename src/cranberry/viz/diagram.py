"""Configuration diagrams (D84-3, D85-13): Graphviz DOT for one built configuration, rendered with ``dot``.

Two diagrams per configuration:

* **Signal flow** (:func:`signal_flow_dot`): every pipeline block is a node, every slot an edge from the
  block that writes it to the blocks that read it (``t`` is left out). Nodes are grouped by stage —
  source / optics / detection / control / offline — plus two clusters that are collapsed into one node
  each unless expanded: ``noise`` (environment noise sources, ``blocks.noise`` except the laser) and
  ``bookkeeping`` (metric probes, the true lock error ``lock_error``, CHK blocks). A read that sees the
  previous sample (the writer runs later in the pipeline: a feedback edge, D55-1) is drawn dashed red.
  Blocks of one ``stack`` (nested controller, PLAN §2.5) form a sub-cluster of the control stage.
  Offline blocks (``offline.blocks``) have no slot edges and are listed in their own cluster.
* **Optical layout** (:func:`optical_layout_dot`): ``topology.beams[].path`` drawn as chains source →
  segments → observer, one edge per consecutive pair with one colour per beam; a segment entered from
  several places is a combiner (double border, "⊕"), one left towards several is a splitter ("⑂");
  reverse traversal (``direction: -1``) is drawn dashed. Locks (``topology.locks``) are notes linked
  from their observer, and to the actuators they drive (traced through the signal flow by
  :func:`lock_drives`), which point at the segment driver they write (``drv.<segment>.<driver>``).

Stage, family and lock tracing use only the block specs, the factories' modules and the blocks'
``reads`` / ``writes``, so a new block family lands in a stage by its package
(:data:`FAMILY_STAGE`); blocks outside ``cranberry.blocks`` fall back to a slot heuristic.

Rendering (:func:`render`) needs the Graphviz ``dot`` executable (Homebrew ``graphviz``; looked up on
``PATH``, then ``/opt/homebrew/bin/dot``). CLI: ``python -m cranberry diagram <preset> [--overlay ...]
--out <dir>``.
"""

from __future__ import annotations

import dataclasses
import shutil
import subprocess
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.core.block import BlockSpec

__all__ = [
    "STAGES",
    "GROUPS",
    "FAMILY_STAGE",
    "Leaf",
    "pipeline_leaves",
    "lock_drives",
    "signal_flow_dot",
    "optical_layout_dot",
    "render",
    "write_diagrams",
    "find_dot",
]

#: Stage clusters of the signal-flow diagram, in drawing order.
STAGES: tuple[str, ...] = ("source", "optics", "detection", "control", "offline")
#: Collapsible clusters (one node each unless named in ``expand``).
GROUPS: tuple[str, ...] = ("noise", "bookkeeping")

#: Block family (``cranberry.blocks.<family>`` or framework module) → stage.
FAMILY_STAGE: dict[str, str] = {
    "source": "source",
    "path": "optics", "fib": "optics", "pol": "optics", "spf": "optics", "los": "optics", "loss": "optics",
    "dmx": "optics", "fsr": "optics",
    "err": "detection", "det": "detection",
    "ctrl": "control", "act": "control", "mod": "control", "delay": "control",
    "meas": "detection", "chk": "detection", "noise": "optics",
}
_DETECTOR_FAMILIES = frozenset({"err", "det"})
#: Families traversed when tracing an actuator back to its detector (controller / actuator chain).
_CHAIN_FAMILIES = frozenset({"ctrl", "act", "mod", "delay"})

_STAGE_STYLE = {
    "source": ("#fde2c4", "#c46a00"), "optics": ("#d6e8fb", "#1f5fa8"), "detection": ("#dff3dc", "#2e7d32"),
    "control": ("#efe0f7", "#7b3fa0"), "offline": ("#eeeeee", "#555555"),
    "noise": ("#f6f1d5", "#8a7a1d"), "bookkeeping": ("#f2f2f2", "#777777"),
}
_BEAM_COLORS = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b", "#e377c2",
                "#7f7f7f", "#bcbd22")
_MAX_EDGE_SLOTS = 3


@dataclass(frozen=True)
class Leaf:
    """One pipeline block as the diagrams see it.

    Attributes:
        index: execution order (0-based).
        path: ``iter_block_specs`` path (``pipeline[22].stack.stages[0].seq[1]``).
        spec: the block spec (kind, name, params).
        family: ``blocks.<family>`` of the factory (``path`` / ``delay`` for framework blocks, ``?``
            for blocks outside the package).
        module: last component of the factory's module (``laser``, ``locks``, ``supervisor``, ...).
        reads / writes: slot names of the built block.
        stage: one of :data:`STAGES` (never ``offline`` for a pipeline block).
        group: ``"noise"`` / ``"bookkeeping"`` / None.
        stack: path of the enclosing top-level ``stack`` (``pipeline[22]``) or None.
    """

    index: int
    path: str
    spec: BlockSpec
    family: str
    module: str
    reads: tuple[str, ...]
    writes: tuple[str, ...]
    stage: str
    group: str | None
    stack: str | None

    @property
    def name(self) -> str:
        return self.spec.name


def _family(kind: str) -> tuple[str, str]:
    from cranberry.core.registry import get_factory

    mod = getattr(get_factory(kind), "__module__", "") or ""
    parts = mod.split(".")
    if len(parts) >= 3 and parts[:2] == ["cranberry", "blocks"]:
        return parts[2], parts[-1]
    if mod == "cranberry.core.delay":
        return "delay", "delay"
    if len(parts) >= 3 and parts[:2] == ["cranberry", "path"]:
        return "path", parts[-1]
    return "?", parts[-1] if parts else "?"


def _classify(family: str, module: str, kind: str, reads: Sequence[str], writes: Sequence[str]) -> tuple[str, str | None]:
    if family == "noise":
        return ("source", None) if module == "laser" else ("optics", "noise")
    if family in ("meas", "chk") or (family == "path" and module == "locks"):
        return "detection", "bookkeeping"
    if kind == "probe_pulser":
        return "optics", None
    if family in FAMILY_STAGE:
        return FAMILY_STAGE[family], None
    # outside cranberry.blocks (test / user blocks): judge by the slots
    if any(s.startswith("beam.") for s in reads):
        return "detection", None
    if any(s.startswith("drv.") for s in writes) and any(s != "t" for s in reads):
        return "control", None
    return "optics", None


def _stack_of(path: str) -> str | None:
    i = path.find(".stack.")
    return path[:i] if i >= 0 else None


def _resolved(spec: RunSpec) -> RunSpec:
    """The spec itself, or its first sweep point when sweep markers are left (diagrams are structural)."""
    from cranberry import config

    if config.sweep_markers(spec):
        from cranberry.runner.sweep import expand

        rs = spec if spec.run.seed_mode else dataclasses.replace(spec, run=dataclasses.replace(spec.run, seed_mode="common"))
        return expand(rs, limit=1)[0].spec
    return spec


def pipeline_leaves(spec: RunSpec) -> list[Leaf]:
    """Every pipeline block of ``spec`` (topology injected, toggles applied, as ``runner.run`` builds them)."""
    from cranberry.config.spec import apply_toggles, iter_block_specs
    from cranberry.core.registry import ensure_loaded, make_block
    from cranberry.runner.run import inject_topology, offline_block_specs

    ensure_loaded()
    spec = _resolved(spec)
    off = {b.name for b in offline_block_specs(spec)}
    toggles = {k: v for k, v in spec.toggles.items() if k.rpartition(".")[0] not in off}
    pipe = apply_toggles(inject_topology(spec.pipeline, spec.topology), toggles)
    out = []
    for i, (path, b) in enumerate(iter_block_specs(pipe)):
        blk = make_block(b)
        fam, mod = _family(b.kind)
        reads, writes = tuple(blk.reads), tuple(blk.writes)
        stage, group = _classify(fam, mod, b.kind, reads, writes)
        stack = _stack_of(path)
        if stack is not None and group is None:
            stage = "control"
        out.append(Leaf(i, path, b, fam, mod, reads, writes, stage, group, stack))
    return out


def _slot_edges(leaves: Sequence[Leaf]) -> list[tuple[int, int, str, bool]]:
    """(writer index, reader index, slot, back edge) for every read slot with a writer (``t`` skipped).
    The writer is the last one before the reader, else the last one overall (previous sample)."""
    writers: dict[str, list[int]] = defaultdict(list)
    for lf in leaves:
        for s in lf.writes:
            writers[s].append(lf.index)
    out = []
    for lf in leaves:
        for s in lf.reads:
            ws = writers.get(s)
            if s == "t" or not ws:
                continue
            before = [w for w in ws if w < lf.index]
            w = before[-1] if before else ws[-1]
            if w != lf.index:
                out.append((w, lf.index, s, w >= lf.index))
    return out


def _q(s: str) -> str:
    return '"' + str(s).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _lbl(*lines: str) -> str:
    return _q("\n".join(lines)).replace("\n", "\\n")


def _setpoint(lock: Mapping[str, Any]) -> float | None:
    v = lock.get("setpoint_rad", lock.get("setpoint"))
    return float(v) if isinstance(v, (int, float)) else None


# ------------------------------------------------------------------------------------------------
# lock tracing


def lock_drives(spec: RunSpec, leaves: Sequence[Leaf] | None = None) -> list[dict[str, Any]]:
    """Every actuator writing a segment driver, with the lock that drives it.

    An actuator is a pipeline block outside the noise cluster that writes ``drv.<segment>.<driver>``.
    Its inputs are traced upstream through controller / actuator blocks (ctrl, act, mod, delay
    families, and err / det blocks that read no beam, e.g. demodulators) to detector blocks (err / det
    blocks reading ``beam.*``). The lock is the ``topology.locks`` entry whose beam combo
    shares the most beams with the beams the detectors read (``beam.<b>.*``), or the one whose true
    error ``lock.<id>.err`` a block of the chain reads (e.g. a lock supervisor); with a single lock and
    a detector found, that lock.

    Returns [{actuator, kind, segment, driver, lock (id or None), detectors [names]}], one row per
    written driver slot.
    """
    leaves = pipeline_leaves(spec) if leaves is None else leaves
    locks = list((_resolved(spec).topology or {}).get("locks") or [])
    lock_ids = {str(lk["id"]) for lk in locks}
    by_slot: dict[str, list[Leaf]] = defaultdict(list)
    for lf in leaves:
        for s in lf.writes:
            by_slot[s].append(lf)
    out = []
    for act in leaves:
        drv = [s for s in act.writes if s.startswith("drv.") and s.count(".") >= 2]
        if not drv or act.group == "noise":
            continue
        detectors: list[Leaf] = []
        named: set[str] = set()
        seen = {act.index}
        todo = [act]
        while todo:
            cur = todo.pop()
            for s in cur.reads:
                parts = s.split(".")
                if len(parts) == 3 and parts[0] == "lock" and parts[2] == "err" and parts[1] in lock_ids:
                    named.add(parts[1])
                for w in by_slot.get(s, ()):
                    if w.index in seen:
                        continue
                    seen.add(w.index)
                    if (w.family in _DETECTOR_FAMILIES or (w.family == "?" and w.stage == "detection")) and any(
                            r.startswith("beam.") for r in w.reads):
                        detectors.append(w)
                    elif w.family in _CHAIN_FAMILIES | _DETECTOR_FAMILIES or (w.family == "?" and w.group is None):
                        todo.append(w)
        beams = {s.split(".")[1] for d in detectors for s in d.reads if s.startswith("beam.")}
        lock = None
        if detectors and locks:
            score = {str(lk["id"]): len(beams & set(lk.get("combo") or {})) for lk in locks}
            best = max(score.values())
            if best > 0:
                lock = next(i for i, v in score.items() if v == best)
            elif len(locks) == 1:
                lock = str(locks[0]["id"])
        if lock is None and len(named) == 1:
            lock = next(iter(named))
        for s in drv:
            _, seg, driver = s.split(".", 2)
            out.append({"actuator": act.name, "kind": act.spec.kind, "segment": seg, "driver": driver,
                        "lock": lock, "detectors": [d.name for d in detectors]})
    return out


# ------------------------------------------------------------------------------------------------
# signal flow


def signal_flow_dot(spec: RunSpec, *, expand: Iterable[str] = (), title: str | None = None) -> str:
    """DOT source of the signal-flow diagram (module docstring). ``expand``: names of :data:`GROUPS`
    to draw block by block (``"all"`` expands every group)."""
    from cranberry.runner.run import offline_block_specs

    expand = set(expand)
    if "all" in expand:
        expand = set(GROUPS)
    unknown = expand - set(GROUPS)
    if unknown:
        raise ValueError(f"expand: unknown group(s) {sorted(unknown)}; known {GROUPS} or 'all'")
    rs = _resolved(spec)
    leaves = pipeline_leaves(rs)
    node_of: dict[int, str] = {}
    members: dict[str, list[Leaf]] = defaultdict(list)
    for lf in leaves:
        if lf.group is not None and lf.group not in expand:
            node_of[lf.index] = f"grp:{lf.group}"
            members[lf.group].append(lf)
        else:
            node_of[lf.index] = f"b:{lf.name}"

    edges: dict[tuple[str, str, bool], list[str]] = defaultdict(list)
    for w, r, s, back in _slot_edges(leaves):
        u, v = node_of[w], node_of[r]
        if u != v:
            edges[(u, v, back)].append(s)

    title = title or rs.preset or "configuration"
    lines = [
        "digraph signal_flow {",
        f"  graph [rankdir=LR, compound=true, newrank=true, fontname=Helvetica, fontsize=14, labelloc=t, "
        f"nodesep=0.25, ranksep=0.6, label={_lbl(f'{title}: signal flow (blocks = nodes, slots = edges)')}];",
        "  node [shape=box, style=\"rounded,filled\", fontname=Helvetica, fontsize=10];",
        "  edge [fontname=Helvetica, fontsize=7, color=\"#555555\", arrowsize=0.6];",
    ]

    def node_line(lf: Leaf, indent: str) -> str:
        fill, line = _STAGE_STYLE[lf.group or lf.stage]
        rate = f" (÷{lf.spec.rate_div})" if lf.spec.rate_div != 1 else ""
        return f"{indent}{_q(node_of[lf.index])} [label={_lbl(lf.name, lf.spec.kind + rate)}, fillcolor={_q(fill)}, color={_q(line)}];"

    clusters: dict[str, list[Leaf]] = defaultdict(list)
    for lf in leaves:
        clusters[lf.group if (lf.group and lf.group in expand) else (None if lf.group else lf.stage)].append(lf)
    for name in (*STAGES, *GROUPS):
        if name == "offline":
            continue
        items = clusters.get(name, [])
        collapsed = name in GROUPS and name not in expand and members.get(name)
        if not items and not collapsed:
            continue
        fill, line = _STAGE_STYLE[name]
        lines.append(f"  subgraph {_q('cluster_' + name)} {{")
        lines.append(f"    label={_q(name)}; style=\"rounded,dashed\"; color={_q(line)}; fontcolor={_q(line)};")
        if collapsed:
            ms = members[name]
            kinds = sorted({m.spec.kind for m in ms})
            names = ", ".join(m.name for m in ms)
            lines.append(f"    {_q('grp:' + name)} [shape=box3d, style=filled, fillcolor={_q(fill)}, color={_q(line)}, "
                         f"label={_lbl(f'{name}: {len(ms)} blocks (collapsed)', *kinds)}, tooltip={_q(names)}];")
        stacks: dict[str, list[Leaf]] = defaultdict(list)
        for lf in items:
            if lf.stack is not None:
                stacks[lf.stack].append(lf)
            else:
                lines.append(node_line(lf, "    "))
        for sp, lfs in stacks.items():
            lines.append(f"    subgraph {_q('cluster_stack_' + sp)} {{")
            lines.append(f"      label={_q('stack ' + sp)}; style=\"rounded\"; color=\"#b48ad0\"; fontsize=9;")
            lines.extend(node_line(lf, "      ") for lf in lfs)
            lines.append("    }")
        lines.append("  }")

    off = offline_block_specs(rs)
    if off:
        fill, line = _STAGE_STYLE["offline"]
        lines.append("  subgraph \"cluster_offline\" {")
        lines.append(f"    label=\"offline (no slots)\"; style=\"rounded,dashed\"; color={_q(line)}; fontcolor={_q(line)};")
        lines.append(f"    \"offline\" [shape=note, style=filled, fillcolor={_q(fill)}, color={_q(line)}, "
                     f"label={_lbl(*(f'{b.name}  [{b.kind}]' for b in off))}];")
        lines.append("  }")

    for (u, v, back), slots in edges.items():
        shown = slots[:_MAX_EDGE_SLOTS] + ([f"+{len(slots) - _MAX_EDGE_SLOTS} more"] if len(slots) > _MAX_EDGE_SLOTS else [])
        attrs = [f"label={_lbl(*shown)}", f"tooltip={_q(', '.join(slots))}"]
        if back:
            attrs += ["style=dashed", "color=\"#c62828\"", "fontcolor=\"#c62828\"", "constraint=false"]
        lines.append(f"  {_q(u)} -> {_q(v)} [{', '.join(attrs)}];")
    lines.append("}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------
# optical layout


def optical_layout_dot(spec: RunSpec, *, title: str | None = None) -> str:
    """DOT source of the optical-layout diagram (module docstring)."""
    rs = _resolved(spec)
    topo = rs.topology or {}
    segs = {str(s["id"]): s for s in topo.get("segments") or []}
    beams = list(topo.get("beams") or [])
    observers = {str(o["id"]): o for o in topo.get("observers") or []}
    locks = list(topo.get("locks") or [])
    color = {str(b["id"]): _BEAM_COLORS[i % len(_BEAM_COLORS)] for i, b in enumerate(beams)}

    pair: dict[tuple[str, str], list[tuple[str, bool]]] = defaultdict(list)
    sources: set[str] = set()
    for b in beams:
        bid = str(b["id"])
        src = f"src:{b['source']}"
        sources.add(str(b["source"]))
        chain = [(src, True)] + [(f"seg:{p['segment']}", int(p.get("direction", 1)) >= 0) for p in b.get("path") or []]
        if b.get("observer") is not None:
            chain.append((f"obs:{b['observer']}", True))
        for (u, _), (v, fwd) in zip(chain, chain[1:]):
            pair[(u, v)].append((bid, fwd))
    preds: dict[str, set[str]] = defaultdict(set)
    succs: dict[str, set[str]] = defaultdict(set)
    for u, v in pair:
        preds[v].add(u)
        succs[u].add(v)

    title = title or rs.preset or "configuration"
    lines = [
        "digraph optical_layout {",
        f"  graph [rankdir=LR, fontname=Helvetica, fontsize=14, labelloc=t, nodesep=0.3, ranksep=0.5, "
        f"label={_lbl(f'{title}: optical layout (beams on segments, locks)')}];",
        "  node [fontname=Helvetica, fontsize=10];",
        "  edge [fontname=Helvetica, fontsize=8, arrowsize=0.6];",
    ]
    for s in sorted(sources):
        lines.append(f"  {_q('src:' + s)} [shape=cds, style=filled, fillcolor=\"#fde2c4\", color=\"#c46a00\", label={_lbl(s, 'source')}];")
    zones: dict[str, list[str]] = defaultdict(list)
    for sid, s in segs.items():
        zones[str(s.get("zone", ""))].append(sid)
    for zone, ids in zones.items():
        if zone:
            lines.append(f"  subgraph {_q('cluster_zone_' + zone)} {{ label={_q('zone ' + zone)}; style=dotted; color=\"#999999\"; fontsize=9;")
        for sid in ids:
            s = segs[sid]
            nid = f"seg:{sid}"
            mark = ("⊕ " if len(preds[nid]) > 1 else "") + ("⑂ " if len(succs[nid]) > 1 else "")
            drivers = ", ".join(sorted((s.get("drivers") or {}).keys()))
            extra = [f"{s.get('kind')}, L = {float(s.get('L_m', 0.0)):g} m"] + ([f"drv: {drivers}"] if drivers else [])
            per = ", peripheries=2" if len(preds[nid]) > 1 else ""
            lines.append(f"    {_q(nid)} [shape=box, style=\"rounded,filled\", fillcolor=\"#d6e8fb\", color=\"#1f5fa8\"{per}, "
                         f"label={_lbl(mark + sid, *extra)}];")
        if zone:
            lines.append("  }")
    for oid, o in observers.items():
        lines.append(f"  {_q('obs:' + oid)} [shape=doublecircle, style=filled, fillcolor=\"#dff3dc\", color=\"#2e7d32\", "
                     f"fontsize=9, label={_lbl(oid, str(o.get('kind', '')))}];")
    for (u, v), bl in pair.items():
        ids = [b for b, _ in bl]
        cols = ":".join(color[b] for b in ids)
        style = "solid" if all(f for _, f in bl) else "dashed"
        lines.append(f"  {_q(u)} -> {_q(v)} [color={_q(cols)}, style={style}, label={_lbl(', '.join(ids))}, fontcolor=\"#333333\"];")

    drives = lock_drives(rs) if rs.pipeline.children else []
    for lk in locks:
        lid = str(lk["id"])
        combo = " ".join(f"{float(c):+g}·{b}" for b, c in (lk.get("combo") or {}).items())
        sp = _setpoint(lk)
        lines.append(f"  {_q('lock:' + lid)} [shape=note, style=filled, fillcolor=\"#ffe0e0\", color=\"#c62828\", "
                     f"label={_lbl('lock ' + lid, combo, *([f'set point {sp:.4g} rad'] if sp is not None else []))}];")
        det = sorted({d for r in drives if r["lock"] == lid for d in r["detectors"]})
        if lk.get("observer") is not None:
            lines.append(f"  {_q('obs:' + str(lk['observer']))} -> {_q('lock:' + lid)} [color=\"#c62828\", style=dashed, "
                         f"fontcolor=\"#c62828\", label={_lbl(*(det or ['error signal']))}, constraint=false];")
    acts: dict[str, dict[str, Any]] = {}
    for r in drives:
        acts.setdefault(r["actuator"], {"kind": r["kind"], "lock": r["lock"], "rows": []})["rows"].append(r)
    for name, a in acts.items():
        red = a["lock"] is not None
        col = "#c62828" if red else "#888888"
        lines.append(f"  {_q('act:' + name)} [shape=component, style=filled, fillcolor={_q('#ffe0e0' if red else '#f2f2f2')}, "
                     f"color={_q(col)}, fontsize=9, label={_lbl(name, a['kind'])}];")
        if red:
            lines.append(f"  {_q('lock:' + a['lock'])} -> {_q('act:' + name)} [color=\"#c62828\"];")
        for r in a["rows"]:
            if f"seg:{r['segment']}" not in {f"seg:{s}" for s in segs}:
                continue
            lines.append(f"  {_q('act:' + name)} -> {_q('seg:' + r['segment'])} [color={_q(col)}, fontcolor={_q(col)}, "
                         f"style=bold, label={_q(r['driver'])}];")
    lines.append("}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------
# rendering


def find_dot() -> str:
    """Path of the Graphviz ``dot`` executable (``PATH``, then ``/opt/homebrew/bin/dot``)."""
    exe = shutil.which("dot") or ("/opt/homebrew/bin/dot" if Path("/opt/homebrew/bin/dot").is_file() else None)
    if exe is None:
        raise FileNotFoundError("Graphviz 'dot' not found (brew install graphviz); the .dot source can still be written")
    return exe


def render(dot_src: str, stem: str | Path, formats: Sequence[str] = ("svg", "png")) -> list[Path]:
    """Write ``<stem>.<fmt>`` for each format: ``dot`` writes the source, others run Graphviz."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    out = []
    for fmt in formats:
        p = stem.with_name(f"{stem.name}.{fmt}")
        if fmt == "dot":
            p.write_text(dot_src)
        else:
            subprocess.run([find_dot(), f"-T{fmt}", "-o", str(p)], input=dot_src, text=True, check=True,
                           capture_output=True)
        out.append(p)
    return out


def write_diagrams(spec: RunSpec, out_dir: str | Path, name: str | None = None, *, expand: Iterable[str] = (),
                   formats: Sequence[str] = ("svg", "png"), title: str | None = None) -> list[Path]:
    """Both diagrams of ``spec`` as ``<out_dir>/<name>.signal_flow.<fmt>`` and ``<name>.optical_layout.<fmt>``."""
    name = name or spec.preset or "config"
    out = Path(out_dir)
    return [*render(signal_flow_dot(spec, expand=expand, title=title or name), out / f"{name}.signal_flow", formats),
            *render(optical_layout_dot(spec, title=title or name), out / f"{name}.optical_layout", formats)]
