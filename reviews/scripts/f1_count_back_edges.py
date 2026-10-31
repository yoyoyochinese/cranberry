"""Phase 5 M1 framework review, F1-1 companion: number of one-sample back edges on every feedback
cycle of the reference presets (D55-1 assumes exactly one per loop).

For every back edge j -> i (reader i runs before writer j) we take the shortest forward/back path
i -> ... -> j through the dataflow graph of core.pipeline.check_pipeline and count the back edges on
the closed cycle. Run: uv run python reviews/scripts/f1_count_back_edges.py
"""

from collections import deque

from cranberry import config
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, from_spec, iter_leaves
from cranberry.core.registry import ensure_loaded
from cranberry.runner.run import inject_topology


def graph(leaves):
    writers = {}
    for i, leaf in enumerate(leaves):
        for w in dict.fromkeys(leaf.writes):
            writers.setdefault(w, []).append(i)
    edges = {i: {} for i in range(len(leaves))}  # u -> {v: is_back}
    back = []
    for i, leaf in enumerate(leaves):
        for r in dict.fromkeys(leaf.reads):
            others = [j for j in writers.get(r, []) if j != i]
            if not others:
                continue
            prior = [j for j in others if j < i]
            if prior:
                edges[max(prior)].setdefault(i, False)
            else:
                j = max(writers[r])
                edges[j][i] = True
                back.append((j, i, r))
    return edges, back


def cycle_back_edges(edges, j, i):
    """Fewest back edges on a path i -> j (0-1 BFS), +1 for the closing edge j -> i."""
    dist, dq = {i: 0}, deque([i])
    while dq:
        u = dq.popleft()
        for v, is_back in edges[u].items():
            nd = dist[u] + int(is_back)
            if v not in dist or nd < dist[v]:
                dist[v] = nd
                (dq.append if is_back else dq.appendleft)(v)
    return None if j not in dist else dist[j] + 1


if __name__ == "__main__":
    ensure_loaded()
    for preset in ("presets/reference/refD_nominal.yaml", "presets/reference/refA.yaml"):
        rs = config.load(preset)
        root = from_spec(inject_topology(rs.pipeline, rs.topology))
        build(root, SlotLayout(rs.run.f_sim_hz), rs.run.seed)
        leaves = list(iter_leaves(root))
        edges, back = graph(leaves)
        print(f"{preset}:")
        for j, i, r in back:
            n = cycle_back_edges(edges, j, i)
            print(f"  back edge {leaves[j].name} -> {leaves[i].name} via {r}: back edges on its cycle = {n}")
