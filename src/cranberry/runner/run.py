"""``run(spec) -> RunResult``: offline stage → build → time loop (chunks) → measurement → record
(protected list "experiment runner"; PLAN §2.8). Owner: WP-F1 (records written via WP-F2's
``runner.record``).

The runner is a pure function of the spec: same spec (incl. seed) → bit-identical results,
independent of ``chunk_samples`` (VT-F2).

Minimal path implemented by WP-F1 (wave 1):

1. ``spec`` is a resolved ``RunSpec`` (a path / mapping is passed to ``config.load``, WP-F2) and is
   checked with ``config.validate`` (D34: no sweep axes left, unique names, every block builds).
2. ``RunSpec.toggles`` overrides are written into the leaf specs (``config.apply_toggles``).
3. With ``record_root``, the run directory is created and the config / seeds / git / env parts are
   written before the loop (``runner.record.record_run``; WP-F2b).
4. ``build`` on a fresh ``SlotLayout(run.f_sim_hz)`` with root seed ``run.seed``; when recording is
   requested, a :class:`Recorder` leaf is appended after the pipeline.
5. The time loop runs ``run.n_samples`` samples in chunks of ``run.chunk_samples`` (compiled loop
   when every leaf has a kernel, else the pure-Python reference; ``mode`` selects explicitly).
   After each chunk the recorder buffer is harvested.
6. ``RunResult`` with the recorded series, seeds and the final state; with ``record_root`` the
   results (``results.npz``: every recorded series + ``t_s``) are added to the run directory.

Not yet wired (later WPs): the offline stage and tables (WP-F3), metrics / PSDs (WP-F4). Their record
parts have hooks in :func:`_write_results` (``metrics``, ``tables``, ``loss_registry``,
``band_ledger``).

Recording (``RunSpec.measure["record"]``, optional; without it no series are recorded)::

    measure:
      record: {slots: [<slot name>, ...], decim: <int >= 1>}

Slot values are recorded at the end of every sample n with ``n % decim == 0`` (after all blocks
ran), at time t = n / f_sim. float64, int64 and complex128 slots of any (non-ring) shape are
supported; ``series[name]`` has shape ``(n_records, *slot.shape)``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock
from cranberry.core.rates import max_steps_in

if TYPE_CHECKING:
    from cranberry.blocks.meas.metrics import MetricsResult
    from cranberry.config.spec import RunSpec
    from cranberry.core.state import State

__all__ = ["RunResult", "run", "Recorder", "RUN_MODES"]

#: ``run(..., mode=)`` values.
RUN_MODES: tuple[str, ...] = ("auto", "compiled", "reference")


@dataclass
class RunResult:
    """Result of one run.

    Attributes:
        spec: the resolved spec that was run.
        run_dir: ``runs/<UTC time>_<preset>_<spechash8>/`` or None when recording is off.
        t_s: time axis of the recorded (decimated) series [s].
        series: {slot name: decimated time series} (true path phases, errors, u, δψ, ...).
        psd: {name: (f [Hz], one-sided PSD)}.
        metrics: metric (i)/(ii), overlap, stray light, η (``blocks.meas.metrics.MetricsResult``).
        events: lock / CHK events as JSON-serialisable dicts.
        table_hashes: {table name: hash} of offline tables used.
        seeds: {"root": entropy, "<stream name>": key, ...} (G11).
        final_state: the State after the last sample (WP-F1 addition; None if not kept).
        mode: "compiled" or "reference" (the loop actually used; WP-F1 addition).
    """

    spec: RunSpec
    run_dir: Path | None
    t_s: np.ndarray
    series: Mapping[str, np.ndarray] = field(default_factory=dict)
    psd: Mapping[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    metrics: MetricsResult | None = None
    events: list[dict[str, Any]] = field(default_factory=list)
    table_hashes: Mapping[str, str] = field(default_factory=dict)
    seeds: Mapping[str, Any] = field(default_factory=dict)
    final_state: State | None = None
    mode: str = ""


class Recorder(BaseBlock):
    """Runner-owned leaf that copies slots into a buffer every ``decim`` samples (``rate_div = decim``).

    Slots: ``rec.<name>.buf`` float64 ``(capacity, n_columns)`` and ``rec.<name>.idx`` int64 (rows
    written since the last :meth:`harvest`). ``capacity = ceil(chunk_samples / decim)``; the runner
    harvests after every chunk.

    Args:
        slots: slot names to record.
        decim: record every ``decim``-th base sample.
        chunk_samples: maximum samples between two harvests.
        name: instance name.
    """

    def __init__(self, slots: Sequence[str], decim: int, chunk_samples: int, name: str = "recorder"):
        if isinstance(decim, bool) or int(decim) != decim or decim < 1:
            raise ValueError(f"measure.record.decim must be an int >= 1, got {decim!r}")
        if not slots or len(set(slots)) != len(slots):
            raise ValueError(f"measure.record.slots must be a non-empty list of distinct names, got {slots!r}")
        self.buf_name, self.idx_name = f"rec.{name}.buf", f"rec.{name}.idx"
        super().__init__("recorder", name, {"slots": list(slots), "decim": int(decim)},
                         reads=tuple(slots), writes=(self.buf_name, self.idx_name), rate_div=int(decim))
        self.capacity = max_steps_in(int(chunk_samples), int(decim))
        self.chunk_samples = int(chunk_samples)

    def init(self, layout, rng):
        super().init(layout, rng)
        cols = 0
        for s in self.reads:
            if s not in layout:
                continue  # reported by check_pipeline
            sl = layout.slot(s)
            if sl.ring_len:
                raise ValueError(f"recorder: ring slot {s!r} cannot be recorded")
            cols += sl.size * (2 if sl.array == "Z" else 1)
        self.n_cols = cols
        layout.add(self.buf_name, (self.capacity, max(cols, 1)), np.float64, "1")
        layout.add(self.idx_name, (), np.int64, "1")

    def bind(self, layout):
        src = {"S": [], "I": [], "Z": []}
        self.columns: list[tuple[str, int, int]] = []  # (slot, first column, n columns)
        col = 0
        # column order: all S elements, then I, then Z (re, im interleaved); remembered per slot
        for arr in ("S", "I", "Z"):
            for s in self.reads:
                sl = layout.slot(s)
                if sl.array != arr:
                    continue
                ncols = sl.size * (2 if arr == "Z" else 1)
                self.columns.append((s, col, ncols))
                src[arr].extend(range(sl.offset, sl.stop))
                col += ncols
        s_src, i_src, z_src = (np.array(src[a], dtype=np.int64) for a in ("S", "I", "Z"))
        self._src = (s_src, i_src, z_src)
        buf, idx, width, cap = layout.slot(self.buf_name).offset, layout.slot(self.idx_name).offset, max(self.n_cols, 1), self.capacity
        ns, ni, nz = len(s_src), len(i_src), len(z_src)
        self._off = (buf, idx, width, cap)

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            row = I[idx]
            if row >= cap:
                raise IndexError("recorder buffer overflow (harvest after every chunk)")
            base = buf + row * width
            for j in range(ns):
                S[base + j] = S[s_src[j]]
            for j in range(ni):
                S[base + ns + j] = I[i_src[j]]
            for j in range(nz):
                S[base + ns + ni + 2 * j] = Z[z_src[j]].real
                S[base + ns + ni + 2 * j + 1] = Z[z_src[j]].imag
            I[idx] = row + 1

        self.kernel = kernel

    def __call__(self, state, params):
        buf, idx, width, cap = self._off
        s_src, i_src, z_src = self._src
        row = int(state.I[idx])
        if row >= cap:
            raise IndexError("recorder buffer overflow (harvest after every chunk)")
        base = buf + row * width
        S = state.S
        ns, ni = len(s_src), len(i_src)
        S[base : base + ns] = S[s_src]
        S[base + ns : base + ns + ni] = state.I[i_src]
        z = state.Z[z_src]
        S[base + ns + ni : base + ns + ni + 2 * len(z_src) : 2] = z.real
        S[base + ns + ni + 1 : base + ns + ni + 2 * len(z_src) : 2] = z.imag
        state.I[idx] = row + 1
        return state

    def harvest(self, state) -> np.ndarray:
        """Return the rows written since the last harvest (copy) and reset the row index."""
        rows = int(state[self.idx_name])
        out = state[self.buf_name][:rows, : self.n_cols].copy()
        state[self.idx_name] = 0
        return out

    def split(self, table: np.ndarray, layout) -> dict[str, np.ndarray]:
        """Columns of harvested rows → {slot: array (n_rows, *shape)}."""
        out: dict[str, np.ndarray] = {}
        for s, c0, nc in self.columns:
            sl = layout.slot(s)
            block = table[:, c0 : c0 + nc]
            if sl.array == "Z":
                arr = block[:, 0::2] + 1j * block[:, 1::2]
            elif sl.array == "I":
                arr = block.astype(np.int64)
            else:
                arr = block
            out[s] = arr.reshape((table.shape[0], *sl.shape))
        return out


def _as_runspec(spec: Any) -> RunSpec:
    from cranberry.config.spec import RunSpec

    if isinstance(spec, RunSpec):
        return spec
    from cranberry import config

    return config.load(spec)  # WP-F2 (raises NotImplementedError until it lands)


def _write_results(run_dir: Path, result: RunResult) -> None:
    """Add the results of ``result`` to its run directory (after the loop)."""
    from cranberry.runner.record import write_record

    parts: dict[str, Any] = {"results": {**result.series, "t_s": result.t_s}}
    # Hooks for later WPs: they fill RunResult (or pass the objects here) and the parts follow.
    if result.metrics is not None:  # WP-F4: metrics (i)/(ii), PSDs
        parts["metrics"] = result.metrics
    if result.table_hashes:  # WP-F3: offline tables used
        parts["tables"] = dict(result.table_hashes)
    # WP-F3: parts["loss_registry"] = <LossRegistry>, parts["band_ledger"] = <BandLedger> (to_json())
    write_record(run_dir, **parts)


def run(
    spec: RunSpec | str | Path | Mapping[str, Any],
    *,
    mode: str = "auto",
    keep_state: bool = True,
    record_root: str | Path | None = None,
    config_input: Any = None,
    invocation: Mapping[str, Any] | None = None,
) -> RunResult:
    """Run one fully resolved configuration (no sweep markers left).

    Args:
        spec: a ``RunSpec``; a path or mapping is loaded with ``config.load`` (WP-F2).
        mode: "auto" (compiled if every leaf has a kernel, else reference), "compiled", "reference".
        keep_state: store the final State in ``RunResult.final_state``.
        record_root: write a run record under this directory (``runs/<UTC>_<preset>_<hash8>/``,
            ``runner.record``); None (default) writes nothing (D34: no silent default directory).
        config_input: the input document for ``config.input.yaml`` (``config.read_input(...)``);
            None → read from ``spec`` when it is a path / mapping, else the spec itself.
        invocation: extra context for ``invocation.json`` (CLI argv, overrides, sweep point).
    """
    from cranberry import config
    from cranberry.config.spec import RunSpec, apply_toggles
    from cranberry.core.compile import compile_loop, is_compilable, run_compiled, run_reference
    from cranberry.core.layout import SlotLayout
    from cranberry.core.pipeline import Seq, build, from_spec, iter_leaves
    from cranberry.core.registry import ensure_loaded
    from cranberry.core.rng import root_seed_sequence, stream_key
    from cranberry.core.state import Params, State

    if mode not in RUN_MODES:
        raise ValueError(f"mode must be one of {RUN_MODES}, got {mode!r}")
    rs = _as_runspec(spec)
    settings = rs.run
    ensure_loaded()
    config.validate(rs)

    run_dir: Path | None = None
    if record_root is not None:
        from cranberry.runner.record import record_run

        if config_input is None and not isinstance(spec, RunSpec):
            config_input = config.read_input(spec)
        run_dir = record_run(rs, Path(record_root), config_input=config_input,
                             **({"invocation": dict(invocation)} if invocation else {}))

    pipeline_spec = apply_toggles(rs.pipeline, rs.toggles)
    root = from_spec(pipeline_spec)
    rec_cfg = (rs.measure or {}).get("record")
    recorder = None
    if rec_cfg is not None:
        missing = [k for k in ("slots", "decim") if k not in rec_cfg]
        if missing:
            raise KeyError(f"measure.record needs {missing} (D34: no defaults)")
        recorder = Recorder(list(rec_cfg["slots"]), rec_cfg["decim"], settings.chunk_samples)
        root = Seq([root, recorder])

    layout = SlotLayout(settings.f_sim_hz)
    block = build(root, layout, settings.seed)
    state, params = State(layout), Params(layout)

    compiled = mode == "compiled" or (mode == "auto" and is_compilable(block))
    loop = compile_loop(block, layout) if compiled else None

    n_total, chunk = settings.n_samples, settings.chunk_samples
    tables: list[np.ndarray] = []
    done = 0
    while done < n_total:
        m = min(chunk, n_total - done)
        if compiled:
            run_compiled(block, state, params, m, loop)
        else:
            run_reference(block, state, params, m)
        if recorder is not None:
            tables.append(recorder.harvest(state))
        done += m

    series: dict[str, np.ndarray] = {}
    t_s = np.zeros(0)
    if recorder is not None:
        table = np.concatenate(tables, axis=0) if tables else np.zeros((0, recorder.n_cols))
        series = recorder.split(table, layout)
        n_rec = np.arange(0, n_total, recorder.rate_div, dtype=np.int64)
        t_s = n_rec / settings.f_sim_hz

    root_ss = root_seed_sequence(settings.seed)
    seeds: dict[str, Any] = {"root": root_ss.entropy}
    seeds.update({leaf.name: list(stream_key(leaf.name)) for leaf in iter_leaves(block) if leaf is not recorder})
    result = RunResult(
        spec=rs,
        run_dir=run_dir,
        t_s=t_s,
        series=series,
        seeds=seeds,
        final_state=state if keep_state else None,
        mode="compiled" if compiled else "reference",
    )
    if run_dir is not None:
        _write_results(run_dir, result)
    return result
