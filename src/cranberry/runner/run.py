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

Offline stage and measurement (WP-F4)
    * **Topology injection** (F6 proposal): a leaf of a kind in :data:`TOPOLOGY_KINDS` without a
      ``topology`` param gets ``RunSpec.topology`` at build time, so presets need not repeat it (YAML
      anchors). An explicit ``topology`` param is kept as given. The recorded config is the spec as
      written (injection is not stored; it is redone on reload).
    * **Offline-only blocks**: ``RunSpec.offline["blocks"]`` is a list of block specs (same form as
      pipeline leaves: ``{block, name, params, rate_div?}``) built for the offline stage only and never
      stepped (source models such as ``single_mode_reference`` / ``opa_waveguide``: no kernel, so
      keeping them out of the loop keeps it compilable). ``RunSpec.toggles`` may address them by name
      (``ablate``). Pipeline leaves take part in the offline stage too.
    * **Offline stage** (:func:`build_offline_stage`, before the time loop): ``source_output`` of every
      block that has one → {block name: ``SourceOutput``}; ``LossRegistry`` with one path per topology
      beam (``Topology.register_loss_paths``) and every block's ``losses()``; ``BandLedger`` with the
      topology beams (``register_band_beams``), every block's ``bands()`` and one metric per metric
      probe port with ``f_max_hz`` / ``weight`` from ``measure.metrics`` (``meas.metrics.define_metrics``).
      Recorded as ``loss_registry.json`` / ``band_ledger.json`` (when non-empty) and ``tables.json``
      (hashes of the source outputs).
    * **Measurement** (after the loop): when the pipeline holds ``metric_i`` / ``metric_ii`` probes,
      their slots are added to ``measure.record`` (if recording) and ``meas.metrics.evaluate_metrics``
      fills ``RunResult.metrics`` (``metrics.json``); the out-of-band ledger variance enters there, in
      ``phase_to_quadrature`` only (G05).

Measurement-stage hooks (WP-I1, D58-6; :func:`evaluate_measurements`, after the metrics)
    Blocks of the measurement stage (``esa_zero_span``, ``esa_swept_spectrum``, ``pulse_sample``;
    WP-B5) have no-op kernels; the runner calls their ``evaluate`` on the results of the loop. They are
    named in ``offline.blocks`` (or the pipeline). Every entry key is required (D34)::

      measure:
        esa:
          <entry>:                       # esa_zero_span: V′(t) from the recorded metric series
            {block, metric, quadrature_rad}
          <entry>:                       # esa_swept_spectrum: R through a bhd_frontend
            {block, metric, quadrature_rad, frontend, exclude_losses}
        pulse_sample:
          <entry>: {block, kernel, r_slots}   # r_slots: K items, each a number R (shot-noise units)
                                              # or {metric, quadrature_rad}

    ``metric`` names a metric with noise levels (``measure.metrics``), e.g. ``i.port0``;
    ``quadrature_rad`` is added to the measured angle Δ = Δ_set + δ − θ_off (0: the locked quadrature,
    π/2: the same jitter around the orthogonal one, prompt §7-1). The variance uses the metric's loss-
    applied V∓ and its χ_HF, χ_dither (``phase_to_quadrature``; nothing is added twice):

    * zero span: V′(t) = (V− + V+)/2 − (V+ − V−)/2 · cos 2Δ(t) · χ_HF χ_dither, Δ(t) from the recorded
      (decimated) metric series, sampled every ``decim/f_sim``; the record must cover ``sweep_time_s``
      (T04 §5.1). A phase-scanned trace is the same entry on a run with the lock off and a ramp
      (Ref A Fig. 3).
    * swept spectrum: R = V′/V_SN with ⟨cos 2Δ⟩ from the full-rate metric moments and η recomputed without
      the LossRegistry elements listed in ``exclude_losses`` (the detector's electronic / imbalance
      elements: ``bhd_frontend`` adds them itself, D58-1); ``frontend.psd(f, R)`` gives the signal /
      LO-only / dark PSDs (T04 §5.5, T13).
    * pulse sample: ``evaluate(r_slots, c=kernel.c, clearance_lin=kernel.clearance_pulse)``.

    Random numbers: ``core.rng.stream(run.seed, "measure.<hook>.<entry>")`` (G11). Results go to
    ``RunResult.measurements[<hook>][<entry>]`` (the block's ``evaluate`` dict plus ``v_prime`` /
    ``r_lin``) and, when recording, to ``results.npz`` as ``measure.<hook>.<entry>.<field>`` (numeric
    fields).

Recording (``RunSpec.measure["record"]``, optional; without it no series are recorded)::

    measure:
      record: {slots: [<slot name>, ...], decim: <int >= 1>}

Slot values are recorded at the end of every sample n with ``n % decim == 0`` (after all blocks
ran), at time t = n / f_sim. float64, int64 and complex128 slots of any (non-ring) shape are
supported; ``series[name]`` has shape ``(n_records, *slot.shape)``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, BlockSpec
from cranberry.core.rates import max_steps_in

if TYPE_CHECKING:
    from cranberry.blocks.meas.metrics import MetricsResult
    from cranberry.config.spec import RunSpec
    from cranberry.core.pipeline import SeqSpec
    from cranberry.core.state import State
    from cranberry.offline.band_ledger import BandLedger
    from cranberry.offline.loss_registry import LossRegistry
    from cranberry.offline.source_output import SourceOutput
    from cranberry.path.topology import Topology

__all__ = [
    "RunResult",
    "run",
    "Recorder",
    "RUN_MODES",
    "TOPOLOGY_KINDS",
    "OfflineStage",
    "inject_topology",
    "offline_block_specs",
    "build_offline_stage",
    "MEASUREMENT_HOOKS",
    "evaluate_measurements",
]

#: ``run(..., mode=)`` values.
RUN_MODES: tuple[str, ...] = ("auto", "compiled", "reference")

#: Block kinds that take the ``topology`` param and receive ``RunSpec.topology`` when it is absent
#: (WP-F4, F6 proposal). WPs adding such kinds extend this set.
TOPOLOGY_KINDS: frozenset[str] = frozenset(
    {"segment_phase", "beam_assembler", "lock_error", "metric_i", "metric_ii", "segment_jones", "opo_angle_map",
     # segment noise sources (WP-B2 kinds, added by WP-I1: they take the topology param too)
     "segment_mech_noise", "segment_thermodynamic", "segment_thermal"}
)


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
        offline: the offline stage (source outputs, LossRegistry, BandLedger; WP-F4 addition).
        measurements: {hook: {entry: result}} of the measurement-stage hooks (``measure.esa``,
            ``measure.pulse_sample``; WP-I1 addition, :func:`evaluate_measurements`).
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
    offline: OfflineStage | None = None
    measurements: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)

    def summary(self) -> dict[str, float]:
        """Flat {name: scalar} for ablation / convergence tables (WP-F4): ``MetricsResult.summary()``
        plus ``<slot>.mean`` / ``<slot>.std`` / ``<slot>.rms`` of every recorded real scalar series."""
        out: dict[str, float] = {}
        for name, x in self.series.items():
            x = np.asarray(x)
            if x.ndim == 1 and x.size and np.isrealobj(x):
                xf = x.astype(float)
                out[f"{name}.mean"] = float(xf.mean())
                out[f"{name}.std"] = float(xf.std())
                out[f"{name}.rms"] = float(np.sqrt(np.mean(xf * xf)))
        if self.metrics is not None:
            out.update(self.metrics.summary())
        return out


@dataclass
class OfflineStage:
    """Products of the offline stage of one run (PLAN §2.1, §2.2; WP-F4).

    Attributes:
        topology: parsed ``RunSpec.topology`` (None when the spec has none).
        blocks: offline-only blocks built from ``RunSpec.offline["blocks"]``.
        source_outputs: {block name: ``SourceOutput``} of every block exposing ``source_output``.
        loss_registry: η per (beam path, target) (G06).
        band_ledger: out-of-band variance per (metric, beam combination) (G05).
        metric_names: metrics defined in the ledger (``i.port<k>`` / ``ii.port<k>``).
    """

    topology: Topology | None
    blocks: list[Any]
    source_outputs: dict[str, SourceOutput]
    loss_registry: LossRegistry
    band_ledger: BandLedger
    metric_names: list[str] = field(default_factory=list)

    def table_hashes(self) -> dict[str, str]:
        """{"source.<block>": sha256} of each SourceOutput (basis, grid, cov, modes, φ_OPA,0)."""
        import hashlib

        out = {}
        for name, so in self.source_outputs.items():
            h = hashlib.sha256(so.basis.encode())
            for arr in (so.grid, so.cov, so.modes):
                if arr is not None:
                    a = np.ascontiguousarray(arr)
                    h.update(str((a.dtype.str, a.shape)).encode())
                    h.update(a.tobytes())
            h.update(float(so.phi_opa0_rad).hex().encode())
            out[f"source.{name}"] = h.hexdigest()
        return out

    def record_parts(self) -> dict[str, Any]:
        """Non-empty registry / ledger parts for the run record."""
        parts: dict[str, Any] = {}
        if self.loss_registry.paths or self.loss_registry.elements:
            parts["loss_registry"] = self.loss_registry
        led = self.band_ledger
        if led.decls or led.metrics or led.beams:
            parts["band_ledger"] = led
        return parts


def _map_leaf_specs(spec: Any, fn) -> Any:
    import dataclasses

    from cranberry.core.pipeline import ParallelSpec, SeqSpec, StackSpec

    if isinstance(spec, BlockSpec):
        return fn(spec)
    if isinstance(spec, (SeqSpec, ParallelSpec)):
        return dataclasses.replace(spec, children=tuple(_map_leaf_specs(c, fn) for c in spec.children))
    if isinstance(spec, StackSpec):
        return dataclasses.replace(spec, stages=tuple(_map_leaf_specs(c, fn) for c in spec.stages),
                                   actuators=tuple(_map_leaf_specs(c, fn) for c in spec.actuators))
    raise TypeError(f"not a pipeline spec: {spec!r}")


def _inject(b: BlockSpec, topology: Mapping[str, Any]) -> BlockSpec:
    if b.kind in TOPOLOGY_KINDS and "topology" not in b.params and topology:
        return BlockSpec(b.kind, b.name, {**b.params, "topology": topology}, b.rate_div)
    return b


def inject_topology(pipeline: SeqSpec, topology: Mapping[str, Any]) -> SeqSpec:
    """``pipeline`` with ``topology`` added to every :data:`TOPOLOGY_KINDS` leaf that lacks it."""
    return _map_leaf_specs(pipeline, lambda b: _inject(b, topology))


def offline_block_specs(spec: RunSpec) -> list[BlockSpec]:
    """Block specs of ``RunSpec.offline["blocks"]`` (topology injected as for pipeline leaves)."""
    from cranberry.config.spec import ConfigError

    raw = (spec.offline or {}).get("blocks", [])
    if not isinstance(raw, (list, tuple)):
        raise ConfigError("offline.blocks must be a list of block specs {block, name, params}")
    out = []
    for i, obj in enumerate(raw):
        if not isinstance(obj, Mapping):
            raise ConfigError(f"offline.blocks[{i}] must be a block spec mapping, got {obj!r}")
        try:
            out.append(_inject(BlockSpec.from_obj(obj), spec.topology))
        except (ValueError, TypeError) as e:
            raise ConfigError(f"offline.blocks[{i}]: {e}") from None
    return out


def build_offline_stage(
    topology_spec: Mapping[str, Any],
    blocks: Sequence[Any],
    offline_blocks: Sequence[Any],
    measure: Mapping[str, Any],
) -> OfflineStage:
    """Offline stage of one run (module docstring): collect source outputs, build the LossRegistry
    and the BandLedger from the topology and every block (``blocks``: pipeline leaves,
    ``offline_blocks``: offline-only blocks) and define the metrics of the metric probes."""
    from cranberry.blocks.meas.metrics import define_metrics, metric_blocks
    from cranberry.offline.band_ledger import BandLedger
    from cranberry.offline.loss_registry import LossRegistry
    from cranberry.path.topology import Topology

    # A topology section without the path-model keys (e.g. parameter trees used only by expressions)
    # is not a Topology: no beam paths then. A complete one is parsed (and must be valid).
    full = bool(topology_spec) and all(k in topology_spec for k in ("segments", "beams", "observers", "locks",
                                                                      "ports", "coef"))
    topo = Topology.from_spec(topology_spec) if full else None
    every = [*offline_blocks, *blocks]
    sources: dict[str, Any] = {}
    for b in every:
        so = getattr(b, "source_output", None)
        if so is not None:
            sources[b.name] = so
    reg = LossRegistry()
    led = BandLedger()
    if topo is not None:
        topo.register_loss_paths(reg)
        topo.register_band_beams(led)
    for b in every:
        reg.add_all(b.losses())
        led.declare_all(b.bands())
    names = define_metrics(led, metric_blocks(blocks), measure)
    return OfflineStage(topo, list(offline_blocks), sources, reg, led, names)


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


#: Measurement-stage hooks of ``RunSpec.measure`` and the block kinds they drive (WP-I1, D58-6).
MEASUREMENT_HOOKS: dict[str, tuple[str, ...]] = {
    "esa": ("esa_zero_span", "esa_swept_spectrum"),
    "pulse_sample": ("pulse_sample",),
}


def _metric_entry(metrics: Any, name: str, what: str) -> tuple[Mapping[str, Any], Any]:
    """(noise levels, metric result object) of metric ``name`` (``i.port<k>`` / ``ii.port<k>``)."""
    if metrics is None or name not in metrics.noise_levels:
        known = [] if metrics is None else sorted(metrics.noise_levels)
        raise KeyError(f"{what}: metric {name!r} has no noise levels (needs a metric probe and measure.metrics; "
                       f"known: {known})")
    short, _, port = name.partition(".port")
    table = metrics.metric_i if short == "i" else metrics.metric_ii
    return metrics.noise_levels[name], table[int(port)]


def _quadrature_variance(nl: Mapping[str, Any], stats: Any, quadrature_rad: float, eta: float | None = None) -> float:
    """V′ [conv units] of a metric at Δ_set + ``quadrature_rad`` from its full-rate moments, with the
    loss-applied V∓ of its noise levels (η replaced by ``eta`` when given) and its χ_HF χ_dither (same
    algebra as ``phase_to_quadrature``: ⟨cos 2Δ⟩ = cos 2a⟨cos 2δ′⟩ − sin 2a⟨sin 2δ′⟩, δ′ = δ − θ_off)."""
    from cranberry.core.units import V_VAC

    vm, vp = float(nl["v_minus"]), float(nl["v_plus"])
    if eta is not None:
        e0 = float(nl["eta"])
        vm, vp = ((v - (1.0 - e0) * V_VAC) / e0 * eta + (1.0 - eta) * V_VAC for v in (vm, vp))
    a = float(nl["setpoint_rad"]) + float(quadrature_rad)
    st = stats.shifted(float(nl["theta_off_rad"]))
    c = math.cos(2 * a) * float(st.mean_cos2) - math.sin(2 * a) * float(st.mean_sin2)
    chi = c * float(nl["chi_hf"]) * float(nl["chi_dither"])
    return 0.5 * (vm + vp) - 0.5 * (vp - vm) * chi


def _need(entry: Any, keys: Sequence[str], what: str) -> None:
    if not isinstance(entry, Mapping):
        raise TypeError(f"{what} must be a mapping, got {entry!r}")
    missing = [k for k in keys if k not in entry]
    if missing:
        raise KeyError(f"{what} needs {missing} (D34: no defaults)")


def _zero_span(blk: Any, cfg: Mapping[str, Any], metrics: Any, f_record_hz: float | None, rng, what: str) -> dict:
    from cranberry.core.units import V_SN

    _need(cfg, ("block", "metric", "quadrature_rad"), what)
    nl, m = _metric_entry(metrics, cfg["metric"], what)
    delta = np.asarray(getattr(m, "delta_rad", np.zeros(0)), dtype=float)
    if f_record_hz is None or delta.size < 2:
        raise ValueError(f"{what}: the zero span reads the recorded metric series (needs measure.record)")
    dt = 1.0 / f_record_hz
    if (delta.size - 1) * dt < float(blk.params["sweep_time_s"]):
        raise ValueError(f"{what}: record of {(delta.size - 1) * dt:.4g} s is shorter than sweep_time_s "
                         f"{blk.params['sweep_time_s']} s")
    big = float(nl["setpoint_rad"]) + delta - float(nl["theta_off_rad"]) + float(cfg["quadrature_rad"])
    vm, vp = float(nl["v_minus"]), float(nl["v_plus"])
    v_t = 0.5 * (vm + vp) - 0.5 * (vp - vm) * np.cos(2.0 * big) * float(nl["chi_hf"]) * float(nl["chi_dither"])
    r = dict(blk.evaluate(v_t, dt_s=dt, rng=rng))
    r["v_prime"] = float(np.mean(v_t))
    r["r_lin"] = r["v_prime"] / V_SN
    return r


def _swept(blk: Any, cfg: Mapping[str, Any], blocks: Mapping[str, Any], metrics: Any, loss_registry: Any, rng,
           what: str) -> dict:
    from cranberry.core.units import V_SN

    _need(cfg, ("block", "metric", "quadrature_rad", "frontend", "exclude_losses"), what)
    nl, m = _metric_entry(metrics, cfg["metric"], what)
    fe = blocks.get(cfg["frontend"])
    if fe is None or fe.kind != "bhd_frontend":
        raise KeyError(f"{what}.frontend: no bhd_frontend block {cfg['frontend']!r}")
    eta = float(nl["eta"])
    excl = list(cfg["exclude_losses"])
    if excl:
        if loss_registry is None or nl["loss_path"] is None:
            raise ValueError(f"{what}: exclude_losses needs the metric's LossRegistry path")
        on_path = loss_registry.elements_on(nl["loss_path"], nl["loss_target"], include_disabled=True)
        unknown = sorted(set(excl) - {e.id for e in on_path})
        if unknown:
            raise KeyError(f"{what}.exclude_losses: {unknown} not on path {nl['loss_path']!r} ({nl['loss_target']})")
        kept = loss_registry.elements_on(nl["loss_path"], nl["loss_target"])
        eta = float(np.prod([e.eta for e in kept if e.id not in excl]))
    v = _quadrature_variance(nl, m.stats, float(cfg["quadrature_rad"]), eta)
    r_lin = v / V_SN

    def part(key):
        return lambda f: fe.psd(f, r_lin)[key]

    r = dict(blk.evaluate(part("signal"), part("shot"), part("dark"), rng=rng))
    r.update({"v_prime": v, "r_lin": r_lin, "eta": eta})
    return r


def _pulse(blk: Any, cfg: Mapping[str, Any], blocks: Mapping[str, Any], metrics: Any, rng, what: str) -> dict:
    from cranberry.core.units import V_SN

    _need(cfg, ("block", "kernel", "r_slots"), what)
    ker = blocks.get(cfg["kernel"])
    if ker is None or ker.kind != "pulsed_bhd_kernel":
        raise KeyError(f"{what}.kernel: no pulsed_bhd_kernel block {cfg['kernel']!r}")
    r_slots = []
    for i, item in enumerate(cfg["r_slots"]):
        if isinstance(item, Mapping):
            _need(item, ("metric", "quadrature_rad"), f"{what}.r_slots[{i}]")
            nl, m = _metric_entry(metrics, item["metric"], f"{what}.r_slots[{i}]")
            r_slots.append(_quadrature_variance(nl, m.stats, float(item["quadrature_rad"])) / V_SN)
        else:
            r_slots.append(float(item))
    r = dict(blk.evaluate(np.array(r_slots), c=ker.c, clearance_lin=ker.clearance_pulse, rng=rng))
    r["r_slots"] = np.array(r_slots)
    return r


def evaluate_measurements(
    measure: Mapping[str, Any],
    blocks: Mapping[str, Any],
    *,
    metrics: Any,
    loss_registry: Any,
    f_record_hz: float | None,
    seed: int,
) -> dict[str, dict[str, Any]]:
    """Run the measurement-stage hooks of ``measure`` after the loop (module docstring).

    Args:
        measure: ``RunSpec.measure``.
        blocks: {name: block} of the offline blocks and pipeline leaves.
        metrics: the run's ``MetricsResult`` (noise levels, recorded metric series, full-rate moments).
        loss_registry: the offline stage's LossRegistry (``exclude_losses``).
        f_record_hz: sample rate of the recorded series (f_sim / decim); None without recording.
        seed: run seed (RNG streams ``measure.<hook>.<entry>``).
    """
    from cranberry.core.rng import stream

    out: dict[str, dict[str, Any]] = {}
    for hook, kinds in MEASUREMENT_HOOKS.items():
        entries = (measure or {}).get(hook)
        if not entries:
            continue
        if not isinstance(entries, Mapping):
            raise TypeError(f"measure.{hook} must be a mapping {{entry: {{...}}}}")
        res: dict[str, Any] = {}
        for entry, cfg in entries.items():
            what = f"measure.{hook}.{entry}"
            _need(cfg, ("block",), what)
            blk = blocks.get(cfg["block"])
            if blk is None or blk.kind not in kinds:
                raise KeyError(f"{what}.block: no block {cfg['block']!r} of kind {kinds} in offline.blocks / pipeline")
            rng = stream(seed, what)
            if blk.kind == "esa_zero_span":
                res[str(entry)] = _zero_span(blk, cfg, metrics, f_record_hz, rng, what)
            elif blk.kind == "esa_swept_spectrum":
                res[str(entry)] = _swept(blk, cfg, blocks, metrics, loss_registry, rng, what)
            else:
                res[str(entry)] = _pulse(blk, cfg, blocks, metrics, rng, what)
        out[hook] = res
    return out


def _measurement_arrays(measurements: Mapping[str, Mapping[str, Any]]) -> dict[str, np.ndarray]:
    """Numeric fields of the hook results as ``measure.<hook>.<entry>.<field>`` arrays (results.npz)."""
    arrays: dict[str, np.ndarray] = {}
    for hook, entries in measurements.items():
        for entry, r in entries.items():
            for k, v in r.items():
                if v is None or isinstance(v, str):
                    continue
                a = np.asarray(v)
                if a.dtype.kind in "biuf":
                    arrays[f"measure.{hook}.{entry}.{k}"] = a
    return arrays


def _as_runspec(spec: Any) -> RunSpec:
    from cranberry.config.spec import RunSpec

    if isinstance(spec, RunSpec):
        return spec
    from cranberry import config

    return config.load(spec)  # WP-F2 (raises NotImplementedError until it lands)


def _write_results(run_dir: Path, result: RunResult) -> None:
    """Add the results of ``result`` to its run directory (after the loop)."""
    from cranberry.runner.record import write_record

    parts: dict[str, Any] = {"results": {**result.series, "t_s": result.t_s,
                                         **_measurement_arrays(result.measurements)}}
    # Hooks for later WPs: they fill RunResult (or pass the objects here) and the parts follow.
    if result.metrics is not None:  # WP-F4: metrics (i)/(ii), PSDs
        parts["metrics"] = result.metrics
    if result.table_hashes:  # WP-F3: offline tables used
        parts["tables"] = dict(result.table_hashes)
    if result.offline is not None:  # LossRegistry / BandLedger (to_json(); WP-F4 offline stage)
        parts.update(result.offline.record_parts())
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
    import dataclasses

    from cranberry import config
    from cranberry.blocks.meas.metrics import evaluate_metrics, metric_blocks
    from cranberry.config.spec import ConfigError, RunSpec, apply_toggles, iter_block_specs
    from cranberry.core.compile import compile_loop, is_compilable, run_compiled, run_reference
    from cranberry.core.layout import SlotLayout
    from cranberry.core.pipeline import Seq, SeqSpec, build, from_spec, iter_leaves
    from cranberry.core.registry import ensure_loaded, make_block
    from cranberry.core.rng import root_seed_sequence, stream_key
    from cranberry.core.state import Params, State

    if mode not in RUN_MODES:
        raise ValueError(f"mode must be one of {RUN_MODES}, got {mode!r}")
    rs = _as_runspec(spec)
    settings = rs.run
    ensure_loaded()

    # topology injection, offline-only blocks and their toggle overrides (WP-F4)
    off_specs = offline_block_specs(rs)
    names = [b.name for b in off_specs]
    off_names = set(names)
    dups = sorted({n for n in names if names.count(n) > 1} | (off_names & {b.name for _, b in iter_block_specs(rs.pipeline)}))
    if dups:
        raise ConfigError(f"block names must be unique across pipeline and offline.blocks: {dups}")
    pipe_toggles = {k: v for k, v in rs.toggles.items() if k.rpartition(".")[0] not in off_names}
    off_toggles = {k: v for k, v in rs.toggles.items() if k not in pipe_toggles}
    pipeline = inject_topology(rs.pipeline, rs.topology)
    config.validate(dataclasses.replace(rs, pipeline=pipeline, toggles=pipe_toggles))
    if off_toggles:
        off_specs = list(apply_toggles(SeqSpec(tuple(off_specs)), off_toggles).children)
    offline_blocks = []
    for i, b in enumerate(off_specs):
        try:
            offline_blocks.append(make_block(b))
        except (KeyError, ValueError, TypeError) as e:
            raise ConfigError(f"offline.blocks[{i}]: block {b.name!r} ({b.kind}): {e}") from None
    pipeline_spec = apply_toggles(pipeline, pipe_toggles)
    root = from_spec(pipeline_spec)
    leaves = list(iter_leaves(root))
    probes = metric_blocks(leaves)
    rec_cfg = (rs.measure or {}).get("record")
    recorder = None
    if rec_cfg is not None:
        missing = [k for k in ("slots", "decim") if k not in rec_cfg]
        if missing:
            raise KeyError(f"measure.record needs {missing} (D34: no defaults)")
        slots = list(rec_cfg["slots"])
        slots += [s for b in probes for p in b.pairs if (s := b.slot(p.port)) not in slots]
        recorder = Recorder(slots, rec_cfg["decim"], settings.chunk_samples)
        root = Seq([root, recorder])

    layout = SlotLayout(settings.f_sim_hz)
    block = build(root, layout, settings.seed)
    # offline stage after build (blocks' bands() need init, e.g. the block rate) and before anything
    # is written (double registrations etc. fail early)
    offline = build_offline_stage(rs.topology, leaves, offline_blocks, rs.measure)

    run_dir: Path | None = None
    if record_root is not None:
        from cranberry.runner.record import record_run

        if config_input is None and not isinstance(spec, RunSpec):
            config_input = config.read_input(spec)
        run_dir = record_run(rs, Path(record_root), config_input=config_input,
                             **({"invocation": dict(invocation)} if invocation else {}))

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
    metrics = None
    if probes:
        metrics = evaluate_metrics(
            leaves, state, t_s=t_s, series=series, topology=offline.topology, source_outputs=offline.source_outputs,
            loss_registry=offline.loss_registry, band_ledger=offline.band_ledger, measure=rs.measure)
    measurements: dict[str, dict[str, Any]] = {}
    if any((rs.measure or {}).get(h) for h in MEASUREMENT_HOOKS):
        named = {b.name: b for b in (*offline_blocks, *leaves)}
        f_rec = settings.f_sim_hz / recorder.rate_div if recorder is not None else None
        measurements = evaluate_measurements(rs.measure, named, metrics=metrics, loss_registry=offline.loss_registry,
                                             f_record_hz=f_rec, seed=settings.seed)
    result = RunResult(
        spec=rs,
        run_dir=run_dir,
        t_s=t_s,
        series=series,
        metrics=metrics,
        table_hashes=offline.table_hashes(),
        seeds=seeds,
        final_state=state if keep_state else None,
        mode="compiled" if compiled else "reference",
        offline=offline,
        measurements=measurements,
    )
    if run_dir is not None:
        _write_results(run_dir, result)
    return result
