"""Metric output containers and the metric blocks (conv §10, D12, D31, D46, D47, G13). Owner: WP-F4.

Metrics (conv §10):

* (i)  per port k at its BHD observer: δ_i,k(t) = wrap_π(φ_LO,k − θ_sq,k − Δ_set) → :class:`MetricI`.
* (ii) per port k ≥ 1 at the GBS inputs, pulses of the same GBS shot:
  Θ_k(t) = wrap_π(θ_sq,k − θ_sq,0 − Θ_k,set) → :class:`MetricII`, with the four components
  (:data:`METRIC_II_COMPONENTS`) obtained by reruns + linearity check (G13, ``runner.decompose``):

  1. ``locked_path_residual``     residual of the locked path noise
  2. ``common_points``            motion of the common points (probe split, reference delivery)
  3. ``laser_decorrelation_oob``  laser phase decorrelation φ_L(t − τ_k) − φ_L(t − τ_0) (conv §7,
                                  4 sin²(πfΔτ_k) S_φ, Δτ_k = τ_k − τ_0), in band AND out of band: the
                                  same-laser Liu lock does NOT remove the white-FM in-band part (D77-1 as
                                  corrected by D78-1). With a finite loop L(f) the in-band term is
                                  |1/(1 + L)|²·4 sin²(πfΔτ_k) S_φ + |H_ideal|² S_φ, and for white FM the
                                  locked/unlocked ratio R(f_u) ≥ 1 (≈ 1 + 3πf_uΔτ for short ports;
                                  1.02–1.05 at the future presets' f_u = 1 kHz, 1.2 at 10 kHz): the lock
                                  adds decorrelation. The linewidth requirement is therefore the open-loop
                                  form σ² = R·2π²h0Δτ (D31-3, D78). The old in-band term
                                  16 sin²(πf(T_a − T_b)) sin²(πfΔτ_k) S_φ (D51-5, T10 eq. (12)/(13)) is
                                  wrong (D62-5, D77-1); delay matching T_a − T_b has no use. The lock acts
                                  on drift, 1/f and random-walk FM only (D77-1; with an orthogonal PM probe
                                  + TEC offload a fraction |K_T|/(k_T + |K_T|) ≈ 4.5 % of the offloaded
                                  drift phase leaks into Θ_k, M6-04).
  4. ``squeezed_only``            segments only the squeezed light traverses (incl. demux phase)

* alongside (ii): wave-packet overlap between ports (:class:`Overlap`), stray light in both units of
  D46 (:class:`StrayLight`), η_k at GBS input and at BHD separately (:class:`Efficiency`), and each
  port's own σ (D46).

Port convention (G13): port-based configs always include reference port 0. Port-less configs
(Ref A–D) report metric (i) under port key 0 when they have a squeezed beam and a BHD, and
``metric_ii_status = "not_applicable"`` with an empty ``metric_ii``.

Metric blocks (WP-F4)
    ``metric_i`` / ``metric_ii`` are time-loop probes (T stage) placed after ``beam_assembler``. The
    squeezed beam's ``beam.<b>.phi`` is its squeezing angle at the observer (the source slot wired in
    ``beam_assembler`` holds θ_sq of the source, conv §4; segments rotate it by +ψ, T04 §4), the LO
    beam's is φ_LO. Per selected port k each block writes

    =================================  =====  ==========================================================
    ``metric.<block>.port<k>``         rad    δ (i) or Θ (ii), wrapped to (−π/2, π/2] (quadrature angles
                                              are defined mod π, T04 §4)
    ``metric.<block>.port<k>.acc``     1      running sums [n, Σδ, Σδ², Σcos 2δ, Σsin 2δ] at the block
                                              rate (full-rate χ = ⟨cos 2δ⟩ independent of recording)
    P ``<block>.port<k>.setpoint_rad``  rad   Δ_set (i) / Θ_set (ii), runtime parameter
    =================================  =====  ==========================================================

    The block value is the **time-loop part**: the static source rotation θ_off
    (``SourceOutput.project``, D55 rule 4) is applied in ``phase_to_quadrature`` and reported in
    ``MetricsResult.noise_levels``. ``metric_i`` reads φ_LO − θ_sq; ``metric_ii`` θ_sq,k − θ_sq,0, both
    beams evaluated at the same sample (the beam assembler gives each beam's phase on arrival, so
    pulses arriving together belong to the same GBS shot).

    Params (all required, D34; no toggles — a probe has no physics to switch)::

        topology     : mapping or Topology (injected from RunSpec.topology by the runner if absent)
        ports        : metric_i:  "all" | [k, ...] | {k: {squeezed: <beam>, lo: <beam>}}
                         "all" = topology ports with an LO and a BHD observer; the explicit mapping
                         serves port-less configs (Ref A–D: {0: {squeezed: sqz, lo: lo}}) and ports
                         whose BHD sees a different squeezed beam than the GBS
                       metric_ii: "all" | [k ≥ 1, ...] | {k: {squeezed: <beam>}} (the mapping must
                         include port 0); "all" with no topology ports selects nothing (G13: n/a)
        setpoint_rad : number (every port) or {k: number}
        t_settle_s   : optional (WP-I3 addition; absent = 0, the old behaviour): the running sums start at
                       t ≥ t_settle_s (runtime parameter), so lock the start-up transient (the initial state is deterministic, Δ(0) = 0; WP-I3/M2C-06)
                       is not part of the moments / noise levels nor of ``MetricI/II`` series (the ESA zero
                       span reads those); the δ slot itself, and ``RunResult.series``, cover every sample

    Port keys k of these mappings (and of ``measure.metrics`` / ``overlap`` / ``stray_light``) are
    ints or digit strings; a ``RunSpec`` (config tree) requires string keys, e.g. ``{"0": {...}}``.

Measurement stage (:func:`evaluate_metrics`, called by ``runner.run`` after the loop): fills
:class:`MetricsResult` from the accumulators (and the recorded series of the metric slots, which the
runner adds to ``measure.record`` when recording is on), the offline stage (``SourceOutput``,
``LossRegistry``, ``BandLedger``) and ``RunSpec.measure``::

    measure:
      metrics:                         # per metric kind; absent kind → no noise levels for it
        i:  {f_max_hz, weight, source, lo_mode, average, include_hf_phase_noise, dither_m_rad}
        ii: {...same keys...}
      overlap:     {k: <slot> | <number>}                                  # optional
      stray_light: {k: {clicks_per_gate: <slot>|<number>, photons_per_mode: <slot>|<number>}}

    f_max_hz     upper limit of the out-of-band integral (BandLedger.define_metric; conv §6: flat
                 source range, ~1/τ_p pulses, ~κ_c OPO)
    weight       sideband weighting w(f) (callable) or null (flat source, T04 §3.2)
    source       name of the offline block whose ``source_output`` gives V±, θ_off; or {k: name}
    lo_mode      LO mode for temporal sources (array on the source grid) or null; or {k: ...}; or the
                 string ``optimal`` (WP-I5 addition): the source's optimal LO mode for a temporal source
                 (``temporal_modes.optimal_lo_mode``), null for single-mode / cw_sideband bases — so one
                 config serves single-mode stand-ins and pulsed sources (LO-shape variants are loss
                 elements of ``lo_temporal_shape``)
    average      exact | small_angle
    include_hf_phase_noise   add the ledger σ²_HF (G05)
    dither_m_rad analytic dither peak amplitudes not simulated in the loop ([] for none)

Metric names (ledger, ``noise_levels``): ``i.port<k>``, ``ii.port<k>``. η: metric (i) uses the
squeezed beam's path with target ``bhd``, metric (ii) target ``gbs`` (D46, G06). Each port's own σ
(D46, ``sigma_rad``) is the standard deviation of its metric (i) (open: T04 Q2 asks whether D16 means
Var(Θ_k) or a per-port σ; both are reported).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from typing import TYPE_CHECKING, Any

import numba
import numpy as np

from cranberry.blocks.meas.phase_to_quadrature import PhaseStats, quadrature_for_metric
from cranberry.core.block import BaseBlock, Block, MissingParameterError, require
from cranberry.core.registry import register_block

if TYPE_CHECKING:
    from cranberry.offline.band_ledger import BandLedger
    from cranberry.offline.loss_registry import LossRegistry
    from cranberry.offline.source_output import SourceOutput
    from cranberry.path.topology import Topology

__all__ = [
    "METRIC_II_COMPONENTS",
    "METRIC_II_STATUS",
    "METRIC_KINDS",
    "METRIC_CONFIG_KEYS",
    "wrap_pi",
    "wrap_quadrature",
    "metric_name",
    "MetricI",
    "MetricIIComponents",
    "MetricII",
    "Overlap",
    "StrayLight",
    "Efficiency",
    "MetricsResult",
    "MetricProbe",
    "make_metric_i",
    "make_metric_ii",
    "metric_blocks",
    "define_metrics",
    "evaluate_metrics",
]

METRIC_II_COMPONENTS: tuple[str, ...] = (
    "locked_path_residual",
    "common_points",
    "laser_decorrelation_oob",
    "squeezed_only",
)
#: ``ok``: computed; ``not_applicable``: no ports (G13); ``decomposition_invalid``: saturation / lock
#: loss / linearity check failed; ``not_run``: decomposition not requested.
METRIC_II_STATUS: tuple[str, ...] = ("ok", "not_applicable", "decomposition_invalid", "not_run")
#: Block kinds of the metric probes and their short metric-kind keys (``measure.metrics``).
METRIC_KINDS: dict[str, str] = {"metric_i": "i", "metric_ii": "ii"}
#: Required keys of ``measure.metrics.<i|ii>`` (D34).
METRIC_CONFIG_KEYS: tuple[str, ...] = (
    "f_max_hz", "weight", "source", "lo_mode", "average", "include_hf_phase_noise", "dither_m_rad",
)


def wrap_pi(x):
    """Wrap angles to (−π, π] (conv §10 wrap_π)."""
    return np.pi - np.mod(np.pi - np.asarray(x, dtype=float), 2.0 * np.pi)


def wrap_quadrature(x):
    """Wrap quadrature angles (defined mod π) to (−π/2, π/2] (T04 §4); used by the metric blocks."""
    return 0.5 * np.pi - np.mod(0.5 * np.pi - np.asarray(x, dtype=float), np.pi)


def metric_name(kind: str, port: int) -> str:
    """Ledger / noise-level name of a metric: ``i.port<k>`` or ``ii.port<k>``."""
    short = METRIC_KINDS.get(kind, kind)
    if short not in ("i", "ii"):
        raise ValueError(f"unknown metric kind {kind!r}")
    return f"{short}.port{int(port)}"


def _stats_dict(stats: PhaseStats | None) -> dict[str, Any]:
    if stats is None:
        return {}
    return {"n": int(stats.n), "mean_rad": float(stats.mean), "std_rad": float(stats.std),
            "rms_rad": float(stats.rms), "mean_cos2": float(stats.mean_cos2), "mean_sin2": float(stats.mean_sin2)}


@dataclass
class MetricI:
    """Metric (i) of port k: δ_i,k(t) [rad] at its BHD observer (time axis = ``MetricsResult.t_s``).

    ``delta_rad`` is the recorded (decimated) series of the time-loop part (empty when not recorded);
    ``stats`` the full-rate moments from the block's accumulators (WP-F4 addition); ``beams`` the
    signed beam combination {lo: +1, squeezed: −1} (WP-F4 addition, ledger / loss paths).
    """

    port: int
    observer: str
    delta_rad: np.ndarray
    setpoint_rad: float
    stats: PhaseStats | None = None
    beams: Mapping[str, float] = field(default_factory=dict)

    @property
    def rms_rad(self) -> float:
        if self.stats is not None:
            return float(self.stats.rms)
        return float(np.sqrt(np.mean(np.square(self.delta_rad))))

    @property
    def std_rad(self) -> float:
        return float(self.stats.std) if self.stats is not None else float(np.std(self.delta_rad))


@dataclass
class MetricIIComponents:
    """Variance split of Θ_k [rad²] (G13).

    Attributes:
        variances_rad2: {component name (METRIC_II_COMPONENTS): variance}, time loop + the component's
            ledger rows (Phase 5 M6-06, ``runner.decompose``).
        total_var_rad2: variance of the all-on run (time loop + ledger).
        linearity_ratio: Σ components / total (valid if within 5 %, T04 V12).
        status: "ok" or "decomposition_invalid".
        reason: why invalid (saturation, lock loss, unstable loop, linearity), or None.
    """

    variances_rad2: Mapping[str, float]
    total_var_rad2: float
    linearity_ratio: float
    status: str
    reason: str | None = None


@dataclass
class MetricII:
    """Metric (ii) of port k ≥ 1: Θ_k(t) [rad] at the GBS input, plus its decomposition.

    ``stats`` / ``beams`` ({squeezed_k: +1, squeezed_0: −1}): WP-F4 additions as in :class:`MetricI`.
    """

    port: int
    observer: str
    theta_rad: np.ndarray
    setpoint_rad: float
    components: MetricIIComponents | None = None
    stats: PhaseStats | None = None
    beams: Mapping[str, float] = field(default_factory=dict)

    @property
    def rms_rad(self) -> float:
        if self.stats is not None:
            return float(self.stats.rms)
        return float(np.sqrt(np.mean(np.square(self.theta_rad))))

    @property
    def std_rad(self) -> float:
        return float(self.stats.std) if self.stats is not None else float(np.std(self.theta_rad))


@dataclass
class Overlap:
    """Wave-packet overlap of port k with port 0 at the GBS input (∈ [0, 1]; series or constant)."""

    port: int
    overlap: np.ndarray


@dataclass
class StrayLight:
    """Stray (probe-derived) light in port k (D46): clicks per detection gate (both polarizations) and
    photons in the squeezed mode (N of conv §2)."""

    port: int
    clicks_per_gate: float
    photons_per_mode: float


@dataclass
class Efficiency:
    """Total transmission of port k's squeezed light to the GBS input and to its BHD (D46, G06)."""

    port: int
    eta_gbs: float | None
    eta_bhd: float | None


def _jsonable_levels(levels: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in levels.items():
        if k == "v_prime_t":
            continue  # series: results.npz
        out[k] = v.tolist() if isinstance(v, np.ndarray) else v
    return out


@dataclass
class MetricsResult:
    """All metrics of one run.

    ``noise_levels`` (WP-F4 addition): {metric name: ``phase_to_quadrature`` output} — V′, R′∓ [dB re
    shot noise], χ and σ per band, η, θ_off (``quadrature_for_metric``).
    """

    t_s: np.ndarray
    metric_i: dict[int, MetricI] = field(default_factory=dict)
    metric_ii: dict[int, MetricII] = field(default_factory=dict)
    metric_ii_status: str = "not_run"
    sigma_rad: dict[int, float] = field(default_factory=dict)
    overlap: dict[int, Overlap] = field(default_factory=dict)
    stray_light: dict[int, StrayLight] = field(default_factory=dict)
    efficiency: dict[int, Efficiency] = field(default_factory=dict)
    noise_levels: dict[str, dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.metric_ii_status not in METRIC_II_STATUS:
            raise ValueError(f"metric_ii_status {self.metric_ii_status!r} not in {METRIC_II_STATUS}")

    def to_json(self) -> dict[str, Any]:
        """Summary for ``metrics.json`` (series go to results.npz)."""
        def m_json(m: MetricI | MetricII) -> dict[str, Any]:
            d = {"port": m.port, "observer": m.observer, "setpoint_rad": float(m.setpoint_rad),
                 "beams": dict(m.beams), "n_recorded": int(np.size(m.delta_rad if isinstance(m, MetricI) else m.theta_rad)),
                 **_stats_dict(m.stats)}
            if isinstance(m, MetricII) and m.components is not None:
                d["components"] = asdict(m.components)
            return d

        return {
            "n_t": int(np.size(self.t_s)),
            "metric_i": {str(k): m_json(m) for k, m in self.metric_i.items()},
            "metric_ii": {str(k): m_json(m) for k, m in self.metric_ii.items()},
            "metric_ii_status": self.metric_ii_status,
            "sigma_rad": {str(k): float(v) for k, v in self.sigma_rad.items()},
            "overlap": {str(k): {"mean": float(np.mean(o.overlap)), "min": float(np.min(o.overlap)),
                                 "n": int(np.size(o.overlap))} for k, o in self.overlap.items()},
            "stray_light": {str(k): asdict(s) for k, s in self.stray_light.items()},
            "efficiency": {str(k): asdict(e) for k, e in self.efficiency.items()},
            "noise_levels": {n: _jsonable_levels(v) for n, v in self.noise_levels.items()},
        }

    def summary(self) -> dict[str, float]:
        """Flat {name: scalar} for ablation / convergence tables (``i.port0.std_rad``,
        ``i.port0.r_minus_db``, ``ii.port1.rms_rad``, ``eff.port1.eta_gbs``, ...)."""
        out: dict[str, float] = {}
        for short, group in (("i", self.metric_i), ("ii", self.metric_ii)):
            for k, m in group.items():
                if m.stats is not None:
                    for key, v in _stats_dict(m.stats).items():
                        out[f"{short}.port{k}.{key}"] = float(v)
        for name, lv in self.noise_levels.items():
            for key, v in lv.items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    out[f"{name}.{key}"] = float(v)
        for k, e in self.efficiency.items():
            for key in ("eta_gbs", "eta_bhd"):
                if getattr(e, key) is not None:
                    out[f"eff.port{k}.{key}"] = float(getattr(e, key))
        for k, o in self.overlap.items():
            out[f"overlap.port{k}.mean"] = float(np.mean(o.overlap))
        return out


# ------------------------------------------------------------------------------------------------
# metric probe blocks


@dataclass(frozen=True)
class _Pair:
    port: int
    plus: str  # beam id entering with +1
    minus: str  # beam id entering with −1
    observer: str


def _port_key(k: Any, where: str) -> int:
    if isinstance(k, bool):
        raise ValueError(f"{where}: port key must be an int, got {k!r}")
    if isinstance(k, (int, np.integer)):
        return int(k)
    if isinstance(k, str) and k.strip().lstrip("-").isdigit():
        return int(k)
    raise ValueError(f"{where}: port key must be an int, got {k!r}")


class MetricProbe(BaseBlock):
    """Common implementation of ``metric_i`` / ``metric_ii`` (module docstring)."""

    toggle_names = ()

    def __init__(self, kind: str, params: Mapping[str, Any], name: str | None):
        from cranberry.path.topology import as_topology

        name = name or kind
        require(params, "topology", "ports", "setpoint_rad", block=name)
        self.topo = as_topology(params["topology"])
        pairs = self._pairs_metric_i(params["ports"], name) if kind == "metric_i" else self._pairs_metric_ii(
            params["ports"], name)
        sp = params["setpoint_rad"]
        self.setpoints: dict[int, float] = {}
        for p in pairs:
            if isinstance(sp, Mapping):
                hit = [v for kk, v in sp.items() if _port_key(kk, f"block {name!r}: setpoint_rad") == p.port]
                if not hit:
                    raise MissingParameterError(f"block {name!r}: setpoint_rad has no entry for port {p.port} (D34)")
                self.setpoints[p.port] = float(hit[0])
            else:
                self.setpoints[p.port] = float(sp)
        self.pairs = pairs
        self.has_settle = "t_settle_s" in params  # absent: no settle param, no read of t (pre-WP-I3 block)
        self.t_settle_s = float(params.get("t_settle_s", 0.0))
        if not self.t_settle_s >= 0.0:
            raise ValueError(f"block {name!r}: t_settle_s must be >= 0")
        self.short = METRIC_KINDS[kind]
        self.name = name  # slot names are needed before BaseBlock.__init__ (which sets it again)
        reads = [*(["t"] if self.has_settle else []),
                 *dict.fromkeys(f"beam.{b}.phi" for p in pairs for b in (p.plus, p.minus))]
        writes = [w for p in pairs for w in (self.slot(p.port), self.acc_slot(p.port))]
        super().__init__(kind, name, params, reads=reads, writes=writes)
        self.slot_units = {s: "rad" for s in reads if s != "t"}
        self.slot_units.update({self.slot(p.port): "rad" for p in pairs})

    # names --------------------------------------------------------------------------------------

    def slot(self, port: int) -> str:
        """δ / Θ slot of ``port``."""
        return f"metric.{self.name}.port{port}"

    def acc_slot(self, port: int) -> str:
        """Accumulator slot of ``port`` ([n, Σδ, Σδ², Σcos 2δ, Σsin 2δ])."""
        return f"metric.{self.name}.port{port}.acc"

    def setpoint_param(self, port: int) -> str:
        return f"{self.name}.port{port}.setpoint_rad"

    @property
    def settle_param(self) -> str:
        return f"{self.name}.t_settle_s"

    def beams_of(self, pair: _Pair) -> dict[str, float]:
        """Signed beam combination of the metric (ledger / losses)."""
        return {pair.plus: 1.0, pair.minus: -1.0}

    def squeezed_beam(self, pair: _Pair) -> str:
        """The squeezed beam whose losses set η of this metric."""
        return pair.minus if self.short == "i" else pair.plus

    # pair selection ------------------------------------------------------------------------------

    def _check_beams(self, name: str, k: int, beams: Mapping[str, str | None], observer: str | None) -> str:
        topo = self.topo
        for role, b in beams.items():
            if b is None or b not in topo.beams:
                raise ValueError(f"block {name!r}: port {k}: {role} beam {b!r} not in the topology")
        obs = {role: topo.beams[b].observer for role, b in beams.items()}
        if observer is None:
            observer = next(iter(obs.values()))
        off = {r: o for r, o in obs.items() if o != observer}
        if off:
            raise ValueError(f"block {name!r}: port {k}: beams must be evaluated at observer {observer!r}, but "
                             f"{off} (give the beams explicitly: ports: {{{k}: {{...}}}})")
        return observer

    def _pairs_metric_i(self, spec: Any, name: str) -> list[_Pair]:
        topo = self.topo
        out: list[_Pair] = []
        if isinstance(spec, Mapping):
            for kk, v in spec.items():
                k = _port_key(kk, f"block {name!r}: ports")
                if not isinstance(v, Mapping):
                    raise ValueError(f"block {name!r}: ports[{k}] must be {{squeezed: <beam>, lo: <beam>}}")
                require(v, "squeezed", "lo", block=f"{name}.ports[{k}]")
                obs = self._check_beams(name, k, {"squeezed": v["squeezed"], "lo": v["lo"]}, None)
                out.append(_Pair(k, str(v["lo"]), str(v["squeezed"]), obs))
        else:
            if spec == "all":
                ks = [k for k, p in topo.ports.items() if p.lo is not None and p.bhd_observer is not None]
            elif isinstance(spec, (list, tuple)):
                ks = [_port_key(k, f"block {name!r}: ports") for k in spec]
            else:
                raise ValueError(f"block {name!r}: ports must be 'all', a list of port indices or a mapping")
            for k in ks:
                p = topo.ports.get(k)
                if p is None or p.lo is None or p.bhd_observer is None:
                    raise ValueError(f"block {name!r}: port {k} needs a topology port with an LO and a BHD observer")
                obs = self._check_beams(name, k, {"squeezed": p.squeezed, "lo": p.lo}, p.bhd_observer)
                out.append(_Pair(k, p.lo, p.squeezed, obs))
        if not out:
            raise ValueError(f"block {name!r}: no port selected for metric (i)")
        if len({p.port for p in out}) != len(out):
            raise ValueError(f"block {name!r}: a port is selected twice")
        return sorted(out, key=lambda p: p.port)

    def _pairs_metric_ii(self, spec: Any, name: str) -> list[_Pair]:
        topo = self.topo
        if isinstance(spec, Mapping):
            beams = {}
            for kk, v in spec.items():
                k = _port_key(kk, f"block {name!r}: ports")
                if not isinstance(v, Mapping):
                    raise ValueError(f"block {name!r}: ports[{k}] must be {{squeezed: <beam>}}")
                require(v, "squeezed", block=f"{name}.ports[{k}]")
                beams[k] = str(v["squeezed"])
            if 0 not in beams:
                raise ValueError(f"block {name!r}: the reference port 0 must be in ports (G13)")
            out = []
            for k in sorted(b for b in beams if b != 0):
                obs = self._check_beams(name, k, {"squeezed": beams[k]}, None)
                self._check_beams(name, 0, {"squeezed": beams[0]}, None)
                out.append(_Pair(k, beams[k], beams[0], obs))
            return out
        if spec == "all":
            ks = sorted(k for k in topo.ports if k != 0)
        elif isinstance(spec, (list, tuple)):
            ks = [_port_key(k, f"block {name!r}: ports") for k in spec]
        else:
            raise ValueError(f"block {name!r}: ports must be 'all', a list of port indices or a mapping")
        if ks and 0 not in topo.ports:
            raise ValueError(f"block {name!r}: metric (ii) needs the reference port 0 (G13)")
        out = []
        for k in ks:
            if k == 0 or k not in topo.ports:
                raise ValueError(f"block {name!r}: port {k} is not a non-reference topology port")
            p, p0 = topo.ports[k], topo.ports[0]
            for q in (p, p0):
                if q.gbs_observer is None:
                    raise ValueError(f"block {name!r}: port {q.index} has no GBS observer")
                self._check_beams(name, q.index, {"squeezed": q.squeezed}, q.gbs_observer)
            out.append(_Pair(k, p.squeezed, p0.squeezed, p.gbs_observer))
        return out

    # lifecycle -----------------------------------------------------------------------------------

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        for p in self.pairs:
            layout.add(self.slot(p.port), (), np.float64, "rad")
            layout.add(self.acc_slot(p.port), (5,), np.float64, "1")
            layout.add_param(self.setpoint_param(p.port), self.setpoints[p.port], (), "rad")
        if self.has_settle:
            layout.add_param(self.settle_param, self.t_settle_s, (), "s")

    def bind(self, layout) -> None:
        cols = [(layout.slot(f"beam.{p.plus}.phi").offset, layout.slot(f"beam.{p.minus}.phi").offset,
                 layout.param(self.setpoint_param(p.port)).offset, layout.slot(self.slot(p.port)).offset,
                 layout.slot(self.acc_slot(p.port)).offset) for p in self.pairs]
        arr = np.array(cols if cols else [[0, 0, 0, 0, 0]], dtype=np.int64)
        PLUS, MINUS, SP, OUT, ACC = (arr[:, j].copy() for j in range(5))
        self._cols = cols
        n = len(self.pairs)
        half_pi, pi = 0.5 * math.pi, math.pi
        T, TS = (layout.slot("t").offset, layout.param(self.settle_param).offset) if self.has_settle else (-1, -1)
        self._t_cols = (T, TS)

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            for j in range(n):
                x = S[PLUS[j]] - S[MINUS[j]] - P[SP[j]]
                d = half_pi - ((half_pi - x) % pi)
                S[OUT[j]] = d
                if T >= 0 and S[T] < P[TS]:
                    continue
                a = ACC[j]
                S[a] += 1.0
                S[a + 1] += d
                S[a + 2] += d * d
                S[a + 3] += math.cos(2.0 * d)
                S[a + 4] += math.sin(2.0 * d)

        self.kernel = kernel

    def __call__(self, state, params):
        S, P = state.S, params.P
        for plus, minus, sp, out, a in self._cols:
            x = S[plus] - S[minus] - P[sp]
            d = 0.5 * math.pi - ((0.5 * math.pi - x) % math.pi)
            S[out] = d
            if self._t_cols[0] >= 0 and S[self._t_cols[0]] < P[self._t_cols[1]]:
                continue
            S[a] += 1.0
            S[a + 1] += d
            S[a + 2] += d * d
            S[a + 3] += math.cos(2.0 * d)
            S[a + 4] += math.sin(2.0 * d)
        return state


@register_block("metric_i")
def make_metric_i(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the metric (i) probe block (module docstring)."""
    return MetricProbe("metric_i", params, name)


@register_block("metric_ii")
def make_metric_ii(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the metric (ii) probe block (module docstring)."""
    return MetricProbe("metric_ii", params, name)


# ------------------------------------------------------------------------------------------------
# measurement stage (runner)


def metric_blocks(leaves) -> list[MetricProbe]:
    """The metric probes among ``leaves`` (execution order); a port measured twice is an error."""
    blocks = [b for b in leaves if isinstance(b, MetricProbe)]
    seen: dict[str, str] = {}
    for b in blocks:
        for p in b.pairs:
            n = metric_name(b.kind, p.port)
            if n in seen:
                raise ValueError(f"metric {n!r} computed by both {seen[n]!r} and {b.name!r}")
            seen[n] = b.name
    return blocks


def _metric_config(measure: Mapping[str, Any], short: str) -> Mapping[str, Any] | None:
    cfg = (measure or {}).get("metrics")
    if cfg is None:
        return None
    if not isinstance(cfg, Mapping):
        raise ValueError("measure.metrics must be a mapping {i: {...}, ii: {...}}")
    unknown = set(cfg) - {"i", "ii"}
    if unknown:
        raise ValueError(f"measure.metrics: unknown metric kind(s) {sorted(unknown)} (known: i, ii)")
    c = cfg.get(short)
    if c is None:
        return None
    require(c, *METRIC_CONFIG_KEYS, block=f"measure.metrics.{short}")
    unknown = set(c) - set(METRIC_CONFIG_KEYS)
    if unknown:
        raise ValueError(f"measure.metrics.{short}: unknown key(s) {sorted(unknown)}")
    return c


def _per_port(value: Any, port: int, what: str) -> Any:
    if isinstance(value, Mapping):
        for k, v in value.items():
            if _port_key(k, what) == port:
                return v
        raise MissingParameterError(f"{what}: no entry for port {port} (D34)")
    return value


def _lo_mode(value: Any, source: SourceOutput) -> Any:
    """Resolve ``lo_mode: optimal`` (module docstring) against ``source``; other values unchanged."""
    if isinstance(value, str):
        if value != "optimal":
            raise ValueError(f"measure.metrics.*.lo_mode: unknown string {value!r} (array, null or 'optimal')")
        if source.basis != "temporal":
            return None
        from cranberry.blocks.source.temporal_modes import optimal_lo_mode

        return optimal_lo_mode(source)[0]
    return value


def define_metrics(ledger: BandLedger, blocks: list[MetricProbe], measure: Mapping[str, Any]) -> list[str]:
    """``ledger.define_metric`` for every (configured) metric of ``blocks`` with its ``f_max_hz`` /
    ``weight`` from ``measure.metrics`` (offline stage). Returns the defined names."""
    names = []
    for b in blocks:
        cfg = _metric_config(measure, b.short)
        if cfg is None:
            continue
        w = cfg["weight"]
        if w is not None and not callable(w):
            raise ValueError(f"measure.metrics.{b.short}.weight must be null or a callable w(f)")
        for p in b.pairs:
            n = metric_name(b.kind, p.port)
            ledger.define_metric(n, f_max_hz=float(_per_port(cfg["f_max_hz"], p.port, f"measure.metrics.{b.short}.f_max_hz")),
                                 weight=w)
            names.append(n)
    return names


def _values(slot_or_value: Any, state, series: Mapping[str, np.ndarray], what: str) -> np.ndarray:
    if isinstance(slot_or_value, str):
        if slot_or_value in series:
            return np.asarray(series[slot_or_value], dtype=float).ravel()
        if state is not None and slot_or_value in state.layout:
            return np.atleast_1d(np.asarray(state[slot_or_value], dtype=float)).ravel()
        raise KeyError(f"{what}: slot {slot_or_value!r} is neither recorded nor in the state")
    return np.atleast_1d(np.asarray(slot_or_value, dtype=float))


def evaluate_metrics(
    leaves,
    state,
    *,
    t_s: np.ndarray,
    series: Mapping[str, np.ndarray],
    topology: Topology | None,
    source_outputs: Mapping[str, SourceOutput],
    loss_registry: LossRegistry | None,
    band_ledger: BandLedger | None,
    measure: Mapping[str, Any],
) -> MetricsResult:
    """Measurement stage: metric (i)/(ii) series and moments, noise levels via ``phase_to_quadrature``,
    η per port, overlap and stray light (module docstring)."""
    blocks = metric_blocks(leaves)
    res = MetricsResult(t_s=np.asarray(t_s))
    sqz_bhd: dict[int, str] = {}
    sqz_gbs: dict[int, str] = {}
    for b in blocks:
        cfg = _metric_config(measure, b.short)
        for p in b.pairs:
            stats = PhaseStats.from_acc(state[b.acc_slot(p.port)])
            ser = np.asarray(series.get(b.slot(p.port), np.zeros(0)), dtype=float)
            if b.t_settle_s > 0.0 and ser.size:  # the measured series starts with the moments (t_settle_s)
                ser = ser[np.asarray(t_s, dtype=float)[: ser.size] >= b.t_settle_s]
            beams = b.beams_of(p)
            sp = b.setpoints[p.port]
            if b.short == "i":
                res.metric_i[p.port] = MetricI(p.port, p.observer, ser, sp, stats, beams)
                res.sigma_rad[p.port] = float(stats.std)
                sqz_bhd[p.port] = p.minus
            else:
                res.metric_ii[p.port] = MetricII(p.port, p.observer, ser, sp, None, stats, beams)
                sqz_gbs[p.port], sqz_gbs[0] = p.plus, p.minus
            if cfg is None:
                continue
            name = metric_name(b.kind, p.port)
            what = f"measure.metrics.{b.short}"
            src_name = _per_port(cfg["source"], p.port, f"{what}.source")
            if src_name not in source_outputs:
                raise KeyError(f"{what}.source: no offline source block {src_name!r} (known: {sorted(source_outputs)})")
            lo_mode = _lo_mode(_per_port(cfg["lo_mode"], p.port, f"{what}.lo_mode"), source_outputs[src_name])
            ref_off = 0.0
            if b.short == "ii":
                ref_src = source_outputs[_per_port(cfg["source"], 0, f"{what}.source")]
                ref_off = ref_src.project(_lo_mode(_per_port(cfg["lo_mode"], 0, f"{what}.lo_mode"), ref_src))[2]
            path = b.squeezed_beam(p)
            has_path = loss_registry is not None and path in loss_registry.paths
            res.noise_levels[name] = quadrature_for_metric(
                stats,  # full-rate moments (the decimated series would alias χ only statistically)
                metric=name, source=source_outputs[src_name], lo_mode=lo_mode,
                loss_registry=loss_registry if has_path else None, path=path if has_path else None,
                target="bhd" if b.short == "i" else "gbs", band_ledger=band_ledger, beams=beams,
                average=str(cfg["average"]), setpoint_rad=sp if b.short == "i" else 0.0,
                include_hf_phase_noise=bool(cfg["include_hf_phase_noise"]),
                dither_m_rad=list(cfg["dither_m_rad"] or []), theta_off_ref_rad=ref_off,
            )
    # η per port (D46): squeezed beam of the BHD metric / of the GBS metric (topology port otherwise)
    ports: set[int] = set(sqz_bhd) | set(sqz_gbs)
    if topology is not None:
        ports |= set(topology.ports)
    for k in sorted(ports):
        tp = topology.ports.get(k) if topology is not None else None

        def eta(beam: str | None, target: str) -> float | None:
            if beam is None or loss_registry is None or beam not in loss_registry.paths:
                return None
            return loss_registry.eta(beam, target)

        b_bhd = sqz_bhd.get(k, tp.squeezed if tp is not None and tp.bhd_observer is not None else None)
        b_gbs = sqz_gbs.get(k, tp.squeezed if tp is not None and tp.gbs_observer is not None else None)
        e = Efficiency(k, eta(b_gbs, "gbs"), eta(b_bhd, "bhd"))
        if e.eta_gbs is not None or e.eta_bhd is not None:
            res.efficiency[k] = e
    res.metric_ii_status = "not_run" if topology is not None and topology.has_ports else "not_applicable"
    for kk, v in ((measure or {}).get("overlap") or {}).items():
        k = _port_key(kk, "measure.overlap")
        res.overlap[k] = Overlap(k, _values(v, state, series, f"measure.overlap[{k}]"))
    for kk, v in ((measure or {}).get("stray_light") or {}).items():
        k = _port_key(kk, "measure.stray_light")
        require(v, "clicks_per_gate", "photons_per_mode", block=f"measure.stray_light[{k}]")
        res.stray_light[k] = StrayLight(
            k, float(np.mean(_values(v["clicks_per_gate"], state, series, f"measure.stray_light[{k}]"))),
            float(np.mean(_values(v["photons_per_mode"], state, series, f"measure.stray_light[{k}]"))))
    return res
