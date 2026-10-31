"""Static ``Topology``: segments, beams, coefficients, locks, ports, observers (protected list "per-path
phase model"; PLAN §2.3.1, conv §9, D11, D47). Owner: WP-F6.

The Topology is fixed at build time; the dynamic quantities live in State slots (``drv.*``, ``psi.*``,
``tau_g.*``, ``J.*``, ``beam.*``; ``core.layout``). Cancellation on shared segments is computed, never
assumed (D11): beam phases are sums over their paths evaluated at the times each beam passes each
subsegment (``path.beams``), and segment × beam coefficients may differ (VT-23).

Spec form (``RunSpec.topology``; every key is required unless marked optional, D34)::

    segments:
      - id: arm1                       # slot names drv.arm1.*, psi.arm1.<beam>, tau_g.arm1.<beam>
        kind: pm                       # SEGMENT_KINDS
        L_m: 8.4
        n_g: 1.468                     # group index (reference / slow axis for PM)
        beta2_s2_per_m: -2.2e-26
        beta3_s3_per_m: 1.3e-40
        zone: lab                      # thermal zone id
        n_sub: 20                      # subsegments N_sub (counter-propagation timing, T08 §2.7, T10)
        spatial: lumped                # lumped | distributed | mixed (checked against `drivers`)
        pol_model: pm_diagonal         # jones | pm_diagonal | off (G04)
        pol_params: {...}              # see below ({} for off)
        actuator_sub: 0                # optional: subsegment of an actuator (x_a), null = none
        drivers:                       # driver -> placement (weights per subsegment, path.coef)
          dT: uniform                  #   scalar slot, weight 1/N on every subsegment
          strain: per_sub              #   slot of shape (N_sub,), weight 1/N each
          stretcher: actuator          #   scalar slot, weight 1 at actuator_sub
          eom: {point: 3}              #   scalar slot, weight 1 at subsegment 3
          tec: {subs: [0, 1, 2]}       #   scalar slot, weight 1/N on the listed subsegments (G08)
    beams:
      - id: pr
        source: dx1                    # source id (its phase slot is wired in beam_assembler params)
        path:                          # ordered, in the beam's own propagation order
          - {segment: arm1, direction: 1, axis: slow}   # sub_in / sub_out optional (whole segment);
                                                         # axis required on pm_diagonal segments only
        freq_offset_hz: 0.0            # conv §4 sign
        jones_in: [1, 0]               # complex entries as numbers or [re, im] pairs; unit norm
        jones_out: [1, 0]              # analyzer at the observer; required iff a jones segment is on
                                       # the path (optional otherwise)
        power_w: 1.0e-3
        observer: det
    observers:
      - {id: det, kind: detector}      # segment / sub optional (then they must match the beam ends)
    locks:
      - {id: lockD, combo: {lo: 1, pr: -1}, setpoint_rad: 1.5707963, observer: det}  # or `setpoint`
    ports: []                          # [{index, squeezed, lo, probe, gbs_observer, bhd_observer}]
    coef:                              # list of entries (or, from Python, {CoefKey: value})
      - {segment: arm1, driver: dT, beam: "*", effect: phase, value: 369.6, coords: {}}

Rules checked by :meth:`Topology.from_spec`:

* Every (segment, declared driver, beam traversing the segment, effect) needs a coefficient (explicit 0
  is fine): a forgotten probe coefficient would otherwise fake a non-cancellation (D11, D34). ``beam:
  "*"`` expands to every beam traversing the segment; an entry naming the beam overrides it.
* ``coords`` (optional) are the coordinates a callable / Table coefficient is evaluated at (G16-a).
* A beam traverses a segment at most once (repeated traversals: split the segment); direction d = ±1;
  d = +1 runs sub_in ≤ sub_out, d = −1 runs sub_in ≥ sub_out.
* ``spatial: lumped`` forbids ``per_sub`` placements; ``distributed`` requires ``per_sub`` for the
  environmental drivers ``dT`` / ``strain`` present; ``mixed`` allows both (the noise block mixes
  common and independent parts, T08 §2.7).
* pm_diagonal ``pol_params`` (T09 §1.2, conv §9): ``K_T_rad_per_m_K`` (∂Δψ_b/∂T per length, negative
  nominal −2.05), ``K_tau_s_per_m_K`` (∂τ_DGD/∂T per length, Ritari −1.78e−15), ``b_g_s_per_m`` (DGD per
  length), ``axis_drivers`` (drivers that act as temperature on the birefringence, normally
  ``[dT, tec]``; the ``phase_direct`` preset of G08 leaves ``tec`` out). The segment's coefficients
  are the slow-axis values; ``segment_phase`` adds the fast-axis terms (toggle ``pm_axis_drift``).
* jones ``pol_params``: exactly one of ``J`` (static 2×2 polarization matrix U) or ``J_slot`` (slot
  ``J.<seg>`` written by ``segment_jones``, WP-B7). The full segment Jones matrix is e^{iψ}·U with ψ
  = ``psi.<seg>.<beam>``; U must not contain the common (driver × coefficient) phase. Optional
  ``J_nr_slot`` (WP-B7 addition): slot of a non-reciprocal factor N (``faraday_nonreciprocal``) at the
  segment's d = +1 exit end. ``beam_assembler`` applies N·U for d = +1 and Uᵀ·N for d = −1, i.e. N is
  NOT transposed for the backward beam (J_back = Jᵀ holds for the reciprocal part only, T09 §1.8).
* Locks: every combo beam has the lock's observer. Ports: port 0 exists when ports exist (G13).

Nominal timing (:meth:`Topology.timing`): subsegment m of segment i has length L/N and group delay
δ = n_g L/(c N); a beam passes the centre of its j-th traversed subsegment (j = 0 at entry) (n − j − ½)δ
before leaving the segment (n traversed subsegments). τ_after(step) is the nominal group delay from the
segment exit to the observer; τ_b is source → observer (conv §7 emission time).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.path.coef import DRIVERS, EFFECTS, CoefKey, CoefValue, evaluate_coef

if TYPE_CHECKING:
    from cranberry.offline.band_ledger import BandLedger
    from cranberry.offline.loss_registry import BeamPath, LossRegistry

__all__ = [
    "C_M_PER_S",
    "SEGMENT_KINDS",
    "SPATIAL_STRUCTURES",
    "POL_MODELS",
    "PLACEMENTS",
    "AXES",
    "OBSERVER_KINDS",
    "TopologyError",
    "Segment",
    "PathStep",
    "Beam",
    "Lock",
    "Port",
    "Observer",
    "StepTiming",
    "BeamTiming",
    "Topology",
    "as_topology",
]

#: Speed of light in vacuum [m/s].
C_M_PER_S: float = 299_792_458.0

SEGMENT_KINDS: tuple[str, ...] = ("pm", "smf", "free_space", "opo", "crystal", "demux", "component")
#: Spatial structure of in-segment noise (conv §9 一様 / 独立 / 混合 = lumped / distributed / mixed).
SPATIAL_STRUCTURES: tuple[str, ...] = ("lumped", "distributed", "mixed")
#: Polarization model per segment (G04): full Jones, PM diagonal fast path, or off (J = e^{iψ}·1).
POL_MODELS: tuple[str, ...] = ("jones", "pm_diagonal", "off")
#: Driver placements (normalised form: ("uniform",), ("per_sub",), ("point", k), ("subs", (k, ...))).
PLACEMENTS: tuple[str, ...] = ("uniform", "per_sub", "point", "subs")
#: PM axis labels (pm_diagonal fast path).
AXES: tuple[str, ...] = ("slow", "fast")
OBSERVER_KINDS: tuple[str, ...] = ("detector", "gbs_input", "interaction")
_ENV_DRIVERS: tuple[str, ...] = ("dT", "strain")
_PM_KEYS: tuple[str, ...] = ("K_T_rad_per_m_K", "K_tau_s_per_m_K", "b_g_s_per_m", "axis_drivers")


class TopologyError(ValueError):
    """Invalid topology spec (all problems found are listed)."""


@dataclass(frozen=True)
class Segment:
    """One path segment (conv §9).

    Attributes:
        id: segment id (used in slot names ``drv.<id>.*``, ``psi.<id>.<beam>``).
        kind: one of :data:`SEGMENT_KINDS`.
        length_m: length L [m].
        n_g: group index.
        beta2_s2_per_m, beta3_s3_per_m: dispersion [s²/m], [s³/m].
        zone: thermal zone id (segments in one zone share ambient temperature).
        n_sub: number of subsegments N_sub (timing of counter-propagation, T08/T10).
        spatial: one of :data:`SPATIAL_STRUCTURES`.
        pol_model: one of :data:`POL_MODELS`; ``pol_params`` holds its parameters.
        actuator_sub: subsegment index x_a of an actuator on this segment, or None.
        drivers: {driver: normalised placement} (module docstring; WP-F6 addition).
    """

    id: str
    kind: str
    length_m: float
    n_g: float
    beta2_s2_per_m: float
    beta3_s3_per_m: float
    zone: str
    n_sub: int
    spatial: str
    pol_model: str
    pol_params: Mapping[str, Any] = field(default_factory=dict)
    actuator_sub: int | None = None
    drivers: Mapping[str, tuple] = field(default_factory=dict)

    @property
    def tau_g0_s(self) -> float:
        """Nominal group delay n_g·L/c [s]."""
        return self.n_g * self.length_m / C_M_PER_S

    @property
    def sub_tau_s(self) -> float:
        """Nominal group delay of one subsegment δ = τ_g0 / N_sub [s]."""
        return self.tau_g0_s / self.n_sub

    def weights(self, driver: str) -> np.ndarray | None:
        """Per-subsegment weights of a scalar-slot driver (shape (N_sub,)), or None for ``per_sub``."""
        pl = self.drivers[driver]
        n = self.n_sub
        w = np.zeros(n)
        if pl[0] == "uniform":
            w[:] = 1.0 / n
        elif pl[0] == "point":
            w[pl[1]] = 1.0
        elif pl[0] == "subs":
            w[list(pl[1])] = 1.0 / n
        else:
            return None
        return w


@dataclass(frozen=True)
class PathStep:
    """One traversal of a segment by a beam: direction d = +1 (source → detector/GBS) or −1, entering
    and leaving subsegment indices (conv §9). A counter-propagating beam samples the same subsegments
    at different times; it sees Jᵀ (reciprocity, G04). ``axis`` is the PM axis label on pm_diagonal
    segments (WP-F6 addition), else None."""

    segment: str
    direction: int
    sub_in: int
    sub_out: int
    axis: str | None = None

    @property
    def subs(self) -> tuple[int, ...]:
        """Traversed subsegment indices in traversal order."""
        if self.direction == 1:
            return tuple(range(self.sub_in, self.sub_out + 1))
        return tuple(range(self.sub_in, self.sub_out - 1, -1))


@dataclass(frozen=True)
class Beam:
    """A beam: source, ordered path, frequency offset, input polarization, nominal power.

    Attributes:
        id: beam id (``beam.<id>.*`` slots).
        source: laser / pump / source block name it originates from.
        path: ordered :class:`PathStep` tuple.
        freq_offset_hz: optical frequency offset from the source carrier [Hz] (conv §4 sign).
        jones_in: input Jones vector (complex, shape (2,)).
        power_w: nominal power [W].
        observer: observer id where the beam phase is evaluated.
        jones_out: analyzer Jones vector at the observer (complex (2,)) or None (WP-F6 addition;
            required when a jones segment is on the path).
    """

    id: str
    source: str
    path: tuple[PathStep, ...]
    freq_offset_hz: float
    jones_in: np.ndarray
    power_w: float
    observer: str
    jones_out: np.ndarray | None = None

    @property
    def segments(self) -> tuple[str, ...]:
        return tuple(s.segment for s in self.path)


@dataclass(frozen=True)
class Lock:
    """Lock = linear combination of beam phases held at a setpoint (conv §9).

    error e(t) = Σ_b combo[b]·φ_b(t) − setpoint_rad, evaluated at ``observer``. Device-specific
    error signs are converted inside the device's block (conv §9).
    """

    id: str
    combo: Mapping[str, float]
    setpoint_rad: float
    observer: str


@dataclass(frozen=True)
class Port:
    """Port k (D12, D31, G13). Port 0 is the reference port and always exists in port-based configs.

    ``squeezed`` / ``lo`` / ``probe`` are beam ids (None if absent); ``gbs_observer`` / ``bhd_observer``
    are observer ids (None if absent)."""

    index: int
    squeezed: str
    lo: str | None
    probe: str | None
    gbs_observer: str | None
    bhd_observer: str | None


@dataclass(frozen=True)
class Observer:
    """Point where beam phases are evaluated (detector, GBS input, interaction point).

    ``kind`` ∈ {"detector", "gbs_input", "interaction"}; ``segment``/``sub`` locate it."""

    id: str
    kind: str
    segment: str | None = None
    sub: int | None = None


@dataclass(frozen=True)
class StepTiming:
    """Nominal timing of one path step (module docstring).

    Attributes:
        segment, direction: as in the :class:`PathStep`.
        subs: traversed subsegments in traversal order.
        to_exit_s: delay from passing the centre of ``subs[j]`` to leaving the segment [s].
        tau_s: nominal group delay of the traversed part [s].
        after_s: nominal group delay from leaving this segment to the observer [s].
    """

    segment: str
    direction: int
    subs: tuple[int, ...]
    to_exit_s: tuple[float, ...]
    tau_s: float
    after_s: float


@dataclass(frozen=True)
class BeamTiming:
    """Nominal timing of a beam: ``tau_s`` = source → observer group delay; one entry per step."""

    beam: str
    tau_s: float
    steps: tuple[StepTiming, ...]


@dataclass(frozen=True)
class Topology:
    """Static description of the optical layout (PLAN §2.3.1).

    ``coef_coords`` (WP-F6 addition) holds the evaluation coordinates of each coefficient (G16-a).
    """

    segments: Mapping[str, Segment]
    beams: Mapping[str, Beam]
    coef: Mapping[CoefKey, CoefValue]
    locks: Mapping[str, Lock]
    ports: Mapping[int, Port]
    observers: Mapping[str, Observer]
    coef_coords: Mapping[CoefKey, Mapping[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_spec(cls, tree: Mapping[str, Any]) -> Topology:
        """Parse ``RunSpec.topology`` (module docstring); raise :class:`TopologyError` listing every
        problem."""
        return _Parser(tree).parse()

    @property
    def has_ports(self) -> bool:
        """False for port-less configs (Ref A–D): metric (ii) is then "not applicable" (G13)."""
        return bool(self.ports)

    # ------------------------------------------------------------------------------------- queries

    def beams_on(self, segment: str) -> tuple[str, ...]:
        """Ids of the beams traversing ``segment`` (topology order)."""
        return tuple(b.id for b in self.beams.values() if segment in b.segments)

    def step(self, beam: str, segment: str) -> PathStep:
        """The step of ``beam`` on ``segment``."""
        for s in self.beams[beam].path:
            if s.segment == segment:
                return s
        raise KeyError(f"beam {beam!r} does not traverse segment {segment!r}")

    def coef_value(self, key: CoefKey) -> float:
        """Evaluated coefficient ``key`` at its coordinates (``path.coef.evaluate_coef``)."""
        return evaluate_coef(self.coef[key], **dict(self.coef_coords.get(key, {})))

    def timing(self, beam: str) -> BeamTiming:
        """Nominal timing of ``beam`` (module docstring)."""
        b = self.beams[beam]
        raw = []
        for s in b.path:
            seg = self.segments[s.segment]
            subs = s.subs
            n = len(subs)
            d = seg.sub_tau_s
            raw.append((s, subs, tuple((n - j - 0.5) * d for j in range(n)), n * d))
        steps, after = [], 0.0
        for s, subs, to_exit, tau in reversed(raw):
            steps.append(StepTiming(s.segment, s.direction, subs, to_exit, tau, after))
            after += tau
        return BeamTiming(beam, after, tuple(reversed(steps)))

    # ------------------------------------------------------------------- registry helpers (F3)

    def register_loss_paths(self, registry: LossRegistry, beams: Iterable[str] | None = None) -> list[BeamPath]:
        """``LossRegistry.define_path`` for each beam (default all): its segments and its observer."""
        ids = list(self.beams) if beams is None else list(beams)
        return [registry.define_path(i, self.beams[i].segments, (self.beams[i].observer,)) for i in ids]

    def register_band_beams(self, ledger: BandLedger, beams: Iterable[str] | None = None) -> None:
        """``BandLedger.define_beam`` for each beam (default all) with its segments."""
        for i in (list(self.beams) if beams is None else list(beams)):
            ledger.define_beam(i, self.beams[i].segments)


def as_topology(obj: Topology | Mapping[str, Any]) -> Topology:
    """A :class:`Topology` from a Topology or its spec mapping."""
    if isinstance(obj, Topology):
        return obj
    if isinstance(obj, Mapping):
        return Topology.from_spec(obj)
    raise TypeError(f"topology must be a Topology or a spec mapping, got {type(obj).__name__}")


# --------------------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------------------


def _complex(x: Any, where: str) -> complex:
    if isinstance(x, (list, tuple)) and len(x) == 2 and all(isinstance(v, (int, float, np.number)) for v in x):
        return complex(float(x[0]), float(x[1]))
    if isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, complex, np.number)):
        raise TopologyError(f"{where}: expected a number, complex or [re, im], got {x!r}")
    return complex(x)


def _cvec(x: Any, where: str) -> np.ndarray:
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if not isinstance(x, (list, tuple)) or len(x) != 2:
        raise TopologyError(f"{where}: Jones vector needs 2 entries, got {x!r}")
    v = np.array([_complex(e, where) for e in x], dtype=np.complex128)
    if abs(np.vdot(v, v).real - 1.0) > 1e-9:
        raise TopologyError(f"{where}: Jones vector must have unit norm, got |v|² = {np.vdot(v, v).real!r}")
    v.flags.writeable = False
    return v


def _cmat(x: Any, where: str) -> np.ndarray:
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if not isinstance(x, (list, tuple)) or len(x) != 2 or any(not isinstance(r, (list, tuple)) or len(r) != 2 for r in x):
        raise TopologyError(f"{where}: 2×2 matrix expected, got {x!r}")
    m = np.array([[_complex(e, where) for e in r] for r in x], dtype=np.complex128)
    m.flags.writeable = False
    return m


class _Parser:
    REQUIRED_TOP = ("segments", "beams", "observers", "locks", "ports", "coef")

    def __init__(self, tree: Mapping[str, Any]):
        if not isinstance(tree, Mapping):
            raise TopologyError(f"topology spec must be a mapping, got {type(tree).__name__}")
        self.tree = tree
        self.errs: list[str] = []

    def err(self, msg: str) -> None:
        self.errs.append(msg)

    def need(self, d: Mapping[str, Any], keys: Iterable[str], where: str) -> bool:
        missing = [k for k in keys if k not in d]
        if missing:
            self.err(f"{where}: missing required key(s) {missing} (D34: no code defaults)")
        return not missing

    def items(self, key: str) -> list[Mapping[str, Any]]:
        v = self.tree.get(key)
        if v is None:
            return []
        if not isinstance(v, (list, tuple)):
            self.err(f"topology.{key} must be a list, got {type(v).__name__}")
            return []
        out = []
        for i, e in enumerate(v):
            if not isinstance(e, Mapping):
                self.err(f"topology.{key}[{i}] must be a mapping")
            else:
                out.append(e)
        return out

    def parse(self) -> Topology:
        missing = [k for k in self.REQUIRED_TOP if k not in self.tree]
        if missing:
            raise TopologyError(f"topology: missing required key(s) {missing} (D34; empty lists are fine)")
        unknown = set(self.tree) - set(self.REQUIRED_TOP)
        if unknown:
            self.err(f"topology: unknown key(s) {sorted(unknown)}")
        segments = self.segments()
        observers = self.observers(segments)
        beams = self.beams(segments, observers)
        locks = self.locks(beams, observers)
        ports = self.ports(beams, observers)
        coef, coords = self.coefs(segments, beams)
        if self.errs:
            raise TopologyError("invalid topology:\n  - " + "\n  - ".join(self.errs))
        return Topology(segments, beams, coef, locks, ports, observers, coords)

    # segments -------------------------------------------------------------------------------------

    def segments(self) -> dict[str, Segment]:
        out: dict[str, Segment] = {}
        keys = ("id", "kind", "L_m", "n_g", "beta2_s2_per_m", "beta3_s3_per_m", "zone", "n_sub", "spatial",
                "pol_model", "pol_params", "drivers")
        for i, s in enumerate(self.items("segments")):
            where = f"segments[{s.get('id', i)!r}]"
            unknown = set(s) - set(keys) - {"actuator_sub"}
            if unknown:
                self.err(f"{where}: unknown key(s) {sorted(unknown)}")
            if not self.need(s, keys, where):
                continue
            sid = s["id"]
            if not isinstance(sid, str) or not sid.isidentifier():
                self.err(f"{where}: id must be an identifier (slot name component)")
                continue
            if sid in out:
                self.err(f"{where}: duplicate segment id")
                continue
            ok = True
            if s["kind"] not in SEGMENT_KINDS:
                self.err(f"{where}: kind must be one of {SEGMENT_KINDS}")
                ok = False
            if s["spatial"] not in SPATIAL_STRUCTURES:
                self.err(f"{where}: spatial must be one of {SPATIAL_STRUCTURES}")
                ok = False
            if s["pol_model"] not in POL_MODELS:
                self.err(f"{where}: pol_model must be one of {POL_MODELS}")
                ok = False
            n_sub = s["n_sub"]
            if isinstance(n_sub, bool) or not isinstance(n_sub, (int, np.integer)) or n_sub < 1:
                self.err(f"{where}: n_sub must be an int >= 1")
                continue
            n_sub = int(n_sub)
            try:
                L, n_g = float(s["L_m"]), float(s["n_g"])
                b2, b3 = float(s["beta2_s2_per_m"]), float(s["beta3_s3_per_m"])
            except (TypeError, ValueError):
                self.err(f"{where}: L_m, n_g, beta2_s2_per_m, beta3_s3_per_m must be numbers")
                continue
            if not (np.isfinite(L) and L >= 0 and np.isfinite(n_g) and n_g > 0):
                self.err(f"{where}: need L_m >= 0 and n_g > 0")
                ok = False
            act = s.get("actuator_sub")
            if act is not None and (isinstance(act, bool) or not isinstance(act, (int, np.integer)) or not 0 <= act < n_sub):
                self.err(f"{where}: actuator_sub must be null or an int in [0, {n_sub})")
                ok = False
            act = None if act is None else int(act)
            drivers = self.drivers(s["drivers"], n_sub, act, where)
            pol = dict(s["pol_params"] or {})
            if ok:
                ok = self.check_pol(s["pol_model"], pol, drivers, where)
                ok = self.check_spatial(s["spatial"], drivers, n_sub, where) and ok
            if ok and drivers is not None:
                out[sid] = Segment(sid, s["kind"], L, n_g, b2, b3, str(s["zone"]), n_sub, s["spatial"],
                                   s["pol_model"], pol, act, drivers)
        return out

    def drivers(self, spec: Any, n_sub: int, act: int | None, where: str) -> dict[str, tuple] | None:
        if not isinstance(spec, Mapping):
            self.err(f"{where}: drivers must be a mapping {{driver: placement}}")
            return None
        out: dict[str, tuple] = {}
        for drv, pl in spec.items():
            w = f"{where}.drivers.{drv}"
            if drv not in DRIVERS:
                self.err(f"{w}: unknown driver (known {DRIVERS})")
                continue
            if pl == "uniform" or pl == "per_sub":
                out[drv] = (pl,)
            elif pl == "actuator":
                if act is None:
                    self.err(f"{w}: placement 'actuator' needs actuator_sub")
                    continue
                out[drv] = ("point", act)
            elif isinstance(pl, Mapping) and set(pl) == {"point"}:
                k = pl["point"]
                if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 0 <= k < n_sub:
                    self.err(f"{w}: point must be an int in [0, {n_sub})")
                    continue
                out[drv] = ("point", int(k))
            elif isinstance(pl, Mapping) and set(pl) == {"subs"}:
                ks = pl["subs"]
                if (not isinstance(ks, (list, tuple)) or not ks or len(set(ks)) != len(ks)
                        or any(isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 0 <= k < n_sub for k in ks)):
                    self.err(f"{w}: subs must be a non-empty list of distinct ints in [0, {n_sub})")
                    continue
                out[drv] = ("subs", tuple(sorted(int(k) for k in ks)))
            else:
                self.err(f"{w}: placement must be uniform | per_sub | actuator | {{point: k}} | {{subs: [...]}}, got {pl!r}")
        return out

    def check_pol(self, model: str, pol: dict[str, Any], drivers: dict[str, tuple] | None, where: str) -> bool:
        if model == "off":
            if pol:
                self.err(f"{where}: pol_model off takes no pol_params")
                return False
            return True
        if model == "pm_diagonal":
            missing = [k for k in _PM_KEYS if k not in pol]
            extra = set(pol) - set(_PM_KEYS)
            if missing or extra:
                self.err(f"{where}: pm_diagonal pol_params need exactly {_PM_KEYS} (missing {missing}, extra {sorted(extra)})")
                return False
            try:
                for k in _PM_KEYS[:3]:
                    pol[k] = float(pol[k])
            except (TypeError, ValueError):
                self.err(f"{where}: pm_diagonal K_T / K_tau / b_g must be numbers")
                return False
            if pol["K_T_rad_per_m_K"] > 0:  # allowed (D34), but D51-2 / conv §9: K_T < 0 (Phase 5 F6-4)
                from cranberry.core.checks import ChkLog

                ChkLog().warn("pm_K_T_sign", where, f"K_T_rad_per_m_K = {pol['K_T_rad_per_m_K']!r} > 0; the PM birefringence "
                              "temperature coefficient is negative (conv §9, D51-2: nominal -2.05 rad/(m K)); a positive "
                              "value flips the sign of the slow-fast phase drift", K_T_rad_per_m_K=pol["K_T_rad_per_m_K"])
            ad = pol["axis_drivers"]
            if not isinstance(ad, (list, tuple)) or any(d not in DRIVERS for d in ad):
                self.err(f"{where}: axis_drivers must be a list of drivers from {DRIVERS}")
                return False
            pol["axis_drivers"] = tuple(ad)
            return True
        # jones (optional J_nr_slot: non-reciprocal factor, WP-B7 addition)
        nr = pol.pop("J_nr_slot", None)
        ok = self._check_jones(pol, where)
        if nr is not None:
            if not isinstance(nr, str):
                self.err(f"{where}: jones pol_params J_nr_slot must be a slot name")
                return False
            pol["J_nr_slot"] = nr
        return ok

    def _check_jones(self, pol: dict[str, Any], where: str) -> bool:
        if set(pol) == {"J"}:
            try:
                pol["J"] = _cmat(pol["J"], f"{where}.pol_params.J")
            except TopologyError as e:
                self.err(str(e))
                return False
            return True
        if set(pol) == {"J_slot"} and isinstance(pol["J_slot"], str):
            return True
        self.err(f"{where}: jones pol_params need exactly one of J (2×2) or J_slot (slot name)")
        return False

    def check_spatial(self, spatial: str, drivers: dict[str, tuple] | None, n_sub: int, where: str) -> bool:
        if drivers is None or n_sub == 1:
            return True
        if spatial == "lumped":
            bad = [d for d, p in drivers.items() if p[0] == "per_sub"]
            if bad:
                self.err(f"{where}: spatial lumped forbids per_sub placement ({bad})")
                return False
        elif spatial == "distributed":
            bad = [d for d, p in drivers.items() if d in _ENV_DRIVERS and p[0] != "per_sub"]
            if bad:
                self.err(f"{where}: spatial distributed needs per_sub placement for {bad}")
                return False
        return True

    # observers / beams ----------------------------------------------------------------------------

    def observers(self, segments: Mapping[str, Segment]) -> dict[str, Observer]:
        out: dict[str, Observer] = {}
        for i, o in enumerate(self.items("observers")):
            where = f"observers[{o.get('id', i)!r}]"
            if set(o) - {"id", "kind", "segment", "sub"}:
                self.err(f"{where}: unknown key(s) {sorted(set(o) - {'id', 'kind', 'segment', 'sub'})}")
            if not self.need(o, ("id", "kind"), where):
                continue
            if o["id"] in out:
                self.err(f"{where}: duplicate observer id")
                continue
            if o["kind"] not in OBSERVER_KINDS:
                self.err(f"{where}: kind must be one of {OBSERVER_KINDS}")
                continue
            seg, sub = o.get("segment"), o.get("sub")
            if seg is not None and seg not in segments:
                self.err(f"{where}: unknown segment {seg!r}")
                continue
            out[o["id"]] = Observer(o["id"], o["kind"], seg, None if sub is None else int(sub))
        return out

    def beams(self, segments: Mapping[str, Segment], observers: Mapping[str, Observer]) -> dict[str, Beam]:
        out: dict[str, Beam] = {}
        keys = ("id", "source", "path", "freq_offset_hz", "jones_in", "power_w", "observer")
        for i, b in enumerate(self.items("beams")):
            where = f"beams[{b.get('id', i)!r}]"
            unknown = set(b) - set(keys) - {"jones_out"}
            if unknown:
                self.err(f"{where}: unknown key(s) {sorted(unknown)}")
            if not self.need(b, keys, where):
                continue
            bid = b["id"]
            if not isinstance(bid, str) or not bid.isidentifier():
                self.err(f"{where}: id must be an identifier (slot name component)")
                continue
            if bid in out:
                self.err(f"{where}: duplicate beam id")
                continue
            if b["observer"] not in observers:
                self.err(f"{where}: unknown observer {b['observer']!r}")
                continue
            steps = self.path(b["path"], segments, where)
            if steps is None:
                continue
            try:
                jin = _cvec(b["jones_in"], f"{where}.jones_in")
                jout = None if b.get("jones_out") is None else _cvec(b["jones_out"], f"{where}.jones_out")
            except TopologyError as e:
                self.err(str(e))
                continue
            if jout is None and any(segments[s.segment].pol_model == "jones" for s in steps):
                self.err(f"{where}: jones_out (analyzer at the observer) is required: a jones segment is on the path")
                continue
            obs = observers[b["observer"]]
            if obs.segment is not None and steps and (obs.segment != steps[-1].segment
                                                       or (obs.sub is not None and obs.sub != steps[-1].sub_out)):
                self.err(f"{where}: observer {obs.id!r} is at ({obs.segment}, {obs.sub}) but the path ends at "
                         f"({steps[-1].segment}, {steps[-1].sub_out})")
                continue
            out[bid] = Beam(bid, str(b["source"]), steps, float(b["freq_offset_hz"]), jin, float(b["power_w"]),
                            b["observer"], jout)
        return out

    def path(self, spec: Any, segments: Mapping[str, Segment], where: str) -> tuple[PathStep, ...] | None:
        if not isinstance(spec, (list, tuple)) or not spec:
            self.err(f"{where}: path must be a non-empty list of steps")
            return None
        steps, seen, ok = [], set(), True
        for j, st in enumerate(spec):
            w = f"{where}.path[{j}]"
            if not isinstance(st, Mapping) or not self.need(st, ("segment", "direction"), w):
                ok = False
                continue
            unknown = set(st) - {"segment", "direction", "sub_in", "sub_out", "axis"}
            if unknown:
                self.err(f"{w}: unknown key(s) {sorted(unknown)}")
                ok = False
            sid = st["segment"]
            if sid not in segments:
                self.err(f"{w}: unknown segment {sid!r}")
                ok = False
                continue
            if sid in seen:
                self.err(f"{w}: beam traverses segment {sid!r} twice (split the segment)")
                ok = False
                continue
            seen.add(sid)
            seg = segments[sid]
            d = st["direction"]
            if isinstance(d, bool) or d not in (1, -1):
                self.err(f"{w}: direction must be +1 or -1")
                ok = False
                continue
            d = int(d)
            n = seg.n_sub
            si = st.get("sub_in", 0 if d == 1 else n - 1)
            so = st.get("sub_out", n - 1 if d == 1 else 0)
            if any(isinstance(x, bool) or not isinstance(x, (int, np.integer)) or not 0 <= x < n for x in (si, so)):
                self.err(f"{w}: sub_in / sub_out must be ints in [0, {n})")
                ok = False
                continue
            if (d == 1 and si > so) or (d == -1 and si < so):
                self.err(f"{w}: direction {d:+d} needs sub_in {'<=' if d == 1 else '>='} sub_out")
                ok = False
                continue
            axis = st.get("axis")
            if seg.pol_model == "pm_diagonal":
                if axis not in AXES:
                    self.err(f"{w}: pm_diagonal segment needs axis in {AXES}")
                    ok = False
                    continue
            elif axis is not None:
                self.err(f"{w}: axis is only used on pm_diagonal segments")
                ok = False
                continue
            steps.append(PathStep(sid, d, int(si), int(so), axis))
        return tuple(steps) if ok else None

    # locks / ports --------------------------------------------------------------------------------

    def locks(self, beams: Mapping[str, Beam], observers: Mapping[str, Observer]) -> dict[str, Lock]:
        out: dict[str, Lock] = {}
        for i, lk in enumerate(self.items("locks")):
            where = f"locks[{lk.get('id', i)!r}]"
            if "setpoint" in lk and "setpoint_rad" in lk:
                self.err(f"{where}: give setpoint_rad or setpoint, not both")
                continue
            lk = dict(lk)
            if "setpoint" in lk:
                lk["setpoint_rad"] = lk.pop("setpoint")
            if set(lk) - {"id", "combo", "setpoint_rad", "observer"}:
                self.err(f"{where}: unknown key(s) {sorted(set(lk) - {'id', 'combo', 'setpoint_rad', 'observer'})}")
            if not self.need(lk, ("id", "combo", "setpoint_rad", "observer"), where):
                continue
            lid = lk["id"]
            if not isinstance(lid, str) or not lid.isidentifier() or lid in out:
                self.err(f"{where}: id must be a unique identifier")
                continue
            combo = lk["combo"]
            if not isinstance(combo, Mapping) or not combo:
                self.err(f"{where}: combo must be a non-empty mapping {{beam: coefficient}}")
                continue
            bad = [b for b in combo if b not in beams]
            if bad:
                self.err(f"{where}: unknown beam(s) {bad}")
                continue
            if lk["observer"] not in observers:
                self.err(f"{where}: unknown observer {lk['observer']!r}")
                continue
            off = [b for b in combo if beams[b].observer != lk["observer"]]
            if off:
                self.err(f"{where}: beam(s) {off} are not evaluated at the lock observer {lk['observer']!r}")
                continue
            out[lid] = Lock(lid, {str(b): float(c) for b, c in combo.items()}, float(lk["setpoint_rad"]), lk["observer"])
        return out

    def ports(self, beams: Mapping[str, Beam], observers: Mapping[str, Observer]) -> dict[int, Port]:
        out: dict[int, Port] = {}
        keys = ("index", "squeezed", "lo", "probe", "gbs_observer", "bhd_observer")
        for i, p in enumerate(self.items("ports")):
            where = f"ports[{p.get('index', i)!r}]"
            if set(p) - set(keys):
                self.err(f"{where}: unknown key(s) {sorted(set(p) - set(keys))}")
            if not self.need(p, keys, where):
                continue
            k = p["index"]
            if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k < 0 or k in out:
                self.err(f"{where}: index must be a unique int >= 0")
                continue
            bad = [p[r] for r in ("squeezed", "lo", "probe") if p[r] is not None and p[r] not in beams]
            bad += [p[r] for r in ("gbs_observer", "bhd_observer") if p[r] is not None and p[r] not in observers]
            if p["squeezed"] is None:
                bad.append("squeezed=None")
            if bad:
                self.err(f"{where}: unknown / missing reference(s) {bad}")
                continue
            out[int(k)] = Port(int(k), p["squeezed"], p["lo"], p["probe"], p["gbs_observer"], p["bhd_observer"])
        if out and 0 not in out:
            self.err("ports: the reference port 0 must exist when ports are given (G13)")
        return out

    # coefficients ---------------------------------------------------------------------------------

    def coefs(self, segments: Mapping[str, Segment], beams: Mapping[str, Beam]):
        spec = self.tree.get("coef")
        entries: list[tuple[str, str, str, str, Any, Mapping[str, Any]]] = []
        if isinstance(spec, Mapping):  # Python form {CoefKey | tuple: value}
            for k, v in spec.items():
                if not isinstance(k, tuple) or len(k) != 4:
                    self.err(f"coef: mapping keys must be (segment, driver, beam, effect), got {k!r}")
                    continue
                entries.append((*map(str, k), v, {}))
        elif isinstance(spec, (list, tuple)):
            for i, e in enumerate(spec):
                where = f"coef[{i}]"
                if not isinstance(e, Mapping) or not self.need(e, ("segment", "driver", "beam", "effect", "value"), where):
                    continue
                if set(e) - {"segment", "driver", "beam", "effect", "value", "coords"}:
                    self.err(f"{where}: unknown key(s) {sorted(set(e) - {'segment', 'driver', 'beam', 'effect', 'value', 'coords'})}")
                coords = e.get("coords") or {}
                if not isinstance(coords, Mapping):
                    self.err(f"{where}: coords must be a mapping")
                    continue
                entries.append((e["segment"], e["driver"], e["beam"], e["effect"], e["value"], dict(coords)))
        elif spec is not None:
            self.err("coef must be a list of entries or a {(segment, driver, beam, effect): value} mapping")
        on = {s: [b.id for b in beams.values() if s in b.segments] for s in segments}
        coef: dict[CoefKey, Any] = {}
        coords: dict[CoefKey, Mapping[str, Any]] = {}
        explicit: set[CoefKey] = set()
        for seg, drv, beam, eff, val, crd in entries:
            where = f"coef ({seg}, {drv}, {beam}, {eff})"
            if seg not in segments:
                self.err(f"{where}: unknown segment")
                continue
            if drv not in segments[seg].drivers:
                self.err(f"{where}: driver not declared on segment {seg!r} (drivers: {sorted(segments[seg].drivers)})")
                continue
            if eff not in EFFECTS:
                self.err(f"{where}: effect must be one of {EFFECTS}")
                continue
            if beam == "*":
                targets = on[seg]
            elif beam in on[seg]:
                targets = [beam]
            else:
                self.err(f"{where}: beam does not traverse segment {seg!r} (beams on it: {on[seg]})")
                continue
            for b in targets:
                key = CoefKey(seg, drv, b, eff)
                if beam != "*":
                    if key in explicit:
                        self.err(f"{where}: given twice")
                        continue
                    explicit.add(key)
                elif key in coef and key not in explicit:
                    self.err(f"{where}: wildcard given twice")
                    continue
                elif key in explicit:
                    continue
                coef[key], coords[key] = val, crd
        missing = [CoefKey(s, d, b, e) for s, seg in segments.items() for d in seg.drivers for b in on[s]
                   for e in EFFECTS if CoefKey(s, d, b, e) not in coef]
        if missing:
            self.err("missing coefficients (explicit 0 is fine; D11/D34): "
                     + ", ".join(f"({k.segment}, {k.driver}, {k.beam}, {k.effect})" for k in missing))
        for key in list(coef):
            try:
                evaluate_coef(coef[key], **dict(coords[key]))
            except (TypeError, ValueError, KeyError) as e:
                self.err(f"coef ({key.segment}, {key.driver}, {key.beam}, {key.effect}): {e}")
        return coef, coords
