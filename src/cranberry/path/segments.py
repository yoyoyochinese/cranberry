"""``segment_phase`` block: drivers × coefficients → ψ (segment × beam, per subsegment) and τ_g
(PLAN §2.3.2, §4.3 PATH, conv §9, T08, T10). Owner: WP-F6.

Reads ``drv.<seg>.<driver>`` (per subsegment; the history ring is kept by this block), writes
``psi.<seg>.<beam>`` [rad] and ``tau_g.<seg>.<beam>`` [s]. Spatial structure (lumped / distributed / mixed) and zones come from the
Topology. TEC writes temperature, not phase (G08); "write phase directly" is reproduced by a
coefficient preset (equal axis coefficients, zero group-delay effect), not by a separate code path.

Implementation (WP-F6)
    For every segment i handled by the block and every beam b traversing it (path step with
    traversed subsegments m_0 … m_{n−1}, ``Topology.timing``):

    * per-subsegment phase density, at the current sample:
      ψ_{i,m,b}(t) = Σ_drv c_{i,drv,b}^{phase} · x_{drv,m}(t), with x_{drv,m} = w_{drv,m}·drv (scalar
      slot, placement weights ``Segment.weights``) or drv[m]/N_sub (``per_sub`` slot of shape (N_sub,));
      pushed into the ring ``psi.<seg>.<beam>.sub`` (shape (n,), one row per step of this block);
    * the segment phase the beam carries when it **leaves the segment now** (PLAN §2.3.2):
      ``psi.<seg>.<beam>`` = Σ_j ψ_{i,m_j,b}(t − (n − j − ½)·δ_i), read from the ring with linear
      fractional interpolation (``core.delay`` arithmetic; δ_i = n_g L/(c N_sub)). A counter-propagating
      beam reads the same subsegments in the opposite order, i.e. at other times (T08 §2.7, T10 §3.2);
    * the segment phase as it **arrives at the beam's observer now** (Phase 5 F6-2, D65-4):
      ``psi.<seg>.<beam>.obs`` = Σ_j ψ_{i,m_j,b}(t − (n − j − ½)·δ_i − τ_after), read from the same ring
      with **one** linear interpolation at the total delay (τ_after: segment exit → observer, nominal,
      ``Topology.timing``). ``beam_assembler`` uses this slot (with emission time on) instead of delaying
      ``psi.<seg>.<beam>`` a second time, so a subsegment's phase reaches the observer through a single
      linear fractional delay (two cascaded stages attenuate more: white driver, μ₁ = 0.3, μ₂ = 0.4 →
      variance 0.40 instead of 0.58). Both are evaluated at this block's step h = rate_div/f_sim;
    * the group delay of the traversed part ``tau_g.<seg>.<beam>`` = τ_0 + Σ_drv c^{group_delay} Σ_j
      x_{drv,m_j}(t) (current sample; thermal / strain group delay is slow, T03 §8). τ_0 = n·δ_i, minus
      b_g·L·n/N for a fast-axis beam on a pm_diagonal segment when ``pm_axis_drift`` is on.

    PM fast path (pm_diagonal, conv §9, T09 §1.2): the coefficients are the slow-axis responses; for a
    fast-axis beam the block adds −K_T·L to the phase coefficient and −K_τ·L to the group-delay
    coefficient of every driver in ``axis_drivers`` (dψ_fast/dT = L(C_T − K_T) = L(C_T + |K_T|), so
    Δψ_b = ψ_slow − ψ_fast changes by K_T·L·ΔT < 0). Polarization matrices of ``jones`` segments are
    applied by ``beam_assembler``.

    Coefficients are runtime parameters (P) named ``<block>.<seg>.<beam>.<driver>.<effect>`` holding
    the total value (topology coefficient + fast-axis term), so they can be changed with
    ``Params.set`` without recompiling.

Slot contract (WP-F6; ``tau_g.<seg>`` of the WP-F0 stub became per beam, because coefficients and the
PM DGD are per segment × beam):

    reads   drv.<seg>.<driver>       scalar, or (N_sub,) for per_sub placements (written by NOISE /
                                      TEC / ACT blocks, which register them)
    writes  psi.<seg>.<beam>         rad   segment phase of a beam leaving the segment now
            psi.<seg>.<beam>.obs     rad   the same arriving at the beam's observer now (one
                                            interpolation at to_exit + τ_after; Phase 5 addition)
            psi.<seg>.<beam>.sub     rad   ring (ring_len, n): per-subsegment phase history
                                            (ring_len = max integer delay incl. τ_after + 2)
            tau_g.<seg>.<beam>       s     group delay of the traversed part

Params (all required, D34)::

    topology : mapping (RunSpec.topology form) or Topology
    segments : "all" | [segment id, ...]
    toggles  : {phase: bool, group_delay: bool, pm_axis_drift: bool}
      phase off          → psi.* = 0 (identity of the path phase model), ring untouched
      group_delay off    → tau_g.* = nominal n·δ (no driver response, no DGD)
      pm_axis_drift off  → no fast-axis terms: ψ_fast = ψ_slow for equal coefficients (T09 test 8)

Delays are nominal (build time): the group-delay change Δτ_g ≲ ps does not move the subsegment read
times (it would shift a ~µs sample grid by 1e−6 of a sample). The frequency-offset phase drift
2πΔν·Δτ_g of a shifted probe belongs in that beam's phase coefficient (T10 §3.3).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.path.coef import EFFECTS, CoefKey
from cranberry.path.topology import Topology, as_topology

__all__ = ["SegmentPhase", "make_segment_phase", "DRIVER_UNITS", "frac_delay", "select_ids", "subsegment_channels"]

#: Unit tags expected for driver slots (stretcher / eom units depend on the actuator preset).
DRIVER_UNITS: dict[str, str] = {"dT": "K", "tec": "K", "strain": "1", "detuning": "rad/s"}

_INT_TOL = 1e-9


def frac_delay(delay_s: float, h_s: float) -> tuple[int, float]:
    """(k, μ) with delay/h = k + μ, μ ∈ [0, 1): y = (1 − μ)·x[n − k] + μ·x[n − k − 1] (linear
    interpolation, ``core.delay``); integer delays (within 1e−9) give μ = 0 exactly."""
    D = delay_s / h_s
    if D < -_INT_TOL:
        raise ValueError(f"negative delay {delay_s!r} s")
    r = round(D)
    if abs(D - r) <= _INT_TOL:
        return int(r), 0.0
    k = int(math.floor(D))
    return k, D - k


def subsegment_channels(topology: Any, segment: str, driver: str, combo: Mapping[str, float], *,
                        emission_time: bool = True) -> list[list[tuple[float, float]]]:
    """Channels for ``core.delay.interp_variance_deficit``: how ``driver`` on ``segment`` reaches the
    signed beam combination Σ_b combo[b]·φ_b at the beams' observers (the delays ``segment_phase``
    interpolates once: to_exit + τ_after, F6-2).

    Taps are (delay [s], weight = combo[b] · phase coefficient (topology, slow axis) · placement
    weight). A ``per_sub`` driver gives one independent channel per subsegment (weight 1/N_sub, the
    driver slot holds per-subsegment values); any other placement is one coherent channel. With
    ``emission_time`` off only the in-segment delays remain (``beam_assembler`` then uses τ_after = 0).
    """
    from cranberry.path.coef import CoefKey

    topo = as_topology(topology)
    seg = topo.segments[segment]
    w = seg.weights(driver)
    per_sub: dict[int, list[tuple[float, float]]] = {}
    for b, cb in combo.items():
        st = next(s for s in topo.timing(b).steps if s.segment == segment)
        c = float(cb) * topo.coef_value(CoefKey(segment, driver, b, "phase"))
        after = st.after_s if emission_time else 0.0
        for m, d in zip(st.subs, st.to_exit_s):
            wm = 1.0 / seg.n_sub if w is None else float(w[m])
            if wm != 0.0:
                per_sub.setdefault(m if w is None else 0, []).append((d + after, c * wm))
    return list(per_sub.values())


def select_ids(spec: Any, known: Mapping[str, Any], what: str, block: str) -> list[str]:
    """``"all"`` → every id of ``known`` (in order); a list → validated ids."""
    if spec == "all":
        return list(known)
    if not isinstance(spec, (list, tuple)) or not spec:
        raise ValueError(f"block {block!r}: {what} must be 'all' or a non-empty list of ids")
    bad = [s for s in spec if s not in known]
    if bad:
        raise ValueError(f"block {block!r}: unknown {what} {bad}")
    if len(set(spec)) != len(spec):
        raise ValueError(f"block {block!r}: duplicate {what}")
    return list(spec)


def _ring_read(S, base, ntr, L, h, j, k, mu):  # noqa: N803
    """Column j of the (L, ntr) ring at ``base``, k + μ rows older than head row h (linear interpolation)."""
    a = h - k
    if a < 0:
        a += L
    y = S[base + a * ntr + j]
    if mu != 0.0:
        b = a - 1
        if b < 0:
            b += L
        y = (1.0 - mu) * y + mu * S[base + b * ntr + j]
    return y


_ring_read_nb = numba.njit(inline="always")(_ring_read)


class _Pair:
    """Build-time description of one (segment, beam) pair."""

    def __init__(self, topo: Topology, seg_id: str, beam: str, h: float, axis_on: bool):
        seg = topo.segments[seg_id]
        step = topo.step(beam, seg_id)
        st = next(s for s in topo.timing(beam).steps if s.segment == seg_id)
        self.seg, self.beam = seg_id, beam
        self.subs = np.array(st.subs, dtype=np.int64)
        kmu = [frac_delay(d, h) for d in st.to_exit_s]
        self.k = np.array([k for k, _ in kmu], dtype=np.int64)
        self.mu = np.array([m for _, m in kmu], dtype=np.float64)
        kmo = [frac_delay(d + st.after_s, h) for d in st.to_exit_s]  # to the observer: one interpolation
        self.ko = np.array([k for k, _ in kmo], dtype=np.int64)
        self.muo = np.array([m for _, m in kmo], dtype=np.float64)
        self.ring_len = int(max(self.k.max(), self.ko.max())) + 2
        fast = seg.pol_model == "pm_diagonal" and step.axis == "fast" and axis_on
        pol = seg.pol_params
        frac = len(st.subs) / seg.n_sub
        self.tau0 = st.tau_s - (pol["b_g_s_per_m"] * seg.length_m * frac if fast else 0.0)
        self.drivers = list(seg.drivers)
        self.coef: dict[tuple[str, str], float] = {}
        for drv in self.drivers:
            for eff in EFFECTS:
                c = topo.coef_value(CoefKey(seg_id, drv, beam, eff))
                if fast and drv in pol["axis_drivers"]:
                    c += -(pol["K_T_rad_per_m_K"] if eff == "phase" else pol["K_tau_s_per_m_K"]) * seg.length_m
                self.coef[(drv, eff)] = c
        self.weights = {d: seg.weights(d) for d in self.drivers}  # None for per_sub
        self.n_sub = seg.n_sub
        self.psi = f"psi.{seg_id}.{beam}"
        self.obs = f"psi.{seg_id}.{beam}.obs"
        self.ring = f"psi.{seg_id}.{beam}.sub"
        self.tau = f"tau_g.{seg_id}.{beam}"


class SegmentPhase(BaseBlock):
    """The ``segment_phase`` block (module docstring)."""

    toggle_names = ("phase", "group_delay", "pm_axis_drift")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or "segment_phase"
        require(params, "topology", "segments", "toggles", block=name)
        topo = as_topology(params["topology"])
        seg_ids = select_ids(params["segments"], topo.segments, "segments", name)
        self.topo = topo
        self.pairs_spec = [(s, b) for s in seg_ids for b in topo.beams_on(s)]
        reads = [f"drv.{s}.{d}" for s in seg_ids for d in topo.segments[s].drivers
                 if topo.beams_on(s)]
        writes = [w for s, b in self.pairs_spec
                  for w in (f"psi.{s}.{b}", f"psi.{s}.{b}.obs", f"psi.{s}.{b}.sub", f"tau_g.{s}.{b}")]
        super().__init__("segment_phase", name, params, reads=reads, writes=writes)
        self.slot_units = {r: DRIVER_UNITS[r.rsplit(".", 1)[1]] for r in reads if r.rsplit(".", 1)[1] in DRIVER_UNITS}
        self.pairs: list[_Pair] = []

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        h = block_dt_s(layout.f_sim_hz, self.rate_div)
        axis_on = self.toggles["pm_axis_drift"]
        self.pairs = [_Pair(self.topo, s, b, h, axis_on) for s, b in self.pairs_spec]
        for p in self.pairs:
            layout.add(p.psi, (), np.float64, "rad")
            layout.add(p.obs, (), np.float64, "rad")
            layout.add(p.ring, (len(p.subs),), np.float64, "rad", ring_len=p.ring_len)
            layout.add(p.tau, (), np.float64, "s", init=p.tau0)
            for (drv, eff), c in p.coef.items():
                layout.add_param(self._pname(p, drv, eff), c, (), "s" if eff == "group_delay" else "rad")
            self.slot_units.update({p.psi: "rad", p.obs: "rad", p.tau: "s"})

    def _pname(self, p: _Pair, drv: str, eff: str) -> str:
        return f"{self.name}.{p.seg}.{p.beam}.{drv}.{eff}"

    # ---------------------------------------------------------------------------------- binding

    def _check_driver_slots(self, layout) -> None:
        for p in self.pairs:
            for d in p.drivers:
                sl = layout.slot(f"drv.{p.seg}.{d}")
                want = () if p.weights[d] is not None else (p.n_sub,)
                if sl.array != "S" or sl.ring_len or sl.shape not in (want, (1,) if want == () else want):
                    raise TypeError(f"block {self.name!r}: slot {sl.name!r} must be a float64 slot of shape "
                                    f"{want} (placement {self.topo.segments[p.seg].drivers[d][0]}), got {sl.shape}")

    def bind(self, layout) -> None:
        self._check_driver_slots(layout)
        for p in self.pairs:
            p.psi_off = layout.slot(p.psi).offset
            p.obs_off = layout.slot(p.obs).offset
            p.tau_off = layout.slot(p.tau).offset
            r = layout.slot(p.ring)
            p.ring_off, p.head = r.offset, r.head
            p.drv_off = {d: layout.slot(f"drv.{p.seg}.{d}").offset for d in p.drivers}
            p.p_off = {k: layout.param(self._pname(p, *k)).offset for k in p.coef}
        self.kernel = self._build_kernel()

    def _build_kernel(self):
        # flat tables (structure only: offsets, indices, weights, nominal delays; coefficients stay in P)
        pr = {k: [] for k in ("psi", "obs", "tau", "roff", "rlen", "head", "ntr", "tr0", "tm0", "ntm")}
        tau0, tr_sub, tr_k, tr_mu, tr_ko, tr_muo = [], [], [], [], [], []
        tm_kind, tm_off, tm_w0, tm_ninv, tm_pph, tm_pgd, wflat = [], [], [], [], [], [], []
        for p in self.pairs:
            pr["psi"].append(p.psi_off)
            pr["obs"].append(p.obs_off)
            pr["tau"].append(p.tau_off)
            pr["roff"].append(p.ring_off)
            pr["rlen"].append(p.ring_len)
            pr["head"].append(p.head)
            pr["ntr"].append(len(p.subs))
            pr["tr0"].append(len(tr_sub))
            pr["tm0"].append(len(tm_kind))
            pr["ntm"].append(len(p.drivers))
            tau0.append(p.tau0)
            tr_sub += p.subs.tolist()
            tr_k += p.k.tolist()
            tr_mu += p.mu.tolist()
            tr_ko += p.ko.tolist()
            tr_muo += p.muo.tolist()
            for d in p.drivers:
                w = p.weights[d]
                tm_kind.append(0 if w is not None else 1)
                tm_off.append(p.drv_off[d])
                tm_w0.append(len(wflat))
                if w is not None:
                    wflat += w.tolist()
                tm_ninv.append(1.0 / p.n_sub)
                tm_pph.append(p.p_off[(d, "phase")])
                tm_pgd.append(p.p_off[(d, "group_delay")])
        i64 = lambda x: np.array(x, dtype=np.int64)  # noqa: E731
        f64 = lambda x: np.array(x, dtype=np.float64)  # noqa: E731
        P_psi, P_obs, P_tau, P_roff, P_rlen, P_head = (i64(pr[k]) for k in ("psi", "obs", "tau", "roff", "rlen", "head"))
        P_ntr, P_tr0, P_tm0, P_ntm = (i64(pr[k]) for k in ("ntr", "tr0", "tm0", "ntm"))
        T0, TR_sub, TR_k, TR_mu = f64(tau0), i64(tr_sub), i64(tr_k), f64(tr_mu)
        TR_ko, TR_muo = i64(tr_ko), f64(tr_muo)
        TM_kind, TM_off, TM_w0, TM_ninv = i64(tm_kind), i64(tm_off), i64(tm_w0), f64(tm_ninv)
        TM_pph, TM_pgd, W = i64(tm_pph), i64(tm_pgd), f64(wflat if wflat else [0.0])
        n_pairs = len(self.pairs)
        phase_on, gd_on = self.toggles["phase"], self.toggles["group_delay"]
        ring_read = _ring_read_nb

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            for p in range(n_pairs):
                ntr, tr0, tm0, ntm = P_ntr[p], P_tr0[p], P_tm0[p], P_ntm[p]
                L, base = P_rlen[p], P_roff[p]
                gd = 0.0
                h = 0
                row = 0
                if phase_on:
                    h = I[P_head[p]] + 1
                    if h == L:
                        h = 0
                    I[P_head[p]] = h
                    row = base + h * ntr
                for j in range(ntr):
                    m = TR_sub[tr0 + j]
                    v = 0.0
                    g = 0.0
                    for t in range(tm0, tm0 + ntm):
                        if TM_kind[t] == 0:
                            x = W[TM_w0[t] + m] * S[TM_off[t]]
                        else:
                            x = S[TM_off[t] + m] * TM_ninv[t]
                        v += P[TM_pph[t]] * x
                        g += P[TM_pgd[t]] * x
                    if phase_on:
                        S[row + j] = v
                    gd += g
                psi = 0.0
                obs = 0.0
                if phase_on:
                    for j in range(ntr):
                        psi += ring_read(S, base, ntr, L, h, j, TR_k[tr0 + j], TR_mu[tr0 + j])
                        obs += ring_read(S, base, ntr, L, h, j, TR_ko[tr0 + j], TR_muo[tr0 + j])
                S[P_psi[p]] = psi
                S[P_obs[p]] = obs
                S[P_tau[p]] = T0[p] + gd if gd_on else T0[p]

        return kernel

    # ------------------------------------------------------------------------------- reference

    def __call__(self, state, params):
        phase_on, gd_on = self.toggles["phase"], self.toggles["group_delay"]
        for p in self.pairs:
            n = len(p.subs)
            dens = np.zeros(n)
            gdv = np.zeros(n)
            for d in p.drivers:
                val = np.asarray(state[f"drv.{p.seg}.{d}"], dtype=float).reshape(-1)
                x = p.weights[d][p.subs] * val[0] if p.weights[d] is not None else val[p.subs] / p.n_sub
                dens += float(params[self._pname(p, d, "phase")]) * x
                gdv += float(params[self._pname(p, d, "group_delay")]) * x
            psi = obs = 0.0
            if phase_on:
                ring = state[p.ring]  # (ring_len, n) view
                head = state.layout.slot(p.ring).head
                h = (int(state.I[head]) + 1) % p.ring_len
                state.I[head] = h
                ring[h] = dens
                rows0 = (h - p.k) % p.ring_len
                rows1 = (h - p.k - 1) % p.ring_len
                cols = np.arange(n)
                y0, y1 = ring[rows0, cols], ring[rows1, cols]
                psi = float(np.sum(np.where(p.mu != 0.0, (1.0 - p.mu) * y0 + p.mu * y1, y0)))
                y0, y1 = ring[(h - p.ko) % p.ring_len, cols], ring[(h - p.ko - 1) % p.ring_len, cols]
                obs = float(np.sum(np.where(p.muo != 0.0, (1.0 - p.muo) * y0 + p.muo * y1, y0)))
            state[p.psi] = psi
            state[p.obs] = obs
            state[p.tau] = p.tau0 + float(np.sum(gdv)) if gd_on else p.tau0
        return state


@register_block("segment_phase")
def make_segment_phase(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``segment_phase`` block (module docstring)."""
    return SegmentPhase(params, name)
