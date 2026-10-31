"""``beam_assembler`` block: beam phase at its observer with emission-time evaluation (PLAN §2.3.2,
conv §7, §9). Owner: WP-F6.

    φ_b(t) = φ_L(t − τ_b) + Σ_{(i,m) ∈ path b} ψ_{i,m,b}(t − τ_{b,after(i,m)})

τ_b: group delay source → observer; τ_{b,after(i,m)}: from the exit of subsegment (i, m) to the
observer. Counter-propagating beams (d = −1) read the same subsegments at other times. Histories are
kept in rings owned by ``segment_phase`` (per-subsegment ψ) and by this block (source phase, segment
phases), read with fractional interpolation (``core.delay`` arithmetic). The
frequency dependence of the propagation phase is NOT added on top (no double counting, conv §7).
Jones matrices (Jᵀ for d = −1) multiply the input Jones vector (G04). Writes ``beam.<id>.*``.

Implementation (WP-F6)
    The delay τ_{b,after(i,m)} is split at the segment exit: ``segment_phase`` already evaluates
    ``psi.<seg>.<beam>`` = Σ_m ψ_{i,m,b} at the times the beam passes each subsegment relative to the
    moment it leaves the segment; this block delays that value by τ_after(i) (segment exit → observer)
    and the source phase by τ_b, both with linear fractional interpolation on rings it owns (nominal
    delays from ``Topology.timing``; ``core.delay`` arithmetic). Per beam b and sample t:

        φ_b(t) = φ_src(t − τ_b) + 2π·f_off·τ_b + Σ_steps psi.<seg>.<b>(t − τ_after) [+ arg t_pol]

    * ``2π f_off τ_b``: a beam shifted by f_off from its source carrier is written in its own rotating
      frame (the carrier term −2π f_off t of conv §4 belongs to the demodulating error block); the
      emission time leaves the static phase 2π f_off τ_b (conv §4, §7; T05 V9). Its drift with
      temperature is a phase coefficient of the beam (T10 §3.3), not a moving delay.
    * Polarization (G04, T09 §1.6): for a beam with a ``jones`` segment on its path, e_out = M_n⋯M_1
      e_in with M = U (d = +1) or Uᵀ (d = −1), U = static ``pol_params.J`` or the slot ``J.<seg>``
      (current sample; polarization drift is slow); t_pol = jones_out^† e_out; the phase gains arg t_pol
      (wrapped, undefined when |t_pol| → 0) and the power |t_pol|². Segments with pol_model ``off`` /
      ``pm_diagonal`` have U = 1 (their per-axis phases are in ψ, ``segment_phase``).
      A jones segment with ``J_nr_slot`` (non-reciprocal factor N, WP-B7 addition) uses M = N·U for
      d = +1 and M = Uᵀ·N for d = −1 (N is not transposed: Faraday, T09 §1.8).

Slots (WP-F6):

    reads   <source slot>            rad   from params ``sources`` (e.g. laser.dx1.phi, pump.p.phi)
            psi.<seg>.<beam>         rad   for every path step (segment_phase / opo_angle_map / ...)
            J.<seg>                  1     complex (2, 2), jones segments with J_slot (polarization on)
            <J_nr_slot>              1     complex (2, 2), non-reciprocal factor N (polarization on)
    writes  beam.<b>.phi             rad   phase at the observer (unwrapped except arg t_pol)
            beam.<b>.power           W     power_w·|t_pol|² (static losses: LossRegistry, G06)
            beam.<b>.jones           1     complex (2,): output Jones vector e_out
            beam.<b>.t_emit          s     emission time t − τ_b
            beam.<b>.hist_<seg>      rad   ring of psi.<seg>.<b> (τ_after delay)
            beam.<block>.src_<src>   rad   ring of the source phase (emission-time history)

Params (all required, D34)::

    topology : mapping (RunSpec.topology form) or Topology
    beams    : "all" | [beam id, ...]
    sources  : {source id: slot name}  (every source of the selected beams)
    toggles  : {emission_time: bool, polarization: bool}
      emission_time off → τ_b = τ_after = 0 (instantaneous propagation, f_off phase 0, t_emit = t)
      polarization off  → J = e^{iψ}·1 on every segment: no matrix product, power = power_w,
                          jones = jones_in (with segment_phase.pm_axis_drift off: T09 test 8)
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
from cranberry.path.segments import frac_delay, select_ids
from cranberry.path.topology import as_topology

__all__ = ["BeamAssembler", "make_beam_assembler"]


def _interp(S, off, L, h, k, mu):  # noqa: N803
    a = h - k
    if a < 0:
        a += L
    y = S[off + a]
    if mu != 0.0:
        b = a - 1
        if b < 0:
            b += L
        y = (1.0 - mu) * y + mu * S[off + b]
    return y


_interp_nb = numba.njit(inline="always")(_interp)


class BeamAssembler(BaseBlock):
    """The ``beam_assembler`` block (module docstring)."""

    toggle_names = ("emission_time", "polarization")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or "beam_assembler"
        require(params, "topology", "beams", "sources", "toggles", block=name)
        topo = as_topology(params["topology"])
        self.topo = topo
        self.beam_ids = select_ids(params["beams"], topo.beams, "beams", name)
        sources = params["sources"]
        if not isinstance(sources, Mapping):
            raise ValueError(f"block {name!r}: sources must be a mapping {{source id: slot name}}")
        used = list(dict.fromkeys(topo.beams[b].source for b in self.beam_ids))
        missing = [s for s in used if s not in sources]
        if missing:
            raise ValueError(f"block {name!r}: no slot given for source(s) {missing} in params['sources']")
        self.sources = {s: str(sources[s]) for s in used}
        pol_on = bool(dict(params["toggles"]).get("polarization", False))
        reads = list(self.sources.values())
        writes: list[str] = []
        for b in self.beam_ids:
            beam = topo.beams[b]
            reads += [f"psi.{s.segment}.{b}" for s in beam.path]
            if pol_on:
                reads += [topo.segments[s.segment].pol_params[k] for s in beam.path for k in ("J_slot", "J_nr_slot")
                          if k in topo.segments[s.segment].pol_params]
            writes += [f"beam.{b}.{f}" for f in ("phi", "power", "jones", "t_emit")]
            writes += [f"beam.{b}.hist_{s.segment}" for s in beam.path]
        writes += [f"beam.{name}.src_{s}" for s in self.sources]
        super().__init__("beam_assembler", name, params, reads=list(dict.fromkeys(reads)), writes=writes)
        self.slot_units = {s: "rad" for s in self.sources.values()}

    def _src_ring(self, src: str) -> str:
        return f"beam.{self.name}.src_{src}"

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        h = block_dt_s(layout.f_sim_hz, self.rate_div)
        emit = self.toggles["emission_time"]
        self.pol_on = self.toggles["polarization"]
        src_len = dict.fromkeys(self.sources, 2)
        self.bd = []
        for b in self.beam_ids:
            beam = self.topo.beams[b]
            tm = self.topo.timing(b)
            tau = tm.tau_s if emit else 0.0
            k, mu = frac_delay(tau, h)
            src_len[beam.source] = max(src_len[beam.source], k + 2)
            steps = []
            for st, ps in zip(tm.steps, beam.path):
                seg = self.topo.segments[st.segment]
                ka, mua = frac_delay(st.after_s if emit else 0.0, h)
                if seg.pol_model != "jones":
                    kind, mat = 0, None
                elif "J" in seg.pol_params:
                    kind, mat = 1, np.asarray(seg.pol_params["J"], dtype=np.complex128)
                else:
                    kind, mat = 2, seg.pol_params["J_slot"]
                nr = seg.pol_params.get("J_nr_slot") if seg.pol_model == "jones" else None
                steps.append({"seg": st.segment, "k": ka, "mu": mua, "ring_len": ka + 2, "kind": kind, "mat": mat, "nr": nr,
                              "T": ps.direction == -1, "psi": f"psi.{st.segment}.{b}", "ring": f"beam.{b}.hist_{st.segment}"})
            has_pol = self.pol_on and any(s["kind"] != 0 for s in steps)
            self.bd.append({"id": b, "src": beam.source, "k": k, "mu": mu, "tau": tau,
                            "static": 2.0 * math.pi * beam.freq_offset_hz * tau, "steps": steps, "pol": has_pol,
                            "jin": np.asarray(beam.jones_in, dtype=np.complex128),
                            "jout": None if beam.jones_out is None else np.asarray(beam.jones_out, dtype=np.complex128),
                            "p0": beam.power_w})
            layout.add(f"beam.{b}.phi", (), np.float64, "rad")
            layout.add(f"beam.{b}.power", (), np.float64, "W", init=beam.power_w)
            layout.add(f"beam.{b}.jones", (2,), np.complex128, "1")
            layout.add(f"beam.{b}.t_emit", (), np.float64, "s", init=-tau)
            for s in steps:
                layout.add(s["ring"], (), np.float64, "rad", ring_len=s["ring_len"])
        self.src_len = src_len
        for s, L in src_len.items():
            layout.add(self._src_ring(s), (), np.float64, "rad", ring_len=L)

    def bind(self, layout) -> None:
        for s in self.sources.values():
            sl = layout.slot(s)
            if sl.array != "S" or sl.shape != () or sl.ring_len:
                raise TypeError(f"block {self.name!r}: source slot {s!r} must be a scalar float64 slot")
        self.src_off = {s: layout.slot(slot).offset for s, slot in self.sources.items()}
        self.src_ring = {s: layout.slot(self._src_ring(s)) for s in self.sources}
        self.t_off = layout.slot("t").offset
        for bd in self.bd:
            for nm in ("phi", "power", "jones", "t_emit"):
                bd[nm] = layout.slot(f"beam.{bd['id']}.{nm}").offset
            for s in bd["steps"]:
                s["psi_off"] = layout.slot(s["psi"]).offset
                r = layout.slot(s["ring"])
                s["ring_off"], s["head"] = r.offset, r.head
                if s["kind"] == 2 and bd["pol"]:
                    sl = layout.slot(s["mat"])
                    if sl.array != "Z" or sl.shape != (2, 2) or sl.ring_len:
                        raise TypeError(f"block {self.name!r}: {s['mat']!r} must be a complex (2, 2) slot")
                    s["mat_off"] = sl.offset
                s["nr_off"] = -1
                if s["nr"] is not None and bd["pol"]:
                    sl = layout.slot(s["nr"])
                    if sl.array != "Z" or sl.shape != (2, 2) or sl.ring_len:
                        raise TypeError(f"block {self.name!r}: {s['nr']!r} must be a complex (2, 2) slot")
                    s["nr_off"] = sl.offset
        self.kernel = self._build_kernel()

    def _build_kernel(self):
        srcs = list(self.sources)
        i64 = lambda x: np.array(x, dtype=np.int64)  # noqa: E731
        f64 = lambda x: np.array(x, dtype=np.float64)  # noqa: E731
        Q_off = i64([self.src_off[s] for s in srcs])
        Q_roff = i64([self.src_ring[s].offset for s in srcs])
        Q_rlen = i64([self.src_ring[s].ring_len for s in srcs])
        Q_head = i64([self.src_ring[s].head for s in srcs])
        B = {k: [] for k in ("src", "k", "mu", "static", "tau", "phi", "pow", "jones", "temit", "p0", "st0", "nst", "pol")}
        jin, jout = [], []
        st = {k: [] for k in ("psi", "roff", "rlen", "head", "k", "mu", "kind", "moff", "T", "nr")}
        mats = [0j]
        for bd in self.bd:
            B["src"].append(srcs.index(bd["src"]))
            for k in ("k", "mu", "static", "tau", "p0"):
                B[k].append(bd[k])
            B["phi"].append(bd["phi"])
            B["pow"].append(bd["power"])
            B["jones"].append(bd["jones"])
            B["temit"].append(bd["t_emit"])
            B["st0"].append(len(st["psi"]))
            B["nst"].append(len(bd["steps"]))
            B["pol"].append(1 if bd["pol"] else 0)
            jin += list(bd["jin"])
            jout += list(np.conj(bd["jout"])) if bd["jout"] is not None else [0j, 0j]
            for s in bd["steps"]:
                st["psi"].append(s["psi_off"])
                st["roff"].append(s["ring_off"])
                st["rlen"].append(s["ring_len"])
                st["head"].append(s["head"])
                st["k"].append(s["k"])
                st["mu"].append(s["mu"])
                kind = s["kind"] if bd["pol"] else 0
                st["kind"].append(kind)
                if kind == 1:
                    st["moff"].append(len(mats))
                    mats += list(np.asarray(s["mat"]).reshape(-1))
                elif kind == 2:
                    st["moff"].append(s["mat_off"])
                else:
                    st["moff"].append(-1)
                st["T"].append(1 if s["T"] else 0)
                st["nr"].append(s["nr_off"] if kind != 0 else -1)
        B_src, B_k, B_phi, B_pow, B_jones, B_temit, B_st0, B_nst, B_pol = (
            i64(B[k]) for k in ("src", "k", "phi", "pow", "jones", "temit", "st0", "nst", "pol"))
        B_mu, B_static, B_tau, B_p0 = (f64(B[k]) for k in ("mu", "static", "tau", "p0"))
        JIN = np.array(jin if jin else [0j], dtype=np.complex128)
        JOUT = np.array(jout if jout else [0j], dtype=np.complex128)
        ST_psi, ST_roff, ST_rlen, ST_head, ST_k, ST_kind, ST_moff, ST_T, ST_nr = (
            i64(st[k] if st[k] else [0]) for k in ("psi", "roff", "rlen", "head", "k", "kind", "moff", "T", "nr"))
        ST_mu = f64(st["mu"] if st["mu"] else [0.0])
        MATS = np.array(mats, dtype=np.complex128)
        n_src, n_beam, t_off = len(srcs), len(self.bd), self.t_off
        interp = _interp_nb

        @numba.njit
        def push(S, I, off, L, head, x):  # noqa: N803
            h = I[head] + 1
            if h == L:
                h = 0
            I[head] = h
            S[off + h] = x
            return h

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            for q in range(n_src):
                push(S, I, Q_roff[q], Q_rlen[q], Q_head[q], S[Q_off[q]])
            for b in range(n_beam):
                q = B_src[b]
                phi = interp(S, Q_roff[q], Q_rlen[q], I[Q_head[q]], B_k[b], B_mu[b]) + B_static[b]
                e0 = JIN[2 * b]
                e1 = JIN[2 * b + 1]
                for s in range(B_st0[b], B_st0[b] + B_nst[b]):
                    h = push(S, I, ST_roff[s], ST_rlen[s], ST_head[s], S[ST_psi[s]])
                    phi += interp(S, ST_roff[s], ST_rlen[s], h, ST_k[s], ST_mu[s])
                    kind = ST_kind[s]
                    if kind != 0:
                        if kind == 1:
                            o = ST_moff[s]
                            m00, m01, m10, m11 = MATS[o], MATS[o + 1], MATS[o + 2], MATS[o + 3]
                        else:
                            o = ST_moff[s]
                            m00, m01, m10, m11 = Z[o], Z[o + 1], Z[o + 2], Z[o + 3]
                        if ST_T[s] == 1:
                            m01, m10 = m10, m01
                        o = ST_nr[s]
                        if o >= 0 and ST_T[s] == 1:  # backward: N first (not transposed)
                            e0, e1 = Z[o] * e0 + Z[o + 1] * e1, Z[o + 2] * e0 + Z[o + 3] * e1
                        e0, e1 = m00 * e0 + m01 * e1, m10 * e0 + m11 * e1
                        if o >= 0 and ST_T[s] == 0:  # forward: N at the exit
                            e0, e1 = Z[o] * e0 + Z[o + 1] * e1, Z[o + 2] * e0 + Z[o + 3] * e1
                p = B_p0[b]
                if B_pol[b] == 1:
                    t = JOUT[2 * b] * e0 + JOUT[2 * b + 1] * e1
                    phi += math.atan2(t.imag, t.real)
                    p = p * (t.real * t.real + t.imag * t.imag)
                S[B_phi[b]] = phi
                S[B_pow[b]] = p
                Z[B_jones[b]] = e0
                Z[B_jones[b] + 1] = e1
                S[B_temit[b]] = S[t_off] - B_tau[b]

        return kernel

    def __call__(self, state, params):
        S, I = state.S, state.I
        heads = {}
        for s, sl in self.src_ring.items():
            h = (int(I[sl.head]) + 1) % sl.ring_len
            I[sl.head] = h
            S[sl.offset + h] = S[self.src_off[s]]
            heads[s] = h
        for bd in self.bd:
            sl = self.src_ring[bd["src"]]
            phi = _interp(S, sl.offset, sl.ring_len, heads[bd["src"]], bd["k"], bd["mu"]) + bd["static"]
            e = bd["jin"].copy()
            for s in bd["steps"]:
                h = (int(I[s["head"]]) + 1) % s["ring_len"]
                I[s["head"]] = h
                S[s["ring_off"] + h] = S[s["psi_off"]]
                phi += _interp(S, s["ring_off"], s["ring_len"], h, s["k"], s["mu"])
                if bd["pol"] and s["kind"] != 0:
                    U = s["mat"] if s["kind"] == 1 else np.array(state[s["mat"]])
                    N = np.array(state[s["nr"]]) if s["nr_off"] >= 0 else np.eye(2)
                    e = U.T @ (N @ e) if s["T"] else N @ (U @ e)
            p = bd["p0"]
            if bd["pol"]:
                t = complex(np.vdot(bd["jout"], e))
                phi += math.atan2(t.imag, t.real)
                p = p * abs(t) ** 2
            S[bd["phi"]] = phi
            S[bd["power"]] = p
            state.Z[bd["jones"]:bd["jones"] + 2] = e
            S[bd["t_emit"]] = state.t - bd["tau"]
        return state


@register_block("beam_assembler")
def make_beam_assembler(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``beam_assembler`` block (module docstring)."""
    return BeamAssembler(params, name)
