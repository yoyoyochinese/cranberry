"""``segment_jones``: per-segment polarization matrix U written to ``J.<seg>`` (G04, PLAN §4.3 PATH,
conv §9, T09 §1.2–§1.7). Owner: WP-B7 (M1 minimal part, M4 SMF refinements).

Only the polarization matrix U is written. The common phase ψ (driver × coefficient) belongs to
``segment_phase`` (F6 contract: the full segment matrix is e^{iψ}·U). Basis: lab (x, y). For PM
segments x = slow axis, and U is slow-referenced (U = 1 on the slow axis when there is no crosstalk).
This matches the pm_diagonal fast path, where the topology coefficients are the slow-axis values and
ψ_fast = ψ_slow − Δψ_b. ``beam_assembler`` applies U, or Uᵀ for d = −1 (reciprocity).

Models (``segments.<id>.model``; the segment must be ``pol_model: jones`` with a ``J_slot``):

* ``pm``: U = R(θ_out)·diag(1, e^{−iΔψ_b})·R(θ_in), with R(θ) = [[cos θ, −sin θ], [sin θ, cos θ]].
  θ_in / θ_out are the slow-axis misalignments at the segment ends (connector / splice / coupler
  port; κ = sin θ, PER = −10 log10 tan² θ; ``per_budget.angle_from_per_db``). With two coupling
  points this reproduces the two-point interference of T09 §1.5,
  U_xx = cos θ_i cos θ_o − sin θ_i sin θ_o e^{−iΔψ_b}.
  Δψ_b(t) = Δψ_b0 + K_T·L·T_eff(t), where K_T < 0 (D51-2, nominal −2.05 rad/(m·K)) and
  T_eff = Σ_{d ∈ axis_drivers} Σ_m w_{d,m}·drv_d. The weights are the topology placement weights,
  so they match ``segment_phase``. Δψ_b0 is the static birefringent phase mod 2π (β_b·L).
* ``smf``: U = exp(−i δ_T b̂·σ/2)·W(t)·U0 (M4: U0 and the ``drift`` choice).
  U0 (``u0``): the static birefringence of the segment. ``identity``, or ``haar``: one Haar-random
  SU(2) per segment, drawn at init from the order-independent sub-stream ``u0.<seg>`` of the block's
  RNG stream (fixed seed ⇒ fixed fiber; T10 test 8, T09 test 2).
  W (``drift``):
    ``diffusion``: random walk on SU(2), each step W ← exp(−i dΩ·σ/2)·W with dΩ_k ~ N(0, D_pol·h), so
      for small angles the output state diffuses on the Poincaré sphere with ⟨Θ²⟩ = 2·D_pol·t
      (one-sided PSD of each Stokes angle D_pol/(2π² f²)). W is re-unitarised to the SU(2) form
      [[p, −q*], [q, p*]] every step (removes the rounding drift).
    ``ou`` (spectrum option): W = exp(−i Ω·σ/2), each Ω_k an Ornstein–Uhlenbeck process of corner
      f_c with Lorentzian one-sided PSD S_Ω(f) = 4σ_Ω²τ/(1 + (2πfτ)²), τ = 1/(2π f_c),
      σ_Ω² = D_pol/(4π f_c). Above f_c this equals ``diffusion`` with the same D_pol; below it the state
      wanders around its mean with rms σ_Ω per Stokes component instead of walking away. Exact
      update Ω ← rΩ + σ_Ω√(1 − r²)·w, r = e^{−h/τ}. ponytail: one Lorentzian; add a tabulated spectrum
      when measured SOP-drift data exists.
  W starts at 1 (Ω = 0): the drift starts from U0 at t = 0. The temperature term is δ_T = k_T·L·T_eff
  about the fixed Stokes axis b̂ (T09 §3.6).
  σ = (diag(1, −1), [[0, 1], [1, 0]], [[0, −i], [i, 0]]) (Gordon–Kogelnik).
  PMD (DGD τ = D_PMD·√L) acts on pulses, not on this CW matrix: block ``pmd_mode_overlap``.

Toggles (T09 §6.2; all off ⇒ U = 1 on every segment, T09 test 8):
    pol_birefringence_drift   PM Δψ_b (static + K_T drift); off ⇒ Δψ_b = 0
    pol_crosstalk             PM end couplings θ_in, θ_out; off ⇒ θ = 0
    pol_smf_drift             SMF polarization (U0, drift W, temperature rotation); off ⇒ U = 1 (no
                              draws). A static-only SMF is u0: haar with D_pol = 0.
PDL (``pdl_element``, ``pol_pdl``), Faraday (``faraday_nonreciprocal``, ``pol_faraday``) and PMD
(``pmd_mode_overlap``, ``pol_pmd``) are separate blocks.

Slots:
    reads   drv.<seg>.<d>        K     for d in axis_drivers (scalar, or (N_sub,) for per_sub)
    writes  <J_slot>             1     complex (2, 2) U (registered here, init 1)
            J.<block>.<seg>.w    1     complex (2, 2) SMF drift state W
            J.<block>.<seg>.u0   1     complex (2, 2) SMF static U0 (set at init)
            J.<block>.<seg>.om   rad   (3,) SMF OU state Ω (drift: ou only)
            rng.<block>.*              draws (3 per SMF segment and step, toggle pol_smf_drift)

Params (all required, D34)::

    topology : mapping or Topology
    segments : {seg id: {model: pm, axis_drivers: [dT], dpsi_b0_rad, K_T_rad_per_m_K,
                         theta_in_rad, theta_out_rad}
                | {model: smf, axis_drivers: [dT], axis_stokes: [1, 0, 0], u0: identity | haar,
                   drift: diffusion | ou, D_pol_rad2_per_s, k_T_rad_per_m_K, f_c_hz (ou only)}}
    toggles  : {pol_birefringence_drift, pol_crosstalk, pol_smf_drift}

The numeric keys are runtime parameters ``<block>.<seg>.<key>`` (sweepable, ``set_params``).
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
from cranberry.core.rng import DrawBuffer, substream
from cranberry.path.topology import as_topology

__all__ = ["SegmentJones", "make_segment_jones", "MODEL_KEYS", "haar_su2"]

KIND = "segment_jones"
#: Runtime keys per model, in the order of the P-offset tables P0..P3 (smf ``ou`` adds f_c_hz as P2).
MODEL_KEYS: dict[str, tuple[str, ...]] = {
    "pm": ("dpsi_b0_rad", "K_T_rad_per_m_K", "theta_in_rad", "theta_out_rad"),
    "smf": ("D_pol_rad2_per_s", "k_T_rad_per_m_K"),
}
_STRUCT = {"pm": ("model", "axis_drivers"), "smf": ("model", "axis_drivers", "axis_stokes", "u0", "drift")}
_CHOICES = {"u0": ("identity", "haar"), "drift": ("diffusion", "ou")}


def _keys(cfg: Mapping[str, Any]) -> tuple[str, ...]:
    """Runtime keys of one segment config (P0..P3 order)."""
    keys = MODEL_KEYS[cfg["model"]]
    return (*keys, "f_c_hz") if cfg.get("drift") == "ou" else keys


def haar_su2(rng: np.random.Generator) -> np.ndarray:
    """Haar-random SU(2): [[a, −b*], [b, a*]] from a uniform point (a, b) on S³."""
    q = rng.standard_normal(4)
    q /= np.linalg.norm(q)
    a, b = complex(q[0], q[1]), complex(q[2], q[3])
    return np.array([[a, -b.conjugate()], [b, a.conjugate()]])


def _build(jit: bool):
    dec = numba.njit if jit else (lambda f: f)

    @dec
    def su2(a, b, c):
        """exp(−i (a, b, c)·σ/2) as (e00, e01, e10, e11)."""
        th = math.sqrt(a * a + b * b + c * c)
        if th == 0.0:
            return 1.0 + 0j, 0j, 0j, 1.0 + 0j
        co, s = math.cos(0.5 * th), math.sin(0.5 * th) / th
        return (complex(co, -s * a), complex(-s * c, -s * b), complex(s * c, -s * b), complex(co, s * a))

    @dec
    def step(S, Z, I, P, MODEL, JOFF, WOFF, LEN, TT0, NTT, TT_OFF, TT_N, TT_W, P0, P1, P2, P3, AX,  # noqa: N803
             DRIFT, U0OFF, OMOFF, bire, xt, smf, h, buf, idx):
        for s in range(MODEL.shape[0]):
            T = 0.0
            for t in range(TT0[s], TT0[s] + NTT[s]):
                acc = 0.0
                for m in range(TT_N[t]):
                    acc += S[TT_OFF[t] + m]
                T += TT_W[t] * acc
            o = JOFF[s]
            if MODEL[s] == 0:
                dpsi = P[P0[s]] + P[P1[s]] * LEN[s] * T if bire else 0.0
                ti = P[P2[s]] if xt else 0.0
                to = P[P3[s]] if xt else 0.0
                f = complex(math.cos(dpsi), -math.sin(dpsi))
                ci, si, co, so = math.cos(ti), math.sin(ti), math.cos(to), math.sin(to)
                Z[o] = co * ci - so * si * f
                Z[o + 1] = -co * si - so * ci * f
                Z[o + 2] = so * ci + co * si * f
                Z[o + 3] = -so * si + co * ci * f
            elif smf:
                w = WOFF[s]
                a, b, c = S[buf + I[idx]], S[buf + I[idx] + 1], S[buf + I[idx] + 2]
                I[idx] += 3
                if DRIFT[s] == 0:  # diffusion: W ← exp(−i dΩ·σ/2)·W, re-unitarised
                    sd = math.sqrt(P[P0[s]] * h)
                    e00, e01, e10, e11 = su2(sd * a, sd * b, sd * c)
                    w00, w10 = Z[w], Z[w + 2]
                    w00, w10 = e00 * w00 + e01 * w10, e10 * w00 + e11 * w10
                    nrm = math.sqrt(abs(w00) ** 2 + abs(w10) ** 2)
                    w00, w10 = w00 / nrm, w10 / nrm
                    w01, w11 = -w10.conjugate(), w00.conjugate()
                else:  # ou: Ω ← rΩ + σ_Ω√(1 − r²)·w, W = exp(−i Ω·σ/2)
                    fc = P[P2[s]]
                    r = math.exp(-2.0 * math.pi * fc * h)
                    sd = math.sqrt(P[P0[s]] / (4.0 * math.pi * fc) * (1.0 - r * r))
                    m = OMOFF[s]
                    S[m] = r * S[m] + sd * a
                    S[m + 1] = r * S[m + 1] + sd * b
                    S[m + 2] = r * S[m + 2] + sd * c
                    w00, w01, w10, w11 = su2(S[m], S[m + 1], S[m + 2])
                Z[w], Z[w + 1], Z[w + 2], Z[w + 3] = w00, w01, w10, w11
                q = U0OFF[s]  # V = W·U0
                v00 = w00 * Z[q] + w01 * Z[q + 2]
                v01 = w00 * Z[q + 1] + w01 * Z[q + 3]
                v10 = w10 * Z[q] + w11 * Z[q + 2]
                v11 = w10 * Z[q + 1] + w11 * Z[q + 3]
                d = P[P1[s]] * LEN[s] * T
                u00, u01, u10, u11 = su2(d * AX[s, 0], d * AX[s, 1], d * AX[s, 2])
                Z[o] = u00 * v00 + u01 * v10
                Z[o + 1] = u00 * v01 + u01 * v11
                Z[o + 2] = u10 * v00 + u11 * v10
                Z[o + 3] = u10 * v01 + u11 * v11
            else:
                Z[o], Z[o + 1], Z[o + 2], Z[o + 3] = 1.0, 0.0, 0.0, 1.0

    return step


_STEP_PY, _STEP_NB = _build(False), _build(True)


class SegmentJones(BaseBlock):
    """The ``segment_jones`` block (module docstring)."""

    toggle_names = ("pol_birefringence_drift", "pol_crosstalk", "pol_smf_drift")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "topology", "segments", "toggles", block=name)
        topo = as_topology(params["topology"])
        segs = params["segments"]
        if not isinstance(segs, Mapping) or not segs:
            raise ValueError(f"block {name!r}: segments must be a non-empty mapping {{seg id: model params}}")
        self.cfg: dict[str, dict[str, Any]] = {}
        reads: list[str] = []
        for sid, cfg in segs.items():
            where = f"block {name!r}, segment {sid!r}"
            if sid not in topo.segments:
                raise ValueError(f"{where}: unknown segment")
            seg = topo.segments[sid]
            if seg.pol_model != "jones" or "J_slot" not in seg.pol_params:
                raise ValueError(f"{where}: needs pol_model jones with a J_slot (got {seg.pol_model!r})")
            require(cfg, "model", block=f"{name}.{sid}")
            model = cfg["model"]
            if model not in MODEL_KEYS:
                raise ValueError(f"{where}: model must be one of {tuple(MODEL_KEYS)}, got {model!r}")
            require(cfg, *_STRUCT[model], block=f"{name}.{sid}")
            need = (*_STRUCT[model], *_keys(cfg))
            require(cfg, *need, block=f"{name}.{sid}")
            extra = set(cfg) - set(need)
            if extra:
                raise ValueError(f"{where}: unknown key(s) {sorted(extra)} for model {model!r}")
            bad = [d for d in cfg["axis_drivers"] if d not in seg.drivers]
            if bad:
                raise ValueError(f"{where}: axis_drivers {bad} not declared on the segment ({sorted(seg.drivers)})")
            if model == "smf":
                ax = np.asarray(cfg["axis_stokes"], dtype=float)
                if ax.shape != (3,) or abs(np.linalg.norm(ax) - 1.0) > 1e-9:
                    raise ValueError(f"{where}: axis_stokes must be a unit 3-vector")
                for k, ok in _CHOICES.items():
                    if cfg[k] not in ok:
                        raise ValueError(f"{where}: {k} must be one of {ok}, got {cfg[k]!r}")
                if cfg["drift"] == "ou" and not float(cfg["f_c_hz"]) > 0:
                    raise ValueError(f"{where}: f_c_hz must be > 0")
            self.cfg[sid] = dict(cfg)
            reads += [f"drv.{sid}.{d}" for d in cfg["axis_drivers"]]
        self.topo = topo
        self.smf_ids = [s for s, c in self.cfg.items() if c["model"] == "smf"]
        writes = [topo.segments[s].pol_params["J_slot"] for s in self.cfg]
        writes += [f"J.{name}.{s}.{f}" for s in self.smf_ids for f in ("w", "u0")]
        writes += [f"J.{name}.{s}.om" for s in self.smf_ids if self.cfg[s]["drift"] == "ou"]
        super().__init__(KIND, name, params, reads=list(dict.fromkeys(reads)), writes=writes)
        self.slot_units = {r: "K" for r in self.reads}
        self.draws: DrawBuffer | None = None

    def _slot(self, sid: str, field: str) -> str:
        return f"J.{self.name}.{sid}.{field}"

    def _pname(self, sid: str, key: str) -> str:
        return f"{self.name}.{sid}.{key}"

    def init(self, layout, rng) -> None:
        super().init(layout, rng)
        eye = (1.0, 0.0, 0.0, 1.0)
        for sid, cfg in self.cfg.items():
            layout.add(self.topo.segments[sid].pol_params["J_slot"], (2, 2), np.complex128, "1", init=eye)
            for key in _keys(cfg):
                layout.add_param(self._pname(sid, key), float(cfg[key]), (), "1")
        #: static U0 per SMF segment (identity unless u0: haar)
        self.u0: dict[str, np.ndarray] = {}
        for sid in self.smf_ids:
            cfg = self.cfg[sid]
            layout.add(self._slot(sid, "w"), (2, 2), np.complex128, "1", init=eye)
            u0 = haar_su2(substream(rng, f"u0.{sid}")) if cfg["u0"] == "haar" else np.eye(2, dtype=complex)
            self.u0[sid] = u0
            layout.add(self._slot(sid, "u0"), (2, 2), np.complex128, "1", init=tuple(u0.ravel()))
            if cfg["drift"] == "ou":
                layout.add(self._slot(sid, "om"), (3,), np.float64, "rad")
        self.h = block_dt_s(layout.f_sim_hz, self.rate_div)
        if self.smf_ids and self.toggles["pol_smf_drift"]:
            self.draws = DrawBuffer(self.name, 3 * len(self.smf_ids), "standard_normal")
            self.draws.init(layout, rng, self.rate_div)
            self.writes = (*self.writes, *self.draws.slot_names)

    def prepare(self, state, params, n0, n_samples):
        if self.draws is not None:
            self.draws.refill(state, n0, n_samples)

    def bind(self, layout) -> None:
        i64 = lambda x: np.array(x, dtype=np.int64)  # noqa: E731
        model, joff, woff, length, tt0, ntt, tt_off, tt_n, tt_w = [], [], [], [], [], [], [], [], []
        drift, u0off, omoff = [], [], []
        poff: list[list[int]] = [[], [], [], []]
        ax = []
        for sid, cfg in self.cfg.items():
            seg = self.topo.segments[sid]
            smf = cfg["model"] == "smf"
            model.append(1 if smf else 0)
            joff.append(layout.slot(seg.pol_params["J_slot"]).offset)
            woff.append(layout.slot(self._slot(sid, "w")).offset if smf else -1)
            u0off.append(layout.slot(self._slot(sid, "u0")).offset if smf else -1)
            drift.append(1 if smf and cfg["drift"] == "ou" else 0)
            omoff.append(layout.slot(self._slot(sid, "om")).offset if drift[-1] else -1)
            length.append(seg.length_m)
            tt0.append(len(tt_off))
            ntt.append(len(cfg["axis_drivers"]))
            for d in cfg["axis_drivers"]:
                sl = layout.slot(f"drv.{sid}.{d}")
                w = seg.weights(d)
                if sl.array != "S" or sl.ring_len:
                    raise TypeError(f"block {self.name!r}: {sl.name!r} must be a float64 slot without ring")
                tt_off.append(sl.offset)
                if w is None:  # per_sub: Σ_m drv[m] / N
                    tt_n.append(seg.n_sub)
                    tt_w.append(1.0 / seg.n_sub)
                else:  # scalar slot with placement weights
                    tt_n.append(1)
                    tt_w.append(float(np.sum(w)))
            keys = _keys(cfg)
            for j in range(4):
                poff[j].append(layout.param(self._pname(sid, keys[j])).offset if j < len(keys) else 0)
            ax.append(cfg.get("axis_stokes", (0.0, 0.0, 0.0)))
        buf, idx = (-1, -1)
        if self.draws is not None:
            self.draws.bind(layout)
            buf, idx = self.draws.buf_off, self.draws.idx_off
        tg = self.toggles
        self._args = (i64(model), i64(joff), i64(woff), np.array(length, dtype=np.float64), i64(tt0), i64(ntt),
                      i64(tt_off), i64(tt_n), np.array(tt_w, dtype=np.float64), *(i64(p) for p in poff),
                      np.array(ax, dtype=np.float64).reshape(-1, 3), i64(drift), i64(u0off), i64(omoff),
                      tg["pol_birefringence_drift"], tg["pol_crosstalk"], tg["pol_smf_drift"], float(self.h),
                      int(buf), int(idx))
        (A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16,
         f0, f1, f2, h, b0, b1) = self._args
        nb = _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, Z, I, P, A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16,
               f0, f1, f2, h, b0, b1)

        self.kernel = kernel

    def __call__(self, state, params):
        _STEP_PY(state.S, state.Z, state.I, params.P, *self._args)
        return state


@register_block(KIND)
def make_segment_jones(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``segment_jones`` block (module docstring)."""
    return SegmentJones(params, name)
