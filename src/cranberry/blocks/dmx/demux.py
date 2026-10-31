"""Time → space demultiplexing: ``demux_design`` (O), ``demux_switch`` (T), ``delay_line_k`` (O + T).

WP-B8; T11 (``notes/time-demux-switch.md``), D08, D39, D51-6; PLAN §4.3 DMX / CHK ``demux_feasibility``.

``demux_design`` (offline table + loss registration + feasibility CHK)
    Port delays (T11 §3.1, D39; guard slots replace the period N_b·T_rep by (N_b + n_guard)·T_rep):
    L_k = L_0 + k·(N_b + n_guard)·c·T_rep/n_g, τ_k = n_g L_k/c (port k receives the k-th block from the
    end, σ(k) = K − 1 − k, so every target pulse reaches the GBS at the same time).
    Switch network (T11 §3.2–§3.3), per-stage pass transmission t and wrong-port fraction ε = 10^(−ER/10)
    (``crosstalk``) [+ the rise-time leak of the pulse nearest an edge, ``rise_time_clipping``]:

    * ``tree`` (binary, S = ⌈log2 K⌉ stages for every port). All switches of one depth share a drive,
      so light that leaks at a set of stages lands on the port whose address differs in those bits:
      T[k′, k] = t^S ε^{h} (1 − ε)^{S − h}, h = popcount(k ⊕ k′) (addresses ≥ K are dead outputs).
      η_sw,k = (t(1 − ε))^S; the total leaked power t^S (1 − (1 − ε)^S) (T11 test 2(c)).
    * ``chain`` (PC + PBS in series, port k leaves at stage k + 1, port K − 1 passes all K − 1 stages).
      At time t_j only the target stage is set to "out", the others to "pass": a stage sends t(1 − ε)
      (target stage) or t·ε (pass stage) out and the rest on.

    T[k′, k] is the power of the slot targeted at port k that reaches port k′. In port k′ it arrives
    (j − j_target(k′))·T_rep away from the target bin, never in it (T11 §3.3): for gated detection the
    leak is a pure loss of the target mode and the leaked light is a stray-light output in non-target
    bins (``leak_in`` = Σ_{k≠k′} T[k′, k], photons per frame = leak_in·sinh²r for equal slots).
    ``switching: passive`` (Liu / Borealis delay loops, D39 (c)): no active stage, ε = 0, no
    feasibility check; t is the per-mode interferometer transmission. ``ln_fiber_tree`` / ``pockels_*``
    are presets of these values (PLAN §4.3), not code paths.
    Insertion loss η_sw,k is registered per port on ``segments[k]`` (G06; never multiplied in here).
    Feasibility (T11 §3.4, ChkLog, never stops a run): edge t_r + 6σ_j + τ_mode + t_ring must fit in the
    gap between two settings, (1 + n_guard)·T_rep (code ``demux_rise_time``), and the first stage toggles
    at f_rep/(2(N_b + n_guard)) ≤ f_tog,max (``demux_toggle_rate``).

``demux_switch`` (time loop, port k's stages; T11 §3.5)
    ψ_sw = Σ_s κ_pm (π δV_s/V_π + δΓ_b,s), written to ``act.<name>.dpsi`` and the segment driver
    ``drv.<segment>.<driver>`` (the demux segment of the squeezed path; θ_sq → θ_sq + ψ_sw).
    Bias δΓ_b,s = ``bias_offset_rad`` + random walk (σ_rw [rad/√s]) per stage (toggle ``bias_drift``);
    it also writes the pass transmission Π_s cos²(δΓ_s/2) (δΓ_s = π δV_s/V_π + δΓ_b,s; the drive noise
    enters at 2nd order) → ``pow.<name>.trans`` and the extra leak 1 − trans → ``leak.<name>.frac``.
    Drive noise (toggle ``drive_phase_noise``): the port samples the drive once per frame T_F, so the
    part above 1/(2T_F) aliases into pulse-to-pulse white jitter σ² = S·κ²π²⟨δV²⟩/V_π² (S stages).
    ``drive_noise_mode: timeseries`` writes one white sample per block step (run the block at the
    frame rate, ``rate_div``); ``declared`` writes none and declares σ² as out-of-band variance
    (``bands()``, added once by the ledger; G05). The random-walk part above the block Nyquist
    (κ²σ_rw²/(2π² f_N)) is declared in both modes.

``delay_line_k`` (port k's delay fiber; T11 §3.6, D51-6)
    Timing error relative to the design, written to ``timing.<name>.dt`` [s]:
    δt_k = δτ_k (``static_mismatch``) + κ_T L_k ⟨ΔT⟩ (``delay_thermal``; mean of the ``in_dT`` slot) +
    k (N_b + n_guard) T_rep (δT_rep/T_rep) (``rep_rate_drift``, T11 §3.4).
    δL_k (dispersion length mismatch, D51-6) is kept apart from δτ_k: GDD mismatch β2·δL_k
    (:attr:`DelayLineK.gdd_mismatch_s2`) for the overlap tables (T03 / WP-B11). Fiber attenuation of
    L_k is registered on ``segment`` (``fiber_loss``; turn it off when B10's quantum map already
    attenuates that segment, D59-4).

Gaussian overlap placeholder (T11 §3.6, superseded by T03 / WP-B11 tables):
η_mm(δt) = exp(−δt²/(4σ_t²)), σ_t = τ_p/(2√(2 ln 2)).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
from scipy.special import erfc, erfinv

from cranberry.blocks.ctrl._kit import DrvTarget, KitBlock, drv_write, drv_write_nb
from cranberry.blocks.fib.tables import _cached
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer
from cranberry.offline.band_ledger import BandDecl
from cranberry.offline.loss_registry import LossElement
from cranberry.offline.tables import Table, TableCache
from cranberry.path.topology import C_M_PER_S

__all__ = ["LAYOUTS", "SWITCHINGS", "DESIGN_COLUMNS", "port_length_m", "stage_counts", "transfer_matrix", "edge_leak",
           "gaussian_overlap", "DemuxDesign", "DemuxSwitch", "DelayLineK", "make_demux_design", "make_demux_switch",
           "make_delay_line_k"]

LAYOUTS = ("tree", "chain")
SWITCHINGS = ("active", "passive")
DESIGN_COLUMNS = ("L_m", "tau_s", "n_stages", "eta_sw", "eta_fiber", "leak_in", "gdd_s2")
_DELAY_KEYS = ("f_rep_hz", "N_b", "n_guard", "L0_m", "n_g")


# --- design formulas ---------------------------------------------------------------------------------

def port_length_m(k, *, f_rep_hz, N_b, n_guard, L0_m, n_g):  # noqa: N803
    """L_k = L_0 + k (N_b + n_guard) c T_rep / n_g (T11 §3.1, D39)."""
    return L0_m + np.asarray(k, dtype=float) * (N_b + n_guard) * C_M_PER_S / (f_rep_hz * n_g)


def stage_counts(K: int, layout: str) -> np.ndarray:  # noqa: N803
    """Stages passed by each port: tree ⌈log2 K⌉, chain min(k + 1, K − 1)."""
    if layout == "tree":
        return np.full(K, math.ceil(math.log2(K)) if K > 1 else 0, dtype=int)
    return np.minimum(np.arange(K) + 1, K - 1)


def transfer_matrix(K: int, layout: str, t: float, eps: float) -> np.ndarray:  # noqa: N803
    """T[k′, k]: power of the slot targeted at port k that reaches port k′ (module docstring)."""
    T = np.zeros((K, K))
    if layout == "tree":
        S = int(stage_counts(K, "tree")[0])
        for k in range(K):
            for kp in range(K):
                h = bin(k ^ kp).count("1")
                T[kp, k] = t**S * eps**h * (1.0 - eps) ** (S - h)
        return T
    for k in range(K):
        p = 1.0
        for s in range(K - 1):  # stage s + 1 sends to port s
            out = t * (1.0 - eps) if s == k else t * eps
            T[s, k] = p * out
            p *= t - out
        T[K - 1, k] = p
    return T


def edge_leak(dt_s, t_rise_s: float):
    """Leak sin²(δΓ/2) of a pulse dt from the centre of an erf edge of 10–90 % rise time t_r (T11 §3.4):
    δΓ = (π/2) erfc(|dt|/σ_e), σ_e = t_r/(2 erfinv(0.8)). 0.5 at the centre."""
    sig = t_rise_s / (2.0 * erfinv(0.8))
    return np.sin(0.25 * np.pi * erfc(np.abs(np.asarray(dt_s, dtype=float)) / sig)) ** 2


def gaussian_overlap(dt_s, tau_p_s: float):
    """η_mm = exp(−δt²/(4σ_t²)), σ_t = τ_p/(2√(2 ln 2)) (T11 §3.6 placeholder for T03)."""
    sig = tau_p_s / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return np.exp(-np.asarray(dt_s, dtype=float) ** 2 / (4.0 * sig * sig))


# --- demux_design ------------------------------------------------------------------------------------

class DemuxDesign(BaseBlock):
    """``demux_design`` (module docstring). Offline: never stepped.

    Params (required, D34): ``K``, ``layout`` (tree | chain), ``switching`` (active | passive),
    ``f_rep_hz``, ``N_b``, ``n_guard``, ``L0_m``, ``n_g``, ``alpha_db_per_m``, ``t_stage``,
    ``stage_gdd_s2``, ``segments`` (K segment ids for the loss elements), ``loss_applies_to``,
    ``loss_source``; active only: ``er_db``, ``t_rise_s``, ``sigma_j_s``, ``t_ring_s``, ``tau_mode_s``,
    ``f_tog_max_hz``. Toggles: ``enabled`` (off ⇒ identity: η = 1, no leak, elements disabled),
    ``insertion_loss`` (t → 1), ``crosstalk`` (ε → 0), ``rise_time_clipping``.
    """

    toggle_names = ("enabled", "insertion_loss", "crosstalk", "rise_time_clipping")
    _ACTIVE = ("er_db", "t_rise_s", "sigma_j_s", "t_ring_s", "tau_mode_s", "f_tog_max_hz")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or "demux_design"
        require(params, "K", "layout", "switching", *_DELAY_KEYS, "alpha_db_per_m", "t_stage", "stage_gdd_s2",
                "segments", "loss_applies_to", "loss_source", block=name)
        if params["layout"] not in LAYOUTS or params["switching"] not in SWITCHINGS:
            raise ValueError(f"block {name!r}: layout must be one of {LAYOUTS}, switching one of {SWITCHINGS}")
        if params["switching"] == "active":
            require(params, *self._ACTIVE, block=name)
        K = int(params["K"])
        if K < 1 or len(params["segments"]) != K or not 0.0 < float(params["t_stage"]) <= 1.0:
            raise ValueError(f"block {name!r}: need K >= 1, len(segments) == K and 0 < t_stage <= 1")
        super().__init__("demux_design", name, params)
        self.K = K
        self.chk_log = ChkLog()

    def __call__(self, state, params):
        return state  # offline block

    @property
    def t_rep_s(self) -> float:
        return 1.0 / float(self.params["f_rep_hz"])

    def stage_eps(self) -> float:
        """Effective per-stage wrong-port fraction ε (0 when off / passive)."""
        p, tg = self.params, self._toggles
        if not tg["enabled"] or p["switching"] == "passive":
            return 0.0
        eps = 10.0 ** (-float(p["er_db"]) / 10.0) if tg["crosstalk"] else 0.0
        if tg["rise_time_clipping"]:
            # ponytail: worst case (pulse next to an edge, centred gap) applied to every pulse; a per-pulse
            # table if block routing (N_b > 1) needs the interior pulses separately.
            eps += float(edge_leak(0.5 * (1 + int(p["n_guard"])) * self.t_rep_s, float(p["t_rise_s"])))
        return eps

    def stage_t(self) -> float:
        on = self._toggles["enabled"] and self._toggles["insertion_loss"]
        return float(self.params["t_stage"]) if on else 1.0

    def matrix(self) -> np.ndarray:
        """Transfer matrix T[k′, k] (identity when off)."""
        if not self._toggles["enabled"]:
            return np.eye(self.K)
        return transfer_matrix(self.K, self.params["layout"], self.stage_t(), self.stage_eps())

    def design(self) -> dict[str, np.ndarray]:
        """Per-port columns :data:`DESIGN_COLUMNS`."""
        p = self.params
        k = np.arange(self.K)
        L = port_length_m(k, **{x: p[x] for x in _DELAY_KEYS})
        T = self.matrix()
        n = stage_counts(self.K, p["layout"]) if self._toggles["enabled"] else np.zeros(self.K, dtype=int)
        return {"L_m": L, "tau_s": float(p["n_g"]) * L / C_M_PER_S, "n_stages": n.astype(float),
                "eta_sw": np.diag(T).copy(), "eta_fiber": 10.0 ** (-float(p["alpha_db_per_m"]) * L / 10.0),
                "leak_in": T.sum(axis=1) - np.diag(T), "gdd_s2": n * float(p["stage_gdd_s2"])}

    def tables(self, cache: TableCache | None = None) -> dict[str, Table]:
        def build():
            d = self.design()
            vals = np.stack([d[c] for c in DESIGN_COLUMNS], -1)
            return {"demux_design": Table("demux_design", {"k": np.arange(self.K, dtype=float)}, vals, "mixed",
                                          {"columns": list(DESIGN_COLUMNS), "sources": ["T11 §3.1–§3.3", "D39"]})}

        return _cached("demux_design", self.name, self.params, cache, ("demux_design",), build)

    def losses(self) -> list[LossElement]:
        p, eta = self.params, np.diag(self.matrix())
        return [LossElement(id=f"{self.name}.port{k}", eta=float(eta[k]), kind="component",
                            applies_to=frozenset(p["loss_applies_to"]), source=str(p["loss_source"]),
                            enabled=self._toggles["enabled"], segment=str(p["segments"][k])) for k in range(self.K)]

    def check(self, log: ChkLog | None = None) -> dict[str, float]:
        """CHK ``demux_feasibility`` (T11 §3.4): warns, never raises."""
        p, log = self.params, log or self.chk_log
        if p["switching"] == "passive" or not self._toggles["enabled"]:
            return {}
        per = int(p["N_b"]) + int(p["n_guard"])
        edge = float(p["t_rise_s"]) + 6.0 * float(p["sigma_j_s"]) + float(p["tau_mode_s"]) + float(p["t_ring_s"])
        gap = (1 + int(p["n_guard"])) * self.t_rep_s
        f_tog = float(p["f_rep_hz"]) / (2.0 * per)
        if edge > gap:
            log.warn("demux_rise_time", self.name, f"edge + jitter + mode + ringing = {edge:.3g} s > gap "
                     f"(1 + n_guard)·T_rep = {gap:.3g} s: pulses fall on switching edges", edge_s=edge, gap_s=gap)
        if f_tog > float(p["f_tog_max_hz"]):
            log.warn("demux_toggle_rate", self.name, f"first stage toggles at f_rep/(2(N_b + n_guard)) = {f_tog:.3g} Hz"
                     f" > f_tog,max = {p['f_tog_max_hz']:.3g} Hz", f_tog_hz=f_tog)
        return {"edge_s": edge, "gap_s": gap, "f_tog_hz": f_tog}


# --- demux_switch ------------------------------------------------------------------------------------

def _make_switch_step(write):
    def step(S, I, P, i_out, i_bias, i_tr, i_lk, ns, pk, pvpi, pdv, pb0, prw, sqh, on, bias_on, dn_ts,  # noqa: N803, E741
             buf, idx, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale):
        psi, tr = 0.0, 1.0
        if on:
            for s in range(ns):
                g = 0.0
                if bias_on:
                    r = S[i_bias + s] + P[prw] * sqh * S[buf + I[idx]]
                    I[idx] += 1
                    S[i_bias + s] = r
                    g = P[pb0] + r
                if dn_ts:
                    g += math.pi * P[pdv] / P[pvpi] * S[buf + I[idx]]
                    I[idx] += 1
                psi += P[pk] * g  # g = δΓ_s = π δV/V_π + δΓ_b (T11 §3.5)
                c = math.cos(0.5 * g)  # pass transmission cos²(δΓ/2); drive noise enters at 2nd order
                tr *= c * c
        S[i_out] = psi
        S[i_tr] = tr
        S[i_lk] = 1.0 - tr
        write(S, I, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale, psi)

    return step


class DemuxSwitch(KitBlock):
    """``demux_switch`` (module docstring).

    Params (required, D34): ``segment``, ``driver``, ``n_sub``, ``sub``, ``ring_len`` (segment driver),
    ``n_stages`` (port k's stages, from the design table), ``drive_noise_mode`` (timeseries | declared);
    runtime ``kappa_pm``, ``v_pi_v``, ``drive_noise_v_rms``, ``bias_offset_rad``, ``bias_rw_rad_per_rts``.
    Toggles: ``enabled`` (off ⇒ ψ = 0, trans 1, leak 0), ``drive_phase_noise``, ``bias_drift``.
    ``bands()`` uses the spec values (noise magnitudes are not changed while running, D57-5).
    """

    KIND = "demux_switch"
    STEP = staticmethod(_make_switch_step(drv_write))
    STEP_NB = numba.njit(_make_switch_step(drv_write_nb))
    toggle_names = ("enabled", "drive_phase_noise", "bias_drift")
    RUNTIME = ("kappa_pm", "v_pi_v", "drive_noise_v_rms", "bias_offset_rad", "bias_rw_rad_per_rts")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "n_stages", "drive_noise_mode", block=name)
        if params["drive_noise_mode"] not in ("timeseries", "declared") or int(params["n_stages"]) < 1:
            raise ValueError(f"block {name!r}: drive_noise_mode must be 'timeseries' | 'declared', n_stages >= 1")
        self.drv = DrvTarget(params, name, "rad", split=True)
        self.ns = int(params["n_stages"])
        self.a = {k: f"{g}.{name}.{k}" for g, k in (("act", "dpsi"), ("act", "bias"), ("pow", "trans"), ("leak", "frac"))}
        super().__init__(params, name, reads=(), writes=(*self.a.values(), self.drv.slot), runtime=self.RUNTIME)
        tg = self.toggles
        self._bias = tg["enabled"] and tg["bias_drift"]
        self._dn_ts = tg["enabled"] and tg["drive_phase_noise"] and params["drive_noise_mode"] == "timeseries"
        self._h = None

    def init(self, layout, rng):
        self.drv.check_rate(self.rate_div)
        super().init(layout, rng)
        self._h = self.h
        per = self.ns * (int(self._bias) + int(self._dn_ts))
        self.draws = DrawBuffer(self.name, per, "standard_normal") if per else None
        if self.draws is not None:
            self.draws.init(layout, rng, self.rate_div)
            self.writes = (*self.writes, *self.draws.slot_names)

    def slots(self):
        a = self.a
        return [(a["dpsi"], (), np.float64, "rad", 0, 0.0), (a["bias"], (self.ns,), np.float64, "rad", 0, 0.0),
                (a["trans"], (), np.float64, "1", 0, 1.0), (a["frac"], (), np.float64, "1", 0, 0.0), self.drv.slot_def()]

    def make_args(self, layout):
        d = self.draws
        if d is not None:
            d.bind(layout)
        return (self.o(layout, self.a["dpsi"]), self.o(layout, self.a["bias"]), self.o(layout, self.a["trans"]),
                self.o(layout, self.a["frac"]), self.ns, *(self.p(layout, k) for k in self.RUNTIME),
                math.sqrt(self.h), self.toggles["enabled"], self._bias, self._dn_ts,
                -1 if d is None else d.buf_off, -1 if d is None else d.idx_off, *self.drv.args(layout))

    def prepare(self, state, params, n0, n_samples):
        if self.draws is not None:
            self.draws.refill(state, n0, n_samples)

    def sigma_fast(self) -> float:
        """Pulse-to-pulse drive jitter σ = √S κ π δV/V_π (spec values)."""
        p = self.params
        return math.sqrt(self.ns) * abs(float(p["kappa_pm"])) * math.pi * float(p["drive_noise_v_rms"]) / float(p["v_pi_v"])

    def bands(self) -> list[BandDecl]:
        if self._h is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        p, tg = self.params, self.toggles
        f_n = 0.5 / self._h
        var = 0.0
        if tg["enabled"] and tg["drive_phase_noise"] and p["drive_noise_mode"] == "declared":
            var += self.sigma_fast() ** 2
        if self._bias:
            var += self.ns * (float(p["kappa_pm"]) * float(p["bias_rw_rad_per_rts"])) ** 2 / (2.0 * math.pi**2 * f_n)
        return [BandDecl(self.name, "phase", (str(p["segment"]),), 0.0, f_n, "rad", variance_outside=var)]


# --- delay_line_k ------------------------------------------------------------------------------------

def _make_delay_step():
    def step(S, I, P, i_out, i_dt, n_dt, ring, i_head, kt_l, rep_k, pdtau, prel, f_static, f_th, f_rep):  # noqa: N803, E741
        dt = P[pdtau] if f_static else 0.0
        if f_th:
            base = i_dt + (I[i_head] * n_dt if ring > 0 else 0)
            m = 0.0
            for j in range(n_dt):
                m += S[base + j]
            dt += kt_l * m / n_dt
        if f_rep:
            dt += rep_k * P[prel]
        S[i_out] = dt

    return step


class DelayLineK(KitBlock):
    """``delay_line_k`` (module docstring).

    Params (required, D34): ``port`` (k), ``f_rep_hz``, ``N_b``, ``n_guard``, ``L0_m``, ``n_g``,
    ``kappa_t_s_per_m_k`` (κ_T), ``dL_m`` (δL_k), ``beta2_s2_per_m``, ``alpha_db_per_m``, ``in_dT``
    (temperature slot [K], scalar or per subsegment, or null), ``segment``, ``loss_applies_to``,
    ``loss_source``; runtime ``dtau_s`` (δτ_k), ``rep_rel_drift`` (δT_rep/T_rep).
    Toggles: ``static_mismatch``, ``delay_thermal``, ``rep_rate_drift``, ``fiber_loss``.
    """

    KIND = "delay_line_k"
    STEP = staticmethod(_make_delay_step())
    STEP_NB = numba.njit(_make_delay_step())
    toggle_names = ("static_mismatch", "delay_thermal", "rep_rate_drift", "fiber_loss")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "port", *_DELAY_KEYS, "kappa_t_s_per_m_k", "dL_m", "beta2_s2_per_m", "alpha_db_per_m",
                "in_dT", "segment", "loss_applies_to", "loss_source", block=name)
        self.out = f"timing.{name}.dt"
        reads = () if params["in_dT"] is None else (params["in_dT"],)
        super().__init__(params, name, reads=reads, writes=(self.out,), runtime=("dtau_s", "rep_rel_drift"))
        p = self.params
        self.k = int(p["port"])
        self.L_m = float(port_length_m(self.k, **{x: p[x] for x in _DELAY_KEYS}))
        self.tau0_s = float(p["n_g"]) * self.L_m / C_M_PER_S

    @property
    def gdd_mismatch_s2(self) -> float:
        """β2 δL_k (D51-6: dispersion part of the mismatch, apart from δτ_k); 0 when ``static_mismatch`` off."""
        return float(self.params["beta2_s2_per_m"]) * float(self.params["dL_m"]) if self.toggles["static_mismatch"] else 0.0

    def slots(self):
        return [(self.out, (), np.float64, "s", 0, 0.0)]

    def make_args(self, layout):
        p, tg = self.params, self.toggles
        if p["in_dT"] is None:
            i_dt, n_dt, ring, i_head = -1, 1, 0, -1
        else:
            s = layout.slot(p["in_dT"])
            i_dt, n_dt, ring = s.offset, max(1, int(np.prod(s.shape, dtype=np.int64))), int(s.ring_len)
            i_head = s.head if ring > 0 else -1
        rep_k = self.k * (int(p["N_b"]) + int(p["n_guard"])) / float(p["f_rep_hz"])
        return (self.o(layout, self.out), i_dt, n_dt, ring, i_head, float(p["kappa_t_s_per_m_k"]) * self.L_m, rep_k,
                self.p(layout, "dtau_s"), self.p(layout, "rep_rel_drift"), tg["static_mismatch"],
                tg["delay_thermal"] and i_dt >= 0, tg["rep_rate_drift"])

    def losses(self) -> list[LossElement]:
        p = self.params
        eta = 10.0 ** (-float(p["alpha_db_per_m"]) * self.L_m / 10.0)
        return [LossElement(id=f"{self.name}.fiber", eta=eta, kind="fiber", applies_to=frozenset(p["loss_applies_to"]),
                            source=str(p["loss_source"]), enabled=self.toggles["fiber_loss"], segment=str(p["segment"]))]


@register_block("demux_design")
def make_demux_design(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``demux_design`` (:class:`DemuxDesign`); its :meth:`DemuxDesign.check` is the CHK
    ``demux_feasibility`` of PLAN §4.3."""
    return DemuxDesign(params, name)


@register_block("demux_switch")
def make_demux_switch(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``demux_switch`` (:class:`DemuxSwitch`)."""
    return DemuxSwitch(params, name)


@register_block("delay_line_k")
def make_delay_line_k(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``delay_line_k`` (:class:`DelayLineK`)."""
    return DelayLineK(params, name)
