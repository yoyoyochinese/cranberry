"""Phase 6 offline stray-light tables (D85-4; ``notes/lit-madsen2022.md`` §4.2 N2, ``notes/lit-aghaeerad2025.md``
§4.2 N4). Both blocks are offline (never stepped), build one ``Table`` for the record and check on demand
(the runner calls :meth:`check` once, D80-2 / M6-12). Photon numbers use the units of ``stray_light_table``
(D46): ``*_gate`` = photons arriving at the GBS input of port k in one detection gate (both polarizations,
before the detector QE), ``*_mode`` = photons in the squeezed temporal mode (conv §2).

``window_stray_light`` (P2 / P3 / P8: a counter-propagating probe switched off during GBS measurement windows)
    Time origin: the probe injection point = the GBS input of port k (as ``probe_pulser``). Probe power envelope
    a(t) ∈ [0, 1], period T = ``period_s``: on (lock window) between the 50 % points [φ, φ + D·T) (mod T,
    φ = ``phase_s``, D = ``duty_lock``), linear power ramps of length e = ``edge_s`` centred on the switching
    times. When off the probe keeps ε_off = ``eps_off`` of its on power. Distributed backscatter (Rayleigh,
    Brillouin, Raman) from z returns 2n_g z/c later, weighted by the round-trip attenuation e^{−2αz} (T10 eq. (8)),
    so after the probe switches off it keeps arriving for W = 2 n_g L/c = 2τ_k:

        h(t) = ∫_0^W a(t − s) e^{−κs} ds / ∫_0^W e^{−κs} ds,   κ = α c/n_g   (W → 0: h = a; κ = 0: the plain box),
        factor f(t) = ε_off + (1 − ε_off)·h(t),

    N_mech(t) = f(t)·N_mech,CW(P_on) with N_CW from ``scattering.stray_photons`` at the on-state average power
    P_on = ``p_on_avg_w`` (distributed scattering depends on P̄ only, T10 §3.4.4, so the AOM frame gating inside
    the lock window does not change it). The CW value uses the same weight (L_bs = ∫e^{−2αz}dz), so h → 1 in CW;
    the weight shortens the tail's mean delay (window mean ≈ (1 − αL/3)× the plain box: −5 % at 1 km, −14 % at
    3.2 km for 0.5 dB/km; WS-03). A discrete reflection at z_i returns a(t − 2 n_g z_i/c) instead of h (reflectance
    = its effective value seen from z = 0, as in ``stray_photons``): one at the probe extraction point z ≈ L_k keeps
    the full on-state level for the whole W after switch-off (overlay axis ``win_refl_r``, WS-08).
    GBS slots: every ``slot_period_s`` (D46: once per frame) inside the fully-off part of the measurement window
    [φ + D·T + e/2, φ + T − e/2), at the slot centres (j + ½)·Δ from its start. Window average = mean over slots.
    The phase φ shifts the window and the slots together: it changes the table's time axis, not the numbers.

    ponytail: h is evaluated at the slot centre (exact gate average while the gate does not straddle a knot of
    the piecewise-linear h; gate ≲ 20 ns ≪ W ≥ µs) and the mode count uses the same quasi-static factor (the
    envelope changes on µs, the mode is ns). Gate jitter is not modelled; add it with the gate_schedule block.

``forward_ref_stray_table`` (P1: forward bright reference pulses interleaved with the quantum pulses, D85-11)
    A reference pulse of peak power P_pk = ``p_peak_w`` (at the GBS input of port k), rectangular length
    t_ref = ``t_ref_s``, offset f_ref = ``freq_offset_hz`` from the squeezed carrier. Contributions to the data
    slots of port k (toggles):

    * ``mod_leak``: the modulator that carves the pulses from CW leaks ε_mod·P_pk continuously:
      mode Φ_leak·|ũ(f_ref)|², gate Φ_leak·T_gate·F(f_ref) (Φ = P/hν0; as a line of ``stray_photons``). With
      ``cal_frames`` a slow gate in series (AOM / shutter, extinction ε_gate = ``eps_gate``) closes the reference
      path during the data frames: Φ_leak = ε_mod·ε_gate·P_pk/hν0 (WS-02). The gate is open during the calibration
      frame and switches in the frame gaps, so the edges, crosstalk and afterpulses of the calibration-frame pulses
      are not reduced by it. ``guard_slot`` cannot use it: the reference path must pass one pulse one slot
      (10–170 ns) before the data, every frame, so it is open in every frame, and an AOM / shutter (rise 10 ns – ms)
      cannot separate the guard slot from the data slots; there ``eps_gate`` must be 1.
    * ``edges``: pulse edges / tails outside the rectangle, P(Δ) = s·P_pk·e^{−(|Δ| − t_ref/2)/τ_e}
      (τ_e = ``edge_tau_s``, both sides; s = 1 for the port's own pulse, ε_dmx for a leaked one). Mode: quasi-static
      P(Δ)/hν0·|ũ(f_ref)|²; gate: the exact integral of the exponential over the gate times F(f_ref).
      ponytail: the quasi-static mode count is exact only for τ_e ≫ the mode length, and it is not a bound: the tail
      is coherent with the laser, N = |∫√Φ(t) u(t) dt|², which for a symmetric mode (Gaussian) is ≥ the quasi-static
      value (convexity; × exp(τ_mode²/(16 ln2·τ_e²)), e.g. ×2.7 for F_liu20's τ_mode at τ_e = 3 ns) and for the
      asymmetric OPO mode lies on either side (0.9–3.8×; Phase 5 WS-04). Use the overlap on ``mode_time`` when the
      edge term matters.
    * ``demux_crosstalk``: the reference pulses of the other ports leak into port k with ε_dmx = 1/ER
      (``er_dmx_lin``; a YAML ``er_dmx_db`` loads as it, ``config.db``; every other port at ε_dmx, conservative for
      a tree) and arrive in the time bins next to port k's data block (never in a data slot: the offset
      (k − k')(N_b + n_guard)·T_rep is not a multiple of the frame); their edges count here.
    * ``forward_scatter``: forward GAWBS (T07 §4.1, ``scattering.gawbs_variances``) and forward Raman
      (``scattering.raman_scatter_density``) of the continuous leak over the port fibre. P_pk is the power at the
      fibre output (the GBS input) and light scattered at z co-propagates with the leak, so the output carries
      P_out·L (not L_eff; WS-09). The bright pulse's own forward scattering co-moves with it (walk-off ≪ T_rep) and
      stays in the reference bin.
    * ``afterpulse``: a pulse of N = s·P_pk·t_ref/hν0 photons clicks the GBS detector with p = 1 − e^{−QE·N};
      each click afterpulses with probability p_ap (``p_afterpulse``), delay distributed as e^{−t/τ_ap}/τ_ap
      (``tau_afterpulse_s``): clicks in a gate at Δ > 0 after the pulse = p·p_ap·∫_gate e^{−t/τ_ap}/τ_ap dt, summed
      over the pulses listed below. Clicks only (no photons in the mode). A data slot closer than
      ``det_recovery_s`` after the listed pulses is blind with probability 1 − Π(1 − p_i) (column ``blind_prob``:
      an efficiency loss, reported, not a stray photon).

    Placement (data slot j = 0 … N_b − 1 of port k's block; Δ = t_slot − t_pulse, B = N_b + n_guard,
    T_f = K·B·T_rep):

    * ``cal_frames`` (primary, D85-11): every ``t_cal_s`` one whole frame carries reference pulses instead of
      squeezed light. Only the data frames before and after a calibration frame see it, each with weight T_f/t_cal:
      own pulses at Δ = T_f − (N_b − 1 − j)·T_rep (after) and Δ = −(T_f − j·T_rep) (before), leaked pulses at
      Δ = (1 + n_guard + j)·T_rep and Δ = −(n_guard + N_b − j)·T_rep. The leaked-pulse positions are the WORST
      PORT's (the first / last port in emission order, one side each) applied to every port; a middle port k sees
      its nearest leaked pulses (k + 1) blocks away (WS-06, conservative). L_k is unchanged.
      ponytail: only the nearest pulse of each kind on each side; the other pulses of the calibration frame would
      add afterpulses, frame-weighted by T_f/t_cal ≈ 10⁻⁵–10⁻⁸. A detector recovery longer than T_f (TES 10–20 µs,
      T_f = 5–6 µs) blinds the 2–4 data frames after a calibration frame; only the first is counted, so
      ``blind_prob`` is low by at most ×4 on a weight of T_f/t_cal (WS-05).
    * ``guard_slot`` (kept as the bad configuration, D85-11, D34): one reference pulse in the last guard slot before
      each block, every frame. A pulse reaches port k's detector every B·T_rep: the own pulse at Δ = (j + 1)·T_rep,
      the other ports' leaked guard pulses at Δ = (m·B + j + 1)·T_rep, m = 1 … K − 1, and the own pulse of the
      previous frame at m = K (all counted: each leaked one clicks at µW-class P_pk and afterpulses, WS-05), and
      the next port's leaked pulse at Δ = −(B − 1 − j)·T_rep after the slot. These positions hold for every port and
      every n_guard ≥ 1. Needs n_guard ≥ 1, so L_k grows by (N_b + n_guard)/(N_b + n_guard − 1) (×2 at N_b = 1,
      D62-3; the reference adds one guard slot to the preset's n_guard − 1) and the D78 linewidth requirement
      Δν ≤ σ²/(2πRΔτ), Δτ = (K − 1)(N_b + n_guard)·T_rep, tightens by the same factor; :meth:`check` reports it
      (CHK ``forward_ref_guard_slot``).
      ponytail: pulses more than one frame back are dropped (afterpulse weight ≤ e^{−T_f/τ_ap} < 1 % for
      τ_ap ≤ 1 µs). The switching of the demux during the guard slot (which splits a guard-slot pulse between two
      ports) is not modelled.

Mode: both blocks take ``mode`` like ``stray_light_table`` (parametric, or ``{source: <temporal source spec>}`` for
the OPO's optimal mode, M5-03). The axis overlays use the τ_mode Gaussian of F_common's ``stray1``, which is right on
the waveguide presets. On the OPO presets ``stray1`` is re-pointed at the OPO mode (F_bor6 / F_liu20
``mode: !replace {source: ...}``), and the Gaussian stand-in errs in a direction set by the mechanism and the probe
frequency (Phase 5 WS-01). For a line at the squeezed carrier (same-frequency probe Rayleigh, the P1 leak) it is
2.6–2.8× HIGH (|ũ(0)|² = τ_eff: Gaussian 38.3 / 21.3 ns vs OPO mode 14.2 / 7.8 ns for F_bor6 / F_liu20). For a
frequency-offset probe it is far LOW, because the Gaussian spectrum has no tail: at 80 MHz ×914 (F_liu20) and
×2.9×10⁶ (F_bor6) low, at 100 GHz ×37 / ×165 low, so the CHK would pass wrongly. On the OPO presets apply the OPO-mode
overlays after the axis overlays: ``stray/opo_mode_bor6`` / ``stray/opo_mode_liu20`` (winstray1) and
``ref/opo_mode_bor6`` / ``ref/opo_mode_liu20`` (frefstray1) copy the preset's ``pulsed_opo_source`` spec into
``mode`` (a test keeps them equal to stray1's).

Unknown hardware values (reference power, extinctions ε_mod and ε_gate, edge tail, afterpulse probability and delay,
recovery time, probe switch-off extinction and ramp, extraction-point reflectance) are sweep axes of the overlays ``presets/axes/ref/aurora_*.yaml`` and ``presets/axes/stray/window_tail.yaml``,
marked UNVERIFIED (Aurora SI not available, D85-14).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.blocks.fib.scattering import (
    STRAY_COLUMNS,
    _check_mode,
    _check_raman,
    _filter,
    _gawbs_modes,
    _hnu,
    check_gawbs_group,
    check_probe_groups,
    gawbs_variances,
    mode_average,
    mode_weight,
    raman_scatter_density,
    resolve_mode,
    stray_photons,
)
from cranberry.blocks.fib.tables import _cached
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.offline.tables import Table, TableCache
from cranberry.path.topology import C_M_PER_S

__all__ = ["WINDOW_COLUMNS", "FORWARD_REF_COLUMNS", "PLACEMENTS", "window_factor", "window_stray", "forward_ref_stray",
           "make_window_stray_light", "make_forward_ref_stray_table"]

WINDOW_COLUMNS = ("factor", *STRAY_COLUMNS)
FORWARD_REF_COLUMNS = ("clicks_per_gate", "photons_per_mode", "mod_leak_gate", "mod_leak_mode", "edges_gate",
                       "edges_mode", "crosstalk_gate", "crosstalk_mode", "forward_gate", "forward_mode",
                       "afterpulse_gate", "blind_prob")
PLACEMENTS = ("cal_frames", "guard_slot")
_DIST = ("rayleigh", "brillouin", "raman")


# --- window_stray_light --------------------------------------------------------------------------------------------

def _ramp(y, e):
    """Linear power ramp 0 → 1 over [0, e] (step for e = 0)."""
    y = np.asarray(y, dtype=float)
    return (y >= 0).astype(float) if e <= 0 else np.clip(y / e, 0.0, 1.0)


def _ramp_int(y, e):
    """∫_{−∞}^{y} _ramp(s, e) ds."""
    y = np.asarray(y, dtype=float)
    if e <= 0:
        return np.maximum(y, 0.0)
    return np.where(y <= 0, 0.0, np.where(y >= e, y - e / 2, y * y / (2 * e)))


def _lin_int_exp(y, k):
    """∫_0^y (y − s) e^{−ks} ds for y > 0 (0 otherwise): (ky − 1 + e^{−ky})/k², y²/2 at k = 0."""
    y = np.maximum(np.asarray(y, dtype=float), 0.0)
    if k <= 0:
        return y * y / 2
    x = k * y
    small = x < 1e-3  # series: the closed form cancels catastrophically
    return np.where(small, y * y * (0.5 - x / 6 + x * x / 24), (x + np.expm1(-np.where(small, 1.0, x))) / (k * k))


def _ramp_int_exp(y, e, k):
    """∫_0^∞ _ramp(y − s, e) e^{−ks} ds (= :func:`_ramp_int` at k = 0)."""
    if k <= 0:
        return _ramp_int(y, e)
    if e <= 0:
        return -np.expm1(-k * np.maximum(np.asarray(y, dtype=float), 0.0)) / k
    return (_lin_int_exp(y, k) - _lin_int_exp(np.asarray(y, dtype=float) - e, k)) / e


def window_factor(t_s, *, period_s: float, duty_lock: float, phase_s: float, edge_s: float, roundtrip_s: float,
                  delay_s: float = 0.0, kappa_per_s: float = 0.0) -> np.ndarray:
    """h(t) of the module docstring: the probe envelope seen at the GBS input at times ``t_s``, averaged over the
    round trips s ∈ [0, W] (W = ``roundtrip_s``) with weight e^{−κs} (κ = ``kappa_per_s`` = αc/n_g; 0 = plain box);
    W = 0 gives the envelope itself, delayed by ``delay_s`` (a discrete reflection)."""
    t = np.asarray(t_s, dtype=float) - float(delay_s)
    T, w, e, k = float(period_s), float(roundtrip_s), float(edge_s), float(kappa_per_s)
    on = duty_lock * T
    m = np.arange(math.floor((t.min() - w - phase_s) / T) - 2, math.ceil((t.max() - phase_s) / T) + 2)
    s = (phase_s - e / 2 + m * T)[None, :]  # ramp-up starts; ramp-down starts at s + on
    x = t[:, None] - s
    if w <= 0.0:
        return np.sum(_ramp(x, e) - _ramp(x - on, e), axis=1)
    # ∫_0^W g(x − s) e^{−κs} ds = G(x) − e^{−κW} G(x − W), G = ∫_0^∞ (one pulse) e^{−κs} ds
    decay = math.exp(-k * w)
    a = (_ramp_int_exp(x, e, k) - _ramp_int_exp(x - on, e, k)
         - decay * (_ramp_int_exp(x - w, e, k) - _ramp_int_exp(x - w - on, e, k)))
    a[x - w >= on + e] = 0.0  # a pulse wholly before t − W contributes exactly 0 (not a rounding residue)
    norm = w if k <= 0 else -math.expm1(-k * w) / k
    return np.sum(a, axis=1) / norm


def _window_slots(p: Mapping[str, Any]) -> np.ndarray:
    T, d, e, dt = (float(p[k]) for k in ("period_s", "duty_lock", "edge_s", "slot_period_s"))
    start = float(p["phase_s"]) + d * T + e / 2
    n = int(math.floor(((1.0 - d) * T - e) / dt + 1e-9))
    return start + (np.arange(max(n, 0)) + 0.5) * dt


def window_stray(p: Mapping[str, Any], toggles: Mapping[str, bool]) -> tuple[np.ndarray, np.ndarray]:
    """(slot times [s], rows [n_slots, len(WINDOW_COLUMNS)]) of a ``window_stray_light`` param set."""
    t = _window_slots(p)
    timing = {k: float(p[k]) for k in ("period_s", "duty_lock", "phase_s", "edge_s")}
    eps = float(p["eps_off"])
    w = 2.0 * float(p["n_g"]) * float(p["L_m"]) / C_M_PER_S
    rows = np.zeros((t.size, len(WINDOW_COLUMNS)))
    if t.size == 0:
        return t, rows
    p_on = float(p["p_on_avg_w"])
    col = {c: i for i, c in enumerate(WINDOW_COLUMNS)}
    kappa = float(p["alpha_per_m"]) * C_M_PER_S / float(p["n_g"])  # e^{−2αz}, z = c·s/(2n_g)
    f = eps + (1.0 - eps) * window_factor(t, roundtrip_s=w, kappa_per_s=kappa, **timing)
    rows[:, col["factor"]] = f
    cw = stray_photons(p, {**toggles, "stray_reflections": False}, p_on)
    for mech in _DIST:
        for kind in ("gate", "mode"):
            rows[:, col[f"{mech}_{kind}"]] = f * cw[f"{mech}_{kind}"]
    if toggles["stray_reflections"]:
        only = dict.fromkeys(("stray_rayleigh", "stray_brillouin", "stray_raman"), False) | {"stray_reflections": True}
        for refl, z in zip(p["reflections"], p["reflection_z_m"]):
            one = stray_photons({**p, "reflections": [list(refl)]}, only, p_on)
            fr = eps + (1.0 - eps) * window_factor(t, roundtrip_s=0.0, delay_s=2.0 * float(p["n_g"]) * z / C_M_PER_S,
                                                   **timing)
            for kind in ("gate", "mode"):
                rows[:, col[f"reflections_{kind}"]] += fr * one[f"reflections_{kind}"]
    for kind, total in (("gate", "clicks_per_gate"), ("mode", "photons_per_mode")):
        rows[:, col[total]] = sum(rows[:, col[f"{m}_{kind}"]] for m in (*_DIST, "reflections"))
    return t, rows


class _OfflineTable(BaseBlock):
    """Offline table block with an on-demand CHK (shared by both kinds)."""

    kind_name = ""
    table_name = ""

    def __call__(self, state, params):
        return state  # offline block, not a time-loop stage

    def _limit_check(self, log: ChkLog, value: float, what: str) -> None:
        lim = self.params.get("limit")
        if lim is not None and value > float(lim["photons_per_mode_max"]):
            log.warn("stray_light_limit", self.name, f"{what} stray photons per squeezed mode {value:.3g} > limit "
                     f"{float(lim['photons_per_mode_max']):.3g} ({lim['source']})", photons_per_mode=value,
                     limit=float(lim["photons_per_mode_max"]))

    def tables(self, cache: TableCache | None = None) -> dict[str, Table]:
        p = self.params
        if not isinstance(p["mode"], Mapping) or "source" in p["mode"]:
            p = {**p, "mode": dict(resolve_mode(p["mode"]))}  # build a nested source once (arrays hash)
        return _cached(self.kind_name, self.name, p, cache, (self.table_name,), lambda: self._build(p))


class WindowStrayLight(_OfflineTable):
    """``window_stray_light`` (module docstring). Table ``window_stray`` (axis ``t_s`` = slot centre time; columns
    :data:`WINDOW_COLUMNS`).

    Params (required, D34): every ``stray_light_table`` key except ``p_peak_w``, ``duty`` and ``limit`` (``toggles``
    {stray_rayleigh, stray_brillouin, stray_raman, stray_reflections}, ``probe``, ``brillouin``, ``rayleigh``,
    ``L_m``, ``alpha_per_m``, ``gamma_per_w_m``, ``f_R``, ``raman_model``, ``lambda0_m``, ``temperature_k``,
    ``t_gate_s``, ``reflections``, ``gbs_filter``, ``wdm_isolation_lin``, ``mode``); ``reflection_z_m`` (one
    distance per reflection); ``n_g``; ``p_on_avg_w``; ``eps_off``; the gate schedule ``period_s``, ``duty_lock``
    (0 < D < 1), ``phase_s``, ``edge_s`` (≤ both windows); ``slot_period_s``. Optional ``limit``
    {photons_per_mode_max, source, statistic: mean | max}: CHK ``stray_light_limit`` on the window mean or the
    worst slot. :meth:`check` also warns ``window_no_slot`` when no GBS slot fits in the measurement window.
    """

    kind_name, table_name = "window_stray_light", "window_stray"
    toggle_names = ("stray_rayleigh", "stray_brillouin", "stray_raman", "stray_reflections")

    def __init__(self, params, name):
        name = name or self.kind_name
        _check_raman(params, name)
        check_probe_groups(params, name)
        require(params, "L_m", "alpha_per_m", "t_gate_s", "reflections", "reflection_z_m", "gbs_filter",
                "wdm_isolation_lin", "mode", "n_g", "p_on_avg_w", "eps_off", "period_s", "duty_lock", "phase_s",
                "edge_s", "slot_period_s", block=name)
        _check_mode(params["mode"], name)
        T, d, e = float(params["period_s"]), float(params["duty_lock"]), float(params["edge_s"])
        if not (T > 0 and 0 < d < 1 and 0 <= e <= min(d, 1 - d) * T and float(params["slot_period_s"]) > 0):
            raise ValueError(f"block {name!r}: need period_s > 0, 0 < duty_lock < 1, 0 <= edge_s <= both windows, "
                             "slot_period_s > 0")
        if len(params["reflection_z_m"]) != len(params["reflections"]):
            raise ValueError(f"block {name!r}: reflection_z_m must give one distance per reflection")
        if "limit" in params:
            require(params["limit"], "photons_per_mode_max", "source", "statistic", block=f"{name}.limit")
            if params["limit"]["statistic"] not in ("mean", "max"):
                raise ValueError(f"block {name!r}: limit.statistic must be 'mean' or 'max'")
        super().__init__(self.kind_name, name, params)
        self.chk_log = ChkLog()
        self.result: dict[str, float] | None = None

    def _build(self, p):
        t, rows = window_stray(p, self._toggles)
        if t.size == 0:
            raise ValueError(f"block {self.name!r}: no GBS slot fits in the measurement window")
        return {self.table_name: Table(self.table_name, {"t_s": t}, rows, "photons",
                                       {"columns": list(WINDOW_COLUMNS),
                                        "sources": ["T10 eq. (8)", "lit-madsen2022 §4.2 N2", "D85-4"]})}

    def check(self, log: ChkLog | None = None) -> dict[str, float]:
        """Window summary (:attr:`result`) and CHK ``stray_light_limit`` / ``window_no_slot`` (class docstring)."""
        log = log or self.chk_log
        p = self.params
        t, rows = window_stray(p, self._toggles)
        T, d, e = float(p["period_s"]), float(p["duty_lock"]), float(p["edge_s"])
        w = 2.0 * float(p["n_g"]) * float(p["L_m"]) / C_M_PER_S
        t_meas = (1.0 - d) * T - e
        cw = stray_photons(p, self._toggles, float(p["p_on_avg_w"]))
        r = {"n_slots": float(t.size), "tail_s": w, "t_meas_s": t_meas, "meas_fraction": t_meas / T,
             "tail_time_fraction": min(w, t_meas) / t_meas if t_meas > 0 else float("nan"),
             "photons_per_mode_on": cw["photons_per_mode"], "clicks_per_gate_on": cw["clicks_per_gate"]}
        if t.size == 0:
            log.warn("window_no_slot", self.name, f"no GBS slot (every {float(p['slot_period_s']):.3g} s) fits in the "
                     f"measurement window {t_meas:.3g} s", t_meas_s=t_meas)
            r.update(dict.fromkeys(("mean_factor", "photons_per_mode_mean", "photons_per_mode_max",
                                    "clicks_per_gate_mean", "clicks_per_gate_max"), float("nan")))
        else:
            ci, cm = WINDOW_COLUMNS.index("clicks_per_gate"), WINDOW_COLUMNS.index("photons_per_mode")
            r.update(mean_factor=float(rows[:, 0].mean()), photons_per_mode_mean=float(rows[:, cm].mean()),
                     photons_per_mode_max=float(rows[:, cm].max()), clicks_per_gate_mean=float(rows[:, ci].mean()),
                     clicks_per_gate_max=float(rows[:, ci].max()))
            if "limit" in p:
                stat = p["limit"]["statistic"]
                self._limit_check(log, r[f"photons_per_mode_{stat}"], f"window {stat}")
        self.result = r
        return r


# --- forward_ref_stray_table ---------------------------------------------------------------------------------------

def _exp_gate(x, gate, tau):
    """∫ e^{−t/τ}/τ dt over [x − g/2, x + g/2] ∩ [0, ∞) (0 for τ = 0)."""
    if tau <= 0:
        return 0.0
    lo, hi = max(x - gate / 2, 0.0), max(x + gate / 2, 0.0)
    return math.exp(-lo / tau) - math.exp(-hi / tau)


def _events(p: Mapping[str, Any], j: int) -> list[tuple[str, float, float]]:
    """[(kind 'own' | 'xt', Δ = t_slot − t_pulse, frame weight)] of data slot j (module docstring)."""
    tr, nb, ng, K = float(p["t_rep_s"]), int(p["N_b"]), int(p["n_guard"]), int(p["K"])
    b = nb + ng
    if p["placement"] == "guard_slot":
        return [("own", (j + 1) * tr, 1.0), ("xt", -(b - 1 - j) * tr, 1.0),
                *(("own" if m == K else "xt", (m * b + j + 1) * tr, 1.0) for m in range(1, K + 1))]
    tf = K * b * tr
    wf = tf / float(p["t_cal_s"])
    return [("own", tf - (nb - 1 - j) * tr, wf), ("own", -(tf - j * tr), wf),
            ("xt", (1 + ng + j) * tr, wf), ("xt", -(ng + nb - j) * tr, wf)]


def _forward_scatter(p: Mapping[str, Any], power_w: float) -> tuple[float, float]:
    """(gate, mode) photons of forward GAWBS + Raman from a co-propagating field of ``power_w`` at f_ref."""
    hnu, f0, mode, filt = _hnu(p["lambda0_m"]), float(p["freq_offset_hz"]), resolve_mode(p["mode"]), p["gbs_filter"]
    g = p["gawbs"]
    fm = _gawbs_modes(g)[0]
    var = gawbs_variances(g, temperature_k=p["temperature_k"], lambda0_m=p["lambda0_m"], length_m=p["L_m"])
    ph = power_w * var / (2 * hnu)  # photons/s in each sideband pair (T17 §6)
    mode_n = float(np.sum(ph * (mode_weight(f0 + fm, mode) + mode_weight(f0 - fm, mode))))
    gate_n = float(p["t_gate_s"]) * float(np.sum(ph * (_filter(f0 + fm, filt) + _filter(f0 - fm, filt))))
    length = float(p["L_m"])  # power_w is the output power: the co-propagating scattered light carries P_out·L

    def n_raman(f):
        w = 2 * np.pi * (f0 - np.asarray(f, dtype=float))  # Ω = bright − scattered
        return power_w * length * raman_scatter_density(w, gamma_per_w_m=p["gamma_per_w_m"], f_R=p["f_R"],
                                                      raman_model=p["raman_model"], temperature_k=p["temperature_k"])

    fb = np.linspace(-3, 3, 6001) * filt["fwhm_hz"]
    gate_n += float(p["t_gate_s"]) * float(np.trapezoid(n_raman(fb) * _filter(fb, filt), fb))
    return gate_n, mode_n + mode_average(n_raman, mode)


def forward_ref_stray(p: Mapping[str, Any], toggles: Mapping[str, bool]) -> np.ndarray:
    """Rows [N_b, len(FORWARD_REF_COLUMNS)] per data slot j (frame-weighted) of a ``forward_ref_stray_table``."""
    hnu, f0, gate = _hnu(p["lambda0_m"]), float(p["freq_offset_hz"]), float(p["t_gate_s"])
    mode = resolve_mode(p["mode"])
    w0, F0 = float(mode_weight(f0, mode)), float(_filter(f0, p["gbs_filter"]))
    phi = float(p["p_peak_w"]) / hnu
    eps_x = 1.0 / float(p["er_dmx_lin"])
    t_ref, tau_e, tau_ap = float(p["t_ref_s"]), float(p["edge_tau_s"]), float(p["tau_afterpulse_s"])
    col = {c: i for i, c in enumerate(FORWARD_REF_COLUMNS)}
    rows = np.zeros((int(p["N_b"]), len(FORWARD_REF_COLUMNS)))
    leak = float(p["eps_mod"]) * float(p["eps_gate"]) * phi  # eps_gate = 1 for guard_slot (checked)
    fwd = _forward_scatter(p, leak * hnu) if toggles["forward_scatter"] else (0.0, 0.0)
    rec = float(p["det_recovery_s"])
    for j in range(rows.shape[0]):
        r = rows[j]
        surv, w_blind = 1.0, 0.0  # P(no click within the recovery time before the slot); every Δ > 0 event of a
        # placement has the same frame weight (guard 1, calibration frames: the frame after)
        if toggles["mod_leak"]:
            r[col["mod_leak_gate"]], r[col["mod_leak_mode"]] = leak * gate * F0, leak * w0
        r[col["forward_gate"]], r[col["forward_mode"]] = fwd
        for kind, dt, wf in _events(p, j):
            if kind == "xt" and not toggles["demux_crosstalk"]:
                continue
            s = 1.0 if kind == "own" else eps_x
            x = abs(dt) - t_ref / 2  # distance from the rectangle's edge
            if toggles["edges"] and tau_e > 0:
                ph_slot = s * phi * math.exp(-max(x, 0.0) / tau_e)
                gn = s * phi * tau_e * _exp_gate(x, gate, tau_e) * F0
                name = "edges" if kind == "own" else "crosstalk"
                r[col[f"{name}_gate"]] += wf * gn
                r[col[f"{name}_mode"]] += wf * ph_slot * w0
            if dt > 0:  # detector after a pulse: afterpulses and recovery
                p_click = -math.expm1(-float(p["qe_det"]) * s * phi * t_ref)
                if toggles["afterpulse"]:
                    r[col["afterpulse_gate"]] += wf * p_click * float(p["p_afterpulse"]) * _exp_gate(dt, gate, tau_ap)
                if dt < rec:
                    surv, w_blind = surv * (1.0 - p_click), wf
        r[col["blind_prob"]] = w_blind * (1.0 - surv)
        r[col["clicks_per_gate"]] = sum(r[col[c]] for c in FORWARD_REF_COLUMNS if c.endswith("_gate"))
        r[col["photons_per_mode"]] = sum(r[col[c]] for c in FORWARD_REF_COLUMNS if c.endswith("_mode"))
    return rows


class ForwardRefStrayTable(_OfflineTable):
    """``forward_ref_stray_table`` (module docstring). Table ``forward_ref_stray`` (axis ``slot`` = data slot j of the
    block, 0 … N_b − 1; columns :data:`FORWARD_REF_COLUMNS`, frame-weighted).

    Params (required, D34): ``placement`` (``cal_frames`` | ``guard_slot``); reference ``p_peak_w``, ``t_ref_s``,
    ``freq_offset_hz``, ``eps_mod``, ``eps_gate`` (slow gate of the reference path during data frames,
    ``cal_frames`` only; 1 = none, must be 1 for ``guard_slot``), ``edge_tau_s``; demux ``er_dmx_lin`` (ER as a
    power ratio > 1; YAML ``er_dmx_db``); detector ``qe_det``, ``p_afterpulse``, ``tau_afterpulse_s``,
    ``det_recovery_s``; timing ``t_rep_s``, ``N_b``, ``n_guard`` (the run's value, ≥ 1 for ``guard_slot``), ``K``,
    ``t_cal_s``; GBS side ``t_gate_s``, ``gbs_filter``, ``mode``; forward scattering ``L_m``, ``gawbs`` group,
    ``gamma_per_w_m``, ``f_R``, ``raman_model``, ``lambda0_m``, ``temperature_k``; D78 linewidth
    ``sigma_target_rad``, ``lock_ratio_R`` (R ≥ 1). Toggles: ``mod_leak``, ``edges``, ``demux_crosstalk``, ``forward_scatter``, ``afterpulse``. Optional ``limit``
    {photons_per_mode_max, source}: CHK ``stray_light_limit`` on the slot mean.
    """

    kind_name, table_name = "forward_ref_stray_table", "forward_ref_stray"
    toggle_names = ("mod_leak", "edges", "demux_crosstalk", "forward_scatter", "afterpulse")

    def __init__(self, params, name):
        name = name or self.kind_name
        _check_raman(params, name)
        check_gawbs_group(params, name)
        require(params, "placement", "p_peak_w", "t_ref_s", "freq_offset_hz", "eps_mod", "eps_gate", "edge_tau_s",
                "er_dmx_lin", "qe_det", "p_afterpulse", "tau_afterpulse_s", "det_recovery_s", "t_rep_s", "N_b",
                "n_guard", "K", "t_cal_s", "t_gate_s", "gbs_filter", "mode", "L_m", "sigma_target_rad",
                "lock_ratio_R", block=name)
        require(params["gbs_filter"], "shape", "fwhm_hz", block=f"{name}.gbs_filter")
        _check_mode(params["mode"], name)
        if params["placement"] not in PLACEMENTS:
            raise ValueError(f"block {name!r}: placement must be one of {PLACEMENTS}, got {params['placement']!r}")
        if params["placement"] == "guard_slot" and int(params["n_guard"]) < 1:
            raise ValueError(f"block {name!r}: placement guard_slot needs n_guard >= 1 (overlay demux/guard_slot)")
        if params["placement"] == "guard_slot" and float(params["eps_gate"]) != 1.0:
            raise ValueError(f"block {name!r}: eps_gate must be 1 for guard_slot (the reference path is open every "
                             "frame; a slow gate cannot separate the guard slot from the data slots)")
        if not 0.0 <= float(params["eps_gate"]) <= 1.0:
            raise ValueError(f"block {name!r}: need 0 <= eps_gate <= 1")
        if int(params["N_b"]) < 1 or int(params["K"]) < 2:
            raise ValueError(f"block {name!r}: need N_b >= 1 and K >= 2")
        if "limit" in params:
            require(params["limit"], "photons_per_mode_max", "source", block=f"{name}.limit")
        super().__init__(self.kind_name, name, params)
        self.chk_log = ChkLog()
        self.result: dict[str, float] | None = None

    def _build(self, p):
        rows = forward_ref_stray(p, self._toggles)
        return {self.table_name: Table(self.table_name, {"slot": np.arange(rows.shape[0], dtype=float)}, rows,
                                       "photons", {"columns": list(FORWARD_REF_COLUMNS),
                                                   "sources": ["lit-aghaeerad2025 §4.2 N4", "D85-4", "D85-11", "D78"]})}

    def check(self, log: ChkLog | None = None) -> dict[str, float]:
        """Slot means of every column, the L_k / linewidth consequence of the placement (:attr:`result`), CHK
        ``stray_light_limit`` and, for ``guard_slot``, ``forward_ref_guard_slot`` (class docstring)."""
        log = log or self.chk_log
        p = self.params
        rows = forward_ref_stray(p, self._toggles)
        r = {f"{c}_mean": float(rows[:, i].mean()) for i, c in enumerate(FORWARD_REF_COLUMNS)}
        r["photons_per_mode_max"] = float(rows[:, FORWARD_REF_COLUMNS.index("photons_per_mode")].max())
        K, nb, ng, tr = int(p["K"]), int(p["N_b"]), int(p["n_guard"]), float(p["t_rep_s"])
        sig2r = float(p["sigma_target_rad"]) ** 2 / (2 * math.pi * float(p["lock_ratio_R"]))
        dtau = (K - 1) * (nb + ng) * tr  # worst port (D78)
        guard = p["placement"] == "guard_slot"
        base_ng = ng - 1 if guard else ng
        dtau0 = (K - 1) * (nb + base_ng) * tr
        r.update(dtau_s=dtau, dnu_max_hz=sig2r / dtau, lk_factor=dtau / dtau0, dnu_max_without_ref_hz=sig2r / dtau0,
                 data_slot_fraction=(nb / (nb + ng)) if guard else (nb / (nb + ng)) * (1 - K * (nb + ng) * tr /
                                                                                      float(p["t_cal_s"])))
        if guard:
            log.warn("forward_ref_guard_slot", self.name,
                     f"reference in the guard slot: L_k × {r['lk_factor']:.3g} (Δτ_{K - 1} = {dtau:.3g} s); D78 linewidth "
                     f"requirement at σ = {float(p['sigma_target_rad']):.3g} rad: {r['dnu_max_hz']:.3g} Hz instead of "
                     f"{r['dnu_max_without_ref_hz']:.3g} Hz", lk_factor=r["lk_factor"], dnu_max_hz=r["dnu_max_hz"],
                     dnu_max_without_ref_hz=r["dnu_max_without_ref_hz"])
        self._limit_check(log, r["photons_per_mode_mean"], "forward reference")
        self.result = r
        return r


@register_block("window_stray_light")
def make_window_stray_light(params, *, name=None) -> Block:
    """Factory of ``window_stray_light`` (:class:`WindowStrayLight`)."""
    return WindowStrayLight(params, name)


@register_block("forward_ref_stray_table")
def make_forward_ref_stray_table(params, *, name=None) -> Block:
    """Factory of ``forward_ref_stray_table`` (:class:`ForwardRefStrayTable`)."""
    return ForwardRefStrayTable(params, name)
