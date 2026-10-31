"""``probe_pulser``: AOM gating of the Liu-type probe and per-port stray-light gating (WP-B8; T10 §3.4.2,
§3.4.4, §6.2; T16 §3.5; PLAN §4.3 MOD, §2.3.2 ``leak``; D40, D46, D59-7).

Time loop (T): the probe is a rectangular pulse train of period ``period_s`` (normally the GBS frame
T_f = K·N_b·T_rep), duty D, start ``t_offset_s`` relative to the squeezed slot at the GBS input, and AOM
off-state power leak ε_AOM (relative to the peak). Average power in the pulses P̄ = ``p_avg_w`` (the
same P̄ as ``liu_probe_interferometer.p_pr_avg_w``), so P_pk = P̄/D. Every step writes the stray-light
table inputs (C15 names) for MEAS:

    leak.<name>.p_avg_w     total average power P_pk (D + ε_AOM (1 − D))       [W]
    leak.<name>.p_peak_w    P_pk                                              [W]
    leak.<name>.duty        D
    leak.<name>.t_offset_s  probe start relative to the squeezed slot          [s]
    leak.<name>.eps_aom     ε_AOM (0 when ``aom_leak`` is off)

The phase / frequency shift of the AOM is ``aom_shift`` (WP-B4); this block only gates.
Toggle ``gating`` off ⇒ CW probe of the same average power (D = 1, ε = 0; T10 §3.4.4 compares at equal
P̄), ``aom_leak`` off ⇒ ε_AOM = 0.

Offline, per port (:func:`gated_stray_photons`, D59-7: B10's ``scattering.stray_photons`` gives the
numbers, this adds only the timing). With the injection point at the GBS input (time origin = the
squeezed slot there, gate [−T_gate/2, T_gate/2) modulo the period):

* distributed backscatter (Rayleigh, Brillouin, Raman) of port k's fiber returns as the probe envelope
  convolved with a uniform box of length W = 2 n_g L/c (T10 eq. (8), αL ≪ 1), so the gate collects
  N = N_CW(P_pk)·[ε + (1 − ε)·I/(W·T_gate)], I = ∫_gate ∫_0^W rect(t − t′) dt′ dt (exact: the inner
  integral is a trapezoid in t). It vanishes only if 2n_gL_k/c + τ_pr + T_gate fits in the period (k ≤ 15 at 32 ports,
  T10 test 7); otherwise it equals the CW value on average;
* a discrete reflection at distance z_i returns the envelope delayed by 2 n_g z_i/c:
  N_i = N_i,CW(P_pk)·[ε + (1 − ε)·|gate ∩ pulse|/T_gate] — timing removes it whenever the delayed pulse
  misses the gate;
* photons per squeezed **mode** are weighted by the mode, not by the gate (conv §2, N = ∫n|ũ|²df; M5-02,
  :func:`mode_factors`). The mode u(t) is ``p['mode']`` (parametric or the source's optimal mode,
  ``scattering.mode_time``) with its intensity centroid on the squeezed slot (t = 0, the gate centre).
  With the probe amplitude train a(t) (rect pulses, peak 1): a line at offset f0 (Rayleigh, Brillouin,
  reflections; coherent over the mode) gives N = Φ_pk·|∫u(t) a(t − τ) e^{+i2πf0t} dt|², averaged over the
  scatterer delays τ ∈ [0, W] (distributed: independent scatterers add in power) or at τ = 2n_g z_i/c
  (reflection i); a ≡ 1 recovers the CW Φ|ũ(f0)|². Broadband Raman: N = N_CW(P_pk)·∫|u|²⟨a(t − τ)⟩_τ dt;
  a line beyond the mode grid's Nyquist frequency uses Φ|ũ(f0)|²·∫|u|²⟨a⟩ (intensity envelope).
  The AOM leak adds ε·N_CW(P_pk) to both columns (incoherent floor, ε ≪ 1).
N_CW(P) is linear in P (``stray_photons`` has no SBS gain), so the scaling is exact. Not included:
scattering in the probe distribution path (pdist_k) — a separate stray table of that segment.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.blocks.fib.scattering import (
    STRAY_COLUMNS,
    _hnu,
    backscatter_sources,
    mode_time,
    mode_weight,
    resolve_mode,
    stray_photons,
)
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.path.topology import C_M_PER_S

__all__ = ["ProbePulser", "make_probe_pulser", "gate_factors", "mode_factors", "gated_stray_photons"]

_LEAK = ("p_avg_w", "p_peak_w", "duty", "t_offset_s", "eps_aom")


def _trapezoid_in_gate(c: float, d: float, a: float, b: float, w: float) -> float:
    """∫_c^d |[t − b, t − a] ∩ [0, w]| dt: the integrand is piecewise linear with knots a, b, a + w, b + w, so
    the trapezoid rule on those knots is exact (and exactly 0 for disjoint supports)."""
    def g(t):
        return max(0.0, min(t - a, w) - max(t - b, 0.0))

    ts = sorted({c, d, *(x for x in (a, b, a + w, b + w) if c < x < d)})
    return sum(0.5 * (g(t0) + g(t1)) * (t1 - t0) for t0, t1 in zip(ts[:-1], ts[1:]))


def gate_factors(*, period_s: float, tau_pr_s: float, t_offset_s: float, t_gate_s: float, roundtrip_s: float,
                 delays_s: Sequence[float] = ()) -> tuple[float, list[float]]:
    """(f_dist, [f_i]): fraction of the CW-equivalent count that a rectangular probe train puts into the gate,
    for distributed backscatter of round trip W = ``roundtrip_s`` and for reflections delayed by ``delays_s``
    (module docstring; envelope peak = 1, so the CW value with the peak power is the reference)."""
    c, d, W = -0.5 * t_gate_s, 0.5 * t_gate_s, float(roundtrip_s)
    span = tau_pr_s + max([W, *delays_s])
    m0 = math.floor((c - span - t_offset_s) / period_s) - 1
    m1 = math.ceil((d - t_offset_s) / period_s) + 1
    dist, refl = 0.0, [0.0] * len(delays_s)
    for m in range(m0, m1 + 1):
        a = t_offset_s + m * period_s
        b = a + tau_pr_s
        dist += _trapezoid_in_gate(c, d, a, b, W)
        for i, di in enumerate(delays_s):
            refl[i] += max(0.0, min(d, b + di) - max(c, a + di))
    return dist / (W * t_gate_s), [r / t_gate_s for r in refl]


def mode_factors(t: np.ndarray, u: np.ndarray, *, period_s: float, tau_pr_s: float, t_offset_s: float,
                 roundtrip_s: float, delays_s: Sequence[float] = (), f0_hz: float | None = 0.0) -> tuple[float, list[float]]:
    """(g_dist, [g_i]): the mode-weighted counterpart of :func:`gate_factors` (module docstring). Line
    (``f0_hz`` given): g = |∫u a(t − τ) e^{i2πf0t} dt|² [s] (CW: |ũ(f0)|²); broadband (``f0_hz`` None):
    g = ∫|u|² a(t − τ) dt [1] (CW: 1). ``g_dist`` averages over τ ∈ [0, W = ``roundtrip_s``], g_i is at
    τ = ``delays_s[i]``. u on the uniform grid ``t`` (gate coordinates), treated as piecewise constant."""
    dt = float(t[1] - t[0])
    if f0_hz is None:
        x = np.abs(u) ** 2 * dt
    else:
        x = u * np.exp(2j * np.pi * float(f0_hz) * t) * dt
    edges = np.concatenate([[t[0] - 0.5 * dt], t + 0.5 * dt])
    cum = np.concatenate([[0.0], np.cumsum(x)])
    e0, e1, P, tp = edges[0], edges[-1], float(period_s), float(tau_pr_s)

    def cu(y):  # ∫_{−∞}^{y} x (linear between the bin edges, constant outside)
        if f0_hz is None:
            return np.interp(y, edges, cum)
        return np.interp(y, edges, cum.real) + 1j * np.interp(y, edges, cum.imag)

    def h(s):  # g for the pulse train whose m = 0 pulse starts at s (P-periodic in s)
        s = np.asarray(s, dtype=float)
        m = np.arange(math.floor((e0 - tp - s.max()) / P), math.ceil((e1 - s.min()) / P) + 1)
        st = s[..., None] + m * P
        amp = np.sum(cu(st + tp) - cu(st), axis=-1)
        return amp if f0_hz is None else np.abs(amp) ** 2

    refl = [float(h(np.array([t_offset_s + d]))[0]) for d in delays_s]
    w = float(roundtrip_s)
    if w <= 0.0:
        return float(h(np.array([t_offset_s]))[0]), refl
    # ∫_0^W h(t_off + τ) dτ = ∫_window h(s) c(s) ds over one period, c(s) = #{m : s + mP − t_off ∈ [0, W]};
    # knots: the bin edges seen by the leading and the trailing pulse edge, the jumps of c, the window ends
    s0 = e0 - tp
    wrap = lambda y: s0 + np.mod(np.asarray(y) - s0, P)  # noqa: E731
    knots = np.unique(np.concatenate([wrap(edges), wrap(edges - tp), wrap([t_offset_s, t_offset_s + w]), [s0, s0 + P]]))
    mid = 0.5 * (knots[1:] + knots[:-1])
    c = np.floor((t_offset_s + w - mid) / P) - np.ceil((t_offset_s - mid) / P) + 1
    hk = h(knots)
    return float(np.sum(0.5 * (hk[1:] + hk[:-1]) * np.diff(knots) * c) / w), refl


def gated_stray_photons(p: Mapping[str, Any], toggles: Mapping[str, bool], *, p_avg_w: float, duty: float,
                        period_s: float, t_offset_s: float, eps_aom: float, gating: bool, n_g: float,
                        reflection_z_m: Sequence[float]) -> dict[str, float]:
    """``stray_photons`` (B10) for one port with the probe timing applied (module docstring). ``p`` is B10's
    stray-light parameter set of the port (``L_m`` = port fiber, ``t_gate_s``, ``reflections`` = [[R_i, ε_i]]),
    ``reflection_z_m`` the distances of those reflections from the injection point."""
    if len(reflection_z_m) != len(p["reflections"]):
        raise ValueError("reflection_z_m must give one distance per entry of p['reflections']")
    if not gating:
        return stray_photons(p, toggles, p_avg_w)  # CW of the same average power (T10 §3.4.4)
    pk = p_avg_w / duty
    w = 2.0 * n_g * float(p["L_m"]) / C_M_PER_S
    f_d, f_r = gate_factors(period_s=period_s, tau_pr_s=duty * period_s, t_offset_s=t_offset_s,
                            t_gate_s=float(p["t_gate_s"]), roundtrip_s=w,
                            delays_s=[2.0 * n_g * z / C_M_PER_S for z in reflection_z_m])
    p = {**p, "mode": resolve_mode(p["mode"])}
    t, u = mode_time(p["mode"])
    timing = dict(period_s=period_s, tau_pr_s=duty * period_s, t_offset_s=t_offset_s, roundtrip_s=w)
    delays = [2.0 * n_g * z / C_M_PER_S for z in reflection_z_m]
    cw = stray_photons(p, {**toggles, "stray_reflections": False}, pk)
    out = {k: v * (eps_aom + (1.0 - eps_aom) * f_d) for k, v in cw.items() if k.endswith("_gate")}
    out.update({k: eps_aom * v for k, v in cw.items() if k.endswith("_mode")})
    lines, _ = backscatter_sources(p, {"rayleigh": toggles["stray_rayleigh"], "brillouin": toggles["stray_brillouin"],
                                       "raman": False}, pk)
    iso = float(p["wdm_isolation_lin"])
    g_bb, _ = mode_factors(t, u, f0_hz=None, **timing)
    nyq = 0.5 / float(t[1] - t[0])
    for name, f0, phi, _ in lines:
        # a line beyond the mode grid's Nyquist frequency (Brillouin ±10.8 GHz on an ns-mode grid) is far
        # outside the mode band: gating acts there as an intensity envelope, |ũ(f0)|²·∫|u|²⟨a⟩
        g = mode_factors(t, u, f0_hz=f0, **timing)[0] if abs(f0) < nyq else float(mode_weight(f0, p["mode"])) * g_bb
        out[f"{name}_mode"] += (1.0 - eps_aom) * phi * (iso if name == "rayleigh" else 1.0) * g
    if toggles["stray_raman"]:
        out["raman_mode"] += (1.0 - eps_aom) * cw["raman_mode"] * g_bb
    if toggles["stray_reflections"]:
        only = dict.fromkeys(("stray_rayleigh", "stray_brillouin", "stray_raman"), False) | {"stray_reflections": True}
        _, g_r = mode_factors(t, u, f0_hz=float(p["probe"]["freq_offset_hz"]), delays_s=delays, **timing)
        for refl, f, g in zip(p["reflections"], f_r, g_r):
            one = stray_photons({**p, "reflections": [list(refl)]}, only, pk)
            out["reflections_gate"] += one["reflections_gate"] * (eps_aom + (1.0 - eps_aom) * f)
            phi = pk * float(refl[0]) * float(refl[1]) / _hnu(p["lambda0_m"])
            out["reflections_mode"] += one["reflections_mode"] * eps_aom + (1.0 - eps_aom) * phi * iso * g
    out["clicks_per_gate"] = sum(out[k] for k in STRAY_COLUMNS if k.endswith("_gate") and k != "clicks_per_gate")
    out["photons_per_mode"] = sum(out[k] for k in STRAY_COLUMNS if k.endswith("_mode") and k != "photons_per_mode")
    return out


def _make_step():
    def step(S, I, P, o_avg, o_pk, o_d, o_t, o_e, pa, pd, pt, pe, gating, leak):  # noqa: N803, E741
        d, e = 1.0, 0.0
        if gating:
            d = P[pd]
            e = P[pe] if leak else 0.0
        pk = P[pa] / d
        S[o_avg] = pk * (d + e * (1.0 - d))
        S[o_pk] = pk
        S[o_d] = d
        S[o_t] = P[pt]
        S[o_e] = e

    return step


class ProbePulser(KitBlock):
    """``probe_pulser`` (module docstring).

    Params (required, D34): ``period_s``; runtime ``p_avg_w``, ``duty`` (0, 1], ``t_offset_s``, ``eps_aom``.
    Toggles: ``gating``, ``aom_leak``.
    """

    KIND = "probe_pulser"
    STEP = staticmethod(_make_step())
    STEP_NB = numba.njit(_make_step())
    toggle_names = ("gating", "aom_leak")
    RUNTIME = ("p_avg_w", "duty", "t_offset_s", "eps_aom")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "period_s", *self.RUNTIME, block=name)
        if not 0.0 < float(params["duty"]) <= 1.0 or float(params["period_s"]) <= 0.0:
            raise ValueError(f"block {name!r}: need 0 < duty <= 1 and period_s > 0")
        self.leak = tuple(f"leak.{name}.{k}" for k in _LEAK)
        super().__init__(params, name, reads=(), writes=self.leak, runtime=self.RUNTIME)

    def slots(self):
        units = ("W", "W", "1", "s", "1")
        return [(s, (), np.float64, u, 0, 0.0) for s, u in zip(self.leak, units)]

    def make_args(self, layout):
        return (*(self.o(layout, s) for s in self.leak), *(self.p(layout, k) for k in self.RUNTIME),
                self.toggles["gating"], self.toggles["aom_leak"])

    def stray(self, p: Mapping[str, Any], toggles: Mapping[str, bool], *, n_g: float,
              reflection_z_m: Sequence[float]) -> dict[str, float]:
        """:func:`gated_stray_photons` with this block's spec timing (one port's B10 parameter set ``p``)."""
        q = self.params
        return gated_stray_photons(p, toggles, p_avg_w=float(q["p_avg_w"]), duty=float(q["duty"]),
                                   period_s=float(q["period_s"]), t_offset_s=float(q["t_offset_s"]),
                                   eps_aom=float(q["eps_aom"]) if self.toggles["aom_leak"] else 0.0,
                                   gating=self.toggles["gating"], n_g=n_g, reflection_z_m=reflection_z_m)


@register_block("probe_pulser")
def make_probe_pulser(params: Mapping[str, Any], *, name: str | None = None) -> ProbePulser:
    """Factory of the ``probe_pulser`` block (module docstring)."""
    return ProbePulser(params, name)
