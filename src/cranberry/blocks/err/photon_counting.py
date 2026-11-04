"""``photon_counting_fringe``: SNSPD click counting of a weak in-band fringe (ERR/DET, time loop; Phase 6 P5,
D85-5; ``notes/lit-pittaluga2021.md`` §3.3, §4.3). Owner: WP6-C.

Two weak in-band fields (two seeds, or a seed and a reference) meet on a 50/50 coupler whose two outputs go to
click detectors. X = φ_a − φ_b + ``phi_offset_rad`` (φ_b = 0 when ``in_b`` is null: one phase input). Per
calibration gate, output d ∈ {1 (+), 0 (−)} sees the Poisson mean (photons at the detector, conv §2)::

    λ_d = r · [ qe · ( μ (1 ± V cos X)/2 + b ) + dark ]

μ = ``mu_per_gate`` (both outputs together, at the detector), V = ``visibility``, qe = ``qe`` (Liu SNSPD 0.93,
``notes/pulsed-opo-source.md``), b = background photons per gate and output (``bg_per_gate`` + slot ``bg_in``,
e.g. probe stray light), dark = ``dark_per_gate``. ``in_b``, ``bg_in``, ``gate_in`` are optional slots
(absent or null = not wired: a null inside an overlay is stripped by the merge, D74-7).

Unequal seeds (optional slots ``ia_in`` / ``ib_in``, both or neither; e.g. the output intensities of two
``seed_psa`` blocks, relative, ½ each when balanced): μ(1 ± V cos X)/2 becomes μ(I_a + I_b ± 2V√(I_a I_b) cos X)/2.
The estimator is not told (the fringe contrast 2√(I_a I_b)/(I_a + I_b) is a slope error, no bias at π/2).

Squeezed vacuum (toggle ``squeezed_vacuum``, runtime ``sqvac_per_gate`` = b'): a calibration slot with the pump on
also holds the squeezer's vacuum; b' (detected photons per gate and output, :func:`sqvac_background`) is added to
λ_d and, being measurable in seed-off gates, subtracted with the dark counts (reviews/phase6_wp6c_p.md C-2).
ponytail: b' enters as its equivalent Poisson mean (exact no-click probability); the squeezed vacuum's
super-Poisson click statistics and its interference with the seed are not modelled. Click detector (toggle
``saturation``): one click with p = 1 − e^{−λ}; off: photon-number-resolving Poisson count. r = recovery after the last click of the same
output (toggle ``dead_time``): 0 for Δt < ``dead_time_s``, else 1 − exp(−(Δt − dead)/``recovery_s``).

Gates: ``f_gate_hz`` × ``cal_fraction`` calibration gates per second (structural: they fix the number of random
draws), spread evenly by a fractional accumulator; with ``gate_in`` (a mask slot, e.g. ``gate_schedule``) gates
are counted only while the mask is > 0.5.

Window (structural ``window_mode``): ``time`` (runtime ``window_t_s``) or ``gates`` (runtime
``window_n_gates``). At the end of a window the counts n_1, n_0 of G gates give
λ̂_d = −ln(1 − n_d/G) (clicks; n_d/G ≤ 1 − 1/(2G)) or n_d/G (PNR), minus the known ``dark_per_gate``;
cos X̂ = (λ̂_1 − λ̂_0)/(V (λ̂_1 + λ̂_0)) (normalised: independent of μ·qe drifts), X̂ = acos(·) ∈ [0, π] and
``out`` = X̂ − ``setpoint_rad`` [rad], held until the next window. It is the maximum-likelihood estimate; at
X = π/2 its variance reaches the Cramér-Rao bound 1/(2 G F), F = m²V²/(e^{m+d} − 1) per output, m = qe μ/2
(clicks; tests). The background is not subtracted (unknown in practice). A window with λ̂_1 + λ̂_0 ≤ 0 keeps the
previous ``out`` / ``x_hat`` (counts are still written).

Toggles: ``counting_noise`` (off: expected counts, no draws), ``saturation``, ``dark``, ``background``
(off: b = 0), ``dead_time``, ``squeezed_vacuum`` (off: b' = 0). Slots: ``out`` (rad); ``cnt.<name>.n0`` / ``n1`` / ``n_gates`` (last window),
``cnt.<name>.x_hat`` (rad), ``cnt.<name>.window_index`` (int, completed windows); internal
``cnt.<name>.acc`` = [n0, n1, G, steps, gate fraction], ``cnt.<name>.t_last`` (last click time per output).

Loop delay (D55-1, D61-2): the block models the window explicitly: the estimate is the window average (centroid
T/2 old) and is held for one window, ≈ T of delay; the reader adds 1/f_sim when it runs earlier in the pipeline.

ponytail: dead time uses the last click only and is skipped with ``counting_noise`` off (expected counts);
afterpulses and the GBS data gates between calibration gates are not modelled here; the SNSPD recovery into the
data slot next to a calibration gate (lit-pittaluga2021 §3.3 condition 3) is the offline ``cal_slot_check``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.ctrl._kit import KitBlock
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["PhotonCountingFringe", "make_photon_counting_fringe", "fisher_per_gate", "sqvac_background", "WINDOW_MODES"]

WINDOW_MODES: dict[str, int] = {"time": 0, "gates": 1}
_RUNTIME = ("mu_per_gate", "visibility", "qe", "dark_per_gate", "bg_per_gate", "phi_offset_rad", "setpoint_rad",
            "dead_time_s", "recovery_s", "sqvac_per_gate")
_TOGGLES = ("counting_noise", "saturation", "dark", "background", "dead_time", "squeezed_vacuum")


def sqvac_background(r: float, eta: float, overlap: float = 1.0) -> float:
    """Equivalent Poisson background b' [detected photons per gate and output] of the squeezed vacuum in a pumped
    calibration slot: each output sees half of the port-k and half of the port-0 squeezed vacuum (squeeze r), each
    transmitted by η_1 = ½·eta·overlap (eta = squeezer → detector incl. qe; overlap = fraction of the squeezed
    vacuum's mode the calibration gate counts, the seed / detected-mode overlap). Two lossy single-mode squeezed
    vacua, vacuum probability p0 = 1/sqrt((a + ½)(b + ½)) each (a, b the quadrature variances), so
    b' = −2 ln p0 = ln[(a + ½)(b + ½)] (reviews/phase6_wp6c_p.md C-2, p6_wp6c_counting.py E)."""
    e1 = 0.5 * eta * overlap
    a = 0.5 * (e1 * math.exp(-2.0 * r) + 1.0 - e1)
    b = 0.5 * (e1 * math.exp(2.0 * r) + 1.0 - e1)
    return math.log((a + 0.5) * (b + 0.5))


def fisher_per_gate(m: float, vis: float, x: float, dark: float = 0.0, saturation: bool = True) -> float:
    """Fisher information on X per gate, both outputs (λ_± = m(1 ± V cos X) + d): Σ (∂λ/∂X)²/(e^λ − 1) for clicks,
    Σ (∂λ/∂X)²/λ for PNR (lit-pittaluga2021 §3.3)."""
    out = 0.0
    for s in (1.0, -1.0):
        lam = m * (1.0 + s * vis * math.cos(x)) + dark
        d = m * vis * math.sin(x)
        out += d * d / (math.expm1(lam) if saturation else lam)
    return out


def _poisson(u, lam):
    """Poisson(λ) by inversion of one uniform u ∈ [0, 1)."""
    k = 0
    p = math.exp(-lam)
    c = p
    while u > c and k < 10000:
        k += 1
        p *= lam / k
        c += p
    return float(k)


_poisson_nb = numba.njit(_poisson)


def _make_step(poisson):
    def step(S, I, P, i_t, i_a, i_b, i_bg, i_gate, i_ia, i_ib, i_out, i_n0, i_n1, i_ng, i_x, i_wi, acc, tl,  # noqa: N803, E741
             p0, rate, h, n_max, mode, buf, idx, noise_on, sat_on, dark_on, bg_on, dead_on, sq_on):
        mu, vis, qe = P[p0], P[p0 + 1], P[p0 + 2]
        dark = P[p0 + 3] if dark_on else 0.0
        if sq_on:
            dark += P[p0 + 9]  # squeezed vacuum b': known (seed-off gates), subtracted like the dark counts
        x = S[i_a] + P[p0 + 5]
        if i_b >= 0:
            x -= S[i_b]
        b = 0.0
        if bg_on:
            b = P[p0 + 4] + (S[i_bg] if i_bg >= 0 else 0.0)
        if i_gate < 0 or S[i_gate] > 0.5:
            S[acc + 4] += rate
        n = int(math.floor(S[acc + 4]))
        S[acc + 4] -= n
        m = 0.5 * qe * mu
        vc = vis * math.cos(x)
        m_sum, m_int = m, m  # balanced seeds (no intensity inputs)
        if i_ia >= 0:  # seed intensities I_a, I_b (relative; balanced = ½ each): μ(I_a + I_b ± 2V√(I_a I_b) cos X)/2
            m_sum = qe * mu * 0.5 * (S[i_ia] + S[i_ib])
            m_int = qe * mu * math.sqrt(max(S[i_ia] * S[i_ib], 0.0))
        t = S[i_t]
        for j in range(n_max):
            tg = t + j * h / max(n, 1)
            for d in range(2):
                if i_ia >= 0:
                    lam = qe * b + dark + m_sum + (m_int * vc if d == 1 else -m_int * vc)
                else:
                    lam = qe * b + dark + m * ((1.0 + vc) if d == 1 else (1.0 - vc))
                if dead_on and noise_on:
                    dt = tg - S[tl + d]
                    dead, rec = P[p0 + 7], P[p0 + 8]
                    lam *= 0.0 if dt < dead else -math.expm1(-(dt - dead) / rec)
                if noise_on:
                    u = S[buf + I[idx]]
                    I[idx] += 1
                    if j >= n:
                        continue  # draw consumed and discarded (G11)
                    if sat_on:
                        k = 1.0 if u < -math.expm1(-lam) else 0.0
                    else:
                        k = poisson(u, lam)
                    if k > 0.0:
                        S[tl + d] = tg
                elif j < n:
                    k = -math.expm1(-lam) if sat_on else lam
                else:
                    continue
                S[acc + d] += k
        S[acc + 2] += n
        S[acc + 3] += 1.0
        if mode == 0:
            done = S[acc + 3] >= max(1.0, math.floor(P[p0 + 10] / h + 0.5))
        else:
            done = S[acc + 2] >= P[p0 + 10]
        if not done:
            return
        g = S[acc + 2]
        S[i_n0] = S[acc]
        S[i_n1] = S[acc + 1]
        S[i_ng] = g
        I[i_wi] += 1
        if g > 0.0:
            l0 = S[acc] / g
            l1 = S[acc + 1] / g
            if sat_on:
                top = 1.0 - 0.5 / g
                l0 = -math.log(1.0 - min(l0, top))
                l1 = -math.log(1.0 - min(l1, top))
            l0 -= dark
            l1 -= dark
            if l0 + l1 > 0.0:
                c = (l1 - l0) / (vis * (l1 + l0))
                xh = math.acos(min(1.0, max(-1.0, c)))
                S[i_x] = xh
                S[i_out] = xh - P[p0 + 6]
        for q in range(4):
            S[acc + q] = 0.0

    return step


_step = _make_step(_poisson)
_step_nb = numba.njit(_make_step(_poisson_nb))


class PhotonCountingFringe(KitBlock):
    KIND = "photon_counting_fringe"
    STEP = staticmethod(_step)
    STEP_NB = _step_nb
    toggle_names = _TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in_a", "out", "f_gate_hz", "cal_fraction", "window_mode", block=name)
        params = {"in_b": None, "bg_in": None, "gate_in": None, "ia_in": None, "ib_in": None,
                  **params}   # optional wiring (an overlay null is stripped)
        if (params["ia_in"] is None) != (params["ib_in"] is None):
            raise ValueError(f"block {name!r}: give both ia_in and ib_in or neither")
        mode = params["window_mode"]
        if mode not in WINDOW_MODES:
            raise ValueError(f"block {name!r}: window_mode must be one of {tuple(WINDOW_MODES)}, got {mode!r}")
        self.win_key = "window_t_s" if mode == "time" else "window_n_gates"
        if not float(params["f_gate_hz"]) > 0.0 or not 0.0 < float(params["cal_fraction"]) <= 1.0:
            raise ValueError(f"block {name!r}: need f_gate_hz > 0 and 0 < cal_fraction <= 1")
        c = f"cnt.{name}"
        self.n_names = {k: f"{c}.{k}" for k in ("n0", "n1", "n_gates", "x_hat", "window_index", "acc", "t_last")}
        reads = tuple(params[k] for k in ("in_a", "in_b", "bg_in", "gate_in", "ia_in", "ib_in") if params[k] is not None)
        super().__init__(params, name, reads=reads + ("t",), writes=(params["out"], *self.n_names.values()),
                         runtime=(*_RUNTIME, self.win_key))
        self.draws = DrawBuffer(self.name, 1, "random")
        if self.toggles["counting_noise"]:
            self.writes = (*self.writes, *self.draws.slot_names)
        self.slot_units = {params["out"]: "rad", self.n_names["x_hat"]: "rad"}

    def slots(self):
        n = self.n_names
        return [(self.params["out"], (), np.float64, "rad", 0, 0.0), (n["x_hat"], (), np.float64, "rad", 0, 0.0),
                (n["n0"], (), np.float64, "1", 0, 0.0), (n["n1"], (), np.float64, "1", 0, 0.0),
                (n["n_gates"], (), np.float64, "1", 0, 0.0), (n["window_index"], (), np.int64, "1", 0, 0),
                (n["acc"], (5,), np.float64, "1", 0, 0.0), (n["t_last"], (2,), np.float64, "s", 0, -1.0e30)]

    def init(self, layout, rng):
        super().init(layout, rng)
        self.rate = float(self.params["f_gate_hz"]) * float(self.params["cal_fraction"]) * self.h
        self.n_max = max(1, math.ceil(self.rate - 1e-12))
        if self.toggles["counting_noise"]:
            self.draws = DrawBuffer(self.name, 2 * self.n_max, "random")
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self.toggles["counting_noise"]:
            self.draws.refill(state, n0, n_samples)

    def make_args(self, layout):
        if self.toggles["counting_noise"]:
            self.draws.bind(layout)
        offs = [self.p(layout, k) for k in (*_RUNTIME, self.win_key)]
        if offs != list(range(offs[0], offs[0] + len(offs))):
            raise RuntimeError(f"block {self.name!r}: runtime params not contiguous in P")
        n, t, o = self.n_names, self.toggles, self.o
        pr = self.params
        return (o(layout, "t"), o(layout, pr["in_a"]), o(layout, pr["in_b"]), o(layout, pr["bg_in"]),
                o(layout, pr["gate_in"]), o(layout, pr["ia_in"]), o(layout, pr["ib_in"]), o(layout, pr["out"]), o(layout, n["n0"]), o(layout, n["n1"]),
                o(layout, n["n_gates"]), o(layout, n["x_hat"]), o(layout, n["window_index"]), o(layout, n["acc"]),
                o(layout, n["t_last"]), offs[0], self.rate, self.h, self.n_max, WINDOW_MODES[pr["window_mode"]],
                self.draws.buf_off, self.draws.idx_off, t["counting_noise"], t["saturation"], t["dark"],
                t["background"], t["dead_time"], t["squeezed_vacuum"])


@register_block("photon_counting_fringe")
def make_photon_counting_fringe(params: Mapping[str, Any], *, name: str | None = None) -> PhotonCountingFringe:
    """Factory of the ``photon_counting_fringe`` block (module docstring)."""
    return PhotonCountingFringe(params, name)
