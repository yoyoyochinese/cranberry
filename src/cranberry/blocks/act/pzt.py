"""Phase actuators: ``pzt_fiber_stretcher``, ``pzt_mirror``, ``hv_amp_914`` (PLAN §4.3 ACT; T16 §3.1, §3.3;
T15 §2.7(c)(d); conv §4). The fibre EOM (MOD) reuses :class:`PhaseActuator`.

Phase actuator step (h = rate_div / f_sim), input u [V] from a delay buffer::

    u ← clip(u, drive_min_v, drive_max_v)                          (toggle saturation; act.<name>.sat)
    u ← play operator, width w = hysteresis_frac·(drive_max − drive_min)   (toggle hysteresis, T16 §3.1.2)
    u ← H(u)  linear dynamics, discretised at build time                   (toggle resonance / bandwidth)
    δψ = polarity · K · u     → act.<name>.dpsi [rad] and drv.<segment>.<driver> (G08, F6 names)

Gain K [rad/V] (``gain_mode``, conv §4; runtime values, so a calibration sweep needs no rebuild):

* ``dopl``: K = 2π ΔOPL_full / (λ0 ΔV_full) (``dopl_full_m``, ``v_full_v``, ``lambda0_m``): the spec gives
  the optical-path-length change, no ξ_eo (915B "optical path length measurement"; 8669-1: 16 µm /
  10 V at 1548.8 nm → 6.49 rad/V, T16 §3.1.3).
* ``elongation``: K = 2π n_eff ξ_eo δL_full / (λ0 ΔV_full) (``dl_full_m``, ``v_full_v``, ``lambda0_m``,
  ``n_eff``, ``xi_eo``): δψ = ξ_eo β δL (conv §4; :func:`xi_eo` gives ξ_eo from p11, p12, ν, n).
* ``direct``: K = ``k_rad_per_v``.
* ``mirror`` (``pzt_mirror``): K = 2·(2π/λ0)·cos θ_i·d_x (``d_x_m_per_v``, ``theta_i_rad``, ``lambda0_m``).
* ``vpi`` (``eom_phase``): K = π / V_π (``v_pi_v``).

``polarity`` = σ_act ∈ {+1, −1} (conv §4: positive voltage increases ψ for σ_act = +1).

Dynamics (structural): ``modes`` = [{f_hz, q, r}, ...] → H = Σ r_i ω_i²/(s² + (ω_i/Q_i) s + ω_i²)
(T16 §3.1.2; r_i < 0 gives anti-resonances) for the PZTs, or ``poles_hz`` (series real poles) for
the EOM bandwidth / the 914. ``discretisation``: ``zoh`` (exact for an input held over each step,
e.g. behind a DAC) or ``tustin`` with ``prewarp_hz`` (null = no prewarp; prewarping at a mode
frequency maps |H| and arg H there exactly).

``hv_amp_914`` (Evanescent 914, T16 §3.1.4): out [V, input-referred, DC gain 1] =
poles(slew(u + noise)) after an input clip to [input_min_v, input_max_v] (±5 V full swing). Reference
poles 29.5 kHz and 112 kHz (D51-11); one pole 28 kHz and the slew-rate limit (``slew_v_per_s``,
input-referred) are alternative models: ``poles_hz: []`` + toggle ``slew`` is the pure slew model.
Noise ``noise_v_per_rthz`` (input-referred white density; per-sample std = e_n·√(f_step/2)).
Toggles: ``poles``, ``slew``, ``saturation``, ``noise``.

Losses (G06): ``insertion_loss_db`` (null = none) on ``segment``; ``loss_applies_to``, ``loss_source``;
toggle ``insertion_loss`` (False keeps the record with η = 1).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
from scipy import signal

from cranberry.blocks.ctrl._kit import KitBlock, DrvTarget, drv_write, drv_write_nb, loss_elements
from cranberry.core.block import require
from cranberry.core.registry import register_block
from cranberry.core.rng import DrawBuffer

__all__ = ["GAIN_MODES", "xi_eo", "modes_ss", "poles_ss", "discretise", "PhaseActuator", "HvAmp914",
           "make_pzt_fiber_stretcher", "make_pzt_mirror", "make_hv_amp_914", "lti_step"]

#: gain mode → runtime keys (module docstring)
GAIN_MODES: dict[str, tuple[str, ...]] = {
    "dopl": ("dopl_full_m", "v_full_v", "lambda0_m"),
    "elongation": ("dl_full_m", "v_full_v", "lambda0_m", "n_eff", "xi_eo"),
    "direct": ("k_rad_per_v",),
    "mirror": ("d_x_m_per_v", "theta_i_rad", "lambda0_m"),
    "vpi": ("v_pi_v",),
}
_GAIN_CODE = {m: i for i, m in enumerate(GAIN_MODES)}
_TWO_PI = 2.0 * math.pi


def xi_eo(p11: float, p12: float, nu: float, n: float) -> float:
    """Elasto-optic factor ξ_eo = 1 − (n²/2)[p12 − ν(p11 + p12)] (Butter & Hocker 1978; T16 §3.1.2)."""
    return 1.0 - 0.5 * n * n * (p12 - nu * (p11 + p12))


def modes_ss(modes) -> tuple[np.ndarray, ...]:
    """Continuous (A, B, C, D) of Σ r ω²/(s² + (ω/Q)s + ω²), one [x, ẋ] block per mode."""
    n = len(modes)
    A, B, C = np.zeros((2 * n, 2 * n)), np.zeros((2 * n, 1)), np.zeros((1, 2 * n))
    for i, m in enumerate(modes):
        require(m, "f_hz", "q", "r", block="mode")
        w = _TWO_PI * float(m["f_hz"])
        A[2 * i, 2 * i + 1] = 1.0
        A[2 * i + 1, 2 * i], A[2 * i + 1, 2 * i + 1] = -w * w, -w / float(m["q"])
        B[2 * i + 1, 0] = float(m["r"]) * w * w
        C[0, 2 * i] = 1.0
    return A, B, C, np.zeros((1, 1))


def poles_ss(poles_hz) -> tuple[np.ndarray, ...]:
    """Continuous (A, B, C, D) of Π 1/(1 + s/(2π p)) (DC gain 1)."""
    p = -_TWO_PI * np.asarray(poles_hz, dtype=float)
    if p.size == 0:
        return (np.zeros((0, 0)),)
    return signal.zpk2ss([], p, float(np.prod(-p)))


def discretise(sys_c, h: float, method: str, prewarp_hz: float | None):
    """Discrete (A, B, C, D) as contiguous float arrays (B, C flattened, D scalar)."""
    A = np.asarray(sys_c[0], dtype=float)
    if A.size == 0:
        return np.zeros((0, 0)), np.zeros(0), np.zeros(0), 1.0
    if method == "zoh":
        Ad, Bd, Cd, Dd, _ = signal.cont2discrete(sys_c, h, method="zoh")
    elif method == "tustin":
        hh = h if prewarp_hz is None else 2.0 * math.tan(math.pi * prewarp_hz * h) / (_TWO_PI * prewarp_hz)
        Ad, Bd, Cd, Dd, _ = signal.cont2discrete(sys_c, hh, method="bilinear")
    else:
        raise ValueError(f"discretisation must be 'zoh' or 'tustin', got {method!r}")
    f = lambda a: np.ascontiguousarray(np.asarray(a, dtype=float))  # noqa: E731
    return f(Ad), f(Bd).ravel(), f(Cd).ravel(), float(np.asarray(Dd).ravel()[0])


def lti_step(S, xo, tmp, Ad, Bd, Cd, Dd, u):  # noqa: N803
    """y = C x + D u; x ← A x + B u on S[xo : xo + n] (scratch S[tmp : tmp + n]); returns y."""
    n = Bd.shape[0]
    y = Dd * u
    for i in range(n):
        y += Cd[i] * S[xo + i]
    for i in range(n):
        acc = Bd[i] * u
        for j in range(n):
            acc += Ad[i, j] * S[xo + j]
        S[tmp + i] = acc
    for i in range(n):
        S[xo + i] = S[tmp + i]
    return y


_lti_step_nb = numba.njit(inline="always")(lti_step)


def _make_phase_step(lti, write):
    def step(S, I, P, i_in, i_out, i_sat, i_hy, xo, tmp, Ad, Bd, Cd, Dd, pg, gmode, p_pol, p_lo, p_hi,  # noqa: N803
             p_hf, sat_on, hyst_on, has_drv, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale):
        u = S[i_in]
        lo, hi = P[p_lo], P[p_hi]
        sat = 0
        if sat_on:
            if u > hi:
                u, sat = hi, 1
            elif u < lo:
                u, sat = lo, -1
        if hyst_on:
            w = 0.5 * P[p_hf] * (hi - lo)
            y = S[i_hy]
            y = u - w if y < u - w else (u + w if y > u + w else y)
            S[i_hy] = y
            u = y
        u = lti(S, xo, tmp, Ad, Bd, Cd, Dd, u)
        if gmode == 0:
            k = _TWO_PI * P[pg[0]] / (P[pg[2]] * P[pg[1]])
        elif gmode == 1:
            k = _TWO_PI * P[pg[3]] * P[pg[4]] * P[pg[0]] / (P[pg[2]] * P[pg[1]])
        elif gmode == 2:
            k = P[pg[0]]
        elif gmode == 3:
            k = 2.0 * _TWO_PI / P[pg[2]] * math.cos(P[pg[1]]) * P[pg[0]]
        else:
            k = math.pi / P[pg[0]]
        dpsi = P[p_pol] * k * u
        S[i_out] = dpsi
        I[i_sat] = sat
        if has_drv:
            write(S, I, d_off, d_n, d_lo, d_hi, d_ring, d_head, d_scale, dpsi)

    return step


class PhaseActuator(KitBlock):
    """Generic phase actuator (module docstring); subclasses fix kind, gain modes, dynamics, toggles."""

    STEP = staticmethod(_make_phase_step(lti_step, drv_write))
    STEP_NB = numba.njit(_make_phase_step(_lti_step_nb, drv_write_nb))
    ALLOWED: tuple[str, ...] = ()
    DYN = "modes"  # param key of the dynamics: "modes" or "poles_hz"
    DYN_TOGGLE = "resonance"
    HYSTERESIS = True

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "gain_mode", self.DYN, "discretisation", "prewarp_hz", "segment", block=name)
        if params["gain_mode"] not in self.ALLOWED:
            raise ValueError(f"block {name!r}: gain_mode must be one of {self.ALLOWED}, got {params['gain_mode']!r}")
        self.drv = None if params["segment"] is None else DrvTarget(params, name, "rad", split=True)
        self.a = {k: f"act.{name}.{k}" for k in ("dpsi", "sat", "x", "tmp", "hyst")}
        runtime = (*GAIN_MODES[params["gain_mode"]], "polarity", "drive_min_v", "drive_max_v")
        if self.HYSTERESIS:
            runtime += ("hysteresis_frac",)
        writes = (*self.a.values(), *(() if self.drv is None else (self.drv.slot,)))
        super().__init__(params, name, reads=(params["in"],), writes=writes, runtime=runtime)
        self.slot_units = {params["in"]: "V"}
        self.sys_d = None

    def dynamics(self):
        return modes_ss(self.params["modes"]) if self.DYN == "modes" else poles_ss(self.params["poles_hz"])

    def init(self, layout, rng):
        if self.drv is not None:
            self.drv.check_rate(self.rate_div)
        h = self.rate_div / layout.f_sim_hz
        sys_c = self.dynamics() if self.toggles[self.DYN_TOGGLE] else (np.zeros((0, 0)),)
        self.sys_d = discretise(sys_c, h, self.params["discretisation"], self.params["prewarp_hz"])
        super().init(layout, rng)

    def slots(self):
        n = max(self.sys_d[1].shape[0], 1)
        a = self.a
        out = [(a["dpsi"], (), np.float64, "rad", 0, 0.0), (a["sat"], (), np.int64, "1", 0, 0),
               (a["x"], (n,), np.float64, "1", 0, 0.0), (a["tmp"], (n,), np.float64, "1", 0, 0.0),
               (a["hyst"], (), np.float64, "V", 0, 0.0)]
        return out + ([] if self.drv is None else [self.drv.slot_def()])

    def make_args(self, layout):
        keys = GAIN_MODES[self.params["gain_mode"]]
        pg = np.array([self.p(layout, k) for k in keys] + [-1] * (5 - len(keys)), np.int64)
        t = self.toggles
        drv = (0, 1, 0, 0, 0, -1, 1.0) if self.drv is None else self.drv.args(layout)
        Ad, Bd, Cd, Dd = self.sys_d
        return (self.o(layout, self.params["in"]), self.o(layout, self.a["dpsi"]), self.o(layout, self.a["sat"]),
                self.o(layout, self.a["hyst"]), self.o(layout, self.a["x"]), self.o(layout, self.a["tmp"]),
                Ad, Bd, Cd, Dd, pg, _GAIN_CODE[self.params["gain_mode"]], self.p(layout, "polarity"),
                self.p(layout, "drive_min_v"), self.p(layout, "drive_max_v"),
                self.p(layout, "hysteresis_frac") if self.HYSTERESIS else -1,
                t["saturation"], self.HYSTERESIS and t["hysteresis"], self.drv is not None, *drv)

    def gain_rad_per_v(self, params) -> float:
        """K (without polarity) from the current runtime values."""
        v = {k: float(params[self.param_name(k)]) for k in GAIN_MODES[self.params["gain_mode"]]}
        m = self.params["gain_mode"]
        if m == "dopl":
            return _TWO_PI * v["dopl_full_m"] / (v["lambda0_m"] * v["v_full_v"])
        if m == "elongation":
            return _TWO_PI * v["n_eff"] * v["xi_eo"] * v["dl_full_m"] / (v["lambda0_m"] * v["v_full_v"])
        if m == "direct":
            return v["k_rad_per_v"]
        if m == "mirror":
            return 2.0 * _TWO_PI / v["lambda0_m"] * math.cos(v["theta_i_rad"]) * v["d_x_m_per_v"]
        return math.pi / v["v_pi_v"]

    def freq_response(self, f_hz: Any) -> np.ndarray:
        """Discrete response of the linear dynamics (DC-normalised, no gain), z = e^{+i2πfh} (after init)."""
        Ad, Bd, Cd, Dd = self.sys_d
        z = np.exp(2j * np.pi * np.atleast_1d(np.asarray(f_hz, float)) * self.h)
        n = Bd.shape[0]
        return np.array([Dd + (Cd @ np.linalg.solve(zz * np.eye(n) - Ad, Bd) if n else 0.0) for zz in z])

    def losses(self):
        return loss_elements(self, self.toggles["insertion_loss"])


class PztFiberStretcher(PhaseActuator):
    KIND = "pzt_fiber_stretcher"
    ALLOWED = ("dopl", "elongation", "direct")
    toggle_names = ("resonance", "hysteresis", "saturation", "insertion_loss")


class PztMirror(PhaseActuator):
    KIND = "pzt_mirror"
    ALLOWED = ("mirror", "direct")
    toggle_names = ("resonance", "hysteresis", "saturation", "insertion_loss")


def _make_amp_step(lti):
    def step(S, I, P, i_in, i_out, i_sat, i_prev, xo, tmp, Ad, Bd, Cd, Dd, p_lo, p_hi, p_en, p_sr, h,  # noqa: N803
             buf, idx, sat_on, noise_on, slew_on):
        u = S[i_in]
        sat = 0
        if sat_on:
            if u > P[p_hi]:
                u, sat = P[p_hi], 1
            elif u < P[p_lo]:
                u, sat = P[p_lo], -1
        if noise_on:
            u += P[p_en] * math.sqrt(0.5 / h) * S[buf + I[idx]]
            I[idx] += 1
        if slew_on:
            d = u - S[i_prev]
            lim = P[p_sr] * h
            d = -lim if d < -lim else (lim if d > lim else d)
            u = S[i_prev] + d
            S[i_prev] = u
        S[i_out] = lti(S, xo, tmp, Ad, Bd, Cd, Dd, u)
        I[i_sat] = sat

    return step


class HvAmp914(KitBlock):
    """``hv_amp_914`` (module docstring)."""

    KIND = "hv_amp_914"
    STEP = staticmethod(_make_amp_step(lti_step))
    STEP_NB = numba.njit(_make_amp_step(_lti_step_nb))
    toggle_names = ("poles", "slew", "saturation", "noise")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "in", "out", "poles_hz", "discretisation", "prewarp_hz", block=name)
        self.a = {k: f"act.{name}.{k}" for k in ("sat", "x", "tmp", "prev")}
        self.draws = DrawBuffer(name, 1, "standard_normal")
        super().__init__(params, name, reads=(params["in"],), writes=(params["out"], *self.a.values()),
                         runtime=("input_min_v", "input_max_v", "noise_v_per_rthz", "slew_v_per_s"))
        if self.toggles["noise"]:
            self.writes = (*self.writes, *self.draws.slot_names)
        self.slot_units = {params["in"]: "V", params["out"]: "V"}

    def init(self, layout, rng):
        h = self.rate_div / layout.f_sim_hz
        sys_c = poles_ss(self.params["poles_hz"]) if self.toggles["poles"] else (np.zeros((0, 0)),)
        self.sys_d = discretise(sys_c, h, self.params["discretisation"], self.params["prewarp_hz"])
        super().init(layout, rng)
        if self.toggles["noise"]:
            self.draws.init(layout, rng, self.rate_div)

    def prepare(self, state, params, n0, n_samples):
        if self.toggles["noise"]:
            self.draws.refill(state, n0, n_samples)

    def slots(self):
        n = max(self.sys_d[1].shape[0], 1)
        a = self.a
        return [(self.params["out"], (), np.float64, "V", 0, 0.0), (a["sat"], (), np.int64, "1", 0, 0),
                (a["x"], (n,), np.float64, "1", 0, 0.0), (a["tmp"], (n,), np.float64, "1", 0, 0.0),
                (a["prev"], (), np.float64, "V", 0, 0.0)]

    def make_args(self, layout):
        if self.toggles["noise"]:
            self.draws.bind(layout)
        t = self.toggles
        return (self.o(layout, self.params["in"]), self.o(layout, self.params["out"]), self.o(layout, self.a["sat"]),
                self.o(layout, self.a["prev"]), self.o(layout, self.a["x"]), self.o(layout, self.a["tmp"]),
                *self.sys_d, self.p(layout, "input_min_v"), self.p(layout, "input_max_v"),
                self.p(layout, "noise_v_per_rthz"), self.p(layout, "slew_v_per_s"), self.h,
                self.draws.buf_off, self.draws.idx_off, t["saturation"], t["noise"], t["slew"])

    freq_response = PhaseActuator.freq_response


@register_block("pzt_fiber_stretcher")
def make_pzt_fiber_stretcher(params: Mapping[str, Any], *, name: str | None = None) -> PztFiberStretcher:
    """Factory of the ``pzt_fiber_stretcher`` block (module docstring)."""
    return PztFiberStretcher(params, name)


@register_block("pzt_mirror")
def make_pzt_mirror(params: Mapping[str, Any], *, name: str | None = None) -> PztMirror:
    """Factory of the ``pzt_mirror`` block (module docstring)."""
    return PztMirror(params, name)


@register_block("hv_amp_914")
def make_hv_amp_914(params: Mapping[str, Any], *, name: str | None = None) -> HvAmp914:
    """Factory of the ``hv_amp_914`` block (module docstring)."""
    return HvAmp914(params, name)
