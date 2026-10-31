"""``filter_cavity``: LO filter cavity / mode-cleaning cavity (T05 §5–§6, E17–E18). Owner: WP-B2.

Near resonance (conv §4: E ∝ A e^{−iω0t}, propagation A → A e^{+iβL}), with detuning δ = ν_L − ν_res:

    T_FC(δ) = T0 / (1 − i 2δ/Δν_FC),   arg T = +arctan(2δ/Δν_FC),   |T|² = T0²/(1 + (2δ/Δν_FC)²)
    H_FC(f) = 1/(1 + i 2f/Δν_FC)        (transfer of laser phase fluctuations through the cavity)

At low f, H_FC ≈ e^{−i2πf τ_eff} with τ_eff = 1/(πΔν_FC): for the LO the cavity is a delay line of that
length (T05 §6), so the LO–probe phase difference sees |1 − H_FC|² S_φ ≈ (2/Δν_FC)² S_ν.

Time-loop block (the cavity is a path segment of the LO beam; it writes that segment × beam phase,
so the segment must not also be handled by ``segment_phase``):

    ψ = arctan(2δ(t)/Δν_FC)  [lock_residual]  +  (H_FC φ_L)(t) − φ_L(t)  [laser_lowpass]

``beam_assembler`` adds φ_L(t − τ_b), so the LO carries H_FC φ_L. The low pass has corner Δν_FC/2
(y_n = a y_{n−1} + (1 − a)(x_n + x_{n−1})/2, a = e^{−πΔν_FC h}: DC group delay τ_eff to O(h²)). Its input
is φ_L(t − ``in_delay_s``): set in_delay_s = τ_b of the beam (emission → leaving the cavity segment) so
the LO is exactly H_FC φ_L(t − τ_b). Without the key the input is φ_L "now" (the original behaviour); that
mismatch leaves a delay-interferometer term (1 − H_FC)[φ_L(t) − φ_L(t − τ_b)] of variance ≈ 4π² h0 τ_b
for white FM, which is not small against π h0/Δν_FC when 2π τ_b Δν_FC is not ≪ 1 (Ref A: τ_b = 25 ns,
Δν_FC = 2.5 MHz → +36 %; WP-I1 addition).

Slots: reads ``in_detuning`` [Hz] (cavity lock residual, e.g. a T16 lock block) and ``in_phase`` [rad]
(laser.<s>.phi); writes ``psi.<segment>.<beam>`` [rad], ``psi.<segment>.<beam>.trans`` (|T|²) and the
low-pass state ``psi.<segment>.<beam>.lp``. Static transmission T0² is registered with ``losses()``
(LO path, G06); the block never multiplies powers itself.

Offline helpers (E17): :func:`ase_suppression` (flat ASE through the Airy comb, → π/(2F) for B_sq ≫
FSR), :func:`ase_epsilon`, :func:`ase_mixed_variance` (V = (1 − ε)V_θ + ε(V₊ + V₋)/2, T05 §5.1).

Params (D34)::

    segment, beam, dnu_fc_hz, t0, source (provenance of T0), toggles {lock_residual, laser_lowpass}
    lock_residual: in_detuning      laser_lowpass: in_phase, optional in_delay_s [s] (≥ 0, see above)
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.offline.loss_registry import LossElement

__all__ = ["FilterCavity", "make_filter_cavity", "transmission", "phase_transfer", "tau_eff_s", "ase_suppression",
           "ase_epsilon", "ase_mixed_variance"]


def transmission(delta_hz, dnu_fc_hz: float, t0: float = 1.0):
    """Complex amplitude transmission T_FC(δ) (module docstring)."""
    return t0 / (1.0 - 2j * np.asarray(delta_hz, dtype=float) / dnu_fc_hz)


def phase_transfer(f_hz, dnu_fc_hz: float):
    """H_FC(f) = 1/(1 + i 2f/Δν_FC)."""
    return 1.0 / (1.0 + 2j * np.asarray(f_hz, dtype=float) / dnu_fc_hz)


def tau_eff_s(dnu_fc_hz: float) -> float:
    """Equivalent delay τ_eff = 1/(πΔν_FC)."""
    return 1.0 / (math.pi * dnu_fc_hz)


def ase_suppression(fsr_hz: float, dnu_fc_hz: float, b_sq_hz: float, n_per_line: int = 64) -> float:
    """⟨|T(ν)/T0|²⟩ over |ν − ν_res| ≤ B_sq for flat ASE (Airy comb 1/(1 + (2F/π)² sin²(πν/FSR)), F = FSR/Δν_FC):
    P_ASE,eff with the cavity / without. Whole periods use the exact period average 1/√(1 + (2F/π)²)."""
    m = (2.0 * fsr_hz / (math.pi * dnu_fc_hz)) ** 2
    n_full = math.floor(b_sq_hz / fsr_hz)
    rem = b_sq_hz - n_full * fsr_hz
    part = 0.0
    if rem > 0.0:
        nu = np.linspace(0.0, rem, int(math.ceil(rem / dnu_fc_hz * n_per_line)) + 2)
        part = float(np.trapezoid(1.0 / (1.0 + m * np.sin(np.pi * nu / fsr_hz) ** 2), nu))
    return (n_full * fsr_hz / math.sqrt(1.0 + m) + part) / b_sq_hz


def ase_epsilon(p_ase_eff_w: float, p_lo_w: float) -> float:
    """Equivalent loss ε_ASE = P_ASE,eff/(P_LO + P_ASE,eff) (T05 §5.1)."""
    return p_ase_eff_w / (p_lo_w + p_ase_eff_w)


def ase_mixed_variance(v_theta, v_plus, v_minus, eps):
    """Measured variance / V_SN with ASE beating against broadband squeezed light (T05 §5.1, [AE])."""
    return (1.0 - eps) * v_theta + eps * 0.5 * (v_plus + v_minus)


FC = namedtuple("FC", "f_lock f_lp det ph o tr lp t02 inv a hist nh k mu")


def _build(jit: bool):
    def step(S, I, c):  # noqa: N803
        psi = 0.0
        tr = c.t02
        if c.f_lock:
            x = 2.0 * S[c.det] * c.inv
            psi += math.atan(x)
            tr = c.t02 / (1.0 + x * x)
        if c.f_lp:
            # φ_L at the emission time of the light leaving the cavity now, t − in_delay_s (history
            # shift register, linear interpolation as core.delay); in_delay_s absent → 0 ("now").
            x = S[c.ph]
            if c.nh > 0:
                for j in range(c.nh - 1, 0, -1):
                    S[c.hist + j] = S[c.hist + j - 1]
                S[c.hist] = x
                x = (1.0 - c.mu) * S[c.hist + c.k] + c.mu * S[c.hist + c.k + 1]
            y = c.a * S[c.lp] + (1.0 - c.a) * 0.5 * (x + S[c.lp + 1])
            S[c.lp] = y
            S[c.lp + 1] = x
            psi += y - x
        S[c.o] = psi
        S[c.tr] = tr

    return numba.njit(step) if jit else step


_STEP_PY, _STEP_NB = _build(False), _build(True)


class FilterCavity(BaseBlock):
    """The ``filter_cavity`` block (module docstring)."""

    toggle_names = ("lock_residual", "laser_lowpass")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or "filter_cavity"
        require(params, "segment", "beam", "dnu_fc_hz", "t0", "source", "toggles", block=name)
        tg = params["toggles"]
        if tg.get("lock_residual"):
            require(params, "in_detuning", block=name)
        if tg.get("laser_lowpass"):
            require(params, "in_phase", block=name)
        if not 0.0 <= float(params["t0"]) <= 1.0:
            raise ValueError(f"block {name!r}: t0 must be in [0, 1]")
        base = f"psi.{params['segment']}.{params['beam']}"
        self.o = {"psi": base, "trans": f"{base}.trans", "lp": f"{base}.lp", "hist": f"{base}.hist"}
        reads = [params["in_detuning"]] if tg.get("lock_residual") else []
        reads += [params["in_phase"]] if tg.get("laser_lowpass") else []
        super().__init__("filter_cavity", name, params, reads=reads,
                         writes=tuple(v for k, v in self.o.items() if k != "hist"))
        self.slot_units = {self.o["psi"]: "rad", self.o["trans"]: "1"}
        if tg.get("lock_residual"):
            self.slot_units[params["in_detuning"]] = "Hz"
        if tg.get("laser_lowpass"):
            self.slot_units[params["in_phase"]] = "rad"

    def init(self, layout, rng):
        super().init(layout, rng)
        t02 = float(self.params["t0"]) ** 2
        layout.add(self.o["psi"], (), np.float64, "rad")
        layout.add(self.o["trans"], (), np.float64, "1", init=t02)
        layout.add(self.o["lp"], (2,), np.float64, "rad")
        self._h = block_dt_s(layout.f_sim_hz, self.rate_div)
        self._k, self._mu, self._nh = 0, 0.0, 0
        tau_in = self.params.get("in_delay_s")
        if self.toggles["laser_lowpass"] and tau_in is not None and float(tau_in) > 0.0:
            from cranberry.path.segments import frac_delay

            self._k, self._mu = frac_delay(float(tau_in), self._h)
            self._nh = self._k + 2
            layout.add(self.o["hist"], (self._nh,), np.float64, "rad")
            self.writes = (*self.writes, self.o["hist"])

    def bind(self, layout):
        tg, p = self.toggles, self.params
        off = lambda s: layout.slot(s).offset  # noqa: E731
        dnu = float(p["dnu_fc_hz"])
        self._c = FC(tg["lock_residual"], tg["laser_lowpass"],
                     off(p["in_detuning"]) if tg["lock_residual"] else -1,
                     off(p["in_phase"]) if tg["laser_lowpass"] else -1,
                     off(self.o["psi"]), off(self.o["trans"]), off(self.o["lp"]), float(p["t0"]) ** 2, 1.0 / dnu,
                     math.exp(-math.pi * dnu * self._h), off(self.o["hist"]) if self._nh else -1, self._nh, self._k,
                     self._mu)
        c, nb = self._c, _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, c)

        self.kernel = kernel

    def __call__(self, state, params):
        _STEP_PY(state.S, state.I, self._c)
        return state

    def losses(self) -> list[LossElement]:
        return [LossElement(f"{self.name}.t0", float(self.params["t0"]) ** 2, "component", frozenset({"lo"}),
                            str(self.params["source"]), True, segment=str(self.params["segment"]))]


@register_block("filter_cavity")
def make_filter_cavity(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``filter_cavity`` block (module docstring)."""
    return FilterCavity(params, name)
