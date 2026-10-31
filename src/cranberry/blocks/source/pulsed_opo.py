"""``pulsed_opo_source``: pulse-pumped below-threshold OPO, time-domain quantum Langevin (T17 §2).
Owner: WP-B11.

Offline (O) block. Product: :attr:`PulsedOpoSource.source_output` (``temporal`` ``SourceOutput``,
G07 / D55-4) plus the optimal LO mode, supermode squeeze parameters r_n, Schmidt number K, τ_eff,
|ũ(f)|², ε_tail(T) and the detuning coefficient tables read by ``opo_angle_map``. Never stepped.

Model (T17 §2.1, symbols of conv §4 / T17 §1)::

    db/dt = −(κ_p + iΔ_p) b + κ_p s(t)                       (pump cavity; b = s if non-resonant)
    dX/dt = M(t) X + √(2κ_T) X_in + √(2κ_L) X_L,   X = (q, p),  â_out = √(2κ_T) â − â_in
    M = [[−κ_c + ε_r, Δ_s + ε_i], [−Δ_s + ε_i, −κ_c − ε_r]],  ε(t) = κ_c x_ss b(t) e^{iφ}

κ_c = π Δν_c, κ_T = ρκ_c, κ_L = (1 − ρ)κ_c, κ_p = π Δν_p. With φ = 0, Δ = 0 the p quadrature is
squeezed (θ_sq = π/2, φ_OPA,0 = 0 of this Langevin form, T17 §1); φ = ``phi_opa0_rad`` gives the D55-4
frame (θ_sq,0 = φ/2 + π/2).

Discretisation (T17 §2.2, δt = ``dt_kappa``/κ_c, PLAN §6.4: 1/(40–160 κ_c)): bins of width δt from the
pump turn-on t = 0 to ``window_kappa``/κ_c. Inside a bin M is constant (``box``: so are the projected
inputs on the bin boxes), the Langevin equation is solved exactly
(augmented matrix exponential for e^{Mτ}, ∫e^{Ms}ds, ∫∫), and the output is bin-averaged:
Y_k = δt^{−1/2}∫_bin X_out dt. Y = T·w with w ~ N(0, 𝟙) (input, loss and the initial intracavity
vacuum), so the shot-noise-normalised 2N×2N kernel is C = T Tᵀ (xxpp, vacuum = 𝟙; ``cov`` = C/2).
Errors scale with δt² (V-P1). ``discretization`` (D34, no default; M5-07):

* ``box``: the white inputs are projected on the bin boxes (4 noise columns per bin). C obeys the
  uncertainty relation only to O((κ_c δt)²) (vacuum eigenvalues within 2(κ_c δt)², V-P2; min symplectic
  eigenvalue 0.99979 at δt = 1/(40κ_c) — 0.002 dB, but a thermal photon number (ν − 1)/2 of −1e−4
  downstream).
* ``van_loan``: the in-bin noise of (X(t_k+1), ∫_bin X dt, box part of the input) is integrated exactly
  (Van Loan's augmented exponential, 6 columns per bin), so C is the exact covariance of the bin
  averages and physical to rounding (use when a consumer needs a physical covariance, e.g. GBS sampling). The initial-cavity-vacuum columns ("the term 2κ_T G(t,t₀)G(t′,t₀)ᵀ")
are **always** included; dropping them gives the spurious perfect squeezing of V-P3 — the toggle
``include_initial_vacuum`` exists only for that regression test. The basis modes are the bin boxes
f_k = δt^{−1/2}·1_bin (orthonormal), ``grid`` = bin centres.

Frame of ``cov`` and pump-path phase (conv §4, PLAN G07): ``cov`` uses b·e^{+i arctan(Δ_p/κ_p)}, i.e.
the static phase of the resonant pump cavity is taken out and reported as
``pump_path_phase_rad[segment] = −arctan(Δ_p/κ_p)`` (enters θ_sq as ½ψ in the time loop); θ_off of
``project`` then holds the signal-detuning rotation and the pump chirp. When ``opo_angle_map`` drives
the OPO segment in the time loop it carries the *total* detuning rotation, so set the static
detuning here to 0 (or the toggle ``detuning`` off) to avoid counting it twice.

Detuning coefficient tables (T17 §2.5, V-P9): for each x_ss of ``angle_table_x_ss`` (None → the
operating point) the optimal mode u_0 at Δ = 0 is the fixed LO and the minimum angle is found at
Δ_s = d·κ_c (Δ_p = 0) and Δ_p = d·κ_c (Δ_s = 0), d = ``angle_probe_d``; with the full b (no phase
taken out). Tables (rad per rad/s, axis ``x_ss``):

* ``<segment>.detuning.squeezed`` columns (∂θ_sq/∂Δ_s, ∂θ_sq/∂Δ_p)
* ``<segment>.detuning.pump``     columns (0, ∂ arg∫b dt/∂Δ_p): phase of the transmitted pump
  (Liu-type lock beam, T17 §2.5: locking it removes ½ of it from θ_sq)

The rotation is linear in d up to d ≈ 0.05 (V-P9); the map is first order. Borealis values at
x_ss = 1.3, ρ = 0.95 (per Δ/κ_c): −0.638 (signal), −0.425 (pump, Δ_p = 2Δ_s), −1.063 (length),
−0.543 (transmitted-pump lock).

``pulsed`` off (T17 §8.1, VT-03): the block returns the CW ``opo_cw_spectrum`` result (``cw_sideband``,
linewidth model, x = x_ss, Δ = Δ_s) bit for bit; needs ``f_s_hz``.

Params (D34: all required where they apply)
    ``cavity_model`` (linewidth: ``fwhm_hz`` | fsr_finesse: ``fsr_hz``, ``finesse``),
    ``escape_efficiency`` ρ, ``x_ss``, ``phi_opa0_rad``, ``delta_s_rad_s``, ``delta_p_rad_s``,
    ``pump_shape`` (rect: ``pump_duration_s`` | envelope: ``pump_envelope`` callable of ``t_s``, peak 1),
    ``dt_kappa``, ``window_kappa``, ``segment``, ``rep_rate_hz`` (→ ``tail_into_next_slot`` =
    ε_tail(1/f_rep), T from the pump turn-on; 1/f_rep beyond the window → analytic 8.8e^{−2κ_cT} + CHK
    ``source_window_short``), ``angle_table_x_ss`` (list or None),
    ``angle_probe_d``, ``discretization`` (box | van_loan); pump_resonant on: ``pump_fwhm_hz``; pulsed off:
    ``f_s_hz``.
    ``toggles`` {pulsed, pump_resonant, detuning, include_initial_vacuum}.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np
from scipy.linalg import expm

from cranberry.blocks.chk.rep_rate_tail import tail_analytic
from cranberry.blocks.source import temporal_modes as tm
from cranberry.blocks.source.opo_cw import make_opo_cw_spectrum
from cranberry.config.resolve import as_callable
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput, wrap_half_pi
from cranberry.offline.tables import Table

__all__ = ["PulsedOpoSource", "make_pulsed_opo_source", "pump_cavity", "transfer_matrix", "schmidt_number",
           "DISCRETIZATIONS"]

KIND = "pulsed_opo_source"
TOGGLES = ("pulsed", "pump_resonant", "detuning", "include_initial_vacuum")
CAVITY_MODELS = ("linewidth", "fsr_finesse")
PUMP_SHAPES = ("rect", "envelope")
DISCRETIZATIONS = ("box", "van_loan")
#: window cut (mode energy beyond the window) above which ``tail`` warns ``source_window_short``
WINDOW_CUT_WARN = 1e-9


def pump_cavity(s: np.ndarray, dt: float, kappa_p: float | None, delta_p: float) -> np.ndarray:
    """Bin averages of b(t) for a piecewise-constant input envelope s (bins of width dt from b = 0).
    ``kappa_p`` None: non-resonant pump, b = s."""
    s = np.asarray(s, dtype=complex)
    if kappa_p is None:
        return s.copy()
    z = kappa_p + 1j * delta_p
    ez, avg = np.exp(-z * dt), (1.0 - np.exp(-z * dt)) / (z * dt)
    out = np.empty_like(s)
    b = 0j
    for k, sk in enumerate(s):
        bss = kappa_p * sk / z
        out[k] = bss + (b - bss) * avg
        b = bss + (b - bss) * ez
    return out


def _bin_noise_van_loan(m: np.ndarray, dt: float, kt: float, kl: float) -> np.ndarray:
    """L (6×6) with L Lᵀ = exact covariance of (ΔX, ΔI, Ξ) over one bin (``discretization: van_loan``):
    ΔX = X(t_k+1) − e^{Mδt}X(t_k), ΔI = δt^{−1}[∫_bin X dt − (∫e^{Ms}ds)X(t_k)] (bin average, so all
    entries are O(κδt) and the factorisation keeps full precision), Ξ = δt^{−1/2}∫_bin dW_in, driven by the
    white input and loss noise. Augmented S = (X, I, Ξ), dS = A S dt + G dW; Van Loan (1978):
    expm([[−A, GGᵀ], [0, Aᵀ]]δt) = [[·, F12], [0, F22]] ⇒ Q = F22ᵀ F12."""
    a = np.zeros((6, 6))
    a[:2, :2], a[2:4, :2] = m, np.eye(2) / dt
    g = np.zeros((6, 4))
    g[:2, :2], g[:2, 2:] = np.sqrt(2.0 * kt) * np.eye(2), np.sqrt(2.0 * kl) * np.eye(2)
    g[4:, :2] = np.eye(2) / np.sqrt(dt)
    big = np.zeros((12, 12))
    big[:6, :6], big[:6, 6:], big[6:, 6:] = -a, g @ g.T, a.T
    f = expm(big * dt)
    q = f[6:, 6:].T @ f[:6, 6:]
    lam, vec = np.linalg.eigh(0.5 * (q + q.T))
    return vec * np.sqrt(np.clip(lam, 0.0, None))


def transfer_matrix(eps: np.ndarray, dt: float, kappa_c: float, rho: float, delta_s: float,
                    include_initial_vacuum: bool, discretization: str = "box") -> np.ndarray:
    """T (2N × ncol): bin-averaged output quadratures (xxpp) = T·w, w standard normal (module
    docstring). ``box``: columns 4k..4k+3 = (q_in, p_in, q_L, p_L) of bin k (box projection of the
    inputs, physical to O((κδt)²)); ``van_loan``: columns 6k..6k+5 = the exact in-bin noise of
    (ΔX, ΔI, Ξ) (:func:`_bin_noise_van_loan`), so C = TTᵀ is physical to rounding. Last two columns:
    initial cavity vacuum."""
    if discretization not in DISCRETIZATIONS:
        raise ValueError(f"discretization must be one of {DISCRETIZATIONS}, got {discretization!r}")
    vl = discretization == "van_loan"
    n = eps.size
    kt, kl = rho * kappa_c, (1.0 - rho) * kappa_c
    per = 6 if vl else 4
    ncol = per * n + 2
    t = np.zeros((2 * n, ncol))
    a = np.zeros((2, ncol))  # intracavity X(t_k) as a linear form in w
    if include_initial_vacuum:
        a[0, per * n] = a[1, per * n + 1] = 1.0
    sq, cy = np.sqrt(dt), np.sqrt(2.0 * kt / dt)
    fin, fl = np.sqrt(2.0 * kt) / sq, np.sqrt(2.0 * kl) / sq
    big = np.zeros((6, 6))
    big[:2, 2:4] = big[2:4, 4:6] = np.eye(2) * dt
    cache: dict[complex, tuple[np.ndarray, ...]] = {}
    for k in range(n):
        e = complex(eps[k])
        if e not in cache:
            mk = np.array([[-kappa_c + e.real, delta_s + e.imag], [-delta_s + e.imag, -kappa_c - e.real]])
            big[:2, :2] = mk * dt
            x = expm(big)
            cache[e] = (x[:2, :2].copy(), x[:2, 2:4].copy(), x[:2, 4:6].copy(),  # e^{Mδt}, ∫e^{Ms}, ∫∫
                        _bin_noise_van_loan(mk, dt, kt, kl) if vl else None)
        em, p1, p2, lk = cache[e]
        c = slice(per * k, per * k + per)
        if vl:
            y = cy * (p1 @ a)
            y[:, c] += cy * dt * lk[2:4] - lk[4:6]
            t[k], t[n + k] = y[0], y[1]
            a = em @ a
            a[:, c] += lk[:2]
            continue
        y = cy * (p1 @ a)
        f = np.array([[fin, 0.0, fl, 0.0], [0.0, fin, 0.0, fl]])  # forcing per unit w of this bin
        y[:, c] += cy * (p2 @ f)
        y[0, 4 * k] -= 1.0
        y[1, 4 * k + 1] -= 1.0
        t[k], t[n + k] = y[0], y[1]
        a = em @ a
        a[:, c] += p1 @ f
    return t


def schmidt_number(r: np.ndarray) -> float:
    """K = (Σ sinh² r_n)² / Σ sinh⁴ r_n (T17 §1; 1/K = purity of Liu SI Eq. S10)."""
    s2 = np.sinh(np.asarray(r, dtype=float)) ** 2
    return float(s2.sum() ** 2 / np.sum(s2**2)) if s2.sum() > 0 else 1.0


class PulsedOpoSource(BaseBlock):
    """The ``pulsed_opo_source`` block (module docstring)."""

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "cavity_model", "escape_efficiency", "x_ss", "phi_opa0_rad", "delta_s_rad_s", "delta_p_rad_s",
                "segment", block=name)
        super().__init__(KIND, name, params)
        self.chk_log = ChkLog()
        tg = self._toggles
        cm = params["cavity_model"]
        if cm not in CAVITY_MODELS:
            raise ValueError(f"block {name!r}: cavity_model must be one of {CAVITY_MODELS}, got {cm!r}")
        require(params, *{"linewidth": ("fwhm_hz",), "fsr_finesse": ("fsr_hz", "finesse")}[cm], block=name)
        if not 0.0 < float(params["escape_efficiency"]) <= 1.0 or float(params["x_ss"]) < 0.0:
            raise ValueError(f"block {name!r}: need 0 < escape_efficiency <= 1 and x_ss >= 0")
        if not tg["pulsed"]:
            require(params, "f_s_hz", block=name)
            return
        require(params, "pump_shape", "dt_kappa", "window_kappa", "rep_rate_hz", "angle_table_x_ss", "angle_probe_d",
                "discretization", block=name)
        if params["discretization"] not in DISCRETIZATIONS:
            raise ValueError(f"block {name!r}: discretization must be one of {DISCRETIZATIONS}")
        if params["pump_shape"] not in PUMP_SHAPES:
            raise ValueError(f"block {name!r}: pump_shape must be one of {PUMP_SHAPES}")
        require(params, {"rect": "pump_duration_s", "envelope": "pump_envelope"}[params["pump_shape"]], block=name)
        if tg["pump_resonant"]:
            require(params, "pump_fwhm_hz", block=name)
        if not 0.0 < float(params["dt_kappa"]) <= 0.05:
            raise ValueError(f"block {name!r}: dt_kappa must be in (0, 0.05] (δt ≤ 1/(20κ_c); PLAN §6.4 uses 1/40–1/160)")
        if not float(params["window_kappa"]) > 0 or not float(params["rep_rate_hz"]) > 0:
            raise ValueError(f"block {name!r}: need window_kappa > 0 and rep_rate_hz > 0")

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    # --------------------------------------------------------------------------------------------- rates

    @property
    def kappa_c(self) -> float:
        """Signal field decay rate κ_c (HWHM) [rad/s]."""
        p = self.params
        fwhm = float(p["fwhm_hz"]) if p["cavity_model"] == "linewidth" else float(p["fsr_hz"]) / float(p["finesse"])
        return np.pi * fwhm

    @property
    def kappa_p(self) -> float | None:
        """Pump field decay rate κ_p [rad/s] (None: non-resonant pump)."""
        return np.pi * float(self.params["pump_fwhm_hz"]) if self._toggles["pump_resonant"] else None

    @property
    def rho(self) -> float:
        return float(self.params["escape_efficiency"])

    def _detunings(self) -> tuple[float, float]:
        if not self._toggles["detuning"]:
            return 0.0, 0.0
        return float(self.params["delta_s_rad_s"]), float(self.params["delta_p_rad_s"])

    # --------------------------------------------------------------------------------------------- grid

    @cached_property
    def dt_s(self) -> float:
        return float(self.params["dt_kappa"]) / self.kappa_c

    @cached_property
    def grid(self) -> np.ndarray:
        """Bin centres [s] (pump turn-on at t = 0)."""
        n = int(round(float(self.params["window_kappa"]) / float(self.params["dt_kappa"])))
        return (np.arange(n) + 0.5) * self.dt_s

    @cached_property
    def envelope(self) -> np.ndarray:
        """Input pump envelope s per bin (bin average for ``rect``; centre sample for ``envelope``)."""
        p, dt, t = self.params, self.dt_s, self.grid
        if p["pump_shape"] == "rect":
            lo = t - 0.5 * dt
            return np.clip((np.minimum(lo + dt, float(p["pump_duration_s"])) - lo) / dt, 0.0, 1.0)
        return np.asarray(as_callable(p["pump_envelope"], ("t_s",))(t), dtype=float) * np.ones_like(t)

    def pump_b(self, delta_p: float) -> np.ndarray:
        return pump_cavity(self.envelope, self.dt_s, self.kappa_p, delta_p)

    def static_pump_phase(self, delta_p: float) -> float:
        """Steady-state pump-cavity phase −arctan(Δ_p/κ_p) (0 if non-resonant)."""
        return 0.0 if self.kappa_p is None else float(-np.arctan(delta_p / self.kappa_p))

    def transfer(self, *, x_ss: float | None = None, delta_s: float = 0.0, delta_p: float = 0.0,
                 rho: float | None = None, remove_static_pump_phase: bool = False) -> np.ndarray:
        """T matrix of :func:`transfer_matrix` for the given operating point (defaults: block values)."""
        x = float(self.params["x_ss"]) if x_ss is None else float(x_ss)
        b = self.pump_b(delta_p)
        if remove_static_pump_phase:
            b = b * np.exp(-1j * self.static_pump_phase(delta_p))
        eps = self.kappa_c * x * b * np.exp(1j * float(self.params["phi_opa0_rad"]))
        return transfer_matrix(eps, self.dt_s, self.kappa_c, self.rho if rho is None else float(rho), delta_s,
                               self._toggles["include_initial_vacuum"], str(self.params["discretization"]))

    # --------------------------------------------------------------------------------------------- products

    @cached_property
    def cov_kernel(self) -> np.ndarray:
        """Shot-noise-normalised covariance kernel C (2N × 2N, xxpp, vacuum = 𝟙) at the operating point."""
        ds, dp = self._detunings()
        t = self.transfer(delta_s=ds, delta_p=dp, remove_static_pump_phase=True)
        c = t @ t.T
        return 0.5 * (c + c.T)

    @cached_property
    def source_output(self) -> SourceOutput:
        p = self.params
        if not self._toggles["pulsed"]:
            return self._cw_block.source_output
        n = self.grid.size
        modes = np.eye(n) / np.sqrt(self.dt_s)
        seg = str(p["segment"])
        _, dp = self._detunings()
        so = SourceOutput(basis="temporal", grid=self.grid, cov=0.5 * self.cov_kernel, phi_opa0_rad=float(p["phi_opa0_rad"]),
                          modes=modes, pump_path_phase_rad={seg: self.static_pump_phase(dp)},
                          coef_tables=self._coef_tables(), provenance={})
        u, vmin = tm.optimal_lo_mode(so)
        trep = 1.0 / float(p["rep_rate_hz"])
        tail = tm.tail_energy(u, self.grid, trep) if trep < self.window_s else self.tail(trep)
        prov = {"block": self.name, "model": "pulsed_langevin", "kappa_c_rad_s": self.kappa_c, "kappa_p_rad_s": self.kappa_p,
                "escape_efficiency": self.rho, "x_ss": float(p["x_ss"]), "dt_s": self.dt_s, "n_bins": n,
                "R_minus": 2.0 * vmin, "tau_eff_s": tm.tau_eff_s(u, self.grid), "toggles": dict(self._toggles),
                "sources": ["T17 §2", "D55-4"]}
        return SourceOutput(basis="temporal", grid=self.grid, cov=so.cov, phi_opa0_rad=so.phi_opa0_rad, modes=modes,
                            pump_path_phase_rad=so.pump_path_phase_rad, coef_tables=so.coef_tables,
                            tail_into_next_slot=min(max(tail, 0.0), 1.0), provenance=prov)

    @cached_property
    def _cw_block(self):
        p = self.params
        return make_opo_cw_spectrum({
            "cavity_model": "linewidth", "fwhm_hz": self.kappa_c / np.pi, "escape_efficiency": self.rho,
            "pump_mode": "x", "x": float(p["x_ss"]), "detuning_rad_s": float(p["delta_s_rad_s"]), "f_s_hz": p["f_s_hz"],
            "phi_opa0_rad": float(p["phi_opa0_rad"]), "segment": str(p["segment"]), "gain_table_grid": None,
            "coef_table_x": None,
            "toggles": {"frequency_dependence": True, "pump_dependent_loss": False, "detuning": self._toggles["detuning"]},
        })

    @cached_property
    def optimal_mode(self) -> tuple[np.ndarray, float]:
        """(u_opt(t) on ``grid``, R_−): LO mode with the best squeezing and its variance / shot noise."""
        u, vmin = tm.optimal_lo_mode(self.source_output)
        return u, 2.0 * vmin

    @property
    def r_minus(self) -> float:
        return float(np.linalg.eigvalsh(self.cov_kernel)[0])

    @property
    def r_plus(self) -> float:
        return float(np.linalg.eigvalsh(self.cov_kernel)[-1])

    @cached_property
    def supermode_r(self) -> np.ndarray:
        """r_n (descending) from the eigenvalues e^{±2r_n} of C at ρ = 1 (T17 §1)."""
        if self.rho == 1.0:
            c = self.cov_kernel
        else:
            ds, dp = self._detunings()
            t = self.transfer(delta_s=ds, delta_p=dp, rho=1.0, remove_static_pump_phase=True)
            c = t @ t.T
        lam = np.linalg.eigvalsh(0.5 * (c + c.T))
        return np.sort(0.5 * np.log(lam[lam > 1.0]))[::-1]

    @property
    def schmidt_number(self) -> float:
        return schmidt_number(self.supermode_r)

    @property
    def tau_eff_s(self) -> float:
        return tm.tau_eff_s(self.optimal_mode[0], self.grid)

    @property
    def window_s(self) -> float:
        """End of the computation window [s] (``window_kappa``/κ_c from the pump turn-on)."""
        return float(self.grid[-1] + 0.5 * self.dt_s)

    @property
    def window_cut(self) -> float:
        """Mode energy the window cuts off, 8.8 e^{−2κ_c·window} (T17 §5.2); > :data:`WINDOW_CUT_WARN` ⇒
        the window is short (the computed mode itself is truncated)."""
        return tail_analytic(self.kappa_c, self.window_s)

    def tail(self, t_s: float) -> float:
        """ε_tail(T): energy of the optimal mode after T from the pump turn-on. Beyond the window (T ≥
        ``window_s``) the mode is not computed, so the analytic 8.8 e^{−2κ_c T} (T17 §5.2,
        ``chk.rep_rate_tail.tail_analytic``) is returned instead of 0 (M5-06: a short window must not hide
        a tail); if the window also cuts a non-negligible part of the mode (:attr:`window_cut`) a
        ``source_window_short`` CHK warning goes to ``chk_log``."""
        if t_s < self.window_s:
            return tm.tail_energy(self.optimal_mode[0], self.grid, t_s)
        ext = tail_analytic(self.kappa_c, t_s)
        if self.window_cut > WINDOW_CUT_WARN:
            self.chk_log.warn("source_window_short", self.name,
                              f"T = {t_s:.3g} s is beyond the source window {self.window_s:.3g} s "
                              f"(window_kappa = {self.params['window_kappa']}, cuts {self.window_cut:.2g} of the "
                              f"mode); tail extrapolated 8.8 e^(-2 κ_c T) = {ext:.3g} — enlarge window_kappa",
                              t_s=float(t_s), window_s=self.window_s, tail_extrapolated=ext)
        return ext

    # --------------------------------------------------------------------------------------------- detuning

    def _lo_angle(self, t: np.ndarray, u: np.ndarray) -> float:
        n = self.grid.size
        a = np.zeros((2, 2 * n))
        c = np.conj(u) * np.sqrt(self.dt_s)  # c_k = ∫u* f_k
        a[0, :n], a[0, n:], a[1, :n], a[1, n:] = c.real, -c.imag, c.imag, c.real
        m = a @ t
        return gaussian.noise_levels_from_cov(m @ m.T)[2]

    def angle_coefficients(self, x_ss: float) -> tuple[float, float, float]:
        """(∂θ_sq/∂Δ_s, ∂θ_sq/∂Δ_p, ∂ arg∫b dt/∂Δ_p) [rad per rad/s] at x_ss (module docstring)."""
        d, kc = float(self.params["angle_probe_d"]), self.kappa_c
        t0 = self.transfer(x_ss=x_ss)
        c0 = t0 @ t0.T
        so = SourceOutput(basis="temporal", grid=self.grid, cov=0.25 * (c0 + c0.T), modes=np.eye(self.grid.size) / np.sqrt(self.dt_s),
                          phi_opa0_rad=float(self.params["phi_opa0_rad"]))
        u0, _ = tm.optimal_lo_mode(so)
        th0 = self._lo_angle(t0, u0)
        cs = float(wrap_half_pi(self._lo_angle(self.transfer(x_ss=x_ss, delta_s=d * kc), u0) - th0)) / (d * kc)
        cp = float(wrap_half_pi(self._lo_angle(self.transfer(x_ss=x_ss, delta_p=d * kc), u0) - th0)) / (d * kc)
        cb = float(np.angle(self.pump_b(d * kc).sum() / self.pump_b(0.0).sum())) / (d * kc)
        return cs, cp, cb

    def _coef_tables(self) -> dict[str, Table]:
        grid = self.params["angle_table_x_ss"]
        xs = np.array([float(self.params["x_ss"])]) if grid is None else np.asarray(grid, dtype=float)
        vals = np.array([self.angle_coefficients(x) for x in xs])
        seg = str(self.params["segment"])
        prov = {"block": self.name, "driver": "detuning", "columns": ["per_delta_s", "per_delta_p"],
                "probe_d": float(self.params["angle_probe_d"]), "sources": ["T17 §2.5"]}
        sq = np.column_stack([vals[:, 0], vals[:, 1]])
        pu = np.column_stack([np.zeros(xs.size), vals[:, 2]])
        return {f"{seg}.detuning.squeezed": Table(f"{seg}.detuning.squeezed", {"x_ss": xs}, sq, "rad/(rad/s)", prov),
                f"{seg}.detuning.pump": Table(f"{seg}.detuning.pump", {"x_ss": xs}, pu, "rad/(rad/s)", prov)}

    def angle_map_params(self, beams: Mapping[str, str]) -> dict[str, Any]:
        """``x_grid`` and ``coef`` params of ``opo_angle_map`` for {beam: role}, role ∈ {squeezed, pump}.

        ``pulsed`` off (M5-08): the CW tables of ``opo_cw_spectrum`` (axis ``x``, one Δ_s coefficient per
        role ``squeezed`` | ``probe``) are used with c_p = 0 — the CW limit has no pump cavity, so Δ_p does
        not act and role ``pump`` (transmitted pump phase) is not available (ValueError)."""
        tabs, seg = self.source_output.coef_tables, str(self.params["segment"])
        if self._toggles["pulsed"]:
            x = tabs[f"{seg}.detuning.squeezed"].grid["x_ss"]
            return {"segment": seg, "x_grid": x.tolist(),
                    "coef": {b: tabs[f"{seg}.detuning.{role}"].values.tolist() for b, role in beams.items()}}
        bad = {r for r in beams.values() if r not in ("squeezed", "probe")}
        if bad:
            raise ValueError(f"block {self.name!r}: pulsed off (CW limit) has roles squeezed/probe only, got {sorted(bad)}")
        x = tabs[f"{seg}.detuning.squeezed"].grid["x"]
        return {"segment": seg, "x_grid": x.tolist(),
                "coef": {b: [[float(c), 0.0] for c in np.ravel(tabs[f"{seg}.detuning.{role}"].values)]
                         for b, role in beams.items()}}


@register_block(KIND)
def make_pulsed_opo_source(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``pulsed_opo_source`` block (module docstring)."""
    return PulsedOpoSource(params, name)
