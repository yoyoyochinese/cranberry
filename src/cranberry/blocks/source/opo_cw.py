"""``opo_cw_spectrum``: CW below-threshold degenerate OPO (T02). Owner: WP-B6.

Offline (O) block. Product: :attr:`OpoCwSpectrum.source_output` (``cw_sideband`` ``SourceOutput``: one
2×2 covariance per sideband frequency f_s) with an A± gain table for the error signals (T14) and the
OPO segment × beam detuning coefficient tables (D47, G16-a, VT-23). The time loop never steps it.

Model (T02 §2.1; conv §2, §4; symbols of conv §4 / T02 §1)
    Quantum Langevin equation of the intracavity mode in the frame rotating at ω0, H_0 = Δ â†â,
    H_int = i(κ_c/2)(β_p â†² − h.c.) at φ_p = 0 (ε = κ x, x = √(P/P_th))::

        dâ/dt = −(κ + iΔ) â + ε â† + √(2κ_T) â_in + √(2κ_L) b̂_in,    â_out = √(2κ_T) â − â_in

    Fourier (conv §8, d/dt → −iΩ): M(Ω) = [[κ − ε − iΩ, −Δ], [Δ, κ + ε − iΩ]] on (q, p),
    A = 2κ_T M⁻¹ − 𝟙, B = 2√(κ_T κ_L) M⁻¹, shot-noise-normalised spectral matrix S = A A† + B B†
    (:func:`langevin_spectrum`). The homodyne spectrum at angle θ is u_θᵀ Re(S) u_θ, so
    ``cov`` = Re(S)/2 (vacuum 1/2, conv §2). Δ = 0 reduces to the closed form (Ref B Eq. 2 corrected,
    T02 §2.2, :func:`closed_form`): R± = 1 ± ρ·4x/((1 ∓ x)² + (Ω/κ)²).

    Literature → conv (T02 §1, §3.2, D48): Ref B's γ = c(T + L)/l is the FWHM in rad/s = 2κ and its
    Ω = 2πf/γ (the 2π is required, V-B4), i.e. 4Ω_B² = (Ω/κ)². Ref B's L is ``loss_ic`` (L_ic).

Frame (D55-4): ``cov`` is at φ_p = 0 in the lab frame; the amplified quadrature of the H_int above is
θ_amp,0 = φ_OPA,0/2, so the Langevin result (amplified q) is rotated by φ_OPA,0/2 and θ_sq,0 =
φ_OPA,0/2 + π/2. A static detuning Δ0 (toggle ``detuning``) rotates the ellipse: θ_off(f_s) of
:meth:`SourceOutput.project` is δθ_sq(Δ0, x, f_s) ≈ −(Δ0/κ)/(1 + x²) for Ω ≪ κ (T02 §2.4).

Cavity (``cavity_model``, no default, D34)
    * ``mirror``   : κ_T = c·T/(2l), κ_L = c·L_ic/(2l) (good-cavity limit, T02 §4; ``T_oc``,
                     ``round_trip_m``, ``loss_ic`` = number or callable of ``x``, e.g. Ref B Eq. (3)
                     L_ic(x) = 0.00249 + 0.00222x). ρ = T/(T + L_ic), FWHM = κ/π.
    * ``linewidth``: κ = π·FWHM, κ_T = ρκ (``fwhm_hz``, ``escape_efficiency`` ρ; Ref C gives the fitted
                     1.75 GHz FWHM directly, T02 §4).

Detuning coefficients (T02 §2.4; first order in Δ at Ω → 0, independent of ρ)
    squeezing angle  dθ_sq/dΔ  = −1/(κ(1 + x²))       (:func:`squeezing_angle_coef`)
    transmitted probe dφ_pr/dΔ = −1/(2κ)  (at the lock-1 extremum, either sign)
    so after lock 2 (LO on the probe) the residual φ_LO − θ_sq = (Δ/κ)(1 − x²)/(2(1 + x²)). Unequal
    segment × beam coefficients on the OPO segment: the path model computes this residual, nothing
    is assumed (D11, D47, VT-23). ``coef_tables`` = Tables ``<segment>.detuning.squeezed`` and
    ``<segment>.detuning.probe`` [rad per rad/s] on axis ``x``; :meth:`coef_entries` turns them into
    topology ``coef`` entries (driver ``detuning``, phase; group delay 0).

Gain table (T14 §1, D51-10): ``gain_table`` = Table ``<name>.gain`` with values (A+, A−) =
(1/(1 − x), 1/(1 + x)) (quasi-static probe gains, ``err.amplitude_gains_opo``) on axis ``pump_w``
(``pump_mode: power``) or ``x`` (``pump_mode: x``).

Losses (G06): the intracavity loss is inside ``cov`` through ρ (not a lumped element downstream), so
``losses()`` is empty. Propagation, mode matching and detector efficiencies (Ref B ζ, ξ², η) belong to
the blocks on the beam path and are applied once by the LossRegistry.

Params (D34: all required where they apply)
    ``cavity_model``; mirror: ``T_oc``, ``round_trip_m``, ``loss_ic``; linewidth: ``fwhm_hz``,
    ``escape_efficiency``. ``pump_mode`` (power | x); power: ``pump_w``, ``p_th_w``; x: ``x``.
    ``detuning_rad_s`` Δ0, ``f_s_hz`` (ascending sideband grid), ``phi_opa0_rad``, ``segment`` (OPO
    segment id), ``gain_table_grid`` (ascending pump-axis grid or None), ``coef_table_x`` (ascending x
    grid or None = the operating point only), ``toggles``:
      frequency_dependence off → every f_s evaluated at Ω = 0 (T02 §8)
      pump_dependent_loss  off → L_ic evaluated at x = 0 (pump-off loss)
      detuning             off → Δ0 = 0
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from functools import cached_property
from typing import Any

import numpy as np

from cranberry.blocks.err._common import amplitude_gains_opo
from cranberry.config.resolve import as_callable
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput
from cranberry.offline.tables import Table

__all__ = ["C_M_PER_S", "OpoCwSpectrum", "make_opo_cw_spectrum", "langevin_spectrum", "closed_form",
           "squeezing_angle_coef", "probe_phase_coef", "residual_coef"]

KIND = "opo_cw_spectrum"
C_M_PER_S = 299_792_458.0
CAVITY_MODELS = ("mirror", "linewidth")
PUMP_MODES = ("power", "x")
TOGGLES = ("frequency_dependence", "pump_dependent_loss", "detuning")


def langevin_spectrum(x: float, kappa_t: float, kappa_l: float, delta: float, omega) -> np.ndarray:
    """S(Ω) = A A† + B B† (shot-noise units, φ_p = 0 frame, (q, p)); shape (n_Ω, 2, 2), complex."""
    kappa = kappa_t + kappa_l
    eps = x * kappa
    om = np.atleast_1d(np.asarray(omega, dtype=float))
    m = np.empty((om.size, 2, 2), dtype=complex)
    m[:, 0, 0] = kappa - eps - 1j * om
    m[:, 0, 1] = -delta
    m[:, 1, 0] = delta
    m[:, 1, 1] = kappa + eps - 1j * om
    minv = np.linalg.inv(m)
    a = 2.0 * kappa_t * minv - np.eye(2)
    b = 2.0 * np.sqrt(kappa_t * kappa_l) * minv
    h = lambda z: np.conj(np.swapaxes(z, -1, -2))  # noqa: E731
    return a @ h(a) + b @ h(b)


def closed_form(x, eta, omega_over_kappa) -> tuple[np.ndarray, np.ndarray]:
    """(R−, R+) = 1 ∓ η·4x/((1 ± x)² + (Ω/κ)²) at Δ = 0 (T02 §2.2, Ref B Eq. 2 with Ω_B = Ω/(2κ))."""
    x, eta, w2 = np.asarray(x, float), np.asarray(eta, float), np.asarray(omega_over_kappa, float) ** 2
    return 1.0 - eta * 4.0 * x / ((1.0 + x) ** 2 + w2), 1.0 + eta * 4.0 * x / ((1.0 - x) ** 2 + w2)


def squeezing_angle_coef(x: float, kappa: float) -> float:
    """dθ_sq/dΔ = −1/(κ(1 + x²)) [rad per rad/s] (T02 §2.4, Ω ≪ κ)."""
    return -1.0 / (kappa * (1.0 + x * x))


def probe_phase_coef(kappa: float) -> float:
    """dφ_pr/dΔ = −1/(2κ) of the transmitted co-propagating probe at the lock-1 extremum (T02 §2.4)."""
    return -0.5 / kappa


def residual_coef(x: float, kappa: float) -> float:
    """d(φ_pr − θ_sq)/dΔ = (1 − x²)/(2κ(1 + x²)): residual after lock 2 (T02 §2.4)."""
    return probe_phase_coef(kappa) - squeezing_angle_coef(x, kappa)


def _ascending(v: Any, what: str, block: str) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    if a.ndim != 1 or a.size < 1 or np.any(np.diff(a) <= 0):
        raise ValueError(f"block {block!r}: {what} must be a non-empty ascending list")
    return a


class OpoCwSpectrum(BaseBlock):
    """The ``opo_cw_spectrum`` block (module docstring)."""

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "cavity_model", "pump_mode", "detuning_rad_s", "f_s_hz", "phi_opa0_rad", "segment",
                "gain_table_grid", "coef_table_x", block=name)
        cm, pm = params["cavity_model"], params["pump_mode"]
        if cm not in CAVITY_MODELS:
            raise ValueError(f"block {name!r}: cavity_model must be one of {CAVITY_MODELS}, got {cm!r}")
        if pm not in PUMP_MODES:
            raise ValueError(f"block {name!r}: pump_mode must be one of {PUMP_MODES}, got {pm!r}")
        require(params, *{"mirror": ("T_oc", "round_trip_m", "loss_ic"),
                          "linewidth": ("fwhm_hz", "escape_efficiency")}[cm], block=name)
        require(params, *{"power": ("pump_w", "p_th_w"), "x": ("x",)}[pm], block=name)
        super().__init__(KIND, name, params)
        self.cavity_model, self.pump_mode = cm, pm
        _ascending(params["f_s_hz"], "f_s_hz", name)
        if np.min(params["f_s_hz"]) < 0:
            raise ValueError(f"block {name!r}: f_s_hz must be >= 0")
        for key in ("gain_table_grid", "coef_table_x"):
            if params[key] is not None:
                _ascending(params[key], key, name)
        if cm == "mirror" and not (float(params["T_oc"]) > 0 and float(params["round_trip_m"]) > 0):
            raise ValueError(f"block {name!r}: need T_oc > 0 and round_trip_m > 0")
        if cm == "linewidth" and not (float(params["fwhm_hz"]) > 0 and 0 < float(params["escape_efficiency"]) <= 1):
            raise ValueError(f"block {name!r}: need fwhm_hz > 0 and 0 < escape_efficiency <= 1")
        if pm == "power" and not (float(params["p_th_w"]) > 0 and float(params["pump_w"]) >= 0):
            raise ValueError(f"block {name!r}: need p_th_w > 0 and pump_w >= 0")
        self._check_x(self.x)

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _check_x(self, x: float) -> None:
        if not 0.0 <= x < 1.0:
            raise ValueError(f"block {self.name!r}: OPO must be below threshold, 0 <= x < 1, got x = {x}")

    # ---------------------------------------------------------------------------------------------

    @property
    def x(self) -> float:
        """Normalised pump amplitude x = √(P/P_th) (or the given x)."""
        p = self.params
        return float(p["x"]) if self.pump_mode == "x" else float(np.sqrt(float(p["pump_w"]) / float(p["p_th_w"])))

    def loss_ic(self, x: float) -> float:
        """Intracavity round-trip loss L_ic (mirror model; at x = 0 when pump_dependent_loss is off)."""
        xe = x if self._toggles["pump_dependent_loss"] else 0.0
        val = float(as_callable(self.params["loss_ic"], ("x",))(xe))
        if not 0.0 <= val < 1.0:
            raise ValueError(f"block {self.name!r}: loss_ic must be in [0, 1), got {val}")
        return val

    def rates(self, x: float) -> tuple[float, float]:
        """(κ_T, κ_L) amplitude decay rates [rad/s] at pump x."""
        p = self.params
        if self.cavity_model == "linewidth":
            kappa, rho = np.pi * float(p["fwhm_hz"]), float(p["escape_efficiency"])
            return rho * kappa, (1.0 - rho) * kappa
        c2l = C_M_PER_S / (2.0 * float(p["round_trip_m"]))
        return c2l * float(p["T_oc"]), c2l * self.loss_ic(x)

    def kappa(self, x: float) -> float:
        """Total amplitude decay rate κ = HWHM [rad/s] at pump x."""
        return sum(self.rates(x))

    def spectrum(self, x: float | None = None) -> np.ndarray:
        """Covariances (n_f, 2, 2) in the D55 frame (module docstring) at pump x (default: operating)."""
        x = self.x if x is None else float(x)
        self._check_x(x)
        kt, kl = self.rates(x)
        delta = float(self.params["detuning_rad_s"]) if self._toggles["detuning"] else 0.0
        fs = np.asarray(self.params["f_s_hz"], dtype=float)
        omega = 2.0 * np.pi * fs if self._toggles["frequency_dependence"] else np.zeros_like(fs)
        s = langevin_spectrum(x, kt, kl, delta, omega)
        cov = 0.5 * s.real
        cov = 0.5 * (cov + np.swapaxes(cov, -1, -2))
        return gaussian.apply_symplectic(cov, gaussian.rotation(0.5 * float(self.params["phi_opa0_rad"])))

    def _x_of(self, pump: float) -> float:
        return float(pump) if self.pump_mode == "x" else float(np.sqrt(pump / float(self.params["p_th_w"])))

    def _coef_tables(self) -> dict[str, Table]:
        grid = self.params["coef_table_x"]
        xs = np.array([self.x]) if grid is None else np.asarray(grid, dtype=float)
        seg = str(self.params["segment"])
        prov = {"block": self.name, "driver": "detuning", "effect": "phase", "sources": ["T02 §2.4"]}
        sq = np.array([squeezing_angle_coef(x, self.kappa(x)) for x in xs])
        pr = np.array([probe_phase_coef(self.kappa(x)) for x in xs])
        return {f"{seg}.detuning.squeezed": Table(f"{seg}.detuning.squeezed", {"x": xs}, sq, "rad/(rad/s)", prov),
                f"{seg}.detuning.probe": Table(f"{seg}.detuning.probe", {"x": xs}, pr, "rad/(rad/s)", prov)}

    def coef_entries(self, squeezed: Iterable[str], probes: Iterable[str]) -> list[dict[str, Any]]:
        """Topology ``coef`` entries of the OPO segment for driver ``detuning`` (module docstring)."""
        tabs, seg = self.source_output.coef_tables, str(self.params["segment"])
        out = []
        for role, beams in (("squeezed", squeezed), ("probe", probes)):
            for b in beams:
                out += [{"segment": seg, "driver": "detuning", "beam": b, "effect": "phase",
                         "value": tabs[f"{seg}.detuning.{role}"], "coords": {"x": self.x}},
                        {"segment": seg, "driver": "detuning", "beam": b, "effect": "group_delay", "value": 0.0}]
        return out

    @cached_property
    def source_output(self) -> SourceOutput:
        p = self.params
        x = self.x
        kt, kl = self.rates(x)
        kappa = kt + kl
        table = None
        if p["gain_table_grid"] is not None:
            grid = np.asarray(p["gain_table_grid"], dtype=float)
            for g in grid:
                self._check_x(self._x_of(g))
            axis = "pump_w" if self.pump_mode == "power" else "x"
            table = Table(f"{self.name}.gain", {axis: grid}, np.array([amplitude_gains_opo(self._x_of(g)) for g in grid]),
                          "1", {"block": self.name, "columns": ["A+", "A-"], "sources": ["T02 §2.2", "T14 §1", "D51-10"]})
        prov = {
            "block": self.name, "cavity_model": self.cavity_model, "x": x, "kappa_rad_s": kappa,
            "fwhm_hz": kappa / np.pi, "escape_efficiency": kt / kappa,
            "loss_ic": self.loss_ic(x) if self.cavity_model == "mirror" else None,
            "detuning_rad_s": float(p["detuning_rad_s"]) if self._toggles["detuning"] else 0.0,
            "toggles": dict(self._toggles), "sources": ["T02 §2.1–§2.4", "D48", "D55-4"],
        }
        return SourceOutput(basis="cw_sideband", grid=np.asarray(p["f_s_hz"], dtype=float), cov=self.spectrum(x),
                            phi_opa0_rad=float(p["phi_opa0_rad"]), gain_table=table, coef_tables=self._coef_tables(),
                            provenance=prov)


@register_block(KIND)
def make_opo_cw_spectrum(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opo_cw_spectrum`` block (module docstring)."""
    return OpoCwSpectrum(params, name)
