"""``opa_supermodes``: pulse-pumped single-pass waveguide OPA, frequency-domain supermodes (T03 §2.1–§2.2).
Owner: WP-B11.

Offline (O) block. Product: :attr:`OpaSupermodes.source_output` (``temporal`` ``SourceOutput``) plus the
supermode squeeze parameters r_n, the Schmidt number K and the supermode functions. Never stepped.

Model (T03 §2.1; conv §4, §8; frame moving with the signal group velocity, undepleted pump, type-0)::

    ∂_z b_j = i D_j b_j + Σ_k G_jk(z) b_k†,     G_jk = g (ΔΩ/2π) Ã_p(z, Ω_j + Ω_k),   g = √α / L_c
    Ã_p(z, Ω) = Ã_p(0, Ω) exp[i(Δβ1 Ω + β2,p Ω²/2) z],   D_j = β2,s Ω_j²/2
    Ã_p(0, Ω) = √P_pk √(2π) T0p e^{−Ω²T0p²/2} e^{iφ_OPA,0}    (field ∝ e^{−T²/(2T0p²)}, T0p = τ_p/(2√ln2))

on the grid Ω_j = (j − N/2)ΔΩ, ΔΩ = 2π/W (W = ``time_window_s``), b_j the discrete frequency modes
(b̃(Ω_j)·√(ΔΩ/2π)). The z-propagator of (b, b†) is a product of ``n_z`` matrix exponentials of the
Hamiltonian generator [[iD, G], [G*, −iD]] at the step midpoints, so b(L) = C b(0) + S b(0)† is
symplectic to round-off for any ``n_z`` (V2). CW check (T03 §2.1): Ã_p → 2π√P δ(Ω) gives
b(L) = b cosh r + b† e^{iφ_OPA,0} sinh r, r = √(αP): θ_sq,0 = φ_OPA,0/2 + π/2 (the D55-4 frame).
Grid resolution is the user's job: ΔΩ must resolve the pump spectrum (~1/T0p) and the GVM sinc
(~2π/τ_GVM) and NΔΩ must cover the phase-matching bandwidth (T03 §5.2: ps pumps need N ≈ 2400–3200).
Grid check (M5-05, CHK ``opa_grid_pm_coverage``, at construction): the difference-frequency phase
matching (Ω_j = −Ω_k = Ω: Δk = −β2,s Ω²) has its first sinc zero at Ω_PM = √(2π/(|β2,s| L)); the
coverage ratio (NΔΩ/2)/Ω_PM < :data:`PM_COVER_MIN` warns (ratio 0 when β2,s = 0 or ``gvd`` off: the
band is unbounded). r_0 and the leading supermode converge on a narrow grid, but **K and the total photon
number Σ sinh² r_n grow with the spectral span NΔΩ** until the band is covered (VT-32 fixture: N =
128/256/512 at W = 8 ps → K = 24/48/88), so quote K only with a covered grid (``pm_coverage``).

Output basis: the DFT modes f_j(T_n) = e^{−iΩ_j T_n}/√W on T_n = (n − N/2)·W/N (orthonormal; â_j =
∫f_j* â dT = b̃(Ω_j)·√(ΔΩ/2π), conv §8 sign). ``cov`` = R Rᵀ/2 with the real symplectic matrix
R = [[Re(C + S), −Im(C − S)], [Im(C + S), Re(C − S)]] (xxpp). The state is pure: the eigenvalues of
2·cov are e^{±2r_n} (Bloch–Messiah, T03 §2.2), the eigenvectors give the supermodes ψ_n (through the
LO-mode convention of ``temporal_modes``). K = (Σ sinh² r_n)²/Σ sinh⁴ r_n (low gain: the JSA Schmidt
number, V3). Low-gain JSA (first Born term, :func:`jsa_low_gain`)::

    F_jk = (ΔΩ/2π) Ã_p(0, Ω_j + Ω_k) e^{iD_j L} ∫_0^L exp[i(Δβ1 S + β2,p S²/2 − D_j − D_k) z] dz,  S = Ω_j + Ω_k
    so S ≈ g·F for g → 0.

Toggles (T03 §8.1): ``multimode`` off → the CW / single-mode limit r = √(αP_pk) via
``SourceOutput.single_mode_reference(loss=0)`` (bit-identical to the §7-1 model, VT-03, T03 V5);
``gvm`` off → Δβ1 = 0; ``gvd`` off → β2,s = β2,p = 0.

Pump-path segment: the crystal's own pump phase ψ_p,xtal is reported as ``pump_path_phase_rad[segment]``
(enters θ_sq as ½ψ in the time loop, conv §4) and is not in ``cov``.

Params (D34): ``length_m`` L_c (D44: the OPA length is a parameter), ``alpha_per_w`` α,
``peak_pump_w`` P_pk, ``pump_fwhm_s`` τ_p, ``gvm_s_per_m`` Δβ1, ``beta2_signal_s2_per_m``,
``beta2_pump_s2_per_m``, ``n_points`` N (even), ``time_window_s`` W, ``n_z``, ``phi_opa0_rad``,
``segment``, ``psi_p_xtal_rad``, ``toggles`` {multimode, gvm, gvd}.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np
from scipy.linalg import expm

from cranberry.blocks.source import temporal_modes as tm
from cranberry.blocks.source.pulsed_opo import schmidt_number
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block
from cranberry.offline.source_output import SourceOutput

__all__ = ["OpaSupermodes", "make_opa_supermodes", "jsa_low_gain", "schmidt_number_jsa", "PM_COVER_MIN"]

KIND = "opa_supermodes"
TOGGLES = ("multimode", "gvm", "gvd")
#: minimum (NΔΩ/2)/Ω_PM before ``opa_grid_pm_coverage`` warns (the sinc tail ∝ 1/Ω² still adds to K)
PM_COVER_MIN = 3.0


def schmidt_number_jsa(f: np.ndarray) -> float:
    """K = 1/Σλ_n², λ_n = normalised squared singular values of a (discretised) JSA."""
    s2 = np.linalg.svd(np.asarray(f), compute_uv=False) ** 2
    lam = s2 / s2.sum()
    return float(1.0 / np.sum(lam**2))


def jsa_low_gain(omega, d_omega, length_m, ap0, dbeta1, b2p, d) -> np.ndarray:
    """Low-gain JSA F (module docstring); ``ap0(Ω)`` the input pump spectrum, ``d`` = D_j."""
    s = omega[:, None] + omega[None, :]
    k = dbeta1 * s + 0.5 * b2p * s**2 - d[:, None] - d[None, :]
    kl = k * length_m
    with np.errstate(invalid="ignore", divide="ignore"):
        integ = np.where(np.abs(kl) < 1e-12, length_m, (np.exp(1j * kl) - 1.0) / (1j * np.where(k == 0, 1.0, k)))
    return d_omega / (2 * np.pi) * ap0(s) * np.exp(1j * d * length_m)[:, None] * integ


class OpaSupermodes(BaseBlock):
    """The ``opa_supermodes`` block (module docstring)."""

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "length_m", "alpha_per_w", "peak_pump_w", "phi_opa0_rad", "segment", "psi_p_xtal_rad", block=name)
        super().__init__(KIND, name, params)
        if not float(params["length_m"]) > 0 or float(params["alpha_per_w"]) < 0 or float(params["peak_pump_w"]) < 0:
            raise ValueError(f"block {name!r}: need length_m > 0, alpha_per_w >= 0, peak_pump_w >= 0")
        if self._toggles["multimode"]:
            require(params, "pump_fwhm_s", "gvm_s_per_m", "beta2_signal_s2_per_m", "beta2_pump_s2_per_m", "n_points",
                    "time_window_s", "n_z", block=name)
            n = int(params["n_points"])
            if n != params["n_points"] or n < 2 or n % 2 or int(params["n_z"]) < 1:
                raise ValueError(f"block {name!r}: n_points must be an even int >= 2 and n_z >= 1")
            if not float(params["pump_fwhm_s"]) > 0 or not float(params["time_window_s"]) > 0:
                raise ValueError(f"block {name!r}: need pump_fwhm_s > 0 and time_window_s > 0")
        self.chk_log = ChkLog()
        if self._toggles["multimode"] and self.pm_coverage < PM_COVER_MIN:
            self.chk_log.warn("opa_grid_pm_coverage", name,
                              f"spectral grid half-span NΔΩ/2 covers {self.pm_coverage:.3g}× the phase-matching "
                              f"bandwidth Ω_PM = {self.omega_pm:.3g} rad/s (< {PM_COVER_MIN:g}): the Schmidt number "
                              "K and the photon number depend on the grid span (r_0 does not)",
                              coverage=self.pm_coverage, omega_pm_rad_s=self.omega_pm)

    @property
    def omega_pm(self) -> float:
        """Ω_PM = √(2π/(|β2,s| L)) [rad/s] (module docstring; inf if β2,s = 0 or ``gvd`` off)."""
        b2s = abs(self._coeffs()[1])
        return float(np.sqrt(2 * np.pi / (b2s * float(self.params["length_m"])))) if b2s > 0 else float("inf")

    @property
    def pm_coverage(self) -> float:
        """(NΔΩ/2)/Ω_PM: how far the grid covers the phase-matching band (0 if unbounded)."""
        return float(int(self.params["n_points"]) * self.d_omega / 2 / self.omega_pm)

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    @property
    def g(self) -> float:
        return float(np.sqrt(float(self.params["alpha_per_w"])) / float(self.params["length_m"]))

    @cached_property
    def omega(self) -> np.ndarray:
        n = int(self.params["n_points"])
        return (np.arange(n) - n // 2) * self.d_omega

    @property
    def d_omega(self) -> float:
        return 2 * np.pi / float(self.params["time_window_s"])

    @cached_property
    def grid(self) -> np.ndarray:
        n, w = int(self.params["n_points"]), float(self.params["time_window_s"])
        return (np.arange(n) - n // 2) * (w / n)

    def _coeffs(self) -> tuple[float, float, float]:
        p, tg = self.params, self._toggles
        dbeta1 = float(p["gvm_s_per_m"]) if tg["gvm"] else 0.0
        b2s, b2p = (float(p["beta2_signal_s2_per_m"]), float(p["beta2_pump_s2_per_m"])) if tg["gvd"] else (0.0, 0.0)
        return dbeta1, b2s, b2p

    def pump_spectrum(self, om) -> np.ndarray:
        """Ã_p(0, Ω) [√W·s] (module docstring)."""
        p = self.params
        t0 = float(p["pump_fwhm_s"]) / (2 * np.sqrt(np.log(2)))
        return (np.sqrt(float(p["peak_pump_w"])) * np.sqrt(2 * np.pi) * t0 * np.exp(-0.5 * (om * t0) ** 2)
                * np.exp(1j * float(p["phi_opa0_rad"])))

    @cached_property
    def propagator(self) -> tuple[np.ndarray, np.ndarray]:
        """(C, S): b(L) = C b(0) + S b(0)†."""
        p = self.params
        om, n = self.omega, self.omega.size
        dbeta1, b2s, b2p = self._coeffs()
        d = 0.5 * b2s * om**2
        s = om[:, None] + om[None, :]
        g0 = self.g * self.d_omega / (2 * np.pi) * self.pump_spectrum(s)
        nz, length = int(p["n_z"]), float(p["length_m"])
        dz = length / nz
        u = np.eye(2 * n, dtype=complex)
        h = np.zeros((2 * n, 2 * n), dtype=complex)
        h[:n, :n], h[n:, n:] = np.diag(1j * d), np.diag(-1j * d)
        for i in range(nz):
            gm = g0 * np.exp(1j * (dbeta1 * s + 0.5 * b2p * s**2) * (i + 0.5) * dz)
            h[:n, n:], h[n:, :n] = gm, gm.conj()
            u = expm(h * dz) @ u
        return u[:n, :n], u[:n, n:]

    def jsa(self) -> np.ndarray:
        """Low-gain JSA F on the block grid (S ≈ g F for g → 0)."""
        dbeta1, b2s, b2p = self._coeffs()
        return jsa_low_gain(self.omega, self.d_omega, float(self.params["length_m"]), self.pump_spectrum, dbeta1, b2p,
                            0.5 * b2s * self.omega**2)

    @cached_property
    def source_output(self) -> SourceOutput:
        p = self.params
        seg = {str(p["segment"]): float(p["psi_p_xtal_rad"])}
        if not self._toggles["multimode"]:
            so = SourceOutput.single_mode_reference(loss=0.0, alpha_per_w=float(p["alpha_per_w"]),
                                                    pump_w=float(p["peak_pump_w"]), phi_opa0_rad=float(p["phi_opa0_rad"]))
            return SourceOutput(basis="single_mode", grid=None, cov=so.cov, phi_opa0_rad=so.phi_opa0_rad,
                                pump_path_phase_rad=seg, provenance={**so.provenance, "block": self.name, "multimode": False})
        c, s = self.propagator
        r = np.block([[(c + s).real, -(c - s).imag], [(c + s).imag, (c - s).real]])
        cov = 0.5 * r @ r.T
        n, w = self.omega.size, float(p["time_window_s"])
        modes = np.exp(-1j * np.outer(self.omega, self.grid)) / np.sqrt(w)
        prov = {"block": self.name, "model": "frequency_domain_supermodes", "n_points": n, "n_z": int(p["n_z"]),
                "toggles": dict(self._toggles), "sources": ["T03 §2.1–§2.2", "D55-4"]}
        return SourceOutput(basis="temporal", grid=self.grid, cov=0.5 * (cov + cov.T), phi_opa0_rad=float(p["phi_opa0_rad"]),
                            modes=modes, pump_path_phase_rad=seg, provenance=prov)

    @cached_property
    def supermode_r(self) -> np.ndarray:
        """r_n (descending) from the eigenvalues e^{±2r_n} of 2·cov (pure state)."""
        if not self._toggles["multimode"]:
            return np.array([np.sqrt(float(self.params["alpha_per_w"]) * float(self.params["peak_pump_w"]))])
        lam = np.linalg.eigvalsh(2.0 * np.asarray(self.source_output.cov))
        return np.sort(0.5 * np.log(lam[lam > 1.0]))[::-1]

    @property
    def schmidt_number(self) -> float:
        return schmidt_number(self.supermode_r)

    @property
    def leading_mode(self) -> np.ndarray:
        """ψ_0(T) on ``grid``: the most squeezed supermode as an LO mode (θ_off = 0 when projected)."""
        return tm.optimal_lo_mode(self.source_output)[0]


@register_block(KIND)
def make_opa_supermodes(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opa_supermodes`` block (module docstring)."""
    return OpaSupermodes(params, name)
