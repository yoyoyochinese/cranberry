"""``opa_waveguide``: single-pass waveguide OPA (T01). Owner: WP-B1.

Offline (O) block. Product: :attr:`OpaWaveguide.source_output` (``single_mode`` ``SourceOutput``)
with an A±(P) gain table (T14) and the crystal as a pump-path segment. The time loop never steps it
(no reads / writes, ``__call__`` is the identity).

Models (``model``, no default, D34):

* ``lumped``       lossless parametric gain; the waveguide loss L_wg(P) is a ``LossElement``
                   (``losses()``, G06) applied downstream like every other loss — the prompt §7-1 /
                   Ref A fit reading (D48).
* ``distributed``  gain and propagation loss a_s act together inside the waveguide (Yamashima-type,
                   T01 §2.3: dR±/dz = ±2gR± − a_s(R± − 1)); the loss is inside ``cov`` and *not*
                   registered (it is not equivalent to a lumped element, T01 §2.3).
* ``coupled_wave`` ``distributed`` + the classical coupled-wave equations for pump and a
                   co-propagating seed/probe (T01 §2.2), so ``pump_depletion`` can act.

Equations (power-normalised envelopes |A|² = P [W], κ = √α/ℓ, δ = Δk/2, frame B = A e^{−iΔkz/2}
of T01 §2.5, i.e. angles relative to the exit pump/2)::

    dB_s/dz = −iδ B_s + c B_s* − (a_s/2) B_s,     c(z) = iκ A_p(z)
    dA_p/dz = −(a_p/2) A_p [+ iκ B_s²  if pump_depletion]
    quadratures (conv §2): d(q, p)/dz = G (q, p),  G = [[Re c, Im c], [Im c, −Re c]] + δ[[0, 1], [−1, 0]] − (a_s/2)𝟙
    covariance:            dσ/dz = Gσ + σGᵀ + a_s V_vac 𝟙   (vacuum in, conv §2 pure loss in the limit g → 0)

Literature → conv (T01 §2.2): the coupled-wave form iκA_p puts the amplified quadrature at
(arg A_p + π/2)/2, so the input pump is set to A_p(0) = √P_in e^{i(φ_OPA,0 − π/2)}: G(0) then amplifies
θ_amp,0 = φ_OPA,0/2 and ``cov`` is in the D55 frame (rule 4: φ_p = 0, θ_sq,0 = φ_OPA,0/2 + π/2).
θ_off = θ_min − θ_sq,0 is the phase-mismatch rotation δθ_mm ≈ −Δkℓ/(4r) (sign fixed by V-T01-7).
The crystal's own pump phase ψ_p,xtal is *not* in ``cov``: it is reported as
``pump_path_phase_rad[segment]`` and enters θ_sq as ½ψ_p,xtal through the time loop (conv §4, T01 §2.5).

Gain table (T14): ``gain_table`` = Table ``<name>.gain`` on axis ``pump_w`` with values (A+, A−):
amplitude gains of a seed on the nominal amplified (θ_amp,0) / squeezed (θ_amp,0 + π/2) quadrature,
|T e_±| for the linear models (= e^{±r} for ``lumped`` without mismatch). With ``pump_depletion`` the
seed of power ``seed_power_w`` is propagated through the nonlinear equations instead.

Pump RIN (T01 §2.7): δr = r·δP/(2P). Declared as a ``rin`` band on the crystal segment (``bands()``,
G05) — never added here; ``provenance["dr_d_rin"] = r_eff/2`` is the gain-fluctuation coefficient.

Params (D34: all required where they apply)
    ``model``, ``alpha_per_w`` α [1/W] (η_SHG, or the Ref A fitted 9.06 /W tied to P = output, D48),
    ``pump_w`` P [W], ``pump_ref`` (input | output: which face P refers to), ``length_m`` ℓ,
    ``phi_opa0_rad`` φ_OPA,0, ``segment`` (crystal segment id: loss element, band target, pump-path
    key), ``psi_p_xtal_rad`` (nominal ψ_p,xtal), ``gain_table_pump_w`` (ascending list or None),
    ``toggles`` {pump_loss, phase_mismatch, pump_depletion, pump_rin}.
    lumped: ``waveguide_loss`` L_wg (number or callable of ``P_w``, e.g. pump-induced loss),
    ``loss_source`` (provenance tag of L_wg). distributed / coupled_wave: ``loss_signal_db_per_cm``.
    coupled_wave: ``seed_power_w`` (co-propagating probe on θ_amp,0; 0 for none).
    pump_loss on: ``loss_pump_db_per_cm``. phase_mismatch on: ``delta_k_per_m`` Δk.
    pump_rin on: ``pump_rin`` {f_lo_hz, f_hi_hz, psd_outside (callable of ``f_hz``) | variance_outside}.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cached_property
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cranberry.config.resolve import as_callable
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.core.units import V_VAC
from cranberry.offline import gaussian
from cranberry.offline.band_ledger import BandDecl
from cranberry.offline.loss_registry import LossElement
from cranberry.offline.source_output import SourceOutput, wrap_half_pi
from cranberry.offline.tables import Table

__all__ = ["MODELS", "PUMP_REFS", "OpaWaveguide", "make_opa_waveguide", "propagate", "db_per_cm_to_per_m"]

KIND = "opa_waveguide"
MODELS = ("lumped", "distributed", "coupled_wave")
PUMP_REFS = ("input", "output")
TOGGLES = ("pump_loss", "phase_mismatch", "pump_depletion", "pump_rin")


def db_per_cm_to_per_m(x_db_per_cm: float) -> float:
    """Power loss coefficient a [1/m] from dB/cm (a ℓ = (dB/cm)(ℓ/cm) ln10/10, T01 §1)."""
    return float(x_db_per_cm) * 100.0 * np.log(10.0) / 10.0


def propagate(
    *,
    alpha_per_w: float,
    pump_in_w: float,
    length_m: float,
    a_s_per_m: float,
    a_p_per_m: float,
    delta_k_per_m: float,
    phi_opa0_rad: float,
    seed: complex = 0j,
    depletion: bool = False,
) -> dict[str, Any]:
    """Integrate the module-docstring equations over the crystal.

    Returns {"cov": σ(ℓ) (2×2, vacuum in), "T": mean-field transfer (2×2, (q, p)),
    "seed_out": B_s(ℓ) [√W], "pump_out": A_p(ℓ) [√W]}.
    """
    kappa = np.sqrt(alpha_per_w) / length_m
    delta = 0.5 * delta_k_per_m
    ap0 = np.sqrt(pump_in_w) * np.exp(1j * (phi_opa0_rad - 0.5 * np.pi))
    rot = np.array([[0.0, 1.0], [-1.0, 0.0]])
    eye = np.eye(2)

    def rhs(_z, y):
        b, ap = y[0] + 1j * y[1], y[2] + 1j * y[3]
        c = 1j * kappa * ap
        db = -1j * delta * b + c * np.conj(b) - 0.5 * a_s_per_m * b
        dap = -0.5 * a_p_per_m * ap + (1j * kappa * b * b if depletion else 0.0)
        g = np.array([[c.real, c.imag], [c.imag, -c.real]]) + delta * rot - 0.5 * a_s_per_m * eye
        t = y[4:8].reshape(2, 2)
        s = np.array([[y[8], y[9]], [y[9], y[10]]])
        ds = g @ s + s @ g.T + a_s_per_m * V_VAC * eye
        return np.concatenate([[db.real, db.imag, dap.real, dap.imag], (g @ t).ravel(), [ds[0, 0], ds[0, 1], ds[1, 1]]])

    y0 = np.array([seed.real, seed.imag, ap0.real, ap0.imag, 1.0, 0.0, 0.0, 1.0, V_VAC, 0.0, V_VAC])
    sol = solve_ivp(rhs, (0.0, length_m), y0, method="DOP853", rtol=1e-12, atol=1e-15)
    if not sol.success:
        raise RuntimeError(f"opa_waveguide: integration failed: {sol.message}")
    y = sol.y[:, -1]
    return {
        "cov": np.array([[y[8], y[9]], [y[9], y[10]]]),
        "T": y[4:8].reshape(2, 2),
        "seed_out": complex(y[0], y[1]),
        "pump_out": complex(y[2], y[3]),
    }


class OpaWaveguide(BaseBlock):
    """The ``opa_waveguide`` block (module docstring)."""

    toggle_names = TOGGLES

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "model", "alpha_per_w", "pump_w", "pump_ref", "length_m", "phi_opa0_rad", "segment",
                "psi_p_xtal_rad", "gain_table_pump_w", block=name)
        model = params["model"]
        if model not in MODELS:
            raise ValueError(f"block {name!r}: model must be one of {MODELS}, got {model!r}")
        if params["pump_ref"] not in PUMP_REFS:
            raise ValueError(f"block {name!r}: pump_ref must be one of {PUMP_REFS}, got {params['pump_ref']!r}")
        super().__init__(KIND, name, params)
        tg = self._toggles
        require(params, *{"lumped": ("waveguide_loss", "loss_source"),
                          "distributed": ("loss_signal_db_per_cm",),
                          "coupled_wave": ("loss_signal_db_per_cm", "seed_power_w")}[model], block=name)
        if tg["pump_loss"]:
            require(params, "loss_pump_db_per_cm", block=name)
        if tg["phase_mismatch"]:
            require(params, "delta_k_per_m", block=name)
        if tg["pump_rin"]:
            require(params, "pump_rin", block=name)
            require(params["pump_rin"], "f_lo_hz", "f_hi_hz", block=f"{name}.pump_rin")
        if tg["pump_depletion"] and model != "coupled_wave":
            raise ValueError(f"block {name!r}: pump_depletion needs model 'coupled_wave' (got {model!r})")
        if float(params["alpha_per_w"]) < 0 or float(params["pump_w"]) < 0 or not float(params["length_m"]) > 0:
            raise ValueError(f"block {name!r}: need alpha_per_w >= 0, pump_w >= 0, length_m > 0")
        grid = params["gain_table_pump_w"]
        if grid is not None and (np.ndim(grid) != 1 or np.any(np.diff(grid) <= 0) or np.min(grid) < 0):
            raise ValueError(f"block {name!r}: gain_table_pump_w must be an ascending list of P >= 0, or None")
        self.model = model

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    # ---------------------------------------------------------------------------------------------

    def _coeffs(self) -> dict[str, float]:
        p, tg = self.params, self._toggles
        return {
            "a_s": 0.0 if self.model == "lumped" else db_per_cm_to_per_m(p["loss_signal_db_per_cm"]),
            "a_p": db_per_cm_to_per_m(p["loss_pump_db_per_cm"]) if tg["pump_loss"] else 0.0,
            "dk": float(p["delta_k_per_m"]) if tg["phase_mismatch"] else 0.0,
        }

    def _pump_in(self, pump_w: float, a_p: float) -> float:
        # ponytail: output reference only undoes the linear pump loss, not depletion (negligible for
        # squeezed vacuum, T01 §2.6); solve for P_in iteratively if mW probes deplete the pump.
        ell = float(self.params["length_m"])
        return pump_w * np.exp(a_p * ell) if self.params["pump_ref"] == "output" else pump_w

    def _run(self, pump_w: float, seed: complex = 0j) -> dict[str, Any]:
        c = self._coeffs()
        return propagate(
            alpha_per_w=float(self.params["alpha_per_w"]), pump_in_w=self._pump_in(pump_w, c["a_p"]),
            length_m=float(self.params["length_m"]), a_s_per_m=c["a_s"], a_p_per_m=c["a_p"],
            delta_k_per_m=c["dk"], phi_opa0_rad=float(self.params["phi_opa0_rad"]), seed=seed,
            depletion=self._toggles["pump_depletion"],
        )

    def _seed(self, quadrature_offset_rad: float) -> complex:
        """Seed (co-propagating probe) on θ_amp,0 + offset; 0 unless coupled_wave."""
        if self.model != "coupled_wave":
            return 0j
        theta = 0.5 * float(self.params["phi_opa0_rad"]) + quadrature_offset_rad
        return np.sqrt(float(self.params["seed_power_w"])) * np.exp(1j * theta)

    def gains(self, pump_w: float) -> tuple[float, float]:
        """(A+, A−): amplitude gains on the nominal amplified / squeezed quadratures (module docstring)."""
        th = 0.5 * float(self.params["phi_opa0_rad"])
        if self._toggles["pump_depletion"] and float(self.params["seed_power_w"]) > 0:
            s0 = np.sqrt(float(self.params["seed_power_w"]))
            return tuple(abs(self._run(pump_w, self._seed(o))["seed_out"]) / s0 for o in (0.0, 0.5 * np.pi))
        t = self._run(pump_w, self._seed(0.0))["T"]
        e_amp, e_sq = np.array([np.cos(th), np.sin(th)]), np.array([-np.sin(th), np.cos(th)])
        return float(np.linalg.norm(t @ e_amp)), float(np.linalg.norm(t @ e_sq))

    def waveguide_loss(self) -> float:
        """Lumped L_wg at ``pump_w`` (lumped model only)."""
        lw = float(as_callable(self.params["waveguide_loss"], ("P_w",))(float(self.params["pump_w"])))
        if not 0.0 <= lw <= 1.0:
            raise ValueError(f"block {self.name!r}: waveguide_loss must be in [0, 1], got {lw}")
        return lw

    @cached_property
    def source_output(self) -> SourceOutput:
        p = self.params
        pump_w = float(p["pump_w"])
        res = self._run(pump_w, self._seed(0.0))
        cov = res["cov"]
        sv = np.linalg.svd(res["T"], compute_uv=False)
        r_eff = 0.5 * np.log(sv[0] / sv[1])
        c = self._coeffs()
        ell = float(p["length_m"])
        vm, _, th_min = gaussian.noise_levels_from_cov(cov)
        if self.model == "lumped":
            l_wg = self.waveguide_loss()
        elif r_eff > 1e-12:
            l_wg = (vm / V_VAC - np.exp(-2 * r_eff)) / (1 - np.exp(-2 * r_eff))  # T01 §2.3 effective loss
        else:
            l_wg = 1.0 - np.exp(-c["a_s"] * ell)
        table = None
        if p["gain_table_pump_w"] is not None:
            grid = np.asarray(p["gain_table_pump_w"], dtype=float)
            table = Table(f"{self.name}.gain", {"pump_w": grid}, np.array([self.gains(x) for x in grid]), "1",
                          {"block": self.name, "model": self.model, "columns": ["A+", "A-"], "pump_ref": p["pump_ref"],
                           "sources": ["T01 §2.2", "T14 §1"]})
        prov = {
            "block": self.name, "model": self.model, "r_nominal": float(np.sqrt(float(p["alpha_per_w"]) * pump_w)),
            "r_eff": float(r_eff), "L_wg": float(l_wg), "pump_in_w": self._pump_in(pump_w, c["a_p"]),
            "pump_out_w": abs(res["pump_out"]) ** 2, "a_s_per_m": c["a_s"], "a_p_per_m": c["a_p"],
            "delta_k_per_m": c["dk"], "dr_d_rin": 0.5 * float(r_eff), "toggles": dict(self._toggles),
            "delta_theta_mm_rad": float(wrap_half_pi(th_min - 0.5 * float(p["phi_opa0_rad"]) - 0.5 * np.pi)),
        }
        return SourceOutput(basis="single_mode", grid=None, cov=cov, phi_opa0_rad=float(p["phi_opa0_rad"]),
                            pump_path_phase_rad={str(p["segment"]): float(p["psi_p_xtal_rad"])},
                            gain_table=table, provenance=prov)

    def losses(self) -> list[LossElement]:
        if self.model != "lumped":
            return []  # distributed loss lives inside cov (module docstring)
        return [LossElement(id=f"{self.name}.wg", eta=1.0 - self.waveguide_loss(), kind="component",
                            applies_to={"gbs", "bhd", "probe"}, source=str(self.params["loss_source"]),
                            enabled=True, segment=str(self.params["segment"]))]

    def bands(self) -> list[BandDecl]:
        if not self._toggles["pump_rin"]:
            return []
        rin = self.params["pump_rin"]
        kw = {}  # BandDecl checks that exactly one of the two is given
        if "psd_outside" in rin:
            kw["psd_outside"] = as_callable(rin["psd_outside"], ("f_hz",))
        if "variance_outside" in rin:
            kw["variance_outside"] = float(rin["variance_outside"])
        return [BandDecl(source=self.name, quantity="rin", targets=(str(self.params["segment"]),),
                         f_lo_hz=float(rin["f_lo_hz"]), f_hi_hz=float(rin["f_hi_hz"]), unit="1", **kw)]


@register_block(KIND)
def make_opa_waveguide(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``opa_waveguide`` block (module docstring)."""
    return OpaWaveguide(params, name)
