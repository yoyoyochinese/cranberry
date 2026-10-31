"""Phase error → measured quadrature variance (conv §6, T04 §3, §13, G05). Owner: WP-F4.

V′ = (V_− + V_+)/2 − (V_+ − V_−)/2 · ⟨cos 2Δ⟩, Δ(t) = Δ_set + δ(t) − θ_off (T04 eq. (1)).

χ is the product of independent parts (conv §6):

* **time loop** (in-loop band): measured directly from the residual series / its running sums,
  ⟨cos 2Δ⟩ = cos 2a·⟨cos 2δ′⟩ − sin 2a·⟨sin 2δ′⟩ with a = Δ_set, δ′ = δ − θ_off (any distribution;
  "the implementation basis", T04 §13);
* **out-of-loop band** (HF): χ_HF = e^{−2σ²_HF} with σ²_HF from the ``BandLedger``. This module is the
  **only** place where ledger variance is added (G05): :func:`quadrature_for_metric` calls
  ``BandLedger.hf_variance(metric, beams)`` once per metric and beam combination;
* **analytic dither** (``dither_m_rad``, peak amplitudes m_i): Π J0(2m_i). Use it only for a dither
  that is *not* in the time series (``dither_mode: analytic``); a dither simulated in the loop is
  already inside ⟨cos 2δ′⟩ and must not be listed again (no double counting).

θ_off is the static rotation of the LO-mode ellipse that the lock does not see (``SourceOutput.project``,
D55 rule 4); the time-loop residual δ of the metric blocks does not contain it, so it is applied here.

Variances are in conv §2 units (vacuum 1/2) as returned by ``SourceOutput.project``; outputs in dB
re shot noise (V_SN = 1/2). ``R′∓`` follow prompt §7-1: the same jitter measured around the squeezed
(Δ = 0) and the anti-squeezed (Δ = π/2) quadrature; ``V′`` is the value at the actual set point.

Integration limits (conv §6, T04 §3.4): lower = 1/(record length) (slower noise is a static offset,
already inside the measured mean); upper = the metric's ``f_max_hz`` in the ledger (flat range of the
source spectrum, ~1/τ_p for pulses, ~κ_c for OPO modes).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.core.units import V_VAC, lin_to_db
from cranberry.offline import gaussian

if TYPE_CHECKING:
    from cranberry.offline.band_ledger import BandLedger
    from cranberry.offline.loss_registry import LossRegistry
    from cranberry.offline.source_output import SourceOutput

__all__ = ["AVERAGE_MODES", "PhaseStats", "phase_to_quadrature", "quadrature_for_metric"]

AVERAGE_MODES: tuple[str, ...] = ("exact", "small_angle")


@dataclass(frozen=True)
class PhaseStats:
    """Moments of a phase-error series δ [rad]: n, ⟨δ⟩, ⟨δ²⟩, ⟨cos 2δ⟩, ⟨sin 2δ⟩.

    Filled from a series (:meth:`from_series`) or from the running sums a metric block accumulates at
    the full loop rate (:meth:`from_acc`, acc = [n, Σδ, Σδ², Σcos 2δ, Σsin 2δ]). Fields may be arrays
    after :meth:`shifted` by an array offset (cw_sideband θ_off(f_s)).
    """

    n: int
    mean: Any
    mean_sq: Any
    mean_cos2: Any
    mean_sin2: Any

    @classmethod
    def from_series(cls, delta_rad: Any) -> PhaseStats:
        d = np.asarray(delta_rad, dtype=float).ravel()
        if d.size == 0:
            raise ValueError("PhaseStats: empty series")
        return cls(int(d.size), float(d.mean()), float(np.mean(d * d)), float(np.mean(np.cos(2 * d))),
                   float(np.mean(np.sin(2 * d))))

    @classmethod
    def from_acc(cls, acc: Any) -> PhaseStats:
        a = np.asarray(acc, dtype=float)
        if a.shape != (5,):
            raise ValueError(f"PhaseStats.from_acc needs [n, Σδ, Σδ², Σcos2δ, Σsin2δ], got shape {a.shape}")
        n = int(round(a[0]))
        if n < 1:
            raise ValueError("PhaseStats.from_acc: no samples accumulated")
        return cls(n, a[1] / n, a[2] / n, a[3] / n, a[4] / n)

    @property
    def var(self) -> Any:
        """Var(δ) = ⟨δ²⟩ − ⟨δ⟩² (≥ 0)."""
        return np.maximum(np.asarray(self.mean_sq) - np.asarray(self.mean) ** 2, 0.0)

    @property
    def std(self) -> Any:
        return np.sqrt(self.var)

    @property
    def rms(self) -> Any:
        return np.sqrt(np.asarray(self.mean_sq))

    def shifted(self, offset_rad: Any) -> PhaseStats:
        """Moments of δ − offset (exact trigonometric identities)."""
        a = np.asarray(offset_rad, dtype=float)
        c, s = np.cos(2 * a), np.sin(2 * a)
        return PhaseStats(
            self.n,
            self.mean - a,
            self.mean_sq - 2 * a * self.mean + a * a,
            self.mean_cos2 * c + self.mean_sin2 * s,
            self.mean_sin2 * c - self.mean_cos2 * s,
        )


def _db(v: Any) -> Any:
    out = lin_to_db(np.asarray(v, dtype=float) / V_VAC)
    return float(out) if np.ndim(out) == 0 else np.asarray(out)


def _f(x: Any) -> Any:
    return float(x) if np.ndim(x) == 0 else np.asarray(x, dtype=float)


def phase_to_quadrature(
    delta_rad: np.ndarray | PhaseStats,
    v_minus: Any,
    v_plus: Any,
    *,
    hf_variance_rad2: float,
    average: str,
    setpoint_rad: float,
    theta_off_rad: Any = 0.0,
    dither_m_rad: Sequence[float] = (),
) -> dict[str, Any]:
    """Map the time-loop residual δ (series or :class:`PhaseStats`) and V± to V′, R′∓ [dB] and σ per band.

    Args:
        delta_rad: δ(t) of one metric (time-loop part, relative to its set point) or its moments.
        v_minus, v_plus: squeezed / anti-squeezed variances after loss [conv §2 units, vacuum 1/2];
            arrays over f_s for a cw_sideband source.
        hf_variance_rad2: out-of-loop-band phase variance σ²_HF (from the ledger, 0 for none).
        average: ``exact`` (χ products, T04 eq. (1)) or ``small_angle`` (prompt §7-1 form with
            σ² = ⟨δ′²⟩ + σ²_HF + Σm²/2; set point must be 0 or π/2 mod π).
        setpoint_rad: Δ_set (0 → squeezing, π/2 → anti-squeezing, conv §5).
        theta_off_rad: static source rotation θ_off (``SourceOutput.project``).
        dither_m_rad: peak amplitudes of sinusoidal dithers not present in δ (χ = Π J0(2m)).

    Returns:
        dict with ``v_prime`` / ``r_minus`` / ``r_plus`` (linear, conv units), the same in dB re shot
        noise (``v_prime_db``, ``r_minus_db``, ``r_plus_db``), ``v_minus_db`` / ``v_plus_db`` (no
        jitter), χ factors (``chi_loop``, ``chi_hf``, ``chi_dither``, ``chi``), σ per band
        (``sigma_loop_rad`` rms of δ′ incl. static offset, ``std_loop_rad``, ``mean_loop_rad`` = δ0,
        ``sigma_hf_rad``, ``sigma_dither_rad``, ``sigma_total_rad``) and, for a series input with a
        scalar θ_off, ``v_prime_t`` (V′(t) of the quasi-static Δ(t), HF/dither χ applied; T04 §5.1).
    """
    if average not in AVERAGE_MODES:
        raise ValueError(f"average must be one of {AVERAGE_MODES}, got {average!r}")
    hf = float(hf_variance_rad2)
    if not hf >= 0.0:
        raise ValueError(f"hf_variance_rad2 must be >= 0, got {hf_variance_rad2!r}")
    series = None
    if isinstance(delta_rad, PhaseStats):
        stats = delta_rad
    else:
        series = np.asarray(delta_rad, dtype=float)
        stats = PhaseStats.from_series(series)
    st = stats.shifted(theta_off_rad)
    vm, vp = np.asarray(v_minus, dtype=float), np.asarray(v_plus, dtype=float)
    mean_v, half_dv = 0.5 * (vm + vp), 0.5 * (vp - vm)
    m = np.asarray(list(dither_m_rad), dtype=float)
    chi_hf = float(np.exp(-2.0 * hf))
    chi_d = float(np.prod(gaussian.chi_sinusoid(m))) if m.size else 1.0
    sig_d2 = float(np.sum(m * m) / 2.0)
    a = float(setpoint_rad)

    def v_at(angle: float) -> Any:  # V′ for Δ = angle + δ′ (exact)
        cos2 = np.cos(2 * angle) * st.mean_cos2 - np.sin(2 * angle) * st.mean_sin2
        return mean_v - half_dv * cos2 * chi_hf * chi_d

    if average == "exact":
        r_minus, r_plus, v_prime = v_at(0.0), v_at(0.5 * np.pi), v_at(a)
    else:
        sigma2 = st.mean_sq + hf + sig_d2
        c2, s2 = np.cos(np.sqrt(sigma2)) ** 2, np.sin(np.sqrt(sigma2)) ** 2
        r_minus = vm * c2 + vp * s2
        r_plus = vp * c2 + vm * s2
        k = np.mod(a, np.pi)
        if min(k, np.pi - k) < 1e-12:
            v_prime = r_minus
        elif abs(k - 0.5 * np.pi) < 1e-12:
            v_prime = r_plus
        else:
            raise ValueError("average 'small_angle' needs setpoint_rad = 0 or π/2 (mod π); use 'exact'")
    out: dict[str, Any] = {
        "average": average,
        "setpoint_rad": a,
        "theta_off_rad": _f(theta_off_rad),
        "v_minus": _f(vm),
        "v_plus": _f(vp),
        "v_minus_db": _db(vm),
        "v_plus_db": _db(vp),
        "v_prime": _f(v_prime),
        "r_minus": _f(r_minus),
        "r_plus": _f(r_plus),
        "v_prime_db": _db(v_prime),
        "r_minus_db": _db(r_minus),
        "r_plus_db": _db(r_plus),
        "chi_loop": _f(st.mean_cos2),
        "chi_hf": chi_hf,
        "chi_dither": chi_d,
        "chi": _f(st.mean_cos2 * chi_hf * chi_d),
        "n_samples": int(st.n),
        "mean_loop_rad": _f(st.mean),
        "std_loop_rad": _f(st.std),
        "sigma_loop_rad": _f(st.rms),
        "sigma_hf_rad": float(np.sqrt(hf)),
        "sigma_dither_rad": float(np.sqrt(sig_d2)),
        "sigma_total_rad": _f(np.sqrt(st.mean_sq + hf + sig_d2)),
        "hf_variance_rad2": hf,
    }
    if series is not None and np.ndim(theta_off_rad) == 0 and vm.ndim == 0:
        d = a + series - float(theta_off_rad)
        out["v_prime_t"] = mean_v - half_dv * np.cos(2 * d) * chi_hf * chi_d
    return out


def quadrature_for_metric(
    delta_rad: np.ndarray | PhaseStats,
    *,
    metric: str,
    source: SourceOutput,
    lo_mode: Any,
    loss_registry: LossRegistry | None,
    path: str | None,
    target: str,
    band_ledger: BandLedger | None,
    beams: Mapping[str, float],
    average: str,
    setpoint_rad: float,
    include_hf_phase_noise: bool,
    dither_m_rad: Sequence[float] = (),
    theta_off_ref_rad: Any = 0.0,
) -> dict[str, Any]:
    """Noise levels of one metric: SourceOutput.project → LossRegistry η → ledger σ²_HF → V′, R′∓.

    Args:
        delta_rad: time-loop residual of the metric (series or :class:`PhaseStats`).
        metric: ledger metric name (``"i.port<k>"`` / ``"ii.port<k>"``; defined with its f_max_hz).
        source: the squeezed beam's ``SourceOutput`` (cov at φ_p = 0, D55 rule 4).
        lo_mode: LO mode for a temporal source (None otherwise).
        loss_registry, path, target: η = ``loss_registry.eta(path, target)`` with target ``"bhd"``
            (metric (i)) or ``"gbs"`` (metric (ii)); None registry/path → η = 1 (losses already in
            the source, e.g. the lumped §7-1 reference).
        band_ledger, beams: σ²_HF = ``band_ledger.hf_variance(metric, beams)`` (signed beam
            combination, e.g. {lo: +1, sqz: −1}); used only if ``include_hf_phase_noise``.
        average, setpoint_rad, dither_m_rad: as :func:`phase_to_quadrature`.
        theta_off_ref_rad: θ_off of the reference (port 0) for metric (ii): Θ_k sees θ_off,k − θ_off,0.
    """
    vm, vp, theta_off = source.project(lo_mode)
    eta = 1.0
    if loss_registry is not None and path is not None:
        eta = loss_registry.eta(path, target)
    vm = eta * np.asarray(vm) + (1.0 - eta) * V_VAC  # pure loss, conv §2
    vp = eta * np.asarray(vp) + (1.0 - eta) * V_VAC
    hf = 0.0
    if include_hf_phase_noise:
        if band_ledger is None:
            raise ValueError(f"metric {metric!r}: include_hf_phase_noise needs a BandLedger")
        hf = band_ledger.hf_variance(metric, dict(beams))
    out = phase_to_quadrature(delta_rad, vm, vp, hf_variance_rad2=hf, average=average, setpoint_rad=setpoint_rad,
                              theta_off_rad=np.asarray(theta_off) - np.asarray(theta_off_ref_rad),
                              dither_m_rad=dither_m_rad)
    out.update({"metric": metric, "eta": float(eta), "loss_target": target, "loss_path": path,
                "beams": {str(b): float(c) for b, c in beams.items()}})
    if source.basis == "cw_sideband":
        out["f_s_hz"] = np.asarray(source.grid)
    return out
