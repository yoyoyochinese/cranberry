"""Fiber quantum-noise layer (WP-B10; T07 §2, G03): ``fiber_quantum_linearized`` (reference) and
``fiber_quantum_tw`` (truncated Wigner, validation only, off in presets).

Linearized method (T07 §2.2). Â = A + √(hν0) â on the B9 time grid (N points, Δt), â_n = √Δt â(T_n),
quadrature vector x = (q_1..q_N, p_1..p_N) (xxpp, conv §2/§9) on the *centred* grid t_s. Per segment
[z_k, z_k+1] (``steps.n_save`` segments of ``n_steps/n_save`` RK4IP steps):

    V ← M_k V M_kᵀ + N_k,     M_k = exact tangent of the B9 discrete map (2N real directions δA = e_n, i e_n)

* loss noise (toggle ``loss_noise``): N = (1 − η)/2 · (𝟙 + M Mᵀ)/(1 + η), η = e^{−αh}: trapezoid of
  ∫ M(z+h,s)(α/2)M(z+h,s)ᵀ ds normalised so that it is exact for a passive linear segment
  (M Mᵀ = η𝟙 ⇒ (1 − η)/2 𝟙, §7-3/§7-4 exact) and a convex mix of physical channels otherwise.
* phase noise of the bright field (``raman_noise``, ``gawbs_noise``): ∂_z â ⊃ i A φ/√(hν0) with
  ⟨φ(z,T)φ(z′,T′)⟩ = δ(z−z′) C(T−T′) ⇒ D(z) = 2(Δt/hν0)·(u uᵀ ∘ [[C, C], [C, C]]), u = (−Im A, Re A);
  N = (h/2)(M D(z_k) Mᵀ + D(z_k+1)). Raman: C = ifft_agrawal(S_R(Ω)) on the grid (S_R = hν0 g_R (n_th + ½),
  periodic like the GNLSE); GAWBS: analytic C(τ) on unwrapped lags (its MHz lines are below the grid resolution).
  The Raman gain/absorption (Im h̃_R) is inside M through the B9 tangent (D&C 2001).
* ``kerr_fluct`` off: the tangent is taken with Kerr and Raman removed (fluctuations see dispersion,
  loss and external XPM only); the classical field and the noise D still use the full solution.
* counter-propagating probe backscatter (``brillouin_backscatter``, ``rayleigh_backscatter``,
  ``raman_backscatter``; ``scattering.backscatter_sources``): stationary phase-insensitive noise added at
  the output (scattered light reaches the output end, its spectrum is unchanged by the passive
  remaining fiber, loss is in L_bs). ⟨â_n† â_m⟩ = Δt G(T_n − T_m), G(τ) = ∫ n(f) e^{+i2πfτ} df (lines:
  Φ e^{i2πf0τ} e^{−πδν|τ|}); quadrature block [[Re N, Im N], [−Im N, Re N]] (projected: ∫n|ũ|² df, conv §2).

Output: :attr:`FiberQuantumLinearized.source_output` — a ``temporal`` ``SourceOutput`` whose modes are
the grid bins (f_n = 1/√Δt on bin n) on the centred grid, so ``project(u)`` takes any LO mode u(T);
``result.noise_cov`` holds each noise source's propagated contribution (decomposition, "added N").

CW (``laser_mode = cw`` with ``cw_shortcut`` on; T07 §2.2 "CW"): per sideband f_s a 2×2 map in the
carrier frame, K = [[−α/2, −β_e(Ω)], [β_e(Ω) + 2γ_eff(Ω) P e^{−αz}, −α/2]] (β_e = even dispersion,
γ_eff = γ(c_K + c_R Re h̃_R(Ω))), sub-steps M = expm(K h) at the mid-point power, the same loss-noise
rule, phase noise on p: D = diag(0, 2 P S(Ω)/hν0); output rotated by the carrier phase φ_NL into the
lab frame; ``cw_sideband`` ``SourceOutput``. Backscatter: V += (n(f_s) + n(−f_s))/2 · 𝟙.
ponytail: Im h̃_R(Ω) (≈ Ω·8 fs, 5e−5 at 1 GHz) is dropped in the CW sideband map (keeps the real 2×2 form).

Truncated Wigner (T07 §2.3): Ψ(0) = A + √(hν0/Δt)(q + ip)/√2 with x ~ N(0, V0), full B9 GNLSE per
trajectory (batched), loss noise as additive increments Var(q) = Var(p) = αh/2 per step (Euler);
V = sample covariance of x(L) (symmetric order = conv V directly).
ponytail: TW carries loss noise only (Raman/GAWBS/backscatter are validated in the linear layer).

Losses (G06): the fiber attenuation acts inside the map (toggle ``loss``, distributed and interleaved
with the Kerr term), so ``losses()`` is empty. A preset must not register the same fiber segment in the
``LossRegistry`` as well (or it switches ``loss`` off here and lets the registry apply η once).

Memory/cost: the full 2N×2N map (2N tangents through RK4IP) — fine for N ≲ 1024.
ponytail: an adjoint (Mᵀ w) step would give one LO mode in O(N) memory; add when N ≫ 1024 is needed.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from typing import Any

import numpy as np
from scipy.linalg import expm

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.gnlse import TOGGLES as GNLSE_TOGGLES
from cranberry.blocks.fib.gnlse import (
    _alpha,
    _check_fiber,
    _weights,
    _xpm_rate,
    cw_sideband_rotation,
    l_eff,
    model_from_params,
    pulse_field,
    raman_response,
    time_axis,
)
from cranberry.core import rng as crng
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.fft import fft_agrawal, freq_grid, ifft_agrawal
from cranberry.core.registry import register_block
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput

__all__ = ["QUANTUM_TOGGLES", "NOISE_SOURCES", "LinearResult", "FiberQuantumLinearized", "FiberQuantumTW",
           "fiber_quantum_linear", "fiber_quantum_tw", "lo_matrix", "project_cov", "stationary_noise_cov",
           "make_fiber_quantum_linearized", "make_fiber_quantum_tw"]

QUANTUM_TOGGLES = ("loss_noise", "kerr_fluct", "raman_noise", "gawbs_noise", "brillouin_backscatter",
                   "rayleigh_backscatter", "raman_backscatter")
NOISE_SOURCES = ("loss", "raman", "gawbs", "brillouin", "rayleigh", "raman_bs")
TW_TOGGLES = (*GNLSE_TOGGLES, "loss_noise")


# --------------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------------

def lo_matrix(u: np.ndarray, dt: float) -> np.ndarray:
    """2×2N matrix a with (q_u, p_u) = a x for the mode u(T) on the grid (â_u = Σ c_n â_n, c_n = u_n*√Δt,
    normalised; same as ``SourceOutput.lo_mode_cov``). aᵀ embeds a mode into the grid."""
    u = np.asarray(u, dtype=complex)
    c = u.conj() * np.sqrt(dt)
    c = c / np.linalg.norm(c)
    return np.block([[c.real[None, :], -c.imag[None, :]], [c.imag[None, :], c.real[None, :]]])


def project_cov(cov: np.ndarray, u: np.ndarray, dt: float) -> np.ndarray:
    """a cov aᵀ (2×2) for a mode u — use on a noise contribution (no vacuum fill)."""
    a = lo_matrix(u, dt)
    return a @ cov @ a.T


def _quad_block(n_c: np.ndarray) -> np.ndarray:
    """Quadrature covariance of a phase-insensitive noise with ⟨â_n† â_m⟩ = N_nm: [[Re, Im], [−Im, Re]]."""
    return np.block([[n_c.real, n_c.imag], [-n_c.imag, n_c.real]])


def _lags(n: int, dt: float) -> np.ndarray:
    """Matrix of T_n − T_m (not wrapped: a Toeplitz matrix of a positive-definite function is PSD,
    Bochner; wrapped lags are not for lines off the grid's frequency bins)."""
    return np.subtract.outer(np.arange(n), np.arange(n)) * dt


def stationary_noise_cov(n: int, dt: float, lines, density=None) -> np.ndarray:
    """2N×2N added covariance of stationary noise: ``lines`` [(name, f0, flux, fwhm)] and a density n(f)
    [1/(s Hz)] sampled on the grid frequencies (module docstring)."""
    tau = _lags(n, dt)
    g = np.zeros((n, n), dtype=complex)
    for _, f0, flux, fwhm in lines:
        g += flux * np.exp(2j * np.pi * f0 * tau - np.pi * fwhm * np.abs(tau))
    if density is not None:
        f = freq_grid(n, dt)
        col = fft_agrawal(density(f), dt) / (n * dt**2)  # Σ_k n_k e^{+2πikl/N}/(NΔt): lag l (DFT wrapper)
        g += col[np.subtract.outer(np.arange(n), np.arange(n)) % n]
    return _quad_block(0.5 * dt * (g + g.conj().T))  # Hermitian up to rounding


def _check_common(params: Mapping[str, Any], name: str) -> None:
    _check_fiber(params, name)
    require(params, "laser_mode", "grid", "steps", "input_state", "phi_opa0_rad", block=name)
    if params["laser_mode"] not in ("cw", "pulsed"):
        raise ValueError(f"block {name!r}: laser_mode must be 'cw' or 'pulsed', got {params['laser_mode']!r}")
    require(params, "power_w" if params["laser_mode"] == "cw" else "pulse", block=name)
    if params["laser_mode"] == "pulsed":
        require(params["pulse"], "shape", "peak_power_w", "t0_s", block=f"{name}.pulse")
    require(params["grid"], "n_t", "t_window_s", block=f"{name}.grid")
    require(params["steps"], "n_steps", "n_save", "control", block=f"{name}.steps")
    if params["steps"]["control"] != "fixed" or int(params["steps"]["n_steps"]) % int(params["steps"]["n_save"]):
        raise ValueError(f"block {name!r}: steps.control must be 'fixed' with n_steps a multiple of n_save "
                         "(segments of the covariance propagation)")
    st = params["input_state"]
    require(st, "kind", block=f"{name}.input_state")
    if st["kind"] == "squeezed":
        require(st, "r", "theta_sq_rad", "mode", block=f"{name}.input_state")
    elif st["kind"] != "vacuum":
        raise ValueError(f"block {name!r}: input_state.kind must be 'vacuum' or 'squeezed', got {st['kind']!r}")


def _grid(params):
    g = params["grid"]
    return time_axis(int(g["n_t"]), float(g["t_window_s"]))


def _classical_input(params, t):
    if params["laser_mode"] == "cw":
        return np.full(t.shape, np.sqrt(float(params["power_w"])), dtype=complex)
    q = params["pulse"]
    return pulse_field(q["shape"], float(q["peak_power_w"]), float(q["t0_s"]), t).astype(complex)


def _input_cov(params, t, dt, v0):
    """Initial grid covariance: ``v0`` (array / SourceOutput) or ``params['input_state']``."""
    n = t.size
    vac = 0.5 * np.eye(2 * n)
    if v0 is None:
        st = params["input_state"]
        if st["kind"] == "vacuum":
            return vac
        m = st["mode"]
        u = pulse_field(m["shape"], 1.0, float(m["t0_s"]), t)
        a = lo_matrix(u, dt)
        return vac + a.T @ (gaussian.squeezed_vacuum(float(st["r"]), float(st["theta_sq_rad"])) - 0.5 * np.eye(2)) @ a
    if isinstance(v0, SourceOutput):
        if v0.basis != "temporal" or v0.grid.shape != t.shape or not np.allclose(v0.grid, t, rtol=0, atol=1e-9 * dt):
            raise ValueError("SourceOutput input must be 'temporal' on the fiber's centred grid (no resampling)")
        cmat = v0.modes.conj() * np.sqrt(dt)  # rows: c_jn = ∫ f_j* bin_n
        a = np.block([[cmat.real, -cmat.imag], [cmat.imag, cmat.real]])
        return vac + a.T @ (v0.cov - 0.5 * np.eye(v0.cov.shape[0])) @ a
    v0 = np.asarray(v0, dtype=float)
    if v0.shape != vac.shape:
        raise ValueError(f"v0 must be ({2 * n}, {2 * n}), got {v0.shape}")
    return v0


# --------------------------------------------------------------------------------------------------
# linearized layer
# --------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class LinearResult:
    """Output of :func:`fiber_quantum_linear`.

    ``cov``: output covariance (temporal: 2N×2N on the centred grid ``t_s``; cw: (n_f, 2, 2)).
    ``noise_cov``: {source: propagated added covariance} for :data:`NOISE_SOURCES` (zeros when off).
    ``transfer``: total M (temporal 2N×2N, cw (n_f, 2, 2), lab frame). ``a_out``: classical output
    (centred grid) or CW amplitude. ``t_s`` / ``f_s_hz``: grid; ``dt``: time step (temporal).
    """

    cov: np.ndarray
    noise_cov: dict[str, np.ndarray]
    transfer: np.ndarray
    a_out: np.ndarray | complex
    t_s: np.ndarray | None
    f_s_hz: np.ndarray | None
    dt: float | None

    def lo_variances(self, u: np.ndarray) -> np.ndarray:
        """2×2 covariance of the LO mode u (temporal)."""
        return project_cov(self.cov, u, self.dt)


def _noise_flags(tg):
    return {"rayleigh": tg["rayleigh_backscatter"], "brillouin": tg["brillouin_backscatter"],
            "raman": tg["raman_backscatter"], "sbs_amplification": False}


def _split_lines(lines, density, tg):
    """{source: (lines, density)} for the three backscatter sources."""
    return {"rayleigh": ([ln for ln in lines if ln[0] == "rayleigh"], None),
            "brillouin": ([ln for ln in lines if ln[0] == "brillouin"], None),
            "raman_bs": ([], density if tg["raman_backscatter"] else None)}


def _loss_noise(m: np.ndarray, eta: float) -> np.ndarray:
    eye = np.eye(m.shape[-1])
    return 0.5 * (1.0 - eta) * (eye + m @ np.swapaxes(m, -1, -2)) / (1.0 + eta)


def fiber_quantum_linear(params: Mapping[str, Any], toggles: Mapping[str, bool], v0=None) -> LinearResult:
    """Linearized propagation (module docstring). ``params``: the ``fiber_quantum_linearized`` params;
    ``toggles``: B9 toggles + :data:`QUANTUM_TOGGLES`; ``v0``: optional input (2N×2N array or a
    ``temporal`` ``SourceOutput`` on the same centred grid; cw: (n_f, 2, 2)) overriding ``input_state``."""
    if params["laser_mode"] == "cw" and toggles["cw_shortcut"]:
        return _linear_cw(params, toggles, v0)
    t, dt = _grid(params)
    n, shift = t.size, t.size // 2
    tg_c = {k: toggles[k] for k in GNLSE_TOGGLES}
    model = model_from_params(params, tg_c, n, dt)
    lin = model if toggles["kerr_fluct"] else model_from_params(params, {**tg_c, "kerr": False, "raman": False}, n, dt)
    hnu = sc._hnu(params["lambda0_m"])
    temp = float(params["temperature_k"])

    # photon-number quadratures: â(Ω) ∝ Ã(Ω)/√(1 + Ω/ω0) when self-steepening is on (its conserved
    # quantity is ∫|Ã|²/(ω0 + Ω); with a constant hν0 the map would not be symplectic, ~Ω/ω0)
    wq = np.sqrt(1.0 / model.shock) if toggles["self_steepening"] else np.ones(n)

    def field_spec(x):
        return model.to_w(np.roll(x[..., :n] + 1j * x[..., n:], -shift, axis=-1))

    def to_spec(x):
        return field_spec(x) / wq

    def to_x(at):
        a = np.roll(model.to_t(at * wq), shift, axis=-1)
        return np.concatenate([a.real, a.imag], axis=-1)

    tq = to_x(field_spec(np.eye(2 * n))).T if toggles["self_steepening"] else None  # field → photon quadratures

    # phase-noise correlation per unit length on the grid
    c = np.zeros((n, n))
    if toggles["raman_noise"] and model.c_raman != 0.0:
        s = sc.raman_phase_psd(model.omega, gamma_per_w_m=model.gamma, f_R=float(params["f_R"]),
                               raman_model=params["raman_model"], temperature_k=temp, lambda0_m=params["lambda0_m"])
        col = ifft_agrawal(s, dt).real
        c_r = col[np.subtract.outer(np.arange(n), np.arange(n)) % n]
    else:
        c_r = c
    c_g = (sc.gawbs_correlation(_lags(n, dt), params["gawbs"], temperature_k=temp, lambda0_m=params["lambda0_m"])
           if toggles["gawbs_noise"] else c)

    def diff(at, cc):
        a = np.roll(model.to_t(at), shift)
        u = np.concatenate([-a.imag, a.real])
        d = 2 * dt / hnu * np.outer(u, u) * np.block([[cc, cc], [cc, cc]])  # field-basis noise
        return d if tq is None else tq @ d @ tq.T

    v = _input_cov(params, t, dt, v0)
    noise = {k: np.zeros_like(v) for k in NOISE_SOURCES}
    m_tot = np.eye(2 * n)
    L, s_ = float(params["L_m"]), params["steps"]
    n_seg, per = int(s_["n_save"]), int(s_["n_steps"]) // int(s_["n_save"])
    h = L / int(s_["n_steps"])
    at = model.to_w(np.roll(_classical_input(params, t), -shift))
    basis = to_spec(np.eye(2 * n))
    z = 0.0
    for _ in range(n_seg):
        at0, dat = at, basis
        for _ in range(per):
            if toggles["kerr_fluct"]:
                at, dat = model.step(at, h, z, dat)
            else:
                _, dat = lin.step(at, h, z, dat)
                at, _ = model.step(at, h, z, None)
            z += h
        m = to_x(dat).T
        m_tot = m @ m_tot
        add = {}
        if toggles["loss_noise"] and model.alpha > 0:
            add["loss"] = _loss_noise(m, np.exp(-model.alpha * h * per))
        for key, cc, on in (("raman", c_r, toggles["raman_noise"]), ("gawbs", c_g, toggles["gawbs_noise"])):
            if on and np.any(cc):
                add[key] = 0.5 * h * per * (m @ diff(at0, cc) @ m.T + diff(at, cc))
        v = m @ v @ m.T + sum(add.values(), np.zeros_like(v))
        for k in NOISE_SOURCES:
            noise[k] = m @ noise[k] @ m.T + add.get(k, 0.0)
    lines, dens = sc.backscatter_sources(params, _noise_flags(toggles), _p_avg(params))
    for k, (ln, de) in _split_lines(lines, dens, toggles).items():
        if ln or de is not None:
            add = stationary_noise_cov(n, dt, ln, de)
            noise[k], v = noise[k] + add, v + add
    v = 0.5 * (v + v.T)
    return LinearResult(v, noise, m_tot, np.roll(model.to_t(at), shift), t, None, dt)


def _p_avg(params):
    return float(params["probe"]["p_peak_w"]) * float(params["probe"]["duty"])


def _linear_cw(params, toggles, v0) -> LinearResult:
    """CW sideband path (module docstring)."""
    require(params, "f_s_hz", block="fiber_quantum_linearized")
    f = np.asarray(params["f_s_hz"], dtype=float)
    w = 2 * np.pi * f
    tg = toggles
    alpha = float(params["alpha_per_m"]) if tg["loss"] else 0.0
    betas = [float(params["beta2_s2_per_m"]) * tg["beta2"], float(params["beta3_s3_per_m"]) * tg["beta3"],
             *[float(b) * tg["beta_higher"] for b in params["beta_higher"]]]
    c_k, c_r = _weights(params, tg)
    gam = float(params["gamma_per_w_m"]) * (c_k + c_r * (raman_response(params["raman_model"], w).real if c_r else 0))
    g_fl = gam if tg["kerr_fluct"] else 0.0 * gam
    p0, L = float(params["power_w"]), float(params["L_m"])
    n_steps = int(params["steps"]["n_steps"])
    h = L / n_steps
    be = cw_sideband_rotation(w, 1.0, betas)  # even dispersion per metre
    hnu, temp = sc._hnu(params["lambda0_m"]), float(params["temperature_k"])
    s_ph = np.zeros_like(w)
    if tg["raman_noise"] and c_r:
        s_ph = s_ph + sc.raman_phase_psd(w, gamma_per_w_m=float(params["gamma_per_w_m"]), f_R=float(params["f_R"]),
                                         raman_model=params["raman_model"], temperature_k=temp,
                                         lambda0_m=params["lambda0_m"])
    s_g = np.zeros_like(w)
    if tg["gawbs_noise"]:
        s_g = sc.gawbs_field_psd(f, params["gawbs"], temperature_k=temp, lambda0_m=params["lambda0_m"])
    pw = lambda z: p0 * np.exp(-alpha * z)  # noqa: E731

    def dmat(z, s):
        d = np.zeros((f.size, 2, 2))
        d[:, 1, 1] = 2 * pw(z) * s / hnu
        return d

    st = params["input_state"]
    if v0 is not None:
        v = np.array(v0, dtype=float)
    elif st["kind"] == "vacuum":
        v = np.broadcast_to(0.5 * np.eye(2), (f.size, 2, 2)).copy()
    else:
        v = np.broadcast_to(gaussian.squeezed_vacuum(float(st["r"]), float(st["theta_sq_rad"])), (f.size, 2, 2)).copy()
    noise = {k: np.zeros_like(v) for k in NOISE_SOURCES}
    m_tot = np.broadcast_to(np.eye(2), v.shape).copy()
    eta = np.exp(-alpha * h)
    for i in range(n_steps):
        z = i * h
        k = np.zeros((f.size, 2, 2))
        k[:, 0, 0] = k[:, 1, 1] = -alpha / 2
        k[:, 0, 1], k[:, 1, 0] = -be, be + 2 * g_fl * pw(z + h / 2)
        m = expm(k * h)
        mt = np.swapaxes(m, -1, -2)
        add = {}
        if tg["loss_noise"] and alpha > 0:
            add["loss"] = _loss_noise(m, eta)
        for key, s in (("raman", s_ph), ("gawbs", s_g)):
            if np.any(s):
                add[key] = 0.5 * h * (m @ dmat(z, s) @ mt + dmat(z + h, s))
        v = m @ v @ mt + sum(add.values(), np.zeros_like(v))
        for key in NOISE_SOURCES:
            noise[key] = m @ noise[key] @ mt + add.get(key, 0.0)
        m_tot = m @ m_tot
    lines, dens = sc.backscatter_sources(params, _noise_flags(tg), _p_avg(params))
    for key, (ln, de) in _split_lines(lines, dens, tg).items():
        nf = np.zeros_like(f)
        for _, f0, flux, fw in ln:
            lor = lambda x, f0=f0, fw=fw: (fw / (2 * np.pi)) / ((x - f0) ** 2 + (fw / 2) ** 2)  # noqa: E731
            nf = nf + flux * 0.5 * (lor(f) + lor(-f))
        if de is not None:
            nf = nf + 0.5 * (de(f) + de(-f))
        add = nf[:, None, None] * np.eye(2)
        noise[key] = noise[key] + add
        v = v + add
    phi = _cw_phi_nl(params, tg)
    rot = gaussian.rotation(phi)
    rot_t = rot.T
    v = rot @ v @ rot_t
    noise = {k_: rot @ x @ rot_t for k_, x in noise.items()}
    a_out = np.sqrt(p0) * np.exp(-alpha * L / 2) * np.exp(1j * phi)
    return LinearResult(0.5 * (v + np.swapaxes(v, -1, -2)), noise, rot @ m_tot, a_out, None, f, None)


def _cw_phi_nl(params, tg):
    """Carrier phase γ_eff P L_eff (+ external XPM), as in B9's ``cw_shortcut``."""
    c_k, c_r = _weights(params, tg)
    le = l_eff(_alpha(params, tg), float(params["L_m"]))
    return (float(params["gamma_per_w_m"]) * (c_k + c_r) * float(params["power_w"]) + _xpm_rate(params, tg)) * le


# --------------------------------------------------------------------------------------------------
# blocks
# --------------------------------------------------------------------------------------------------

class FiberQuantumLinearized(BaseBlock):
    """``fiber_quantum_linearized`` offline block (module docstring). Products: :attr:`result`
    (:class:`LinearResult`) and :attr:`source_output`.

    Params (all required, D34; every noise source's parameters are required even when its toggle is
    off, so ablation can switch any toggle):

    * ``toggles``: B9's (loss, beta2, beta3, beta_higher, kerr, raman, self_steepening, xpm_external,
      cw_shortcut) for the classical solution and the fluctuation map, plus :data:`QUANTUM_TOGGLES`
    * fiber: the ``fiber_gnlse`` keys (``L_m``, ``alpha_per_m``, ``beta2_s2_per_m``, ``beta3_s3_per_m``,
      ``beta_higher``, ``gamma_per_w_m``, ``lambda0_m``, ``raman_model``, ``f_R``, ``xpm_b``,
      ``xpm_power_w``), ``temperature_k``
    * bright field: ``laser_mode`` cw → ``power_w`` (+ ``f_s_hz`` sideband grid with ``cw_shortcut``);
      pulsed → ``pulse`` {shape, peak_power_w (0 = squeezed vacuum alone), t0_s}
    * ``grid`` {n_t, t_window_s}, ``steps`` {n_steps, n_save (= covariance segments), control: fixed}
    * ``input_state``: {kind: vacuum} | {kind: squeezed, r, theta_sq_rad, mode: {shape: sech|gaussian, t0_s}}
      (cw: the mode is ignored, flat over f_s); ``phi_opa0_rad`` (frame of the output, D55-4)
    * ``gawbs`` (see ``scattering``), ``probe`` {p_peak_w, duty, freq_offset_hz, linewidth_hz,
      polarization, eps_xt, t_perp}, ``brillouin`` {g_b_m_per_w, a_eff_m2, nu_b_hz, dnu_b_hz},
      ``rayleigh`` {b_r_per_m}
    """

    toggle_names = (*GNLSE_TOGGLES, *QUANTUM_TOGGLES)
    kind_name = "fiber_quantum_linearized"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.kind_name
        _check_common(params, name)
        require(params, "temperature_k", block=name)
        sc.check_gawbs_group(params, name)
        sc.check_probe_groups(params, name)
        require(params["probe"], "p_peak_w", "duty", block=f"{name}.probe")
        super().__init__(self.kind_name, name, params)
        if params["laser_mode"] == "cw" and self._toggles["cw_shortcut"]:
            require(params, "f_s_hz", block=name)

    def __call__(self, state, params):
        return state  # offline block

    @cached_property
    def result(self) -> LinearResult:
        return fiber_quantum_linear(self.params, self._toggles)

    @cached_property
    def source_output(self) -> SourceOutput:
        r = self.result
        prov = {"model": self.kind_name, "block": self.name, "sources": ["T07 §2.2"]}
        if r.t_s is None:
            return SourceOutput("cw_sideband", r.f_s_hz, r.cov, self.params["phi_opa0_rad"], provenance=prov)
        modes = np.eye(r.t_s.size) / np.sqrt(r.dt)
        return SourceOutput("temporal", r.t_s, r.cov, self.params["phi_opa0_rad"], modes=modes, provenance=prov)


def fiber_quantum_tw(params: Mapping[str, Any], toggles: Mapping[str, bool], n_traj: int,
                     rng: np.random.Generator, v0=None) -> np.ndarray:
    """Truncated-Wigner covariance (2N×2N, centred grid) of the fiber output (module docstring)."""
    t, dt = _grid(params)
    n, shift = t.size, t.size // 2
    model = model_from_params(params, {k: toggles[k] for k in GNLSE_TOGGLES}, n, dt)
    scale = np.sqrt(sc._hnu(params["lambda0_m"]) / dt)
    v = _input_cov(params, t, dt, v0)
    x = rng.multivariate_normal(np.zeros(2 * n), v, size=int(n_traj), method="cholesky")
    psi = _classical_input(params, t)[None, :] + scale * (x[:, :n] + 1j * x[:, n:]) / np.sqrt(2)
    at = model.to_w(np.roll(psi, -shift, axis=-1))
    s_ = params["steps"]
    h = float(params["L_m"]) / int(s_["n_steps"])
    sig = np.sqrt(model.alpha * h / 2) if toggles["loss_noise"] else 0.0
    for i in range(int(s_["n_steps"])):
        at, _ = model.step(at, h, i * h, None)
        if sig:
            xi = rng.normal(0.0, sig, size=(int(n_traj), 2, n))
            at = at + model.to_w(scale * (xi[:, 0] + 1j * xi[:, 1]) / np.sqrt(2))
    psi = np.roll(model.to_t(at), shift, axis=-1)
    d = (psi - psi.mean(axis=0)) / scale * np.sqrt(2)
    xs = np.concatenate([d.real, d.imag], axis=-1)
    return np.cov(xs, rowvar=False)


class FiberQuantumTW(BaseBlock):
    """``fiber_quantum_tw`` (optional validation block, off in presets; T07 §2.3, V9). Params: the
    ``fiber_quantum_linearized`` fiber / bright-field / grid / steps / input keys (pulsed or CW on the
    grid), ``n_traj``, ``seed``; ``toggles``: B9's + ``loss_noise``. Product: :attr:`cov`."""

    toggle_names = TW_TOGGLES
    kind_name = "fiber_quantum_tw"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.kind_name
        _check_common(params, name)
        require(params, "n_traj", "seed", block=name)
        super().__init__(self.kind_name, name, params)

    def __call__(self, state, params):
        return state

    @cached_property
    def cov(self) -> np.ndarray:
        return fiber_quantum_tw(self.params, self._toggles, int(self.params["n_traj"]),
                                crng.stream(int(self.params["seed"]), self.name))


@register_block("fiber_quantum_linearized")
def make_fiber_quantum_linearized(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``fiber_quantum_linearized`` (:class:`FiberQuantumLinearized`)."""
    return FiberQuantumLinearized(params, name)


@register_block("fiber_quantum_tw")
def make_fiber_quantum_tw(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``fiber_quantum_tw`` (:class:`FiberQuantumTW`)."""
    return FiberQuantumTW(params, name)
