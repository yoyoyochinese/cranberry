"""M4 review (WP-B7): polarization checks.

  1  T09 test 7 (D62-6): Gaussian amplitude overlap of the crosstalk satellite at τ = 1.36 ps; the τ_p that
     gives 1e−4; arrival sign of the fast-axis satellite from U(Ω) under conv §8.
  2  PMD timing_reference: projected mode f = ½[u(T − τ/2) + u(T + τ/2)] (equal PSP split); η against a
     time-shifted reference u(T − s): psp ⇔ s = τ/2, centered ⇔ s = 0 = the best shift; ``pulse_projection``.
  3  segment_jones SMF statistics: ``diffusion`` ⟨cos Θ⟩ = e^{−D_pol t} (⟨Θ²⟩ = 2·D_pol·t for small angles); ``ou`` stationary
     per-component variance σ_Ω² = D_pol/(4π f_c) (ensemble of independent segments in one block).
  4  Faraday with Jᵀ: Faraday-mirror double pass R(ρ)·R(ρ) = R(2ρ), optical activity cancels.

Run: uv run python reviews/scripts/m4_pol_checks.py
"""
import math

import numpy as np

import cranberry.path.testing as pt
from cranberry.blocks.pol.pmd_mode_overlap import make_pmd_mode_overlap, pulse_projection
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State


def stokes(e):
    a, b = e[..., 0], e[..., 1]
    return np.stack([abs(a) ** 2 - abs(b) ** 2, 2 * (a.conj() * b).real, 2 * (a.conj() * b).imag], -1)


def smf_ensemble(n_seg, drift, D, fc, fs, n_steps, seed=5):
    ids = [f"s{i}" for i in range(n_seg)]
    segs = [pt.seg_spec(s, 1.0, kind="smf", pol_model="jones", pol_params={"J_slot": f"J.{s}"}, drivers={"dT": "uniform"})
            for s in ids]
    topo = pt.topo_spec(segs, [pt.beam_spec("b", [pt.step_spec(ids[0])], "o", jones_out=(1, 0))], ["o"], pt.coef_spec(ids[0], "dT", "*", 0.0))
    cfg = {"model": "smf", "axis_drivers": ["dT"], "axis_stokes": [1, 0, 0], "u0": "identity", "drift": drift,
           "D_pol_rad2_per_s": D, "k_T_rad_per_m_K": 0.0}
    if drift == "ou":
        cfg["f_c_hz"] = fc
    blk_spec = {"block": "segment_jones", "name": "sj", "params": {
        "topology": topo, "segments": {s: dict(cfg) for s in ids},
        "toggles": {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}}}
    lay = SlotLayout(fs)
    for s in ids:
        lay.add(f"drv.{s}.dT", (), np.float64, "K")
    blk = build(spec_from_obj([blk_spec]), lay, seed, external=[f"drv.{s}.dT" for s in ids])
    st, pr = State(lay), Params(lay)
    run_reference(blk, st, pr, n_steps)
    U = np.array([np.array(st[f"J.{s}"]) for s in ids])
    return U


if __name__ == "__main__":
    print("1. T09 test 7")
    tau = 1.36e-12
    sig = lambda tp: tp / (2 * math.sqrt(2 * math.log(2)))  # noqa: E731
    O = math.exp(-tau**2 / (8 * sig(1e-12) ** 2))
    tp4 = tau / math.sqrt(8 * math.log(1e4)) * 2 * math.sqrt(2 * math.log(2))
    print(f"  O(1.36 ps) at τ_p = 1 ps: {O:.4f}; O = 1e−4 needs τ_p ≤ {tp4 * 1e12:.3f} ps")
    assert abs(O - 0.277) < 1e-3 and abs(tp4 - 0.373e-12) < 1e-15
    n, dt = 1 << 13, 5e-15
    t = (np.arange(n) - n // 2) * dt
    u = np.exp(-t**2 / (4 * sig(100e-15) ** 2)).astype(complex)
    u /= math.sqrt(np.sum(abs(u) ** 2) * dt)
    U_fast = lambda om: np.array([[[0, 0], [0, np.exp(-1j * tau * w)]] for w in om])  # noqa: E731
    f = pulse_projection(u, dt, U_fast, (0, 1), (0, 1))
    tc = np.sum(t * abs(f) ** 2) / np.sum(abs(f) ** 2)
    print(f"  fast-axis component with U = e^(−iτΩ) arrives at {tc * 1e12:+.3f} ps (early)")
    assert abs(tc + tau) < dt

    print("2. PMD reference timing (equal PSP split)")
    for tp, tdgd in ((100e-15, 25e-15), (100e-15, 79e-15), (1e-12, 79e-15)):
        s0 = sig(tp)
        u = np.exp(-t**2 / (4 * s0**2)).astype(complex)
        u /= math.sqrt(np.sum(abs(u) ** 2) * dt)

        def U_pmd(om, tdgd=tdgd):  # PSPs = ±45° linear, delays ±τ/2; input x, analyzer x
            P1 = np.array([[1, 1], [1, 1]]) / 2
            P2 = np.array([[1, -1], [-1, 1]]) / 2
            return np.array([P1 * np.exp(1j * w * tdgd / 2) + P2 * np.exp(-1j * w * tdgd / 2) for w in om])
        fm = pulse_projection(u, dt, U_pmd, (1, 0), (1, 0))
        shifts = np.linspace(-tdgd, tdgd, 401)
        nrm = math.sqrt(np.sum(abs(u) ** 2) * dt) / u.real.max()
        eta = [abs(np.sum(np.exp(-(t - s_) ** 2 / (4 * s0**2)) / nrm * fm) * dt) ** 2 for s_ in shifts]
        best = shifts[int(np.argmax(eta))]
        base = {"D_PMD_s_per_sqrt_m": tdgd, "L_m": 1.0, "tau_p_fwhm_s": tp, "segment": "smf", "applies_to": ["gbs"],
                "source": "[review]", "toggles": {"pol_pmd": True}}
        e_psp = make_pmd_mode_overlap({**base, "timing_reference": "psp"}).eta
        e_cen = make_pmd_mode_overlap({**base, "timing_reference": "centered"}).eta
        e_half = float(np.interp(tdgd / 2, shifts, eta))
        print(f"  τ_p {tp * 1e15:5.0f} fs, DGD {tdgd * 1e15:3.0f} fs: best shift {best * 1e15:+.2f} fs, "
              f"η(best) {max(eta):.6f} = centered {e_cen:.6f}; η(s = τ/2) {e_half:.6f} = psp {e_psp:.6f}")
        assert abs(best) <= 2 * (shifts[1] - shifts[0]) and abs(max(eta) - e_cen) < 1e-6 and abs(e_half - e_psp) < 1e-5

    print("3. segment_jones SMF drift statistics")
    D, fs, nst = 2.0, 1e3, 200
    U = smf_ensemble(400, "diffusion", D, None, fs, nst)
    s = stokes(U @ np.array([1, 0]))  # start state x ↔ s = (1, 0, 0)
    c = float(np.mean(s[:, 0]))  # rotational diffusion, per-axis variance D·t: ⟨cos Θ⟩ = e^{−D t} (⟨Θ²⟩ → 2Dt)
    print(f"  diffusion: ⟨cos Θ⟩ = {c:.4f} vs e^(−D t) = {math.exp(-D * nst / fs):.4f}")
    assert abs(c - math.exp(-D * nst / fs)) < 0.04
    fc = 5.0
    U = smf_ensemble(400, "ou", D, fc, fs, 3000)  # 3000 steps = 15 τ
    om = []
    for Ui in U:  # W = exp(−iΩ·σ/2) ⇒ Ω from the SU(2) log
        a, b = Ui[0, 0], Ui[1, 0]
        th = 2 * math.acos(max(-1.0, min(1.0, a.real)))
        k = th / (2 * math.sin(th / 2)) if th > 1e-12 else 1.0
        om.append([-2 * a.imag * k, -2 * b.imag * k, 2 * b.real * k])
    var = np.var(np.array(om), axis=0).mean()
    print(f"  ou: stationary Var(Ω_i) = {var:.4f} vs D/(4π f_c) = {D / (4 * math.pi * fc):.4f}")
    assert abs(var / (D / (4 * math.pi * fc)) - 1) < 0.12

    print("4. Faraday vs reciprocal rotation, double pass with Jᵀ backward")
    R = lambda r: np.array([[math.cos(r), -math.sin(r)], [math.sin(r), math.cos(r)]])  # noqa: E731
    rho = 0.3
    far = R(rho) @ R(rho)  # backward N (not transposed) · forward N
    act = R(rho).T @ R(rho)
    assert np.allclose(far, R(2 * rho)) and np.allclose(act, np.eye(2))
    print("  Faraday mirror R(2ρ), optical activity cancels: OK")
    print("all checks passed")
