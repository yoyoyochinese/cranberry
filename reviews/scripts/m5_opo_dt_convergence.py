"""M5 review: δt convergence of pulsed_opo_source (T17 V-P1, D60-5 errata, PLAN §6.4).

r_0(x_ss), K, R_− (ρ < 1) and the best Gaussian / rect / one-sided-exp LO overlaps on the Borealis and
Liu fixtures for δt = 1/(20, 40, 80, 160 κ_c), with Richardson extrapolation (error ∝ δt²).
Also the physicality of C at the fixture δt: smallest symplectic eigenvalue of C (shot-noise units,
physical ⇔ ≥ 1) and max |λ_i λ_{2N+1−i} − 1| at ρ = 1.

Run: uv run python reviews/scripts/m5_opo_dt_convergence.py
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cranberry.blocks.source.pulsed_testing import BOREALIS, LIU, pulsed_opo


def best_overlaps(b):
    u, g = b.optimal_mode[0], b.grid
    dt = g[1] - g[0]
    kc = b.kappa_c
    sub = (np.arange(64) + 0.5) / 64 - 0.5

    def eta(f):
        v = f(g[:, None] + sub[None, :] * dt)
        w, nrm = v.mean(1), np.sum(v**2) / 64 * dt
        return abs(np.sum(np.conj(w) * u) * dt) ** 2 / nrm

    c = float(np.sum(g * np.abs(u) ** 2) * dt)
    s = float(np.sqrt(np.sum((g - c) ** 2 * np.abs(u) ** 2) * dt))
    opts = {"xatol": 1e-9 * s, "fatol": 1e-14, "maxiter": 8000}
    rg = minimize(lambda x: -eta(lambda t: np.exp(-0.5 * ((t - x[0]) / x[1]) ** 2)), (c, s), method="Nelder-Mead", options=opts)
    rr = minimize(lambda x: -eta(lambda t: ((t >= x[0]) & (t < x[0] + abs(x[1]))).astype(float)), (c - 1.5 * s, 3 * s),
                  method="Nelder-Mead", options=opts)
    re = minimize_scalar(lambda d: -eta(lambda t: np.where(t >= d, np.exp(-kc * np.maximum(t - d, 0)), 0.0)),
                         bounds=(g[0], c + 3 * s), method="bounded", options={"xatol": 1e-9 * s})
    return -rg.fun, -rr.fun, -re.fun


def symplectic_min(c):
    n = c.shape[0] // 2
    om = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    return float(np.min(np.abs(np.linalg.eigvals(1j * om @ c))))


def row(base, x, dtk, rho):
    b1 = pulsed_opo(base, x_ss=x, escape_efficiency=1.0, dt_kappa=dtk)
    b = pulsed_opo(base, x_ss=x, escape_efficiency=rho, dt_kappa=dtk)
    return b1.supermode_r[0], b1.schmidt_number, 10 * np.log10(b.optimal_mode[1]), *best_overlaps(b1)


if __name__ == "__main__":
    names = ["r0", "K", "R-_dB(rho)", "eta_gauss", "eta_rect", "eta_exp"]
    extrap = {}
    for label, base, x, rho in (("Borealis x_ss=2.5 rho=0.95", BOREALIS, 2.5, 0.95),
                                ("Borealis x_ss=2.4 rho=0.95", BOREALIS, 2.4, 0.95),
                                ("Liu x_ss=2.4 rho=0.99", LIU, 2.4, 0.99)):
        res = {n: np.array(row(base, x, 1 / n, rho)) for n in (20, 40, 80, 160)}
        print(f"\n{label}")
        for n, v in res.items():
            print(f"  dt=1/({n:3d} kc): " + "  ".join(f"{k}={a:.5f}" for k, a in zip(names, v)))
        rich = res[160] + (res[160] - res[80]) / 3
        extrap[label] = rich
        print("  Richardson     : " + "  ".join(f"{k}={a:.5f}" for k, a in zip(names, rich)))
        print("  ratio (20-40)/(40-80), (40-80)/(80-160):",
              np.round((res[20] - res[40]) / (res[40] - res[80]), 2), np.round((res[40] - res[80]) / (res[80] - res[160]), 2))
    for label, base, rho in (("Borealis", BOREALIS, 0.95), ("Liu", LIU, 0.99)):
        for dtk in (1 / 40, 1 / 80):
            b = pulsed_opo(base, escape_efficiency=rho, dt_kappa=dtk)
            b1 = pulsed_opo(base, escape_efficiency=1.0, dt_kappa=dtk)
            lam = np.linalg.eigvalsh(b1.cov_kernel)
            print(f"{label} dt=1/{round(1 / dtk)}: min symplectic eig (rho={rho}) = {symplectic_min(b.cov_kernel):.6f};"
                  f" rho=1: max|lam lam'-1| = {np.max(np.abs(lam * lam[::-1] - 1)):.2e},"
                  f" min symplectic eig = {symplectic_min(b1.cov_kernel):.6f}")
    # self-checks (D60-5 errata; T17 §5.2 R_− row)
    bor = extrap["Borealis x_ss=2.5 rho=0.95"]
    assert abs(bor[0] - 1.038) < 1e-3 and abs(bor[0] - 1.029) > 5e-3  # r_0(2.5): D60-5 1.04 confirmed, T17 1.029 refuted
    assert abs(bor[2] + 7.717) < 2e-3                                   # R_− = −7.72 dB, not T17's −7.66 dB
    assert all(abs(v[3] - 0.955) < 1e-3 for v in extrap.values())       # Gaussian LO 0.955 (T17: 0.951)
