"""Phase 5 review of WP6-C (D92-1, D92-3): click statistics, MLE, Cramér-Rao, RNG chunk invariance, and the squeezed
vacuum background of a before-squeezer seed. Independent of the block's own helpers except where the block itself is
the object under test (sections C and F).

  uv run python reviews/scripts/p6_wp6c_counting.py            (≈ 1 min)

A  Fisher information: closed form (fisher_per_gate) vs a numerical derivative of the binomial log-likelihood; the
   CRB with μ·qe known (1 parameter) vs unknown (2 parameters, the normalised estimator's situation).
B  MLE: numerical maximisation of the 2-output binomial likelihood over (m, X) vs the block's closed form.
C  Monte Carlo (numpy, independent) and the block itself at the nominal point (G = T_int/t_frame = 7812 gates,
   μ = 3.43, qe 0.93, V 0.95, dark 1e-6): bias and variance vs CRB, at X = π/2 and off it.
D  Bias from the unsubtracted background b and from dead time off π/2 (both vanish at the π/2 lock point).
E  Squeezed vacuum in the calibration gate (seed before the squeezer: the pump must be on in that slot): exact vacuum
   probability of the two lossy squeezed-vacuum halves at each output → equivalent Poisson background b', Fisher
   information and σ_X per window, re-optimised μ.
F  RNG chunk invariance of the block (different piece boundaries, bitwise).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cranberry.blocks.err.photon_counting import fisher_per_gate
from cranberry.core.compile import run_compiled
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

QE, V, DARK = 0.93, 0.95, 1.0e-6
G_NOM = round(0.05 / 6.4e-6)     # F_liu20: T_int 50 ms, one calibration gate per 6.4 µs frame
MU_NOM = 3.43


def lam(m, x, d=0.0, b=0.0, vis=V):
    return np.array([m * (1 + vis * math.cos(x)) + d + b, m * (1 - vis * math.cos(x)) + d + b])


def loglik(n, g, m, x, d=0.0):
    la = lam(m, x, d)
    p = -np.expm1(-la)
    return float(np.sum(n * np.log(p) + (g - n) * (-la)))


def fisher_matrix(m, x, d=0.0, h=1e-5):
    """2×2 Fisher matrix in (m, X) per gate from the binomial model (numerical derivatives of λ)."""
    def dl(i):
        e = np.zeros(2)
        e[i] = h
        return (lam(m + e[0], x + e[1], d) - lam(m - e[0], x - e[1], d)) / (2 * h)
    la = lam(m, x, d)
    w = 1.0 / np.expm1(la)
    J = np.array([dl(0), dl(1)])          # J[i, out]
    return (J * w) @ J.T


def estimator(n1, n0, g, d=0.0, vis=V):
    """The block's closed form (re-implemented here): λ̂ = −ln(1 − n/G) clamped, minus dark, normalised ratio."""
    top = 1.0 - 0.5 / g
    l1 = -np.log(1.0 - np.minimum(n1 / g, top)) - d
    l0 = -np.log(1.0 - np.minimum(n0 / g, top)) - d
    c = (l1 - l0) / (vis * (l1 + l0))
    return np.arccos(np.clip(c, -1, 1))


def sec_a():
    print("A  Fisher information per gate (both outputs), m = qe μ/2")
    m = 0.5 * QE * MU_NOM
    for x in (math.pi / 2, math.pi / 2 + 0.3, math.pi / 2 + 0.8):
        F = fisher_matrix(m, x, DARK)
        f_closed = fisher_per_gate(m, V, x, DARK)
        crb1 = 1.0 / F[1, 1]
        crb2 = np.linalg.inv(F)[1, 1]
        print(f"   X−π/2 = {x - math.pi/2:+.1f}: F_XX closed {f_closed:.6f} numeric {F[1, 1]:.6f};"
              f"  CRB(μ unknown)/CRB(μ known) = {crb2 / crb1:.4f}")


def sec_b(rng):
    print("B  closed form vs numerical MLE over (m, X) (binomial, 2 outputs)")
    m = 0.5 * QE * MU_NOM
    worst = 0.0
    for x in (math.pi / 2, 1.2, 2.3):
        for _ in range(20):
            la = lam(m, x, DARK)
            n = rng.binomial(G_NOM, -np.expm1(-la))
            xc = float(estimator(n[0], n[1], G_NOM, DARK))
            res = minimize(lambda p: -loglik(n, G_NOM, p[0], p[1], DARK), x0=[m, xc], method="Nelder-Mead",
                           options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 4000})
            worst = max(worst, abs(res.x[1] - xc))
    print(f"   max |X̂_closed − X̂_numerical MLE| = {worst:.2e} rad (60 windows)")


def block_windows(x, mu, g, n_windows, seed=11):
    objs = [{"block": "ramp", "name": "a", "params": {"out": "v.a", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                      "period_s": 1.0, "offset_v": x, "toggles": {"enabled": True}}},
            {"block": "photon_counting_fringe", "name": "pcf",
             "params": {"in_a": "v.a", "in_b": None, "bg_in": None, "gate_in": None, "out": "v.pcf", "f_gate_hz": 1e6,
                        "cal_fraction": 1.0, "window_mode": "gates", "window_n_gates": g, "mu_per_gate": mu,
                        "visibility": V, "qe": QE, "dark_per_gate": DARK, "bg_per_gate": 0.0, "phi_offset_rad": 0.0,
                        "setpoint_rad": math.pi / 2, "dead_time_s": 1e-8, "recovery_s": 1.1e-8, "sqvac_per_gate": 0.0,
                        "toggles": {"counting_noise": True, "saturation": True, "dark": True, "background": True,
                                    "dead_time": False, "squeezed_vacuum": False}}}]
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(objs), lay, seed)
    st, pr = State(lay), Params(lay)
    out = []
    for _ in range(n_windows):
        run_compiled(blk, st, pr, g)
        out.append(float(st["cnt.pcf.x_hat"]))
    return np.array(out)


def sec_c(rng):
    print(f"C  estimator at the nominal point (G = {G_NOM}, μ = {MU_NOM}): bias and variance / CRB(μ unknown)")
    m = 0.5 * QE * MU_NOM
    for x in (math.pi / 2, math.pi / 2 + 0.3):
        crb = np.linalg.inv(fisher_matrix(m, x, DARK))[1, 1] / G_NOM
        la = lam(m, x, DARK)
        n = rng.binomial(G_NOM, -np.expm1(-la), size=(200000, 2))
        e = estimator(n[:, 0], n[:, 1], G_NOM, DARK) - x
        eb = block_windows(x, MU_NOM, G_NOM, 6000) - x
        print(f"   X−π/2 = {x - math.pi/2:+.1f}: numpy bias {e.mean():+.2e} (±{e.std()/math.sqrt(e.size):.1e}) var/CRB "
              f"{e.var()/crb:.4f};  block bias {eb.mean():+.2e} (±{eb.std()/math.sqrt(eb.size):.1e}) var/CRB "
              f"{eb.var()/crb:.3f};  σ_X = {1e3*math.sqrt(crb):.2f} mrad/window")


def sec_d():
    print("D  bias of X̂ (expected counts) from an unsubtracted background b (photons/gate/output at the detector)"
          " and from a dead-time factor r± ")
    m = 0.5 * QE * MU_NOM
    for x in (math.pi / 2, math.pi / 2 + 0.3):
        for b in (0.0, 0.1, 1.0):
            la = lam(m, x, DARK, b)
            p = -np.expm1(-la)
            xh = float(estimator(G_NOM * p[0], G_NOM * p[1], G_NOM, DARK))
            print(f"   X−π/2 = {x - math.pi/2:+.1f}, b = {b:4.1f}: X̂ − X = {xh - x:+.4f} rad;"
                  f"  slope dX̂/dX at π/2 = m/(m+b) = {m/(m+b):.3f}")
    # dead time: r = 0.97 on the brighter output only (asymmetric recovery)
    la = lam(m, math.pi / 2 + 0.3, DARK) * np.array([0.97, 1.0])
    p = -np.expm1(-la)
    print(f"   dead-time r = 0.97 on one output at X−π/2 = +0.3: X̂ − X = "
          f"{float(estimator(G_NOM*p[0], G_NOM*p[1], G_NOM, DARK)) - (math.pi/2 + 0.3):+.4f} rad (0 at π/2 when r+ = r−)")


def p0_sqvac(r, eta):
    """Vacuum probability of a single-mode squeezed vacuum (squeeze r) after transmission eta (Gaussian formula
    p0 = 1/sqrt(det(σ + I/2)), vacuum σ = I/2)."""
    a = 0.5 * (eta * math.exp(-2 * r) + 1 - eta)
    b = 0.5 * (eta * math.exp(2 * r) + 1 - eta)
    return 1.0 / math.sqrt((a + 0.5) * (b + 0.5))


def sec_e():
    print("E  squeezed vacuum in the calibration gate (pump on: seed before the squeezer). Each output sees half of the"
          " port-k and half of the port-0 squeezed vacuum (×η·qe), equivalent Poisson b' = −ln p0")
    f0 = fisher_per_gate(0.5 * QE * MU_NOM, V, math.pi / 2, DARK)
    sig0 = 1.0 / math.sqrt(G_NOM * f0)
    print(f"   no background: F = {f0:.3f}/gate, σ_X = {1e3*sig0:.2f} mrad per 50 ms window")
    for r in (0.5, 1.0, 1.8):
        for eta in (0.59, 0.68):
            bp = -2.0 * math.log(p0_sqvac(r, 0.5 * eta))
            f_nom = fisher_per_gate(0.5 * QE * MU_NOM, V, math.pi / 2, DARK + bp)
            res = minimize_scalar(lambda m: -fisher_per_gate(m, V, math.pi / 2, DARK + bp), bounds=(0.01, 10),
                                  method="bounded")
            f_opt = -res.fun
            print(f"   r = {r:.1f} (sinh² r = {math.sinh(r)**2:5.2f}), η·qe = {eta:.2f}: b' = {bp:5.3f};  μ = 3.43: "
                  f"F = {f_nom:.4f} (σ_X ×{math.sqrt(f0/f_nom):5.2f});  best μ = {2*res.x/QE:4.2f}: σ_X ×{math.sqrt(f0/f_opt):5.2f}")


def sec_f():
    print("F  RNG chunk invariance of the block (piece boundaries 1+999+…, 7, 65536)")
    res = []
    for pieces in ([1, 999, 4000, 15000], [7] * 2857 + [1], [20000]):
        objs = [{"block": "ramp", "name": "a", "params": {"out": "v.a", "waveform": "sawtooth", "amplitude_v": 0.3,
                                                          "period_s": 3e-3, "offset_v": 1.4, "toggles": {"enabled": True}}},
                {"block": "photon_counting_fringe", "name": "pcf",
                 "params": {"in_a": "v.a", "out": "v.pcf", "f_gate_hz": 2.7e6, "cal_fraction": 0.61,
                            "window_mode": "time", "window_t_s": 3.3e-4, "mu_per_gate": 3.4, "visibility": V, "qe": QE,
                            "dark_per_gate": 1e-3, "bg_per_gate": 0.01, "phi_offset_rad": 0.0, "setpoint_rad": 1.5,
                            "dead_time_s": 1e-7, "recovery_s": 2e-7, "sqvac_per_gate": 0.0,
                            "toggles": {"counting_noise": True, "saturation": True, "dark": True, "background": True,
                                        "dead_time": True, "squeezed_vacuum": False}}}]
        lay = SlotLayout(1e6)
        blk = build(spec_from_obj(objs), lay, 5)
        st, pr = State(lay), Params(lay)
        for n in pieces:
            run_compiled(blk, st, pr, n)
        res.append(tuple(float(np.sum(st[f"cnt.pcf.{k}"])) for k in ("n0", "n1", "x_hat", "t_last", "acc"))
                   + (int(st["cnt.pcf.window_index"]),))
    print(f"   identical: {res[0] == res[1] == res[2]}  {res[0]}")


def psa(alpha, r):
    """Coherent seed α through a degenerate phase-sensitive amplifier, amplified quadrature along the real axis (pump
    phase 0): α → α cosh r + α* sinh r."""
    return alpha * math.cosh(r) + np.conj(alpha) * math.sinh(r)


def sec_g():
    print("G  seed before the squeezer: the π/2 seed-preparation offset through a phase-sensitive amplifier (gain e^{±r})")
    for r in (0.5, 1.0, 1.8):
        # (i) port 0 prepared at 0, port k at π/2 (both eigen-quadratures): offset kept, intensities e^{±2r}
        a0, ak = psa(1.0 + 0j, r), psa(1j, r)
        vis = 2 * abs(a0) * abs(ak) / (abs(a0) ** 2 + abs(ak) ** 2)
        eps = 1e-4
        gain_k = (np.angle(psa(1j * np.exp(1j * eps), r)) - np.angle(ak)) / eps
        gain_0 = (np.angle(psa(np.exp(1j * eps), r)) - np.angle(a0)) / eps
        # (ii) ±π/4 preparation (equal output intensities): the offset collapses
        off = np.angle(psa(np.exp(1j * math.pi / 4), r)) - np.angle(psa(np.exp(-1j * math.pi / 4), r))
        print(f"   r = {r:.1f}: (0, π/2) prep → offset {np.angle(ak) - np.angle(a0):.3f} rad, |α_k/α_0|² = "
              f"{abs(ak)**2/abs(a0)**2:.2e}, visibility {vis:.3f} for equal input seeds; output-phase gain on the seed's "
              f"own input phase: port k {gain_k:.1f}, port 0 {gain_0:.3f} (high-gain θ-units: 0);  (±π/4) prep → "
              f"offset {off:.3f} rad instead of π/2")


def main():
    rng = np.random.default_rng(2026)
    sec_a()
    sec_b(rng)
    sec_c(rng)
    sec_d()
    sec_e()
    sec_f()
    sec_g()


if __name__ == "__main__":
    main()
