"""Phase 5 review of the WP6-C fix (D98, reviews/phase6_wp6c_fix.md): independent checks of seed_psa, the squeezed-vacuum
background b' of photon_counting_fringe, cal_slot_check and the slow-loop equivalent delay. Read-only (numpy only for
the physics; the blocks are called only to compare). About 1 minute.

  uv run python reviews/scripts/p6_wp6c_fix_check.py [A B C D E F]

A  seed_psa against S† a S of conv §3–§4 (field map from the operator, the squeeze angle from the variance), dψ/dδ,
   the ±π/4 offset, the pump-phase sensitivity
B  b' of a squeezed vacuum at an on/off detector: exact Gaussian vacuum probability of the two-port 50:50 output mode
   vs sqvac_background (two independent half modes); its dependence on Θ and r
C  the coherent seed displaced on top of the squeezed vacuum: exact click probabilities of the two outputs, the
   estimator X̂ the block would form (b' known and subtracted), its bias vs the seed / squeeze-axis angle β
D  estimator variance / Cramér-Rao with b' (Monte Carlo) vs gates per window: the few-count CHK threshold
E  cal_slot_check numbers (F_liu20) and the detector-side time of the next data slot after a seed
F  slow loop: exact discrete stability limit with latency τ_L vs the check's τ_eq = (π/4)T + τ_L
"""

from __future__ import annotations

import math
import sys

import numpy as np

from cranberry.blocks.chk.cal_slot import cal_slot
from cranberry.blocks.err.photon_counting import fisher_per_gate, sqvac_background
from cranberry.blocks.source.seed_psa import _step as psa_step

R, ETA, QE, MU, V = 1.8, 0.63, 0.93, 4.2, 0.95


def sec(name):
    print(f"\n=== {name} " + "=" * (70 - len(name)))


# ---------------------------------------------------------------- A
def squeeze_matrix(r, phi_xi):
    """Symplectic matrix (xxpp, 1 mode = (q, p)) of S(ξ): from S† a S = a cosh r − a† e^{iφ} sinh r (conv §3)."""
    c, s = math.cosh(r), math.sinh(r)
    # a' = c a − s e^{iφ} a†  →  q' = (a' + a'†)/√2, p' = (a' − a'†)/(i√2)
    return np.array([[c - s * math.cos(phi_xi), -s * math.sin(phi_xi)],
                     [-s * math.sin(phi_xi), c + s * math.cos(phi_xi)]])


def part_a():
    sec("A seed_psa vs the operator (conv §3-§4)")
    rng = np.random.default_rng(3)
    worst_ang = worst_int = worst_blk = 0.0
    for _ in range(2000):
        r, phi_p, phi_s = rng.uniform(0, 3), rng.uniform(-np.pi, np.pi), rng.uniform(-np.pi, np.pi)
        phi_xi = phi_p + math.pi                      # conv §4: ξ = κt|β_p| e^{i(φ_p + π)}
        M = squeeze_matrix(r, phi_xi)
        # squeeze angle from the variance: Var(x_θ) of S|0> minimal at θ_sq
        cov = 0.5 * M @ M.T
        w, vec = np.linalg.eigh(cov)
        th_sq = math.atan2(vec[1, 0], vec[0, 0])
        assert abs(((th_sq - (phi_p / 2 + math.pi / 2)) + math.pi / 2) % math.pi - math.pi / 2) < 1e-9
        d_in = np.array([math.cos(phi_s), math.sin(phi_s)]) * math.sqrt(2)   # |α| = 1, <x_θ> = √2|α|cos(θ − φ_s)
        d_out = M @ d_in
        ang, inten = math.atan2(d_out[1], d_out[0]), 0.5 * d_out @ d_out
        th_amp = phi_p / 2
        delta = phi_s - th_amp
        ang_f = th_amp + math.atan2(math.exp(-r) * math.sin(delta), math.exp(r) * math.cos(delta))
        int_f = math.exp(2 * r) * math.cos(delta) ** 2 + math.exp(-2 * r) * math.sin(delta) ** 2
        worst_ang = max(worst_ang, abs(math.remainder(ang - ang_f, 2 * math.pi)))
        worst_int = max(worst_int, abs(inten / int_f - 1))
        # the block: in_sqz = θ_sq, in_seed = φ_s, delta_rad = 0 → its δ = φ_s − θ_sq = (φ_s − θ_amp) − π/2;
        # with delta_rad = π/2 its δ is φ_s − θ_amp (module docstring: delta_rad is measured against θ_amp)
        S = np.zeros(4)
        S[0], S[1] = phi_s, phi_p / 2 + math.pi / 2
        P = np.array([r, math.pi / 2, 1.0])
        psa_step(S, None, P, 0, 1, 2, 3, 0, 1, 2, True)
        worst_blk = max(worst_blk, abs(math.remainder(S[2] - ang, 2 * math.pi)), abs(S[3] / inten - 1))
    print(f"operator vs closed form: max |Δangle| {worst_ang:.2e} rad, max |Δintensity|/I {worst_int:.2e}; "
          f"block (delta_rad = π/2 ↔ φ_s − θ_amp) vs operator {worst_blk:.2e}")
    print("  => with in_seed − in_sqz = 0 the block's delta_rad IS the angle to θ_amp (beam phases carry only the "
          "deviations): out = θ_sq − π/2 + ψ(δ)")

    def psi(d, r=R):
        return math.atan2(math.exp(-r) * math.sin(d), math.exp(r) * math.cos(d))

    h = 1e-6
    for d, name in ((math.pi / 2, "π/2 (port k, quadrature)"), (0.0, "0 (port 0)"), (math.pi / 4, "+π/4"),
                    (-math.pi / 4, "−π/4")):
        g = (psi(d + h) - psi(d - h)) / (2 * h)
        g_doc = math.exp(-2 * R) / (math.cos(d) ** 2 + math.exp(-4 * R) * math.sin(d) ** 2)
        # pump phase φ_p: θ_sq moves by δφ_p/2 and δ by −δφ_p/2 → d out/dφ_p = ½(1 − ψ')
        print(f"  δ = {name:24s}: dψ/dδ {g:9.4f} (docstring {g_doc:9.4f});  d out/d φ_p = ½(1 − ψ') = "
              f"{0.5 * (1 - g):+8.4f};  d out/d(path common to seed and sqz) = 1")
    off = psi(math.pi / 4) - psi(-math.pi / 4)
    print(f"  ±π/4 offset ψ(π/4) − ψ(−π/4) = {off:.5f} = 2 atan(e^-2r) {2 * math.atan(math.exp(-2 * R)):.5f}; "
          f"gain of out_k on θ_sq,k (not shared with the seed) 1 − ψ' = {1 - (psi(math.pi/4 + h) - psi(math.pi/4 - h)) / (2*h):.4f}")
    print(f"  quadrature: common SEED−PUMP noise ε → X gains ψ'(π/2) − ψ'(0) = {math.exp(2*R) - math.exp(-2*R):.3f} ε "
          "(written into Θ_k by the slow loop)")


# ---------------------------------------------------------------- B, C
def lossy_sv(r, eta, theta):
    """Covariance (q, p) of a squeezed vacuum (squeeze axis θ) after transmission η (conv §2-§3)."""
    a = 0.5 * (eta * math.exp(-2 * r) + 1 - eta)
    b = 0.5 * (eta * math.exp(2 * r) + 1 - eta)
    c, s = math.cos(theta), math.sin(theta)
    Rm = np.array([[c, -s], [s, c]])
    return Rm @ np.diag([a, b]) @ Rm.T


def p0(cov, d):
    """Vacuum probability of a Gaussian state (cov, mean d): exp(−½ dᵀ(σ + I/2)⁻¹d)/√det(σ + I/2)."""
    m = cov + 0.5 * np.eye(2)
    return math.exp(-0.5 * d @ np.linalg.solve(m, d)) / math.sqrt(np.linalg.det(m))


def outputs(r, eta, th_k, th_0, alpha_k, alpha_0):
    """Covariance and mean of the two 50:50 outputs (a_k ± a_0)/√2 of two independent lossy squeezed vacua with the
    coherent seeds α_k, α_0 (amplitudes at the detector)."""
    ck, c0 = lossy_sv(r, eta, th_k), lossy_sv(r, eta, th_0)
    dk = math.sqrt(2) * np.array([alpha_k.real, alpha_k.imag])
    d0 = math.sqrt(2) * np.array([alpha_0.real, alpha_0.imag])
    cov = 0.5 * (ck + c0)                 # both outputs (the sign of a_0 does not change a zero-mean covariance)
    return cov, (dk + d0) / math.sqrt(2), (dk - d0) / math.sqrt(2)


def part_b():
    sec("B squeezed-vacuum background b' (on/off detector, vacuum probability)")
    print(f"r = {R}, η·qe = {ETA} (overlap 1)")
    code = sqvac_background(R, ETA)
    print(f"  sqvac_background (two independent half modes, η/2 each): b' = {code:.4f}")
    for th in (0.0, 0.01, 0.1, math.pi / 4, math.pi / 2):
        cov, _, _ = outputs(R, ETA, th, 0.0, 0j, 0j)
        b = -math.log(p0(cov, np.zeros(2)))
        print(f"  exact output mode, Θ = {th:6.3f} rad: b' = {b:.4f}  ({b / code:.3f} × code)")
    n = ETA * math.sinh(R) ** 2
    print(f"  checks: Θ = 0 single-mode SV ½ln[(a+½)(b+½)] with full η; Θ = π/2 thermal ln(1 + η sinh²r) = "
          f"{math.log(1 + n):.4f}; mean photons per output η sinh²r = {n:.3f} in every case")
    for r in (0.5, 1.0, 1.8, 2.5):
        cov, _, _ = outputs(r, ETA, 0.0, 0.0, 0j, 0j)
        print(f"  r = {r}: code {sqvac_background(r, ETA):.4f}, exact (Θ = 0) {-math.log(p0(cov, np.zeros(2))):.4f}")
    # pump fluctuation of b': common to both outputs
    dr = 1e-4
    db = (sqvac_background(R + dr, ETA) - sqvac_background(R - dr, ETA)) / (2 * dr)
    print(f"  db'/dr = {db:.3f} (code): a 1 % pump-power step (δr ≈ r·0.5 % = {R * 0.005:.4f}) moves b' by "
          f"{db * R * 0.005:.4f}, identical in both outputs → a slope error at X = π/2, no bias")


def x_hat(l1, l0, vis):
    return math.acos(max(-1.0, min(1.0, (l1 - l0) / (vis * (l1 + l0)))))


def part_c():
    sec("C seed displaced on the squeezed vacuum: exact clicks vs the block's model")
    m = 0.5 * QE * MU                     # detected seed photons per output (balanced, block's m)
    amp = math.sqrt(m)                     # |α_k| = |α_0| at the detector with |α_k ± α_0|²/2 = m(1 ± cos X)
    covz, _, _ = outputs(R, ETA, 0.0, 0.0, 0j, 0j)
    bq = -math.log(p0(covz, np.zeros(2)))  # the exact b' (Θ = 0), assumed known and subtracted
    print(f"m = qe·μ/2 = {m:.3f} per output, V = 1, Θ = 0, b' (exact, subtracted) = {bq:.4f}")
    print("  β = seed-0 phase − squeeze axis at the GBS input; X = π/2 (lock point)")
    rows = []
    for beta in np.linspace(0, math.pi / 2, 7):
        X = math.pi / 2
        ak = amp * np.exp(1j * (beta + X))
        a0 = amp * np.exp(1j * beta)
        cov, d1, d2 = outputs(R, ETA, 0.0, 0.0, ak, a0)
        l1 = -math.log(p0(cov, d1)) - bq
        l2 = -math.log(p0(cov, d2)) - bq
        xh = x_hat(l1, l2, 1.0)
        rows.append((beta, xh - math.pi / 2))
        print(f"  β = {beta:5.3f}: λ_+ − b' = {l1:.4f}, λ_− − b' = {l2:.4f} (block model: {m:.4f} each) → "
              f"X̂ − π/2 = {xh - math.pi / 2:+.4f} rad")
    h = 1e-4
    def bias(beta, X=math.pi / 2):
        ak, a0 = amp * np.exp(1j * (beta + X)), amp * np.exp(1j * beta)
        cov, d1, d2 = outputs(R, ETA, 0.0, 0.0, ak, a0)
        return x_hat(-math.log(p0(cov, d1)) - bq, -math.log(p0(cov, d2)) - bq, 1.0) - math.pi / 2
    sl = (bias(h) - bias(-h)) / (2 * h)
    A = 1 / (covz[0, 0] + 0.5)
    B = 1 / (covz[1, 1] + 0.5)
    print(f"  closed form: λ_± − b' = m[(A+B)/2 ∓ (A−B)/2 · sin 2β], A = 1/(a+½) = {A:.3f}, B = 1/(b+½) = {B:.3f}; "
          f"cos X̂ = −(A−B)/(A+B)·sin 2β = −{(A - B) / (A + B):.3f} sin 2β")
    print(f"  dX̂/dβ at β = 0: {sl:+.3f} rad/rad (so β must be held at 0 or π/2 and its noise enters X̂ ×{abs(sl):.2f})")
    # slope of X̂ against the true X at β = 0 (fringe contrast)
    def xh_of(X, beta=0.0):
        ak, a0 = amp * np.exp(1j * (beta + X)), amp * np.exp(1j * beta)
        cov, d1, d2 = outputs(R, ETA, 0.0, 0.0, ak, a0)
        return x_hat(-math.log(p0(cov, d1)) - bq, -math.log(p0(cov, d2)) - bq, 1.0)
    g = (xh_of(math.pi / 2 + h) - xh_of(math.pi / 2 - h)) / (2 * h)
    print(f"  at β = 0: dX̂/dX = {g:.3f} (block model 1): a gain error only")
    # before_squeezer quadrature preparation: port k on the squeezed axis, port 0 on the amplified axis (β = −π/2 + X)
    print("  psa quadrature (seed k on θ_sq, seed 0 on θ_amp): β fixed by the preparation, outputs at ±π/4 to the axes:"
          f" bias {bias(math.pi / 2):+.2e} rad (symmetric), slope dX̂/dβ {(bias(math.pi/2 + h) - bias(math.pi/2 - h)) / (2*h):+.3f}")
    # Fisher information per gate (both outputs, V = 1) of the exact click law vs the block's Poisson-equivalent law
    def fisher_exact(beta, X=math.pi / 2, hh=1e-5):
        def pc(XX):
            ak, a0 = amp * np.exp(1j * (beta + XX)), amp * np.exp(1j * beta)
            cov, d1, d2 = outputs(R, ETA, 0.0, 0.0, ak, a0)
            return np.array([1 - p0(cov, d1), 1 - p0(cov, d2)])
        p, dp = pc(X), (pc(X + hh) - pc(X - hh)) / (2 * hh)
        return float(np.sum(dp ** 2 / (p * (1 - p))))
    f_code = fisher_per_gate(m, 1.0, math.pi / 2, dark=sqvac_background(R, ETA))
    f_ex_b = fisher_per_gate(m, 1.0, math.pi / 2, dark=bq)
    print(f"  Fisher per gate (V = 1): block law with code b' {f_code:.4f}, block law with exact b' {f_ex_b:.4f}, "
          f"exact Gaussian β = 0 {fisher_exact(0.0):.4f}, β = π/8 {fisher_exact(math.pi / 8):.4f}, "
          f"β = π/4 {fisher_exact(math.pi / 4):.4f}")
    # gated pump: no squeezed vacuum
    print("  dbs_cal_pump gated (r = 0 in the calibration slot): A = B = 1, no β dependence")
    # overlap 0.5 as loss
    cz, _, _ = outputs(R, 0.5 * ETA, 0.0, 0.0, 0j, 0j)
    A2, B2 = 1 / (cz[0, 0] + 0.5), 1 / (cz[1, 1] + 0.5)
    print(f"  overlap 0.5 (as loss): (A−B)/(A+B) = {(A2 - B2) / (A2 + B2):.3f}; r = 0.53 (F_bor6): "
          f"{(lambda c: (1/(c[0,0]+.5) - 1/(c[1,1]+.5)) / (1/(c[0,0]+.5) + 1/(c[1,1]+.5)))(outputs(0.53, ETA, 0, 0, 0j, 0j)[0]):.3f}")


# ---------------------------------------------------------------- D
def part_d():
    sec("D estimator vs Cramér-Rao with b' (Monte Carlo, numpy)")
    rng = np.random.default_rng(11)
    b = sqvac_background(R, ETA)
    m = 0.5 * QE * MU
    x = math.pi / 2
    f = fisher_per_gate(m, V, x, dark=b)
    for mm in (0.5 * QE * 3.43, m):
        print(f"  F per gate (both outputs) m = {mm:.3f}: b' = 0 {fisher_per_gate(mm, V, x):.4f}, b' = {b:.3f} "
              f"{fisher_per_gate(mm, V, x, dark=b):.4f}")
    ms = np.linspace(0.5, 4, 3501)
    fs_ = [fisher_per_gate(q, V, x, dark=b) for q in ms]
    print(f"  optimum with b' {b:.3f}: m = {ms[int(np.argmax(fs_))]:.3f} → μ = {2 * ms[int(np.argmax(fs_))] / QE:.2f}")
    lam = np.array([b + m * (1 - V * math.cos(x)), b + m * (1 + V * math.cos(x))])
    pc = -np.expm1(-lam)
    print(f"  click probability per gate and output {pc[0]:.4f} (expected no-click gates per window = G·{1 - pc[0]:.4f})")
    for G in (100, 300, 1000, 3000, 6250):
        n = 40000
        k = rng.binomial(G, pc[:, None], size=(2, n)).astype(float)
        top = 1 - 0.5 / G
        lh = -np.log(1 - np.minimum(k / G, top)) - b
        c = np.clip((lh[1] - lh[0]) / (V * (lh[1] + lh[0])), -1, 1)
        xh = np.arccos(c)
        allclick = np.mean((k[0] == G) | (k[1] == G))
        print(f"  G = {G:5d}: var/CRB {np.var(xh) * G * f:6.3f}, bias {np.mean(xh) - x:+.2e} rad, windows with an "
              f"all-click output {allclick:.3f}, expected no-click gates {G * (1 - pc[0]):6.1f}")


# ---------------------------------------------------------------- E
def part_e():
    sec("E cal_slot_check (F_liu20: K 32, N_b 4, 20 MHz)")
    b = sqvac_background(R, ETA)
    base = dict(placement="cal_frames", N_b=4, n_guard=0, K=32, t_rep_s=5e-8, cal_fraction=0.2, t_int_s=0.05,
                mu_per_gate=MU, qe=QE, bg_per_gate=b + 100 * 2e-9, dead_time_s=1e-8, recovery_s=1.1e-8,
                sigma_target_rad=0.01, lock_ratio_R=1.0, min_gates=100, recovery_loss_max=0.01)
    for name, over in (("cal_frames", {}), ("guard", dict(placement="guard", n_guard=1, cal_fraction=1.0)),
                       ("ideal_free", dict(placement="ideal_free", cal_fraction=1.0))):
        r = cal_slot({**base, **over})
        print(f"  {name:10s}: gates/frame {r['gates_per_frame']:.0f}, G/window {r['gates_per_window']:.0f}, data "
              f"{r['data_fraction']:.2f}, L_k ×{r['lk_factor']:.3f}, Δτ {r['dtau_s'] * 1e6:.2f} µs, Δν_max "
              f"{r['dnu_max_hz']:.3f} Hz (without slot {r['dnu_max_without_slot_hz']:.3f}), p_click {r['p_click']:.4f}, "
              f"gap {r['recovery_gap_s'] * 1e9:.0f} ns, η_next {r['eta_next_slot']:.4f}, weight "
              f"{r['recovery_weight']:.2f}, mean loss {r['recovery_loss_mean']:.4f}")
    # independent: D78 Δν ≤ σ²/(2πRΔτ), Δτ = (K−1)(N_b + n_guard)/f_rep
    for ng in (0, 1):
        dt = 31 * (4 + ng) * 5e-8
        print(f"  D78 independent: n_guard {ng}: Δτ {dt * 1e6:.2f} µs, Δν_max {1e-4 / (2 * math.pi * dt):.3f} Hz")
    # recovery, by hand
    lam = 0.5 * QE * MU + b + 100 * 2e-9
    r50 = -math.expm1(-(50e-9 - 1e-8) / 1.1e-8)
    print(f"  by hand: λ = {lam:.4f}, p = {-math.expm1(-lam):.4f}, r(50 ns) = {r50:.5f}, η_next = "
          f"{1 + math.expm1(-lam) * (1 - r50):.4f}")
    # detector timeline: GBS bursts of B slots every T_f; the next data burst after a calibration frame
    for ng in (0, 1):
        B = 4 + ng
        tf = 32 * B * 5e-8
        gap_det = tf - (4 - 1) * 5e-8
        print(f"  detector side, n_guard {ng}: T_f {tf * 1e6:.2f} µs; own seeds of a calibration frame → first data "
              f"slot of the next frame {gap_det * 1e6:.2f} µs later (r = {-math.expm1(-(gap_det - 1e-8) / 1.1e-8):.6f}); "
              f"seed gates inside one burst T_rep apart (r = {r50:.4f})")


# ---------------------------------------------------------------- F
def part_f():
    sec("F slow loop: exact discrete limit with latency τ_L vs τ_eq = (π/4)T + τ_L")

    def stable(a, phi):
        # e_{n+1} = c1 e_n + c2 e_{n−1} + c3 e_{n−2} (window average → hold → integrator, latency φT; header)
        c1 = 1 - a * (1 - phi) ** 2 / 2
        c2 = -a + a * phi ** 2 / 2 + a * (1 - phi) ** 2 / 2
        c3 = -a * phi ** 2 / 2
        return max(abs(np.roots([1, -c1, -c2, -c3]))) < 1

    def limit(phi):
        lo, hi = 1e-3, 4.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if stable(mid, phi) else (lo, mid)
        return lo / (2 * math.pi)          # f_c·T

    for phi in (0.0, 2e-5, 0.01, 0.1, 0.3, 0.5, 1.0):
        print(f"  τ_L/T = {phi:7.5f}: exact f_c·T < {limit(phi):.4f};  check boundary 1/(π + 4τ_L/T) = "
              f"{1 / (math.pi + 4 * phi):.4f}")
    # time-domain cross-check of the recursion: per-window growth factor vs the largest root
    for a_fc, phi in ((0.30, 0.0), (0.33, 0.0), (0.25, 0.3), (0.2, 1.0)):
        part_f_ts(2 * math.pi * a_fc, phi)


def roots_max(a, phi):
    c1 = 1 - a * (1 - phi) ** 2 / 2
    c2 = -a + a * phi ** 2 / 2 + a * (1 - phi) ** 2 / 2
    c3 = -a * phi ** 2 / 2
    return max(abs(np.roots([1, -c1, -c2, -c3])))


def part_f_ts(a, phi, n_win=80, sub=1000):
    """Time-stepped: window mean of x = −U(t − τ_L), held one window, U ramps at a·ê/T (T = 1)."""
    lag = int(round(phi * sub))
    U = np.zeros(n_win * sub + 1)
    held, es = 1.0, []
    for n in range(n_win):
        acc = 0.0
        for j in range(sub):
            i = n * sub + j
            U[i + 1] = U[i] + a * held / sub
            acc -= U[i - lag] if i - lag >= 0 else 0.0
        held = acc / sub
        es.append(abs(held))
    es = np.array(es)
    g = (np.sqrt(np.mean(es[60:80] ** 2)) / np.sqrt(np.mean(es[40:60] ** 2))) ** (1 / 20)
    print(f"    f_c·T = {a / (2 * math.pi):.3f}, τ_L/T = {phi}: time-stepped growth per window {g:.4f}, "
          f"recursion max|z| {roots_max(a, phi):.4f}")


PARTS = {"A": part_a, "B": part_b, "C": part_c, "D": part_d, "E": part_e, "F": part_f}

if __name__ == "__main__":
    for k in (sys.argv[1:] or PARTS):
        PARTS[k]()
