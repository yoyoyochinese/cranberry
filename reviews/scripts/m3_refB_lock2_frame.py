"""Phase 5 M3 review (Ref B): lock-2 frame under OPO detuning (D71-2), from the fields.

Independent of the WP-B6 / WP-I4 code except for nothing: plain numpy on the intracavity Langevin matrix of T02 §2.1
(φ_p = 0 frame, quadratures (q, p), d/dt → −iΩ):

    M(Ω) = [[κ − ε − iΩ, −Δ], [Δ, κ + ε − iΩ]],  ε = xκ
    squeezed field   X_out = (2κ_T M⁻¹ − 1) X_in + 2√(κ_T κ_L) M⁻¹ X_L       (ellipse: S = A A† + B B†)
    transmitted probe (coherent, injected through the HR mirror): α_out ∝ M(0)⁻¹ u(φ_in),  u = (cos φ, sin φ)

Lock 1 (dithered probe phase, DC intensity of the transmitted probe at the PD after the OPO, Fig. 1) holds φ_in at the
de-amplification extremum of |M⁻¹u|. Lock 2 (the same probe dither, HD output demodulated at f_m) has its zero where the
LO is orthogonal to the dither sideband: Re[l·M(Ω_m)⁻¹ u'(φ_in)] = 0 (b-quadrature, demod phase −π/2) — the HD signal is
l·α_out(φ_in + m sin Ω_m t), its first harmonic is m·l·M⁻¹u'.

Parts
  A  static (Ω_m → 0): output frame, lock-1 extremum shift, ellipse rotation (ρ, f_s), lock-2 zero, residual formula
  B  dither at 295 kHz through the cavity (complex A±(Ω_m)): lock-2 slope and lock-1-error transfer vs the quasi-static
     A± used by the ERR blocks
  C  dynamic: the preset's lock 1 does not sense Δ (pr_pd has no OPO coefficient). Physically the lock-1 extremum moves by
     ΔΦ1 = Δ/κ, so Δ above the lock-1 band but inside the lock-2 band adds −(A−/A+)·S1·T2·Δ/(2κ) to the metric. Linear loop
     model of refB_nominal (integrators, PZT resonances, τ_L) and the preset's OPO-length PSD.
  D  laser frequency noise is a detuning too (Δ = ω_cav − ω_L; the Lock-Beam length lock makes the cavity follow the
     laser only below f_ugf). The preset's Δ has the cavity length only; the laser reaches Δ_metric only through the
     geometric LO − sqz delay (1.67 ns). Variance of both paths for white FM h0 (and the 89 / 177 kHz FM peaks).
  E  D71-5: RIN at the dither frequency as lock-2 / lock-1 sensor noise (LO RIN through the HD's CMRR leak ε, probe RIN
     on the beat and on the lock-1 PD), vs the (squeezed) LO shot noise the preset uses. Integrator loops: in-loop
     residual σ² = S_e·(π/2)·f_u.

Run: uv run python reviews/scripts/m3_refB_lock2_frame.py
"""

import math

import numpy as np

C = 299_792_458.0
T_OC, L_RT, LAM = 0.123, 0.5, 860e-9
NU0 = C / LAM


def l_ic(x):
    return 0.00249 + 0.00222 * x                     # Ref B Eq. (3)


def kappas(x, loss=True):
    kt = C * T_OC / (2 * L_RT)
    kl = C * l_ic(x) / (2 * L_RT) if loss else 0.0
    return kt, kl


def mmat(x, kappa, delta, om=0.0):
    e = x * kappa
    return np.array([[kappa - e - 1j * om, -delta], [delta, kappa + e - 1j * om]])


def ang(v):
    """angle of a real 2-vector, folded to (−π/2, π/2]."""
    a = math.atan2(v[1], v[0])
    while a <= -math.pi / 2:
        a += math.pi
    while a > math.pi / 2:
        a -= math.pi
    return a


def fold(a):
    return (a + math.pi / 2) % math.pi - math.pi / 2


def ellipse_min_angle(x, delta, f_s=0.0, loss=True):
    kt, kl = kappas(x, loss)
    k = kt + kl
    mi = np.linalg.inv(mmat(x, k, delta, 2 * math.pi * f_s))
    a = 2 * kt * mi - np.eye(2)
    b = 2 * math.sqrt(kt * kl) * mi
    s = (a @ a.conj().T + b @ b.conj().T).real
    w, v = np.linalg.eigh(0.5 * (s + s.T))
    return ang(v[:, 0])                                # eigenvector of the smallest eigenvalue


def probe_frames(x, delta):
    """(input angle of the de-amplification extremum, output carrier angle there)."""
    kt, kl = kappas(x)
    mi = np.linalg.inv(mmat(x, kt + kl, delta)).real
    phis = np.linspace(0, math.pi, 200001)
    g = np.einsum("ij,jk->ik", mi, np.vstack([np.cos(phis), np.sin(phis)]))
    n = (g ** 2).sum(0)
    i = int(np.argmin(n))
    # refine (parabola)
    y0, y1, y2 = n[i - 1], n[i], n[i + 1]
    d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    phi = phis[i] + d * (phis[1] - phis[0])
    out = mi @ np.array([math.cos(phi), math.sin(phi)])
    return fold(phi), ang(out)


def lock2_zero(x, delta, phi_in, f_m=0.0):
    """LO angle with Re[l·M(Ω_m)⁻¹u'(φ_in)] = 0 (the dither-synchronous quadrature of the demodulated HD output)."""
    kt, kl = kappas(x)
    mi = np.linalg.inv(mmat(x, kt + kl, delta, 2 * math.pi * f_m))
    up = np.array([-math.sin(phi_in), math.cos(phi_in)])
    s = mi @ up                                        # complex sideband vector
    sr = s.real
    return ang(np.array([-sr[1], sr[0]]))              # LO ⟂ Re(sideband)


def part_a():
    print("A. static frames (Ω_m → 0), Δ/κ = 1e-4; coefficients in units of Δ/κ")
    for x in (0.236, 0.527, 0.745, 0.913):
        kt, kl = kappas(x)
        k = kt + kl
        d = 1e-4 * k
        in0, out0 = probe_frames(x, 0.0)
        in1, out1 = probe_frames(x, d)
        th0 = ellipse_min_angle(x, 0.0)
        th1 = ellipse_min_angle(x, d)
        d0 = 1e-4 * kt                                 # lossless cavity: Δ/κ with its own κ = κ_T
        th1_lossless = ellipse_min_angle(x, d0, loss=False) - ellipse_min_angle(x, 0.0, loss=False)
        th1_1mhz = ellipse_min_angle(x, d, f_s=1e6) - ellipse_min_angle(x, 0.0, f_s=1e6)
        z = lock2_zero(x, d, in1)
        z_wrong = lock2_zero(x, d, in0)               # lock 1 held at the Δ = 0 extremum (preset's lock 1)
        r = 1e-4
        print(f"  x={x:.3f}: Δφ_in,ext={fold(in1 - in0) / r:+.5f} (½)  Δφ_out={fold(out1 - out0) / r:+.5f} (−½)"
              f"  θ_sq: ρ<1 {fold(th1 - th0) / r:+.5f}, ρ=1 {th1_lossless / r:+.5f}, 1 MHz {th1_1mhz / r:+.5f}"
              f"  [−1/(1+x²) = {-1 / (1 + x * x):+.5f}]")
        print(f"           lock-2 zero − output carrier = {fold(z - out1) / r:+.2e};  LO − θ_sq = {fold(z - th1) / r:+.5f}"
              f"  [(1−x²)/(2(1+x²)) = {(1 - x * x) / (2 * (1 + x * x)):+.5f}];"
              f"  lock 1 at the Δ=0 extremum: zero moves {fold(z_wrong - z) / r:+.5f}"
              f" [−(A−/A+)/2 = {-(1 - x) / (1 + x) / 2:+.5f}]")
    # empty cavity: carrier phase through the cavity
    kt, kl = kappas(0.0)
    k = kt + kl
    d = 1e-4 * k
    t = 1 / (k + 1j * d)
    print(f"  x = 0 (no gain): transmitted carrier phase / (Δ/κ) = {math.atan2(t.imag, t.real) / 1e-4:+.5f} (= −1: phase of"
          f" the probe THROUGH the cavity; −½ is the output frame relative to the pump-defined frame at the extremum)")


def part_b():
    print("\nB. 295 kHz dither through the cavity (Δ = 0): quasi-static A± (ERR blocks) vs complex A±(Ω_m)")
    f_m = 2.95e5
    for x in (0.236, 0.527, 0.745, 0.850, 0.913):
        kt, kl = kappas(x)
        k = kt + kl
        w = 2 * math.pi * f_m / k
        ap, am = 1 / (1 - x), 1 / (1 + x)
        apw, amw = 1 / (1 - x - 1j * w), 1 / (1 + x - 1j * w)
        # lock-2 slope ∝ Re A+(Ω_m) (b-quadrature), lock-1 slope ∝ A+·Re A+(Ω_m) − A−·Re A−(Ω_m)
        s2 = apw.real / ap
        s1 = (ap * apw.real - am * amw.real) / (ap * ap - am * am)
        lag = math.degrees(math.atan2(apw.imag, apw.real))
        # lock-1 offset δu → lock-2 zero ε: tan ε = Re A−(Ω)/Re A+(Ω) · tan δu
        du = 0.1
        eps_qs = math.atan(am / ap * math.tan(du))
        z = lock2_zero(x, 0.0, math.pi / 2 + du, f_m)
        eps_dyn = fold(z - math.pi / 2)
        print(f"  x={x:.3f}: Ω_m/κ={w:.4f}  lag of the amplified sideband {lag:5.1f}°  lock-2 slope ×{s2:.3f}"
              f"  lock-1 slope ×{s1:.3f}  ε(δu=0.1): quasi-static {eps_qs * 1e3:.2f} mrad, cavity {eps_dyn * 1e3:.2f} mrad")


def two_pole_psd(f, h, f_lo, f_c):
    """mech.two_pole shape: h/f² between f_lo and f_c, flat below f_lo, f⁻⁴ above f_c (T08 §2.5)."""
    return h / ((f * f + f_lo * f_lo) * (1 + (f / f_c) ** 2))


def part_c():
    print("\nC. dynamic: lock 1 blind to Δ (preset) vs lock 1 at the Δ-shifted extremum (physics); linear loops of refB_nominal")
    x = math.sqrt(0.1 / 0.18)
    kt, kl = kappas(x)
    k = kt + kl
    q, frac, tau = 10.0, 0.3, 1.0e-6
    f = np.logspace(0, 6, 60001)
    s = 2j * math.pi * f

    def loop(f_res):
        fu = frac * f_res / q
        p = 1 / (1 - (f / f_res) ** 2 + 1j * f / (f_res * q))
        return (2 * math.pi * fu / s) * p * np.exp(-s * tau)

    l1, l2 = loop(1.4e4), loop(2.2e4)
    s1, t2 = 1 / (1 + l1), l2 / (1 + l2)
    rr = (1 - x) / (1 + x)
    # metric per unit Δ [rad per rad/s]
    main = -t2 / (2 * k) + 1 / (k * (1 + x * x))       # preset (u = 0)
    extra = -t2 * rr * s1 / (2 * k)                    # physics adds this (u_true = −S1 Δ/(2κ))
    phys = main + extra
    i = np.argmax(np.abs(extra) / np.abs(main))
    print(f"  max |extra/main| = {np.abs(extra[i] / main[i]):.3f} at {f[i]:.0f} Hz; |T2·S1| max = {np.max(np.abs(t2 * s1)):.3f}")
    # OPO round-trip length noise of refB_nominal after the length lock (f_ugf = 1 kHz integrator)
    g = 2 * math.pi * NU0 / L_RT
    sd = g ** 2 * two_pole_psd(f, 1e-18, 1.0, 2e3) * np.abs(1j * f / (1j * f + 1e3)) ** 2
    var = lambda h: np.trapezoid(np.abs(h) ** 2 * sd, f)  # noqa: E731
    print(f"  Δ residual rms {math.sqrt(np.trapezoid(sd, f)) / k:.2e} κ;  metric rms: preset {math.degrees(math.sqrt(var(main))):.4f}°"
          f", physics {math.degrees(math.sqrt(var(phys))):.4f}°  (θ̃ = 1.5°)")
    # worst case of the axes: opo_h = 1e-16, f_ugf = 100 Hz
    sd2 = g ** 2 * two_pole_psd(f, 1e-16, 1.0, 2e3) * np.abs(1j * f / (1j * f + 1e2)) ** 2
    v2 = lambda h: np.trapezoid(np.abs(h) ** 2 * sd2, f)  # noqa: E731
    print(f"  axis corner (opo_h 1e-16, f_ugf 100 Hz): preset {math.degrees(math.sqrt(v2(main))):.3f}°,"
          f" physics {math.degrees(math.sqrt(v2(phys))):.3f}°")


def part_d():
    print("\nD. laser FM → OPO detuning (missing) vs the modelled geometric delay; refB_nominal loops, f_ugf = 1 kHz")
    x = math.sqrt(0.1 / 0.18)
    kt, kl = kappas(x)
    k = kt + kl
    q, frac, tau = 10.0, 0.3, 1.0e-6
    f_hi = k / (2 * math.pi)                            # flat-source range of the metric (conv §6)
    f = np.logspace(0, math.log10(f_hi), 40001)
    s = 2j * math.pi * f
    fu2, fr2 = frac * 2.2e4 / q, 2.2e4
    l2 = (2 * math.pi * fu2 / s) / (1 - (f / fr2) ** 2 + 1j * f / (fr2 * q)) * np.exp(-s * tau)
    t2 = l2 / (1 + l2)
    hp = 1j * f / (1j * f + 1e3)                        # cavity follows the laser below f_ugf
    per_delta = (-t2 / (2 * k) + 1 / (k * (1 + x * x))) * hp
    tau_geo = (2.5 - 2.0) / C                            # |sqz − LO| path mismatch of the preset
    for h0 in (10.0, 100.0, 1.0e3, 1.0e4):
        s_delta = (2 * math.pi) ** 2 * h0                # S_Δ = (2π)² S_ν
        v_det = np.trapezoid(np.abs(per_delta) ** 2 * s_delta, f)
        v_geo = np.trapezoid(4 * np.sin(math.pi * f * tau_geo) ** 2 * h0 / f ** 2, f)
        print(f"  h0 = {h0:7.0f} Hz²/Hz: via detuning {math.degrees(math.sqrt(v_det)):.3f}°, geometric delay"
              f" {math.degrees(math.sqrt(v_geo)):.4f}°; θ̃ = 1.5° → {math.sqrt(1.5 ** 2 + math.degrees(math.sqrt(v_det)) ** 2):.3f}°")
    tau_eff = 1 / (k * (1 + x * x))
    print(f"  equivalent extra delay of θ_sq: 1/(κ(1+x²)) = {tau_eff * 1e9:.1f} ns (vs the modelled OPO segment"
          f" n_g L/c = {0.5 / C * 1e9:.2f} ns)")
    for area in (1e4, 1e6, 1e8):                        # FM peaks: δν rms² = area, above every loop band
        d = 2 * math.pi * math.sqrt(area) / (k * (1 + x * x))
        print(f"  FM peak area {area:.0e} Hz²: θ_sq rotation rms {d * 1e3:.4f} mrad per peak")


def part_e():
    print("\nE. RIN at f_m as sensor noise (refB_nominal: P_LO 5.5 mW, P_pr 10 nW, P_pd 1 µW, m 0.1, f_u 420 / 660 Hz)")
    x = math.sqrt(0.1 / 0.18)
    ap, am = 1 / (1 - x), 1 / (1 + x)
    resp = 0.998 * 1.602176634e-19 * LAM / (6.62607015e-34 * C)
    p_lo, p_pr, p_pd, m = 5.5e-3, 1e-8, 1e-6, 0.1
    j1m, j0m = m / 2 - m ** 3 / 16, 1 - m ** 2 / 4
    j12m, j02m = m - m ** 3 / 2, 1 - m ** 2
    r_fm = 0.13                                          # R−(295 kHz) at 100 mW (r_hd_fm)
    k2 = 2 * resp * math.sqrt(p_lo * p_pr) * 2 * j1m * ap   # A per rad of Δ
    shot = 2 * 1.602176634e-19 * resp * p_lo * r_fm
    for rin_db in (-160.0, -140.0, -120.0, -110.0):
        s_rin = 10 ** (rin_db / 10)
        for eps in (0.01, 0.056):
            leak = eps ** 2 * (resp * p_lo) ** 2 * s_rin
            beat = (2 * resp * math.sqrt(p_lo * p_pr) * am * j0m) ** 2 * s_rin / 4 * 2   # probe + LO RIN on the carrier
            s_d = 2 * (shot + leak + beat) / k2 ** 2
            sig2 = math.sqrt(s_d * math.pi / 2 * 660.0)
            print(f"  lock 2: RIN {rin_db:6.1f} dB/Hz, ε = {eps:.3f}: (leak + beat)/shot = {(leak + beat) / shot:9.2e},"
                  f" residual {sig2 * 1e3:7.3f} mrad")
        gbar, dg = 0.5 * (ap * ap + am * am), 0.5 * (ap * ap - am * am)
        pdc = gbar - dg * j02m                           # DC power at the de-amplification extremum / P0
        k1 = 2 * dg * j12m                               # per rad of Φ1 / P0
        s_phi1 = 2 * s_rin * pdc ** 2 / k1 ** 2
        sig1 = math.sqrt(s_phi1 * math.pi / 2 * 420.0)
        print(f"  lock 1: RIN {rin_db:6.1f} dB/Hz: Φ1 residual {sig1 * 1e3:.4f} mrad → Δ ≤ (A−/A+)Φ1/2 = {am / ap * sig1 / 2 * 1e3:.5f} mrad")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
