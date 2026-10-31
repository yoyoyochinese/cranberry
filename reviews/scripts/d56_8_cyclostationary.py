"""D56-8 adjudication: dither-lock sensor noise is cyclostationary (Phase 5 M1 review).

Independent of the cranberry code (numpy + scipy.special only).

Model (Ref C quantum-noise lock, T14 §6): BHD output x(t) Gaussian, measured quadrature
Δ(t) = Δ0 + m sin θ, θ = 2π f_m t (θ is the DITHER phase). Shot-noise units: variance density
R- cos²Δ + R+ sin²Δ + el. Square-law radiometer y = x², demodulated with 2 sin θ (the error signal
sits in the sin θ component), low-passed. One-sided PSDs throughout (conv §7).

Claim (D56-8): S_e = 2⟨2 sin²θ S_v(Δ0 + m sinθ)⟩ instead of 2 S_v(Δ0).
Part 1: Monte Carlo of the actual square-law radiometer (not of the white-modulated model).
Part 2: Ref C numbers (flat band B = 292 MHz, b_loop = 10 kHz), square law and log law.
Part 3: RAM offset, with and without the electronic (dark) noise being modulated by LO RAM.
"""
import numpy as np
from scipy.special import j1

db = lambda x: 10 ** (x / 10)  # noqa: E731
RP0, RM0, EL, B, M, BLOOP = db(10.0), db(-4.7), db(-11.0), 292e6, 0.0785, 10e3
RP, RM = RP0 + EL, RM0 + EL  # incl. dark noise (as T14 §6.3)
TH = np.linspace(0, 2 * np.pi, 1 << 14, endpoint=False)
S = np.sin(TH)


def P(d, rp=RP, rm=RM):
    """Band power per unit bandwidth: P̄ − ΔP cos 2Δ."""
    return (rp + rm) / 2 - (rp - rm) / 2 * np.cos(2 * d)


# ---------------------------------------------------------------- part 1: radiometer Monte Carlo
def mc(d0, m, rp, rm, el, nper=64, ncyc_blk=64, nblk=4000, seed=1):
    """x = √rm g1 cosΔ + √rp g2 sinΔ + √el g3, white per sample (band fs/2 ≫ f_m = fs/nper), y = x².
    Returns the measured one-sided PSD of the demodulated noise and both predictions (fs = 1;
    a white sequence of variance σ² has one-sided PSD 2σ²/fs; Var(mean over T) = S/(2T))."""
    rng = np.random.default_rng(seed)
    th = 2 * np.pi * np.arange(nper * ncyc_blk) / nper
    d = d0 + m * np.sin(th)
    ref = 2 * np.sin(th)
    means = np.empty(nblk)
    for k in range(nblk):
        g = rng.standard_normal((3, th.size))
        x = np.sqrt(rm) * g[0] * np.cos(d) + np.sqrt(rp) * g[1] * np.sin(d) + np.sqrt(el) * g[2]
        means[k] = np.mean(ref * x * x)
    s_meas = 2 * means.var() * th.size

    def sv(dd):  # Var(x²) = 2V²  ->  one-sided PSD 4V²
        return 4 * (rm * np.cos(dd) ** 2 + rp * np.sin(dd) ** 2 + el) ** 2

    return s_meas, 2 * sv(d0), 2 * np.mean(2 * S**2 * sv(d0 + m * S))


print("Part 1: radiometer Monte Carlo (square law), PSD of demodulated noise; 3σ ≈ %.3f" % (3 * np.sqrt(2 / 4000)))
for d0, m, lab in [(0.0, M, "Ref C m, squeezing pt"), (0.0, 0.3, "m = 0.3 (exaggerated)"),
                   (np.pi / 2, M, "anti-squeezing pt")]:
    sm, ss, sc = mc(d0, m, RP0, RM0, EL)
    print(f"  {lab:24s} meas {sm:9.4f}  stationary {ss:9.4f} (meas/stat {sm/ss:5.3f})  "
          f"cyclo {sc:9.4f} (meas/cyclo {sm/sc:5.3f})")

# ---------------------------------------------------------------- part 2: Ref C residual
K = 4 * (RP - RM) / 2 * j1(2 * M)  # dI/dΔ per unit bandwidth (square law, Δ0 = 0)
w_sq = np.mean(2 * S**2 * P(M * S) ** 2) / P(0.0) ** 2
S_t14 = RM**2 / (B * (RP - RM) ** 2 * j1(2 * M) ** 2)  # T14 §6.3 flat-band form
print("\nPart 2: Ref C sensor-limited residual (continuous first-order loop, b_loop = 10 kHz)")
print(f"  square law: S_Δ(T14) = {S_t14:.4e} rad²/Hz -> {1e3*np.sqrt(S_t14*BLOOP):.3f} mrad;"
      f" weight {w_sq:.4f} -> {1e3*np.sqrt(w_sq*S_t14*BLOOP):.3f} mrad")
eps = 2 * np.pi * (2 * BLOOP / np.pi) / 1e6
print(f"  discrete AR(1) factor of the VT-43 loop (f_sim 1 MHz): 1/(1-eps/2) = {1/(1-eps/2):.4f} ->"
      f" {1e3*np.sqrt(w_sq*S_t14*BLOOP/(1-eps/2)):.3f} mrad (= the '3.62')")
# log law: S_v ∝ κ·2Q/P² = κ·2/B (flat band) -> noise weight exactly 1; slope = c⟨2 sinθ ∂_Δ ln P⟩
k_stat = 4 * (RP - RM) / 2 * j1(2 * M) / P(0.0)  # T14 §6.2 log slope (per c = 10k/ln10)
k_exact = np.mean(2 * S * (RP - RM) * np.sin(2 * M * S) / P(M * S))
S_log = 1.32 * (2 * 2 / B) / k_exact**2
print(f"  log law (κ 1.32): T14 {1e3*np.sqrt(1.32*S_t14*BLOOP):.3f} mrad; exact/stationary slope "
      f"{k_exact/k_stat:.4f}, noise weight 1 -> {1e3*np.sqrt(S_log*BLOOP):.3f} mrad "
      f"(PSD factor {(k_stat/k_exact)**2:.4f})")

# ---------------------------------------------------------------- part 3: RAM offset (a = 1e-2, φ_a = 0)
a = 1e-2
t14 = a * RM / (2 * (RP - RM) * j1(2 * M))
ram_all = a * np.mean(2 * S**2 * P(M * S)) / K  # el modulated too (code, D56-8)
ram_opt = a * np.mean(2 * S**2 * P(M * S, RP0, RM0)) / K  # only LO-proportional noise modulated
print("\nPart 3: RAM offset a = 1e-2")
print(f"  square law: T14 form {1e3*t14:.3f}; exact, el modulated (code) {1e3*ram_all:.3f}; "
      f"exact, el NOT modulated {1e3*ram_opt:.3f} mrad  (cycle weight {np.mean(2*S**2*P(M*S))/P(0.0):.4f})")
lg_all = a / k_exact
lg_opt = a * np.mean(2 * S**2 * P(M * S, RP0, RM0) / P(M * S)) / k_exact
print(f"  log law:    stationary {1e3*a/k_stat:.3f}; exact, el modulated {1e3*lg_all:.3f}; "
      f"el NOT modulated {1e3*lg_opt:.3f} mrad")

# self-check: the Monte Carlo follows the cyclostationary formula, not the stationary one
sm, ss, sc = mc(0.0, 0.3, RP0, RM0, EL, seed=7)
assert abs(sm / sc - 1) < 0.07 and abs(sm / ss - 1) > 0.3
print("\nself-check OK")
