"""Phase 5 review (M1, Ref A): independent check of D58-4 (ESA log-average trace scatter).

Claim: Var(log trace) / Var(power trace, relative) = ∫Li2(|rho|^2)dτ / ∫|rho|^2 dτ  (not pi^2/6), and
B_stat = (∫|H|^2)^2/∫|H|^4 = sqrt(2)·ENBW for a Gaussian |H|^2.

Parts: (1) Li2 identity Cov(ln X, ln Y) = Li2(|rho|^2) by Monte Carlo of correlated complex Gaussians;
(2) ratio and B_stat/ENBW for rect, Gaussian and 4-pole synchronously tuned |H|^2 (quadrature);
(3) time-domain MC: RBW filter -> |z|^2 -> {power, 10log10} -> {boxcar, one-pole VBW}, scatter vs formula;
(4) Ref A numbers (RBW 1 MHz, VBW 100 Hz).
Run: uv run python reviews/scripts/d58_4_log_scatter.py
"""

import math

import numpy as np
from scipy import signal
from scipy.special import spence

DB = 10 / math.log(10)
li2 = lambda x: spence(1 - x)  # noqa: E731  (scipy spence(z) = Li2(1 - z))
rng = np.random.default_rng(20261031)

# (1) Li2 identity --------------------------------------------------------------------------------
print("(1) Cov(ln X, ln Y) vs Li2(|rho|^2), X,Y = |z|^2 of correlated unit complex Gaussians")
n = 4_000_000
for r in (0.3, 0.6, 0.9, 1.0):
    a = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
    b = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
    z2 = r * a + math.sqrt(1 - r * r) * b
    c = np.cov(np.log(np.abs(a) ** 2), np.log(np.abs(z2) ** 2))[0, 1]
    print(f"  |rho|={r:.1f}: MC {c:.4f}   Li2 {li2(r * r):.4f}   (pi^2/6 = {math.pi**2/6:.4f})")

# (2) ratio and B_stat for three RBW shapes (|H|^2 normalized to ENBW = 1) -----------------------
a4 = 16 / (5 * math.pi)  # 4-pole sync-tuned: ∫(1+(f/a)^2)^-4 df = 5πa/16 = 1
shapes = {
    "rect": lambda f: (np.abs(f) <= 0.5).astype(float),
    "gaussian": lambda f: np.exp(-math.pi * f * f),
    "4pole": lambda f: (1 + (f / a4) ** 2) ** -4,
}
print("\n(2) per RBW shape (ENBW = 1): B_stat/ENBW, ratio = ∫Li2(|rho|^2)/∫|rho|^2, -3 dB width/ENBW")
const = {}
fgrid = np.linspace(-40, 40, 2**20 + 1)
dfg = fgrid[1] - fgrid[0]
for name, h2 in shapes.items():
    s = h2(fgrid)
    bstat = (s.sum() * dfg) ** 2 / ((s * s).sum() * dfg)
    if name == "rect":
        tau = np.linspace(-4000, 4000, 8_000_001)
        rho2 = np.sinc(tau) ** 2
    else:  # rho(τ) = FT of |H|^2 / ∫|H|^2 (|H|^2 even -> cosine transform)
        tau = np.linspace(-30, 30, 2401)
        rho = np.array([np.sum(s * np.cos(2 * math.pi * fgrid * t)) for t in tau]) / s.sum()
        rho2 = rho**2
    ratio = li2(np.clip(rho2, 0, 1)).sum() / rho2.sum()
    fs_ = np.linspace(0, 2, 200001)
    w3 = 2 * fs_[np.argmax(h2(fs_) <= 0.5)]
    const[name] = (bstat, ratio)
    print(f"  {name:8s}: B_stat/ENBW = {bstat:.4f}   ratio = {ratio:.4f}   -3dB width/ENBW = {w3:.4f}")
zeta25 = sum(k**-2.5 for k in range(1, 200000)) + 2 / 3 * 200000**-1.5
print(f"  Gaussian closed form: ratio = zeta(5/2) = {zeta25:.4f}, B_stat/ENBW = sqrt(2) = {math.sqrt(2):.4f}")

# (3) time-domain MC ------------------------------------------------------------------------------
fs, vbw, N = 4.0, 2e-3, 2**22
print(f"\n(3) time-domain MC (ENBW = 1, fs = {fs}, VBW = {vbw}): measured var / [pi·VBW/B_stat]")
f = np.fft.fftfreq(N, 1 / fs)
al = -math.expm1(-2 * math.pi * vbw / fs)
m = int(fs / (math.pi * vbw))  # boxcar length with the same noise bandwidth as the one-pole ((π/2)VBW)
for name, h2 in shapes.items():
    bstat, ratio = const[name]
    h = np.sqrt(h2(f))
    h /= math.sqrt(np.mean(h * h))
    acc = {k: [] for k in ("pow_1pole", "log_1pole", "pow_boxcar", "log_boxcar")}
    for _ in range(6):
        w = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / math.sqrt(2)
        p = np.abs(np.fft.ifft(np.fft.fft(w) * h)) ** 2
        for tag, x in (("pow", p / p.mean()), ("log", np.log(p))):
            acc[tag + "_1pole"].append(signal.lfilter([al], [1, al - 1], x - x.mean())[20000:].var())
            nb = x.size // m
            acc[tag + "_boxcar"].append(x[: nb * m].reshape(nb, m).mean(1).var())
    pred = math.pi * vbw / bstat
    for k, v in acc.items():
        exp_ = ratio if k.startswith("log") else 1.0
        print(f"  {name:8s} {k:10s}: {np.mean(v) / pred:.3f} ± {np.std(v) / math.sqrt(len(v)) / pred:.3f}"
              f"   expect {exp_:.3f}")

# (4) Ref A numbers -------------------------------------------------------------------------------
print("\n(4) Ref A trace scatter [dB], RBW 1 MHz, VBW 100 Hz")
for name in shapes:
    bstat, ratio = const[name]
    for k in (1.0, 1.0645, 1.128):
        rel2 = math.pi * 100 / (bstat * k * 1e6)
        print(f"  {name:8s} ENBW={k:.4f}·RBW: power {DB*math.sqrt(rel2):.4f}   log(true) {DB*math.sqrt(rel2*ratio):.4f}"
              f"   log(pi^2/6) {DB*math.sqrt(rel2*math.pi**2/6):.4f}")
