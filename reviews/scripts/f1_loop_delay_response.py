"""Phase 5 M1 framework review: phase / magnitude accuracy of the delay_buffer (linear vs thiran1)
at loop-relevant frequencies, the Thiran small-d regime, and the effective delay of a decimated
(rate_div = r) hold on the f_sim grid.

Run: uv run python reviews/scripts/f1_loop_delay_response.py
"""

import numpy as np

from cranberry.core.delay import DelayLine, design_delay

FS = 1.0e6


def report(tau_s, fs=FS, freqs=(1e4, 6e4, 1e5, 2.5e5)):
    f = np.array(freqs)
    print(f"tau_s = {tau_s*1e6:.3f} us at f_sim = {fs/1e6:g} MHz (D = {tau_s*fs:.3f} samples)")
    for m in ("linear", "thiran1"):
        d = design_delay(tau_s, fs, m)
        h = d.freq_response(f, fs)
        ex = np.exp(-2j * np.pi * f * tau_s)
        dph = np.degrees(np.angle(h / ex))
        print(f"  {m:8s} k={d.k} mu={d.mu:.3f} a={d.a:+.3f} | phase error [deg] "
              + " ".join(f"{x:+7.3f}" for x in dph) + " | |H| " + " ".join(f"{abs(x):.4f}" for x in h))
    print("  at f [kHz]: " + " ".join(f"{x/1e3:g}" for x in f))


def thiran_ringing(D):
    line = DelayLine(D / FS, FS, "thiran1")
    x = np.zeros(4000)
    x[0] = 1.0
    y = line.process(x)
    n99 = int(np.flatnonzero(np.cumsum(y**2) >= 0.99 * np.sum(y**2))[0])
    print(f"thiran1 D = {D}: a = {line.design.a:.4f}, 99 % of impulse energy after {n99} samples, "
          f"y[1..4] = {np.round(y[1:5], 4)}")


def hold_delay(r, f=(1e3, 1e4, 1e5)):
    """y[n] = x[r*floor(n/r)] (decimated block, held). LTI part of the time-varying system:
    H(f) = (1/r) Σ_{j<r} e^{-iωj} (average over phases) → group delay (r-1)/2 samples."""
    f = np.array(f)
    w = 2 * np.pi * f / FS
    h = np.mean([np.exp(-1j * w * j) for j in range(r)], axis=0)
    gd = -np.angle(h) / w
    print(f"hold rate_div={r}: effective delay {np.round(gd, 4)} samples at f = {f/1e3} kHz "
          f"(continuous ZOH of period r/f_sim would be r/2 = {r/2} samples)")


if __name__ == "__main__":
    report(5.67e-6 - 1e-6 - 1.5e-6)   # refD delay buffer (D61-2)
    report(0.3e-6)                    # small fractional delay (< 0.5 sample)
    report(0.3e-6, fs=20e6)
    for D in (0.3, 0.05, 0.01):
        thiran_ringing(D)
    for r in (2, 4):
        hold_delay(r)
