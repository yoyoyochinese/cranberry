"""Phase 5 M1 review (Ref D): linear stability of the D53 gains over the refD 915B-resonance sweep grid.

Continuous L(s) = K_0·C(s)·H_in·H_out·e^{−sτ}·A_914·H_915B with C = K_p + K_i/s + K_d s (D53 via DSC1 p.20),
τ = 5.7 µs, 120/100 kHz poles, 914 poles 29.5/112 kHz, H_915B = ω²/(s² + sω/Q + ω²). Stability by the
Nyquist criterion for an open loop with one integrator and no RHP poles: the closed loop is stable iff
the winding number of 1 + L along s = iω (ω: 0⁺ → ∞, mirrored) around 0 is zero. Reported with the
largest |L| at a −180° crossing (gain margin). V_0 = 2.90 V (refD_nominal, 2 × 0.3 dB coupler excess)
and 3.33 V (T_opt = 1).

Run: uv run python reviews/scripts/m1_refD_stability_grid.py
"""

import math

import numpy as np

K_ACT = 2 * math.pi * 16e-6 / (1.5488e-6 * 10.0)
KP, KI, KD = 2.41e-4, 1.406e-3 / 2e-6, 1.108e-2 * 2e-6


def loop(f, v0, fr, q):
    s = 2j * np.pi * f
    c = KP + KI / s + KD * s
    p = np.exp(-s * 5.7e-6) / ((1 + s / (2 * np.pi * 1.2e5)) * (1 + s / (2 * np.pi * 1.0e5)))
    p = p / ((1 + s / (2 * np.pi * 2.95e4)) * (1 + s / (2 * np.pi * 1.12e5)))
    w = 2 * np.pi * fr
    return v0 * K_ACT * c * p * w * w / (s * s + s * w / q + w * w)


def stable(v0, fr, q):
    f = np.geomspace(1e-2, 1e7, 400001)
    g = 1 + loop(f, v0, fr, q)
    # positive-frequency half; the integrator closes the contour through the right half plane at f → 0
    wind = (np.unwrap(np.angle(g))[-1] - np.unwrap(np.angle(g))[0]) / (2 * np.pi)
    lg = loop(f, v0, fr, q)
    ph = np.unwrap(np.angle(lg))
    cross = np.where(np.diff(np.floor((ph + np.pi) / (2 * np.pi))) != 0)[0]
    gm = max((abs(lg[i]) for i in cross), default=0.0)
    # for a type-1 loop with a stable remainder, 1 + L goes from ∞·(−i) at 0⁺ to 1 at ∞: stable ⇔ Δarg = +π/2
    return abs(wind - 0.25) < 0.05, gm


def main():
    frs = np.geomspace(5e3, 6e4, 4)
    qs = np.geomspace(3.0, 100.0, 4)
    for v0 in (2.90, 3.33):
        print(f"V0 = {v0} V  (rows f_res [kHz], cols Q = {np.round(qs, 2)}): S = stable, U = unstable, (max |L| at −180°)")
        for fr in frs:
            row = []
            for q in qs:
                ok, gm = stable(v0, fr, q)
                row.append(f"{'S' if ok else 'U'}({gm:5.2f})")
            print(f"  {fr / 1e3:6.2f}  " + "  ".join(row))
    ok, gm = stable(2.90, 2.0e4, 10.0)
    print(f"refD_nominal (20 kHz, Q 10, V0 2.90): {'stable' if ok else 'UNSTABLE'}, max |L| at −180° = {gm:.3f} "
          f"(GM {-20 * math.log10(gm):.1f} dB)")


if __name__ == "__main__":
    main()
