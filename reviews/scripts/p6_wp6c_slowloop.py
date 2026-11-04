"""Phase 5 review of WP6-C (D92-3/-4): independent model of the slow loop (window average → hold → continuous
integrator ramp, the block chain pcf1 → dbs1.pid) — exact discrete stability, and the component-2 residual on an f⁻²
common-point disturbance, without the simulator.

  uv run python reviews/scripts/p6_wp6c_slowloop.py           (≈ 20 s)

1  Stability: with a = 2π f_c T_int, consecutive window errors obey e_{n+2} = (1 − a/2) e_{n+1} − (a/2) e_n + Δd̄,
   stable for 0 < a < 2 (f_c T < 1/π = 0.318); chk.dbs_stability's L = (2π f_c/s) e^{−s(T + τ_L)} puts the limit at
   f_c T = 0.25 (conservative). Checked by iterating the time-domain model.
2  Residual: OU (Lorentzian) disturbance S = h/(f² + f_lo²), f_lo = 0.1 Hz (env.mech.dist), 30 s records, mean
   removed, 400 realisations: rms ratio closed/open loop for f_c T = 0.05, 0.1, 0.15, 0.2. Multiplying by the
   author's open-loop 59.9 mrad gives the expected closed-loop component 2 (D92-4: 18.4, analytic 18.7).
"""

from __future__ import annotations

import math

import numpy as np

DT = 1.0e-3           # 1 kHz model rate (f_c ≈ 2 Hz)
T_INT = 0.05
F_LO = 0.1
DUR = 30.0


def ou(rng, n_real, n):
    a = math.exp(-2 * math.pi * F_LO * DT)
    w = rng.standard_normal((n_real, n)) * math.sqrt(1 - a * a)
    x = np.empty((n_real, n))
    x[:, 0] = rng.standard_normal(n_real)          # stationary start
    for i in range(1, n):
        x[:, i] = a * x[:, i - 1] + w[:, i]
    return x


def closed(d, fct):
    """Error X(t) = d(t) − y(t); y ramps during window n+1 with slope ki·ē_n (ē_n: average of X over window n)."""
    m = round(T_INT / DT)
    a = 2 * math.pi * fct                           # k·ki·T
    slope = a / T_INT                               # per second, per rad of error
    n_real, n = d.shape
    x = np.empty_like(d)
    y = np.zeros(n_real)
    e_prev = np.zeros(n_real)
    for w0 in range(0, n, m):
        idx = np.arange(w0, min(w0 + m, n))
        t = (idx - w0 + 1) * DT
        yy = y[:, None] + slope * e_prev[:, None] * t[None, :]
        x[:, idx] = d[:, idx] - yy
        y = yy[:, -1]
        e_prev = x[:, idx].mean(axis=1)
    return x


def main():
    rng = np.random.default_rng(1)
    for fct in (0.25, 0.30, 0.31, 0.33):
        d = np.zeros((1, round(20.0 / DT)))
        d[0, :] = 1.0                               # step
        x = closed(d, fct)
        tail = np.abs(x[0, -round(2.0 / DT):]).max()
        print(f"1  f_c·T = {fct:.2f}: |error| over the last 2 s of a step response {tail:.2e}"
              f"  ({'stable' if tail < 1e-3 else 'not converged / unstable'}; discrete limit 1/π = 0.318)")
    d = ou(rng, 400, round(DUR / DT))
    base = d - d.mean(axis=1, keepdims=True)
    vb = np.mean(base ** 2)
    for fct in (0.05, 0.1, 0.15, 0.2):
        x = closed(d, fct)
        x = x - x.mean(axis=1, keepdims=True)
        r = math.sqrt(np.mean(x ** 2) / vb)
        print(f"2  f_c·T = {fct:.2f}: rms ratio {r:.4f} → 59.9 mrad × ratio = {59.9 * r:.1f} mrad")


if __name__ == "__main__":
    main()
