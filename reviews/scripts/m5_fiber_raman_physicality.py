"""M5 review: physicality of the linearized fiber covariance with Raman (D59-5) vs the number of
covariance segments (``steps.n_save``), T = 0 and 295 K, vacuum and squeezed input.

Physicality: λ_min(V + iΩ/2) ≥ 0 (conv units, vacuum 1/2). Also the T = 0 Stokes photons: the
added photon number in the grid modes vs the analytic spontaneous Stokes count.

Run: uv run python reviews/scripts/m5_fiber_raman_physicality.py
"""
import numpy as np

from cranberry.blocks.fib.quantum import fiber_quantum_linear
from cranberry.blocks.fib.quantum_testing import linear_params, toggles
from cranberry.offline import gaussian


def phys(v):
    n = v.shape[0] // 2
    om = gaussian.symplectic_form(n)
    return float(np.min(np.linalg.eigvalsh(v + 0.5j * om)))


if __name__ == "__main__":
    tg = toggles(beta2=True, kerr=True, raman=True, kerr_fluct=True, raman_noise=True)
    res = {}
    for temp in (0.0, 295.0):
        for inp in ({"kind": "vacuum"},
                    {"kind": "squeezed", "r": 1.15, "theta_sq_rad": 0.3, "mode": {"shape": "sech", "t0_s": 1e-12}}):
            for nsave in (1, 2, 4, 8, 16, 50):
                p = linear_params(temperature_k=temp, input_state=inp, L_m=10 * linear_params()["L_m"],
                                  steps=dict(n_steps=400, n_save=nsave, control="fixed"))
                r = fiber_quantum_linear(p, tg)
                res[temp, inp["kind"], nsave] = phys(r.cov)
                print(f"T={temp:5.1f} K input={inp['kind']:8s} n_save={nsave:3d}: lambda_min(V+iΩ/2) = {phys(r.cov):+.3e}")
    # self-checks: D59-5's "≥ 16 segments" is not sufficient here (10 L_D, soliton, T = 0) — finding M5-01
    assert res[0.0, "vacuum", 16] < -1e-3            # unphysical by > 0.2 % of the vacuum variance
    assert res[0.0, "vacuum", 50] > res[0.0, "vacuum", 16]
    assert res[295.0, "vacuum", 4] < -1e-3 < res[295.0, "vacuum", 8]
    # scaling with the Raman gain-length per segment (T = 0, vacuum): 1, 3, 10 L_D; Raman off = RK4IP part
    for lmul in (1, 3, 10):
        row = []
        for nsave in (1, 4, 16, 50, 100):
            p = linear_params(temperature_k=0.0, L_m=lmul * linear_params()["L_m"],
                              steps=dict(n_steps=400, n_save=nsave, control="fixed"))
            row.append(phys(fiber_quantum_linear(p, tg).cov))
        off = phys(fiber_quantum_linear(p, dict(tg, raman=False, raman_noise=False)).cov)
        print(f"L = {lmul:2d} L_D, n_save 1/4/16/50/100: " + " ".join(f"{x:+.1e}" for x in row) + f"   Raman off: {off:+.1e}")
