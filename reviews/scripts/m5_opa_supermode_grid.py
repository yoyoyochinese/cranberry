"""M5 review: grid / n_z convergence of opa_supermodes (T03 §5.2 notes K depends on the grid for
100 fs pumps). VT-32 fixture (1 mm, α = 1 /W, 0.25 W, GVM 0.3 ps/mm) at τ_p = 1 ps and 100 fs.

Run: uv run python reviews/scripts/m5_opa_supermode_grid.py
"""
import numpy as np

from cranberry.blocks.source.opa_supermodes import make_opa_supermodes

TOG = {"multimode": True, "gvm": True, "gvd": True}
OPA = dict(length_m=1e-3, alpha_per_w=1.0, peak_pump_w=0.25, pump_fwhm_s=1e-12, gvm_s_per_m=3e-10,
           beta2_signal_s2_per_m=1e-25, beta2_pump_s2_per_m=4e-25, n_points=128, time_window_s=8e-12, n_z=20,
           phi_opa0_rad=0.0, segment="xtal", psi_p_xtal_rad=0.0, toggles=TOG)


def run(**kw):
    b = make_opa_supermodes({**OPA, **kw})
    r = b.supermode_r
    u = b.leading_mode
    spec = np.abs(np.fft.fft(u)) ** 2
    edge = float(np.sum(np.fft.fftshift(spec)[: u.size // 16]) + np.sum(np.fft.fftshift(spec)[-u.size // 16:])) / spec.sum()
    return r[0], b.schmidt_number, edge


if __name__ == "__main__":
    ks = {}
    for tp in (1e-12, 100e-15):
        print(f"tau_p = {tp * 1e15:.0f} fs (window W, N points, n_z): r0, K, leading-mode spectral energy in outer 1/8")
        for w, n, nz in ((8e-12, 128, 20), (8e-12, 256, 20), (16e-12, 256, 20), (8e-12, 512, 20), (8e-12, 256, 40)):
            r0, k, e = run(pump_fwhm_s=tp, time_window_s=w, n_points=n, n_z=nz)
            print(f"  W={w * 1e12:4.0f} ps N={n:4d} n_z={nz:3d}: r0={r0:.5f}  K={k:8.3f}  edge={e:.2e}")
            ks[tp, w, n, nz] = (r0, k)
    # finding M5-07: K doubles with the spectral span N/W while r0 is converged
    assert ks[1e-12, 8e-12, 512, 20][1] / ks[1e-12, 8e-12, 256, 20][1] > 1.7
    assert abs(ks[1e-12, 8e-12, 512, 20][0] / ks[1e-12, 8e-12, 256, 20][0] - 1) < 1e-3
