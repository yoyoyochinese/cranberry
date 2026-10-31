"""M5 review: independent checks of fiber_gnlse (WP-B9) with every term on.

(a) step halving with Kerr + Raman + self-steepening + β2 + β3 on an N = 2 soliton (errors vs a
    reference at h/32; RK4IP global error ∝ h⁴ ⇒ ratio ≈ 16);
(b) sign checks in the conv §8 frame: Raman SSFS moves the spectral centroid to Ω < 0 (red);
    self-steepening (no dispersion) delays the intensity peak (T > 0) and extends the blue (Ω > 0) side;
(c) conservation: photon number ∫|Ã|²/(ω0+Ω) conserved with Raman + SS (loss off), energy not.

Run: uv run python reviews/scripts/m5_gnlse_checks.py
"""
import numpy as np

from cranberry.blocks.fib.gnlse import model_from_params, pulse_field, solve_pulsed, time_axis

B2, B3, GAM, LAM = -21.7e-27, 0.127e-39, 1.3e-3, 1545.32e-9
T0 = 100e-15


def params(**over):
    p = dict(L_m=0.0, alpha_per_m=0.0, beta2_s2_per_m=B2, beta3_s3_per_m=B3, beta_higher=[], gamma_per_w_m=GAM,
             lambda0_m=LAM, raman_model="blow_wood", f_R=0.18, xpm_b=0.0, xpm_power_w=0.0)
    p.update(over)
    return p


def run(tg, length, n_steps, n=2048, win=40 * T0, p0=None, t0=T0):
    t, dt = time_axis(n, win)
    p = params(L_m=length)
    m = model_from_params(p, tg, n, dt)
    p0 = p0 if p0 is not None else 4 * abs(B2) / (GAM * t0**2)  # N = 2
    a0 = pulse_field("sech", p0, t0, t)
    r = solve_pulsed(m, a0, t, length, n_steps=n_steps, n_save=1, control="fixed", rtol=None, monitor_edges=False)
    return r, m, t


ALL = dict(loss=False, beta2=True, beta3=True, beta_higher=False, kerr=True, raman=True, self_steepening=True,
           xpm_external=False, cw_shortcut=False)

if __name__ == "__main__":
    ld = T0**2 / abs(B2)
    L = 0.5 * np.pi * ld  # one soliton period
    ref = run(ALL, L, 64 * 32)[0].a_out
    errs = []
    for ns in (64, 128, 256, 512):
        a = run(ALL, L, ns)[0].a_out
        errs.append(np.linalg.norm(a - ref) / np.linalg.norm(ref))
    print("(a) N=2 soliton, all terms: rel. errors", ["%.3e" % e for e in errs],
          "ratios", ["%.1f" % (errs[i] / errs[i + 1]) for i in range(3)])
    assert 13 < errs[2] / errs[3] < 19  # 4th order
    r, m, t = run(ALL, 5 * ld, 2000, p0=abs(B2) / (GAM * T0**2))  # fundamental soliton, Raman on
    print(f"(b) SSFS: centroid after 5 L_D = {m.centroid(r.spectrum(-1)) / (2 * np.pi) / 1e12:+.4f} THz (expect < 0)")
    assert m.centroid(r.spectrum(-1)) < 0
    tg = dict(ALL, beta2=False, beta3=False, raman=False)
    r, m, t = run(tg, 200.0, 4000, p0=50.0, t0=1e-12, win=40e-12)  # φ_NL ≈ 13 rad
    ia = np.abs(r.a_out) ** 2
    spec = np.abs(np.fft.fftshift(r.spectrum(-1))) ** 2
    om = np.fft.fftshift(m.omega)
    w = spec / spec.sum()
    print(f"    SS only: peak time {t[np.argmax(ia)] * 1e15:+.2f} fs (expect > 0); spectral skewness "
          f"{np.sum(w * om**3) / np.sum(w * om**2) ** 1.5:+.3f} (expect > 0, blue tail)")
    assert t[np.argmax(ia)] > 0 and np.sum(w * om**3) > 0
    r, m, t = run(ALL, 2 * ld, 1000)
    print(f"(c) photon number rel. change {r.photon_number[1] / r.photon_number[0] - 1:+.2e},"
          f" energy rel. change {r.energy_j[1] / r.energy_j[0] - 1:+.2e}")
    assert abs(r.photon_number[1] / r.photon_number[0] - 1) < 1e-7 and abs(r.energy_j[1] / r.energy_j[0] - 1) > 1e-3
