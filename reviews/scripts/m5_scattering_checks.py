"""M5 review: scattering formulas of WP-B10 against T07 / D48 / D59-8.

(a) Blow–Wood first moment ∫t h_R dt by direct quadrature vs ``raman_first_moment_s`` and T07 §3.1's
    "⟨T⟩_h = 8.7 fs, g_R′ = 4.1e−18"; g_R peak (T07: 0.68 /(W km) at 13.1 THz).
(b) D48 §7-6: thermal part → 0 as T → 0, Stokes remains; detailed balance aS/S = e^{−ħΩ/kT}.
(c) GAWBS R0,m peak frequencies (125 µm and 80 µm) vs T07 §4.1 table; R0,7 vs 321.75 MHz (London);
    R0,7 σ² per metre vs T07 (1.33e−10 rad²/m at 295 K).
(d) Phase-noise ↔ photon consistency: P·L·S_R/hν0 equals the mean of the Stokes and anti-Stokes
    photon densities (symmetric order).

Run: uv run python reviews/scripts/m5_scattering_checks.py
"""
import numpy as np
from scipy.integrate import quad

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.gnlse import raman_first_moment_s

T1, T2, GAM, FR, LAM = 12.2e-15, 32e-15, 1.3e-3, 0.18, 1545.32e-9

if __name__ == "__main__":
    h = lambda t: (T1**2 + T2**2) / (T1 * T2**2) * np.exp(-t / T2) * np.sin(t / T1)  # noqa: E731
    m1 = quad(lambda x: x * 1e-15 * h(x * 1e-15) * 1e-15, 0, 2000, limit=400, epsabs=0, epsrel=1e-12)[0]  # fs units
    print(f"(a) ∫t h dt = {m1 * 1e15:.4f} fs (code {raman_first_moment_s('blow_wood') * 1e15:.4f} fs; T07 8.7 fs);"
          f" g_R' = {2 * GAM * FR * m1:.3e} /(W m)/(rad/s) (T07 4.1e-18)")
    f = np.linspace(0.1e12, 40e12, 40000)
    g = sc.raman_gain(2 * np.pi * f, gamma_per_w_m=GAM, f_R=FR, raman_model="blow_wood")
    print(f"    g_R peak {g.max() * 1e3:.3f} /(W km) at {f[np.argmax(g)] / 1e12:.2f} THz (T07 0.68 at 13.1)")
    w = 2 * np.pi * np.array([0.0, 1e9, 1e12, 13e12])
    kw = dict(gamma_per_w_m=GAM, f_R=FR, raman_model="blow_wood")
    for temp in (295.0, 1e-3, 0.0):
        th = sc.raman_scatter_density(w, temperature_k=temp, thermal_only=True, **kw)
        st = sc.raman_scatter_density(w, temperature_k=temp, **kw)
        an = sc.raman_scatter_density(-w, temperature_k=temp, **kw)
        print(f"(b) T={temp:g} K: thermal {np.array2string(th, precision=3)}  Stokes {np.array2string(st, precision=3)}"
              f"  anti-Stokes {np.array2string(an, precision=3)}")
    ws = 2 * np.pi * 1e12
    s, a = (sc.raman_scatter_density(x, temperature_k=295.0, **kw) for x in (ws, -ws))
    print(f"    detailed balance aS/S = {a / s:.6f} vs exp(−ħΩ/kT) = {np.exp(-sc.HBAR * ws / (sc.K_B * 295.0)):.6f}")
    for d in (125e-6, 80e-6):
        print(f"(c) d={d * 1e6:.0f} µm R0,m [MHz]:", np.round(sc.gawbs_peak_frequencies(d, 5996.0, 3740.0, 12) / 1e6, 1))
    g7 = dict(cladding_diameter_m=125e-6, v_l_m_s=5996.0, v_s_m_s=3740.0, n_modes=7,
              gamma_om_per_w_m=[0.0] * 6 + [2e-3], linewidth_hz=4.3e6, mode_sum_factor=1.0)
    var = sc.gawbs_variances(g7, temperature_k=295.0, lambda0_m=LAM, length_m=1.0)
    print(f"    R0,7 σ² per metre = {var[-1]:.3e} rad² (T07 1.33e-10)")
    for wd in (2 * np.pi * 1e9, 2 * np.pi * 5e12):
        sr = sc.raman_phase_psd(wd, temperature_k=295.0, lambda0_m=LAM, **kw)
        mean = 0.5 * (sc.raman_scatter_density(wd, temperature_k=295.0, **kw)
                      + sc.raman_scatter_density(-wd, temperature_k=295.0, **kw))
        print(f"(d) Ω/2π={wd / 2 / np.pi:.0e} Hz: S_R/hν0 = {sr / sc._hnu(LAM):.6e}, mean(S, aS) = {mean:.6e}")
        assert abs(sr / sc._hnu(LAM) / mean - 1) < 1e-12
    # self-checks
    assert abs(m1 - raman_first_moment_s("blow_wood")) < 1e-6 * m1 and abs(m1 * 1e15 - 8.7) > 0.5  # code right; T07 "8.7 fs" wrong
    assert abs(g.max() * 1e3 - 0.68) < 0.01
    assert abs(sc.gawbs_peak_frequencies(125e-6, 5996.0, 3740.0, 7)[-1] / 321.75e6 - 1) < 0.006
    assert abs(var[-1] / 1.33e-10 - 1) < 0.01
