"""M4 review (WP-B8 probe_pulser after the D72 / M5-02 fix): per-port gated stray light, both units.

  1  duty → 1 (a ≡ 1) reproduces B10's CW ``stray_photons`` in both columns (gate and mode), for a ps
     Gaussian mode and an ns exponential (OPO-like) mode, Rayleigh + Brillouin + Raman + reflections.
  2  ports k ≤ 15 vs 16 (100 MHz, 32 ports, T_f = 320 ns, τ_pr = 10 ns, T_gate = 1 ns): there is an offset
     with zero distributed photons in the gate AND in the mode for k = 15, none for k = 16, and for
     k = 16 the offset average of the mode column equals the CW value (narrow line).
  3  independent brute-force: random scatterer positions, mode overlap |∫u a(t − τ_i) e^{i2πf0 t}|² summed
     in power, vs ``mode_factors`` (distributed average).

Run: uv run python reviews/scripts/m4_stray_gating_modes.py
"""
import numpy as np

import cranberry.blocks.dmx.testing as dt
from cranberry.blocks.fib.scattering import mode_time, stray_photons
from cranberry.blocks.mod.probe_pulser import gated_stray_photons, mode_factors
from cranberry.path.topology import C_M_PER_S

N_G = dt.N_G
ALL = {"stray_rayleigh": True, "stray_brillouin": True, "stray_raman": True, "stray_reflections": True}
RAY = {"stray_rayleigh": True, "stray_brillouin": False, "stray_raman": False, "stray_reflections": False}


def L_of(k, f_rep=100e6):
    return k * C_M_PER_S / (f_rep * N_G)


if __name__ == "__main__":
    print("1. duty = 1 ⇒ CW (both columns)")
    for mode in ({"shape": "gaussian", "t_fwhm_s": 1e-12}, {"shape": "exponential", "tau_eff_s": 14.5e-9}):
        p = dt.stray_params(L_of(20), mode=mode, reflections=[[1e-6, 1.0], [3e-7, 1.0]], t_gate_s=25e-9,
                            gbs_filter={"shape": "flat", "fwhm_hz": 40e9})
        cw = stray_photons(p, ALL, 1e-6)
        g = gated_stray_photons(p, ALL, p_avg_w=1e-6, duty=1.0, period_s=320e-9, t_offset_s=3e-9, eps_aom=0.0,
                                gating=True, n_g=N_G, reflection_z_m=[0.5, 7.0])
        for col in ("clicks_per_gate", "rayleigh_mode", "raman_mode", "reflections_mode"):
            assert abs(g[col] / cw[col] - 1) < 2e-3, (mode, col, g[col], cw[col])
        rb = g["brillouin_mode"] / cw["brillouin_mode"]
        print(f"  {mode['shape']:11s}: brillouin_mode gated(a ≡ 1)/CW = {rb:.4f}, share of photons/mode "
              f"{cw['brillouin_mode'] / cw['photons_per_mode']:.2e}")
        if mode["shape"] == "exponential":
            # finding: ±10.8 GHz is below the 18 GHz grid Nyquist, so the gated path uses the discrete sum of
            # the sampled exponential (|1 − e^{iωΔt}| vs ωΔt) while CW uses the analytic Lorentzian
            k = 2.0 / 14.5e-9
            t_, _ = mode_time(mode)
            dts = t_[1] - t_[0]
            w = 2 * np.pi * 10.8e9
            disc = abs(dts / (1 - np.exp((-k + 1j * w) * dts))) ** 2 * 2 * k
            assert abs(rb - disc / (2 * k / (k * k + w * w))) < 0.02 and rb > 1.3
        else:
            assert abs(rb - 1) < 2e-3
        print(f"  {mode['shape']:11s}: clicks/gate {g['clicks_per_gate']:.4e} (CW {cw['clicks_per_gate']:.4e}), "
              f"photons/mode {g['photons_per_mode']:.4e} (CW {cw['photons_per_mode']:.4e})")

    print("2. ports 15 vs 16, both columns (Rayleigh, ps mode)")
    kw = dict(p_avg_w=1e-3 * 10 / 320, duty=10 / 320, period_s=320e-9, eps_aom=0.0, gating=True, n_g=N_G,
              reflection_z_m=[])
    offs = np.arange(640) * 0.5e-9
    res = {}
    for k in (15, 16):
        p = dt.stray_params(L_of(k))
        rows = [gated_stray_photons(p, RAY, t_offset_s=o, **kw) for o in offs]
        gate = np.array([r["rayleigh_gate"] for r in rows])
        mode = np.array([r["rayleigh_mode"] for r in rows])
        cw = stray_photons(p, RAY, kw["p_avg_w"])
        res[k] = (gate, mode, cw)
        print(f"  k = {k}: min gate {gate.min():.2e}, min mode {mode.min():.2e}; mean/CW gate {gate.mean() / cw['rayleigh_gate']:.4f},"
              f" mode {mode.mean() / cw['rayleigh_mode']:.4f}")
    g15, m15, _ = res[15]
    assert np.any((g15 == 0) & (m15 == 0))
    g16, m16, cw16 = res[16]
    assert g16.min() > 0 and m16.min() > 0
    assert abs(m16.mean() / cw16["rayleigh_mode"] - 1) < 2e-3

    print("3. brute force: random scatterers vs mode_factors (ns exponential mode, narrow line at f0)")
    rng = np.random.default_rng(4)
    t, u = mode_time({"shape": "exponential", "tau_eff_s": 14.5e-9})
    dts = t[1] - t[0]
    per, tp, off, W = 320e-9, 40e-9, -30e-9, 2 * N_G * L_of(9) / C_M_PER_S
    for f0 in (0.0, 12e6):
        g_mf, _ = mode_factors(t, u, period_s=per, tau_pr_s=tp, t_offset_s=off, roundtrip_s=W, f0_hz=f0)
        taus = rng.uniform(0, W, 4000)
        vals = []
        for tau in taus:
            s = (t - off - tau) % per
            a = (s < tp).astype(float)
            vals.append(abs(np.sum(u * a * np.exp(2j * np.pi * f0 * t)) * dts) ** 2)
        mc = float(np.mean(vals))
        se = float(np.std(vals) / np.sqrt(len(vals)))
        print(f"  f0 = {f0 / 1e6:4.0f} MHz: mode_factors {g_mf:.4e} s, Monte Carlo {mc:.4e} ± {se:.1e} s")
        assert abs(g_mf - mc) < 4 * se + 1e-3 * mc
    print("all checks passed")
