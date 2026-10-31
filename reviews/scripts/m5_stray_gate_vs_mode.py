"""M5 review: does the WP-B8 gate factor (flux averaged uniformly over the detection gate) give the
photons per squeezed *mode* (conv §2, N = ∫ n|ũ|² df, i.e. weighted by the mode, not by the gate)?

A discrete reflection returns the 10 ns probe pulse delayed by 2n_g z/c. ``gate_factors`` scales both
clicks/gate and photons/mode by f = |gate ∩ pulse|/T_gate. The photons in the Borealis-type optimal
mode u(t) are, for a line much narrower than the mode bandwidth (Rayleigh / reflection, coherence
≫ mode), N = Φ_pk |∫u(t)·1_pulse(t) dt|² (= f_mode · Φ_pk τ_eff), and for broadband light
N = Φ_pk ∫|u|²·1_pulse dt. The gate is centred on the mode centroid (25 ns, Borealis evaluation window).

Run: uv run python reviews/scripts/m5_stray_gate_vs_mode.py
"""
import numpy as np

from cranberry.blocks.mod.probe_pulser import gate_factors
from cranberry.blocks.source.pulsed_testing import BOREALIS, pulsed_opo

if __name__ == "__main__":
    b = pulsed_opo(BOREALIS, x_ss=2.4)
    u, g = b.optimal_mode[0].real, b.grid
    dt = g[1] - g[0]
    tc = float(np.sum(g * u**2) * dt)
    tau_eff = float(np.sum(u) * dt) ** 2
    t_gate, tau_pr, period = 25e-9, 10e-9, 1e-6
    print(f"mode centroid {tc * 1e9:.2f} ns after pump turn-on, tau_eff {tau_eff * 1e9:.2f} ns, gate {t_gate * 1e9:.0f} ns")
    print(" pulse start rel. gate centre [ns] | gate factor f | mode factor (narrow line) | mode factor (broadband)")
    rows = []
    for start in (-15e-9, -12.5e-9, -7.5e-9, -2.5e-9, 2.5e-9, 7.5e-9):
        # reflection with zero extra delay: probe pulse occupies [start, start + tau_pr) in gate coordinates
        _, (f,) = gate_factors(period_s=period, tau_pr_s=tau_pr, t_offset_s=start, t_gate_s=t_gate,
                               roundtrip_s=1e-9, delays_s=[0.0])
        ind = ((g - tc >= start) & (g - tc < start + tau_pr)).astype(float)
        fm_line = float(np.sum(u * ind) * dt) ** 2 / tau_eff
        fm_bb = float(np.sum(u**2 * ind) * dt)
        print(f" {start * 1e9:+6.1f} | {f:.3f} | {fm_line:.3f} | {fm_bb:.3f}")
        rows.append((f, fm_line))
    # finding M5-02: the same gate factor 0.4 corresponds to mode factors from 0.035 to 0.47
    same = [m for f, m in rows if abs(f - 0.4) < 1e-9]
    assert max(same) / min(same) > 10
