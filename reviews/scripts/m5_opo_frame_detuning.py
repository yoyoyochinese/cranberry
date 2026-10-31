"""M5 review: frame of pulsed_opo_source (D60-1) and the detuning coefficients.

1. Static Δ_p: the angle measured with the Δ = 0 optimal LO from the *full* pump b(t) must equal the
   angle of ``cov`` (static pump-cavity phase removed) + ½·pump_path_phase_rad (conv §4: ½ψ).
2. Sign/size of the signal-detuning rotation vs the CW block (T02 −d/(1 + x²)) in the long-pump limit.
3. First-order map vs the exact angle at d = 0.05 for Δ_p = 2Δ_s (length change), Borealis x_ss = 1.3.

Run: uv run python reviews/scripts/m5_opo_frame_detuning.py
"""
import numpy as np

from cranberry.blocks.source.pulsed_testing import BOREALIS, pulsed_opo
from cranberry.offline.source_output import wrap_half_pi


def angle(b, t):
    u0 = b.optimal_mode[0]
    return b._lo_angle(t, u0)


if __name__ == "__main__":
    kc = np.pi * BOREALIS["fwhm_hz"]
    kp = np.pi * BOREALIS["pump_fwhm_hz"]
    print("1. static pump detuning: theta(full b) vs theta(cov) + psi/2")
    for dpk in (0.05, 0.2, 0.5):
        b = pulsed_opo(BOREALIS, x_ss=1.3, delta_p_rad_s=dpk * kc, toggles={"detuning": True})
        u0 = pulsed_opo(BOREALIS, x_ss=1.3).optimal_mode[0]
        th_full = b._lo_angle(b.transfer(delta_p=dpk * kc), u0)
        th_cov = b._lo_angle(b.transfer(delta_p=dpk * kc, remove_static_pump_phase=True), u0)
        psi = b.source_output.pump_path_phase_rad["opo"]
        print(f"  Δp/κc={dpk}: full={th_full:.6f} cov={th_cov:.6f} psi={psi:.6f} (−atan(Δp/κp)={-np.arctan(dpk*kc/kp):.6f})"
              f"  full−(cov+psi/2)={float(wrap_half_pi(th_full - th_cov - psi / 2)):.2e}")
        assert abs(float(wrap_half_pi(th_full - th_cov - psi / 2))) < 1e-12  # D60-1 frame is consistent
    print("2. long pump (non-resonant, T_p κ_c = 30, x = 0.5) vs CW −d/(1+x²)")
    b = pulsed_opo(BOREALIS, x_ss=0.5, escape_efficiency=1.0, pump_duration_s=30 / kc, window_kappa=40.0,
                   toggles={"pump_resonant": False})
    cs = b.angle_coefficients(0.5)[0] * kc
    print(f"  pulsed: dθ/d = {cs:.4f}; CW T02: {-1 / 1.25:.4f}")
    assert abs(cs + 0.8) < 0.01
    print("3. length change d = 0.05 (Δp = 2Δs): exact vs first-order map")
    b = pulsed_opo(BOREALIS, x_ss=1.3, angle_probe_d=0.01)
    cs, cp, _ = b.angle_coefficients(1.3)
    u0 = b.optimal_mode[0]
    th0 = b._lo_angle(b.transfer(), u0)
    for d in (0.01, 0.05, 0.1):
        ex = float(wrap_half_pi(b._lo_angle(b.transfer(delta_s=d * kc, delta_p=2 * d * kc), u0) - th0))
        print(f"  d={d}: exact {ex / d:.4f}·d, map {(cs + 2 * cp) * kc:.4f}·d")
        assert abs(ex / d - (cs + 2 * cp) * kc) < 2e-3
