"""Shared parameter sets for the WP-B10 tests (fiber quantum layer, scattering tables). Not a block."""

from __future__ import annotations

from typing import Any

from cranberry.blocks.fib.gnlse import TOGGLES as GNLSE_TOGGLES
from cranberry.blocks.fib.quantum import QUANTUM_TOGGLES

B2, GAMMA, LAMBDA0 = -21.7e-27, 1.3e-3, 1545.32e-9
#: T07 §6 common values (UNVERIFIED typical: g_B, A_eff, ν_B, Δν_B; London et al. R0,7; T10 b_R)
GAWBS_SMF = dict(cladding_diameter_m=125e-6, v_l_m_s=5996.0, v_s_m_s=3740.0, n_modes=12,
                 gamma_om_per_w_m=[0.0] * 6 + [2e-3] + [0.0] * 5, linewidth_hz=4.3e6, mode_sum_factor=1.0)
BRILLOUIN = dict(g_b_m_per_w=5e-11, a_eff_m2=80e-12, nu_b_hz=10.8e9, dnu_b_hz=30e6)
PROBE = dict(p_peak_w=1e-3, duty=1.0, freq_offset_hz=0.0, linewidth_hz=1e5, polarization="same", eps_xt=0.0,
             t_perp=0.0)


def toggles(value: bool = False, **over: bool) -> dict[str, bool]:
    """Every ``fiber_quantum_linearized`` toggle set to ``value``, then ``over``."""
    return {**dict.fromkeys((*GNLSE_TOGGLES, *QUANTUM_TOGGLES), value), **over}


def linear_params(**over: Any) -> dict[str, Any]:
    """A soliton-order-1 pulse (T0 = 1 ps) over one dispersion length on a 64-point grid."""
    t0 = 1e-12
    ld = t0**2 / abs(B2)
    p = dict(L_m=ld, alpha_per_m=0.0, beta2_s2_per_m=B2, beta3_s3_per_m=0.0, beta_higher=[], gamma_per_w_m=GAMMA,
             lambda0_m=LAMBDA0, raman_model="blow_wood", f_R=0.18, xpm_b=0.0, xpm_power_w=0.0, temperature_k=295.0,
             laser_mode="pulsed", pulse=dict(shape="sech", peak_power_w=abs(B2) / (GAMMA * t0**2), t0_s=t0),
             grid=dict(n_t=64, t_window_s=20 * t0), steps=dict(n_steps=200, n_save=4, control="fixed"),
             input_state=dict(kind="vacuum"), phi_opa0_rad=0.0, gawbs=dict(GAWBS_SMF), probe=dict(PROBE),
             brillouin=dict(BRILLOUIN), rayleigh=dict(b_r_per_m=6.2e-8))
    p.update(over)
    return p
