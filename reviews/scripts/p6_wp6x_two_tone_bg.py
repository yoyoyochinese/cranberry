"""Phase 5 review of WP6-X (D94-5, the open point): group birefringence of orthogonal (fast-axis) probe tones in
probe/two_tone.

The tones pr_k / pr_km run on the FAST axis of F_k (two_tone.yaml beam paths: axis: fast), but the static phase
E = 2πδν(τ_pr − τ_ref) uses Topology.timing (n_g·L, slow axis). The fast axis adds b_g·L_k (F_common pol_params
b_g_s_per_m = 1.36e-12 s/m, PM15-U25D, Ritari 2004), so port k sees E_k = E_0 + 2πδν·b_g·(L_k − L_0) beyond the grid
condition, and the bias −atan(r tan E_k) and the slope cos E_k differ per port. The thermal coefficient of b_g
(K_tau_s_per_m_K = −1.78e-15) adds ±2πδν·K_τ·L_k per kelvin to E (not in the ±δν rows either).
Prints, for F_liu20_nominal + probe/two_tone (7 GHz on grid) and the ports of the 32-port design: the extra E mod 2π,
the bias at r = 0.1 / 0.2 (E_0 = 0 trimmed), cos E_k; the corrected grid (probe-axis group index) and the b_g accuracy it
needs for 0.1 rad at the longest port.
Usage: uv run python reviews/scripts/p6_wp6x_two_tone_bg.py
"""
import math
from pathlib import Path

from cranberry import config
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
C = 299792458.0
B_G, K_TAU = 1.36e-12, -1.78e-15

spec = config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=["probe/two_tone"])
spec = expand(spec, {"offline.params.tt_dnu_target_hz": [7.0e9], "offline.params.tt_on_grid": [True],
                     "offline.params.tt_r": [0.2], "offline.params.tt_dL0_detune_rad": [0.0]})[0].spec
p = spec.offline["params"]
dnu, n_g = p["tt_dnu_hz"], p["n_g"]
dl1 = (p["N_b"] + p["n_guard"]) * C / (p["f_rep_hz"] * n_g)       # port spacing in length
print(f"F_liu20 + two_tone: δν = {dnu:.6g} Hz, m = {p['tt_m']}, port spacing ΔL = {dl1:.2f} m, n_g = {n_g}, "
      f"K_ports = {p['K_ports']}")
print(f"{'port':>5s} {'L_k − L_0 [m]':>14s} {'extra E [rad]':>14s} {'mod 2π':>8s} {'cos E':>7s} {'bias r=0.1':>11s}"
      f" {'bias r=0.2':>11s} {'dE/dT [rad/K]':>14s}")
for k in (1, 2, 3, 8, 15, 31):
    dl = k * dl1
    e = 2 * math.pi * dnu * B_G * dl
    em = math.remainder(e, 2 * math.pi)
    b1, b2 = (-math.atan(r * math.tan(em)) for r in (0.1, 0.2))
    print(f"{k:5d} {dl:14.1f} {e:14.1f} {em:+8.3f} {math.cos(em):+7.3f} {b1:+11.3f} {b2:+11.3f}"
          f" {2 * math.pi * dnu * K_TAU * (p['L0_m'] + dl):+14.4f}")
f_fix = 1.0 / (1.0 + C * B_G / n_g)
print(f"\ngrid with the fast-axis group index n_g + c·b_g = {n_g + C * B_G:.5f}: δν × {f_fix:.6f} "
      f"(shift {dnu * (f_fix - 1) / 1e6:+.2f} MHz)")
L_max = p["L0_m"] + (p["K_ports"] - 1) * dl1
tol = 0.1 / (2 * math.pi * dnu * L_max)
print(f"b_g accuracy for |ΔE| < 0.1 rad at the longest port (L = {L_max:.0f} m): {tol:.2e} s/m ({tol / B_G:.1e} relative)")
