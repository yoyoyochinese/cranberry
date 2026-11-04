"""Phase 5 review of WP6-X (D94-5, the open point): group birefringence of orthogonal (fast-axis) probe tones in
probe/two_tone.

The tones pr_k / pr_km run on the FAST axis of F_k (two_tone.yaml beam paths: axis: fast), but the static phase
E = 2πδν(τ_pr − τ_ref) uses Topology.timing (n_g·L, slow axis). The fast axis adds b_g·L_k (F_common pol_params
b_g_s_per_m = 1.36e-12 s/m, PM15-U25D, Ritari 2004), so port k sees E_k = E_0 + 2πδν·b_g·(L_k − L_0) beyond the grid
condition, and the bias −atan(r tan E_k) and the slope cos E_k differ per port. The thermal coefficient of b_g
(K_tau_s_per_m_K = −1.78e-15) adds ±2πδν·K_τ·L_k per kelvin to E (not in the ±δν rows either).
Prints, for F_liu20_nominal + probe/two_tone (7 GHz on grid) and the ports of the 32-port design: the extra E mod 2π,
the bias at r = 0.1 / 0.2 (E_0 = 0 trimmed), cos E_k; the corrected grid (probe-axis group index) and the b_g accuracy it
needs for 0.1 rad at the longest port. (The analytic part uses the pre-fix grid, tt_grid_probe_axis: false. Sign: in
this model τ_fast = τ_slow − b_g·L (path/segments.py), so the extra E is −2πδν·b_g·ΔL and the fast-axis index is
n_g − c·b_g; the magnitudes below are sign-independent.)

After the fix (Phase 6 worker, X-7): E_0, E_2, E_31 from the resolved topology (Topology.timing, b_g as the
reference-arm equivalent delay of probe/two_tone) for the fixed grid, the pre-fix grid and probe/two_tone_copol.
Usage: uv run python reviews/scripts/p6_wp6x_two_tone_bg.py
"""
import math
from pathlib import Path

from cranberry import config
from cranberry.path.topology import as_topology
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
C = 299792458.0
B_G, K_TAU = 1.36e-12, -1.78e-15

spec = config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=["probe/two_tone"])
spec = expand(spec, {"offline.params.tt_dnu_target_hz": [7.0e9], "offline.params.tt_on_grid": [True],
                     "offline.params.tt_r": [0.2], "offline.params.tt_dL0_detune_rad": [0.0],
                     "offline.params.tt_grid_probe_axis": [False]})[0].spec
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


def e_ports(overlays, grid_probe_axis):
    """{0: E_0, 2: E_2, 31: E_31} [rad, mod 2π] of the + tone from the resolved topology, plus δν and n_g,pr."""
    out = {}
    for k in (2, 31):
        s = config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=list(overlays))
        s = expand(s, {"offline.params.tt_dnu_target_hz": [7.0e9], "offline.params.tt_on_grid": [True],
                       "offline.params.tt_r": [0.2], "offline.params.tt_dL0_detune_rad": [0.0],
                       "offline.params.port_k": [k], "offline.params.tt_grid_probe_axis": [grid_probe_axis]})[0].spec
        q, topo = s.offline["params"], as_topology(s.topology)

        def e(pr, ref, q=q, topo=topo):
            d = topo.timing(pr).tau_s - topo.timing(ref).tau_s
            return math.remainder(2 * math.pi * q["tt_dnu_hz"] * d, 2 * math.pi)

        out[0], out[k] = e("pr0", "ref0"), e("pr1", "ref1")
        out["dnu"], out["n_pr"] = q["tt_dnu_hz"], q["tt_n_g_probe"]
    return out


print("\nafter the fix: E from Topology.timing, F_liu20, 7 GHz on grid, E_0 trimmed to 0 (bias at r = 0.2):")
print(f"{'variant':>38s} {'δν [Hz]':>15s} {'n_g,pr':>8s} {'E_0':>8s} {'E_2':>8s} {'E_31':>8s} {'cos E_2':>8s}"
      f" {'bias_2':>8s} {'bias_31':>8s}")
for what, ovs, gp in (("fast axis, probe-axis grid (fix)", ["probe/two_tone"], True),
                      ("fast axis, pre-fix grid (bad, D34)", ["probe/two_tone"], False),
                      ("co-polarised (probe/two_tone_copol)", ["probe/two_tone", "probe/two_tone_copol"], True)):
    o = e_ports(ovs, gp)
    b2, b31 = (-math.atan(0.2 * math.tan(o[k])) for k in (2, 31))
    print(f"{what:>38s} {o['dnu']:15.9e} {o['n_pr']:8.5f} {o[0]:+8.4f} {o[2]:+8.4f} {o[31]:+8.4f} {math.cos(o[2]):+8.4f}"
          f" {b2:+8.4f} {b31:+8.4f}")
