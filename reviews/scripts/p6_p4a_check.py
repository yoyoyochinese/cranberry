"""Phase 6, P4 (a): proxy bias of a different-wavelength probe from the dispersion of the thermo-optic coefficient.

Checks the two formulas of docs/phase6_prompt.md (P4 (a)) against the path model (WP-F6) and the preset coefficients.

* Path-length actuator holding the reference (probe) phase:
      δφ_s = (2πL n_s/λ_s)(β_s/n_s − β_r/n_r) δT ≈ −2πΔν·L·(κ_τ − (n_g/n)·k_T/ω_s)·δT        (first order in Δν)
  (α, the thermal expansion, cancels; k_T = (ω/c)(β + nα), κ_τ = (dn_g/dT + n_g α)/c, n = phase index.)
* Phase-only actuator: residual = δφ_s^free − δφ_r^free. The prompt's (1 − λ_s/λ_r)·δφ_s^free assumes equal optical
  path changes; with dispersion the first-order value is −2πΔν·L·κ_τ·δT (κ_τ·ω_s/k_T ≈ 1.08 times the prompt's).

The presets (probe/freq_*.yaml) give the probe dT / tec coefficient (k_T + 2πΔν κ_τ)·L and the stretcher 1.0 for every
beam: the stretcher acts as a phase-only actuator. A path-length stretcher needs the beam coefficient
g_r/g_s = (ν_r n_r)/(ν_s n_s) ≈ 1 + (Δν/ν_s)(n_g/n), which no preset sets (and the presets have no phase index n).

Part 1 reads the resolved coefficients of F_liu20_nominal + probe/pol_same + probe/freq_*.
Part 2 runs segment_phase → beam_assembler (path.testing helpers) with a 1 K step and an ideal static lock written
on the stretcher (as coded / as a path-length stretcher) or on the TEC, and compares the squeezed-light residual with
the formulas. Part 3 sizes the residual per port length and over the swept k_T, κ_τ and an assumed phase index n.
Usage: uv run python reviews/scripts/p6_p4a_check.py
"""
import math
from pathlib import Path

import numpy as np

import cranberry.path.testing as pt
from cranberry import config
from cranberry.path.topology import as_topology
from cranberry.runner.sweep import expand

ROOT = Path(__file__).resolve().parents[2]
C = 299792458.0
N_PHASE = 1.444   # phase index [AE]: not a preset parameter (only n_g = 1.468 exists)

# ---------------------------------------------------------------- part 1: preset coefficients
rows = {}
for ov in ("freq_80MHz", "freq_7GHz", "freq_100GHz", "freq_1THz"):
    spec = expand(config.load(ROOT / "presets/future/F_liu20_nominal.yaml", overlays=["probe/pol_same", f"probe/{ov}"]))[0].spec
    p = spec.offline["params"]
    topo = as_topology(spec.topology)
    c = {(k.driver, k.beam): v for k, v in topo.coef.items() if k.segment == "F1" and k.effect == "phase"}
    rows[ov] = (p, c)
    dnu, L = p["probe_freq_offset_hz"], p["L_k_m"]
    assert c[("stretcher", "pr1")] == c[("stretcher", "sqz1")] == 1.0          # phase-only actuator
    assert c[("strain", "pr1")] == c[("strain", "sqz1")]                        # strain not scaled with frequency
    assert math.isclose(c[("dT", "pr1")] - c[("dT", "sqz1")], 2 * math.pi * dnu * p["kappa_tau_s_per_m_K"] * L, rel_tol=1e-9)
    assert c[("tec", "pr1")] == c[("dT", "pr1")] and c[("tec", "sqz1")] == c[("dT", "sqz1")]
    print(f"{ov:12s} L_k = {L:7.1f} m  c_s = {c[('dT', 'sqz1')]:.4e}  c_r − c_s = {c[('dT', 'pr1')] - c[('dT', 'sqz1')]:+.4e} rad/K"
          f"  stretcher g_s = g_r = {c[('stretcher', 'sqz1')]}")

# ---------------------------------------------------------------- part 2: time loop through the path model


def run(c_s, c_r, g_s, g_r, t_s, t_r, dnu, actuator):
    """1 K step on dT at n = 50; the ideal static lock writes the actuator at the same sample. Returns (Δφ_s, Δφ_r)."""
    seg = pt.seg_spec("F", 30.0, n_sub=4, drivers={"dT": "uniform", "stretcher": "actuator", "tec": "uniform"}, actuator_sub=0)
    beams = [pt.beam_spec("sq", [pt.step_spec("F", 1)], "o"),
             pt.beam_spec("pr", [pt.step_spec("F", -1)], "o", freq_offset_hz=dnu)]
    coef = (pt.coef_spec("F", "dT", "sq", c_s) + pt.coef_spec("F", "dT", "pr", c_r)
            + pt.coef_spec("F", "stretcher", "sq", g_s) + pt.coef_spec("F", "stretcher", "pr", g_r)
            + pt.coef_spec("F", "tec", "sq", t_s) + pt.coef_spec("F", "tec", "pr", t_r))
    topo = pt.topo_spec([seg], beams, ["o"], coef)
    n = 400
    dT, st, tec = np.zeros(n), np.zeros(n), np.zeros(n)
    dT[50:] = 1.0
    if actuator == "stretcher":
        st[50:] = -c_r / g_r            # holds the probe phase: c_r·1 + g_r·u = 0
    else:
        tec[50:] = -c_r / t_r
    pipe = [pt.pathtest("pathtest_table", "a", "drv.F.dT", (), "K", n_rows=n),
            pt.pathtest("pathtest_table", "b", "drv.F.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "d", "drv.F.tec", (), "K", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "beam.pr.phi"], n)]
    rec, _, _ = pt.run_and_record(pipe, 1e8, n, {"drv.F.dT_table": dT, "drv.F.stretcher_table": st, "drv.F.tec_table": tec})
    return rec["beam.sq.phi"][-1] - rec["beam.sq.phi"][0], rec["beam.pr.phi"][-1] - rec["beam.pr.phi"][0]


print("\nideal static lock, 1 K, port k of F_liu20_nominal (residual on the squeezed light, rad/K)")
print(f"{'Δν':>10s} {'code: stretcher':>16s} {'−2πΔνκ_τL':>12s} {'prompt phase-only':>18s} {'path-length stretcher':>22s}"
      f" {'prompt path-length':>19s} {'TEC':>8s}")
for ov, (p, c) in rows.items():
    dnu, L, kT, ktau = p["probe_freq_offset_hz"], p["L_k_m"], p["k_T_rad_per_m_K"], p["kappa_tau_s_per_m_K"]
    nu0 = C / p["lambda0_m"]
    w0 = 2 * math.pi * nu0
    c_s, c_r = c[("dT", "sqz1")], c[("dT", "pr1")]
    g_pl = (nu0 + dnu) * (N_PHASE + (p["n_g"] - N_PHASE) * dnu / nu0) / (nu0 * N_PHASE)   # (ν_r n_r)/(ν_s n_s), n_r to 1st order
    r_code, pr_code = run(c_s, c_r, 1.0, 1.0, c_s, c_r, dnu, "stretcher")
    r_pl, pr_pl = run(c_s, c_r, 1.0, g_pl, c_s, c_r, dnu, "stretcher")
    r_tec, pr_tec = run(c_s, c_r, 1.0, 1.0, c_s, c_r, dnu, "tec")
    f_code = -2 * math.pi * dnu * ktau * L
    f_phase = (1 - (nu0 + dnu) / nu0) * c_s                     # (1 − λ_s/λ_r)·δφ_s^free
    f_pl = -2 * math.pi * dnu * L * (ktau - (p["n_g"] / N_PHASE) * kT / w0)
    for pr in (pr_code, pr_pl, pr_tec):
        assert abs(pr) < 1e-9 * abs(c_r)                         # the lock holds the probe
    assert math.isclose(r_code, f_code, rel_tol=1e-9)
    assert math.isclose(r_code / f_phase, ktau * w0 / kT, rel_tol=1e-9)
    assert math.isclose(r_pl, f_pl, rel_tol=5 * dnu / nu0 + 1e-9, abs_tol=1e-9)   # agrees to O(Δν/ν)
    assert abs(r_tec) < 1e-9 * abs(c_s)
    print(f"{dnu:10.3g} {r_code:16.4e} {f_code:12.4e} {f_phase:18.4e} {r_pl:22.4e} {f_pl:19.4e} {r_tec:8.1e}")

# ---------------------------------------------------------------- part 3: size per length, parameter spread
print("\nΔν = 1 THz (8 nm), residual per kelvin [rad/K]: code stretcher (phase-only) vs path-length stretcher")
print("path-length coefficient R = κ_τ − (n_g/n)·k_T/ω_s [s/(m·K)] is a small difference of two swept parameters:")
nu0 = C / 1.54532e-6
w0 = 2 * math.pi * nu0
Rs = []
for kT in (40.0, 44.0, 48.0):
    for ktau in (37.4e-15, 39e-15):
        for n in (1.444, 1.468):
            R = ktau - (1.468 / n) * kT / w0
            Rs.append(R)
            print(f"  k_T = {kT:4.1f}  κ_τ = {ktau * 1e15:4.1f} ps/(km·K)  n = {n:5.3f}  R = {R:+.2e}  (R/κ_τ = {R / ktau:+.3f})")
for L in (63.0, 317.0, 1060.0, 1268.0, 2100.0):
    code = 2 * math.pi * 1e12 * 39e-15 * L
    lo, hi = (2 * math.pi * 1e12 * min(Rs) * L, 2 * math.pi * 1e12 * max(Rs) * L)
    print(f"  L = {L:6.0f} m: code (phase-only) {code:7.1f} rad/K; path-length {lo:+7.1f} … {hi:+7.1f} rad/K")
print("\nall checks passed")
