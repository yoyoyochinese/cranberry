"""Phase 5 M3 review (Ref B): θ̃ along the pump axis (D71-4) — closed loops of refB_nominal vs Fig. 3's data.

1. refB_nominal closed loop (0.3 s, preset seed) at the nine Fig. 3 pump powers, gains NORMALISED with A±(x) (preset).
2. Same with the loop gains FIXED at their 100 mW values (ki of pid1 / pid2 frozen; an electronic servo whose gain is
   not re-tuned per pump power): |L(f_res)| scales with ΔG(x)/ΔG(0.745) (lock 1) and A+(x)/A+(0.745) (lock 2).
3. Fig. 3 observed points (read from the 600-dpi rendering of p. 4326; 0 dB at the axis, 37 px/dB at our scale; the 100 mW
   readings −9.00 / +15.11 dB reproduce the quoted −9.01 / +15.12 dB): circuit noise −21.7 dB removed as Ref B §4 does,
   then (x, θ̃) solved jointly per point from (R′−, R′+) with L(x) of Eq. (3), η = ζξ²η, f = 1 MHz, Eq. (1).

Run: uv run python reviews/scripts/m3_refB_theta_vs_pump.py [--no-sim]
"""

import importlib.util
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
PUMPS = [10, 20, 50, 65, 90, 100, 120, 130, 150]

# ---- Fig. 3 readings [dB] (observed, NOT circuit-noise subtracted) -------------------------------------------------
OBS = {  # P [mW]: (R′− obs, R′+ obs)
    10: (-3.16, 2.81), 20: (-5.00, 5.49), 50: (-6.89, 9.27), 65: (-7.84, 11.81), 90: (-8.51, 14.11),
    100: (-9.00, 15.11), 120: (-9.00, 17.30), 130: (-8.78, 18.59), 150: (-8.43, 19.46),
}
C_LIN = 10 ** (-2.17)
ETA_D = 0.99 * 0.988**2 * 0.998
C = 299_792_458.0


def subtract(db):
    return (10 ** (db / 10)) * (1 + C_LIN) - C_LIN


def r_pm(x):
    l = 0.00249 + 0.00222 * x
    k = C * (0.123 + l) / 1.0
    w2 = (2 * math.pi * 1e6 / k) ** 2
    eta = 0.123 / (0.123 + l) * ETA_D
    return 1 - eta * 4 * x / ((1 + x) ** 2 + w2), 1 + eta * 4 * x / ((1 - x) ** 2 + w2)


def eq1(x, th):
    rm, rp = r_pm(x)
    c2, s2 = math.cos(th) ** 2, math.sin(th) ** 2
    return rm * c2 + rp * s2, rp * c2 + rm * s2


def infer(p_mw):
    m, p = (subtract(v) for v in OBS[p_mw])
    th = math.radians(1.5)
    x = None
    for _ in range(50):
        x = brentq(lambda xx: eq1(xx, th)[1] - p, 1e-4, 0.9999)
        rm, rp = r_pm(x)
        s2 = (m - rm) / (rp - rm)
        if s2 < 0:
            return x, float("nan"), m, rm
        th = math.asin(math.sqrt(s2))
    return x, th, m, rm


def fig3_part():
    print("3. Fig. 3 observed points: joint (x, θ̃) (circuit noise removed)")
    print("   P[mW]  x_nom   x_fit   R′−_sub[dB]  R−(x_fit)[dB]  θ̃[°]  θ̃ range for ±0.14 dB on R′−")
    for pw in PUMPS:
        x, th, m, rm = infer(pw)
        rng = []
        for dm in (-0.14, 0.14):
            save = OBS[pw]
            OBS[pw] = (save[0] + dm, save[1])
            rng.append(math.degrees(infer(pw)[1]))
            OBS[pw] = save
        print(f"   {pw:4d}  {math.sqrt(pw / 180):.3f}  {x:.3f}   {10 * math.log10(m):7.2f}     {10 * math.log10(rm):7.2f}"
              f"     {math.degrees(th):5.2f}  [{rng[0]:.2f}, {rng[1]:.2f}]")
    # model curve vs observed at the nominal x, both θ̃
    print("   model (x_nom, circuit noise added) − observed R′− [dB]: θ̃=1.5° / 3.9°")
    for pw in PUMPS:
        x = math.sqrt(pw / 180)
        out = []
        for th in (1.5, 3.9):
            rm = eq1(x, math.radians(th))[0]
            out.append(10 * math.log10((rm + C_LIN) / (1 + C_LIN)) - OBS[pw][0])
        print(f"   {pw:4d}: {out[0]:+.2f} / {out[1]:+.2f}")


def sim_part():
    sp = importlib.util.spec_from_file_location("rs", ROOT / "tests/integration/refB_support.py")
    rs = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(rs)
    os.chdir(ROOT)
    from cranberry.runner.run import run

    def theta(res):
        nl = res.metrics.noise_levels["i.port0"]
        return math.degrees(float(nl["sigma_total_rad"][1])), float(nl["r_minus_db"][1]), float(nl["v_minus_db"][1]), \
            float(nl["v_plus_db"][1])

    print("1. normalised gains (preset)")
    norm = {}
    for pw in PUMPS:
        th, rm, vm, vp = theta(run(rs.load(params={"pump_w": pw * 1e-3}, run={"duration_s": 0.3})))
        rm15 = rs.eq1_db(vm, vp, math.radians(1.5))[0]
        norm[pw] = th
        print(f"   {pw:4d} mW: θ̃ = {th:.3f}°  R′− = {rm:.3f} dB (θ̃ = 1.5°: {rm15:.3f}, diff {rm - rm15:+.3f} dB)")

    p100 = rs.params_of(rs.load(params={"pump_w": 0.1}))
    ki1 = -2 * math.pi * p100["f_u1_hz"] / (p100["k1_v_per_rad"] * 2 * p100["pzt_k_rad_per_v"])
    ki2 = -2 * math.pi * p100["f_u2_hz"] / (p100["k2_v_per_rad"] * p100["pzt_k_rad_per_v"])
    print("2. gains fixed at the 100 mW values")
    for pw in PUMPS:
        p = rs.params_of(rs.load(params={"pump_w": pw * 1e-3}))
        g1 = p["k1_v_per_rad"] / p100["k1_v_per_rad"]
        g2 = p["k2_v_per_rad"] / p100["k2_v_per_rad"]
        lf1, lf2 = 0.3 * g1, 0.3 * g2
        if max(lf1, lf2) >= 1.0:
            print(f"   {pw:4d} mW: |L(f_res)| = {lf1:.2f} / {lf2:.2f} ≥ 1 → linearly unstable (not run)")
            continue
        pipe = [{"name": "pid1", "params": {"ki_per_s": ki1}}, {"name": "pid2", "params": {"ki_per_s": ki2}}]
        th, rm, vm, vp = theta(run(rs.load(params={"pump_w": pw * 1e-3}, run={"duration_s": 0.3}, pipeline=pipe)))
        print(f"   {pw:4d} mW: |L(f_res)| = {lf1:.2f} / {lf2:.2f}  θ̃ = {th:.3f}°  R′− = {rm:.3f} dB")


if __name__ == "__main__":
    if "--no-sim" not in sys.argv:
        sim_part()
    fig3_part()
