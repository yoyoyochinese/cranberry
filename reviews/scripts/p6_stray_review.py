"""Phase 5 review of WP6-S (D88): window_stray_light / forward_ref_stray_table, independent checks.

Parts (each prints a table; asserts at the end of each part are the self-check of the review's claims):
  1. D88-3 numbers: P1 modulator leak (1 µW, −50 dB) per preset, with the overlay's Gaussian τ_mode stand-in and with
     the OPO optimal mode (resolved from the preset's own stray1 spec); the ε_mod·P_pk limit for 5.8e−4 photons/mode.
  2. P2 window tail: window mean factor at D = 0.65, T = 1 ms (and the test point T = 100 µs) — block vs. an independent
     brute-force backscatter integral over z, with and without the fibre attenuation e^{−2αz}.
  3. OPO-mode limitation: stray1 (window tail) photons per mode, Gaussian stand-in vs OPO mode, per mechanism.
  4. Edge tail: quasi-static mode count (block) vs the exact coherent overlap |∫√Φ(t) u(t) dt|².
  5. Placement bookkeeping: brute-force pulse timeline (K ports, emission order, fibre delays) vs _events.
  6. Guard slot: L_k factor and the D78 Δν_max per preset.
Usage: uv run python reviews/scripts/p6_stray_review.py
"""
import math
import tempfile
from pathlib import Path

import numpy as np
import yaml

from cranberry import config
from cranberry.blocks.fib import window_stray as ws
from cranberry.blocks.fib.scattering import mode_time, mode_weight, resolve_mode, stray_photons
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block

ROOT = Path(__file__).resolve().parents[2]
C = 299792458.0
H = 6.62607015e-34
TMP = Path(tempfile.mkdtemp())
PRESETS = ("F_wg6", "F_wg100", "F_bor6", "F_liu20")
OPO = ("F_bor6", "F_liu20")
# fix round (WS-02 / WS-08): the new axes fref_eps_gate and win_refl_r are pinned to their "none" values, so parts 1–6
# evaluate the same configuration as the review did
REF_POINT = dict(fref_p_peak_w=1.0e-6, fref_t_pulse_s=1.0e-9, fref_eps_mod=1.0e-5, fref_edge_tau_s=1.0e-9,
                 det_p_afterpulse=1.0e-3, det_tau_afterpulse_s=1.0e-7, det_recovery_s=4.3e-8, fref_eps_gate=1.0)
WIN_POINT = dict(win_period_s=1.0e-3, win_duty_lock=0.65, win_edge_s=1.0e-7, win_eps_off=1.0e-5, win_limit_statistic="mean",
                 win_refl_r=0.0)
OPO_OVL = {"F_bor6": "bor6", "F_liu20": "liu20"}
WG_POINT = dict(tau_p_s=1.0e-12)  # F_wg100 pulse-width axis (00_summary nominal 1 ps)


def spec_of(base, overlay, point, extra=None, pre=(), post=()):
    f = TMP / f"pt_{base}_{overlay.replace('/', '_')}.yaml"
    f.write_text(yaml.safe_dump({"offline": {"params": {**point, **(extra or {}), **(WG_POINT if base == "F_wg100" else {})}}}))
    s = config.load(ROOT / f"presets/future/{base}_nominal.yaml", overlays=[*pre, overlay, *post, str(f)])
    assert not config.sweep_markers(s), (base, config.sweep_markers(s))
    return s


def blk(spec, name):
    b = next(b for b in spec.offline["blocks"] if b["name"] == name)
    return make_block(BlockSpec(b["block"], name, b["params"]))


def tau_eff(mode):
    t, u = mode_time(mode)
    return abs(np.sum(u) * (t[1] - t[0])) ** 2


# ------------------------------------------------------------------------------------------------ part 1
print("== 1. P1 modulator leak, P_pk = 1 µW, ε_mod = −50 dB (D88-3: F_bor6 1.1, F_liu20 0.6 photons/mode)")
print(f"{'preset':8s} {'τ_eff Gauss':>11s} {'τ_eff OPO':>10s} {'block (overlay)':>15s} {'indep Gauss':>11s} {'OPO mode':>9s}"
      f" {'εP max Gauss':>12s} {'εP max OPO':>11s} {'clicks/gate':>11s}")
leak = {}
for base in PRESETS:
    s = spec_of(base, "ref/aurora_calframe", REF_POINT)
    p = s.offline["params"]
    b = blk(s, "frefstray1")
    r = b.check()
    hnu = H * C / p["lambda0_m"]
    te_g = p["tau_mode_s"] * math.sqrt(math.pi / math.log(2))
    stray_mode = next(x for x in s.offline["blocks"] if x["name"] == "stray1")["params"]["mode"]
    te_o = tau_eff(resolve_mode(stray_mode)) if base in OPO else te_g
    phi = 1e-5 * 1e-6 / hnu
    nmax = p["stray_n_mode_max"]
    leak[base] = (r["mod_leak_mode_mean"], phi * te_g, phi * te_o, nmax * hnu / te_g, nmax * hnu / te_o)
    print(f"{base:8s} {te_g * 1e9:9.2f}ns {te_o * 1e9:8.2f}ns {r['mod_leak_mode_mean']:15.3g} {phi * te_g:11.3g} "
          f"{phi * te_o:9.3g} {nmax * hnu / te_g:12.3g} {nmax * hnu / te_o:11.3g} {r['mod_leak_gate_mean']:11.3g}")
    assert math.isclose(r["mod_leak_mode_mean"], phi * te_g, rel_tol=1e-6)  # block = Φ_leak·τ_eff (Gaussian stand-in)
    # the leak dominates the slot mean at this point (D88-3)
    assert r["mod_leak_mode_mean"] > 0.9 * r["photons_per_mode_mean"]
# D88-3's 1.1 / 0.6 are the OPO-mode values, the overlay itself prints 2.6–2.8× more (Gaussian τ_eff is longer)
assert abs(leak["F_bor6"][2] - 1.1) < 0.06 and abs(leak["F_liu20"][2] - 0.6) < 0.05
assert 2.4 < leak["F_bor6"][0] / leak["F_bor6"][2] < 3.0 and 2.4 < leak["F_liu20"][0] / leak["F_liu20"][2] < 3.0

# ------------------------------------------------------------------------------------------------ part 2
print("\n== 2. P2 window tail mean factor (ε_off = 1e−5, edge 100 ns, one slot per frame)")


def brute_window(p, alpha_per_m, nz=4000):
    """Independent: P_bs(t) ∝ ∫_0^L a(t − 2 n_g z/c) e^{−2αz} dz / ∫_0^L e^{−2αz} dz, a = trapezoidal envelope."""
    T, D, e, ph = p["period_s"], p["duty_lock"], p["edge_s"], p["phase_s"]
    t_off = ph + D * T                                       # 50 % switch-off point
    n = int(math.floor(((1 - D) * T - e) / p["slot_period_s"] + 1e-9))
    ts = t_off + e / 2 + (np.arange(n) + 0.5) * p["slot_period_s"]
    z = (np.arange(nz) + 0.5) * p["L_m"] / nz
    wz = np.exp(-2 * alpha_per_m * z)

    def a(t):
        x = np.mod(t - ph, T)                                # time since the 50 % switch-on point (mod T)
        up = np.clip((x + e / 2) / e, 0, 1) if e > 0 else (x >= 0) * 1.0
        dn = np.clip((D * T + e / 2 - x) / e, 0, 1) if e > 0 else (x < D * T) * 1.0
        up2 = np.clip((x - T + e / 2) / e, 0, 1)             # next period's ramp-up
        return np.maximum(np.minimum(up, dn), up2)

    h = np.array([np.sum(a(t - 2 * p["n_g"] * z / C) * wz) / np.sum(wz) for t in ts])
    return ts, p["eps_off"] + (1 - p["eps_off"]) * h


print(f"{'preset':8s} {'T':>7s} {'L_k [m]':>8s} {'W [µs]':>7s} {'block':>8s} {'brute α=0':>9s} {'brute α':>8s} {'box/α':>6s}"
      f" {'D88-3':>6s}")
d883 = {"F_wg100": 0.046, "F_wg6": 0.015, "F_bor6": 0.015, "F_liu20": 0.018}
win = {}
for base in PRESETS:
    for T in (1e-3, 1e-4):
        s = spec_of(base, "stray/window_tail", {**WIN_POINT, "win_period_s": T})
        b = blk(s, "winstray1")
        p = b.params
        r = b.check()
        _, f0 = brute_window(p, 0.0)
        _, fa = brute_window(p, p["alpha_per_m"])
        win[(base, T)] = (r, f0.mean(), fa.mean())
        print(f"{base:8s} {T * 1e6:5.0f}µs {p['L_m']:8.0f} {r['tail_s'] * 1e6:7.2f} {r['mean_factor']:8.4f} {f0.mean():9.4f} "
              f"{fa.mean():8.4f} {r['mean_factor'] / fa.mean():6.3f} {d883[base] if T == 1e-3 else '':>6}")
        # fix round (WS-03): the block now carries the e^{−2αz} weight, so it matches the attenuated brute force (was:
        # the attenuation-free box f0)
        assert math.isclose(r["mean_factor"], fa.mean(), rel_tol=2e-3)
    # D88-3's attenuation-free numbers are the box (f0) values; the block's (with attenuation) are 5–16 % lower
    assert abs(win[(base, 1e-3)][1] - d883[base]) < 0.0015
# fix round (WS-03): was 1.10 < box/attenuated < 1.16 (box overestimate ~13 % on F_wg100); now block = attenuated
assert 0.998 < win[("F_wg100", 1e-3)][0]["mean_factor"] / win[("F_wg100", 1e-3)][2] < 1.002
assert 1.10 < win[("F_wg100", 1e-3)][1] / win[("F_wg100", 1e-3)][2] < 1.16                 # the old box's overestimate

# ------------------------------------------------------------------------------------------------ part 3
print("\n== 3. OPO-mode limitation: stray1 photons per mode at the run's P̄ (window 'on' value), Gaussian vs OPO mode")
print(f"{'preset':8s} {'probe':>12s} {'P̄':>7s} {'mech':>10s} {'Gauss':>10s} {'OPO':>10s} {'OPO/Gauss':>9s}")
ratio3 = {}
for base in OPO:
    for pre in ((), ("probe/freq_80MHz",), ("probe/freq_7GHz",), ("probe/freq_100GHz",)):
        s = spec_of(base, "stray/window_tail", WIN_POINT, pre=pre)
        p = s.offline["params"]
        b = blk(s, "winstray1")
        pp = dict(b.params)
        g = stray_photons(pp, b._toggles, p["p_pr_avg_w"])
        stray_mode = next(x for x in s.offline["blocks"] if x["name"] == "stray1")["params"]["mode"]
        o = stray_photons({**pp, "mode": resolve_mode(stray_mode)}, b._toggles, p["p_pr_avg_w"])
        lab = pre[0].split("_")[-1] if pre else "same"
        for m in ("rayleigh", "brillouin", "raman", "photons_per"):
            k = "photons_per_mode" if m == "photons_per" else f"{m}_mode"
            rr = o[k] / g[k] if g[k] > 0 else float("inf")
            if m == "photons_per" or not pre:
                print(f"{base:8s} {lab:>12s} {p['p_pr_avg_w']:7.1e} {m:>10s} {g[k]:10.3g} {o[k]:10.3g} {rr:9.3g}")
        ratio3[(base, lab)] = o["photons_per_mode"] / g["photons_per_mode"]
# same-frequency probe: the Gaussian stand-in is 2.7× HIGH (docstring / D88-4 say "2.7× below": sign reversed)
assert 1 / 2.9 < ratio3[("F_liu20", "same")] < 1 / 2.5
# frequency-offset probe (80 MHz): the Gaussian has no spectral tail, the stand-in is orders of magnitude LOW
assert ratio3[("F_liu20", "80MHz")] > 1e2 and ratio3[("F_bor6", "80MHz")] > 1e2

# ------------------------------------------------------------------------------------------------ part 4
print("\n== 4. Edge tail into the first data slot (guard_slot, j = 0, Δ = T_rep): quasi-static vs exact coherent overlap")
print(f"{'preset':8s} {'mode':>6s} {'τ_e':>6s} {'Δ [ns]':>7s} {'quasi-static':>12s} {'exact':>10s} {'exact/qs':>8s}")
edge_ratio = {}
for base in PRESETS:
    s = spec_of(base, "ref/aurora_guardslot", REF_POINT)
    p = s.offline["params"]
    tr = 1.0 / p["f_rep_hz"]
    modes = {"gauss": {"shape": "gaussian", "t_fwhm_s": p["tau_mode_s"]}}
    if base in OPO:
        modes["OPO"] = resolve_mode(next(x for x in s.offline["blocks"] if x["name"] == "stray1")["params"]["mode"])
    for mn, mode in modes.items():
        t, u = mode_time(mode)
        dt = t[1] - t[0]
        for tau_e in (1e-9, 3e-9, 1e-8):
            t_ref = 1e-9
            x = tr - t_ref / 2
            qs = math.exp(-x / tau_e) * float(mode_weight(0.0, mode))       # block: P(Δ)/P_pk·|ũ(0)|²
            # field of the pulse that ends t_ref/2 after its centre at −Δ: √P ∝ e^{−(Δ + t − t_ref/2)/(2τ_e)}
            amp = np.where(tr + t > t_ref / 2, np.exp(-np.maximum(tr + t - t_ref / 2, 0) / (2 * tau_e)), 1.0)
            ex = abs(np.sum(amp * u) * dt) ** 2
            edge_ratio[(base, mn, tau_e)] = ex / qs
            print(f"{base:8s} {mn:>6s} {tau_e * 1e9:4.0f}ns {tr * 1e9:7.1f} {qs:12.3g} {ex:10.3g} {ex / qs:8.3g}")
# docstring: "an upper bound of the coherent overlap otherwise". False: for a symmetric (Gaussian) mode quasi-static is a
# LOWER bound (∫e^{−t/2τ}u ≥ ∫u by convexity); for the asymmetric OPO mode it is neither (0.88–3.8×)
assert all(v >= 1.0 - 1e-6 for (b_, m_, _), v in edge_ratio.items() if m_ == "gauss")
assert edge_ratio[("F_liu20", "gauss", 3e-9)] > 2.5
assert any(v < 0.95 for (b_, m_, _), v in edge_ratio.items() if m_ == "OPO")

# ------------------------------------------------------------------------------------------------ part 5
print("\n== 5. Placement bookkeeping: brute-force pulse timeline vs _events (nearest pulses of each kind)")


def timeline(K, nb, ng, tr, placement, k):
    """Arrival times at port k's GBS input. Emission order: port K−1 first (L_k grows with k, D39); a block = ng guard
    slots then nb data slots. Port k's fibre delay τ_k = k·(nb + ng)·tr (+ const). Frame 0 = calibration frame."""
    B = nb + ng
    tf = K * B * tr
    data = [f * tf + (K - 1) * B * tr + (ng + j) * tr for f in ((-1, 1) if placement == "cal_frames" else (0,))
            for j in range(nb)]
    own, xt = [], []
    for f in ((0,) if placement == "cal_frames" else (-1, 0, 1)):
        for kp in range(K):
            slots = range(ng, B) if placement == "cal_frames" else (ng - 1,)
            for i in slots:
                t = f * tf + (K - 1 - kp) * B * tr + i * tr + k * B * tr
                (own if kp == k else xt).append(t)
    return np.array(data), np.array(own), np.array(xt)


def nearest(ts, pulses):
    d = ts[:, None] - pulses[None, :]
    pos = np.where(d > 0, d, np.inf).min(1)
    neg = np.where(d < 0, -d, np.inf).min(1)
    return pos, neg


mism = []
for (K, nb, ng, tr, pl) in [(32, 4, 0, 50e-9, "cal_frames"), (32, 1, 0, 1 / 6e6, "cal_frames"),
                            (32, 4, 1, 50e-9, "guard_slot"), (32, 1, 1, 1 / 6e6, "guard_slot"), (32, 4, 2, 50e-9, "guard_slot")]:
    p = dict(t_rep_s=tr, N_b=nb, n_guard=ng, K=K, t_cal_s=1.0, placement=pl)
    for k in (0, 1, 15, 31):
        data, own, xt = timeline(K, nb, ng, tr, pl, k)
        if pl == "cal_frames":
            data_after = data[nb:]                                  # frame +1 (after the calibration frame)
            data_before = data[:nb]                                 # frame −1
            op, _ = nearest(data_after, own)
            _, on = nearest(data_before, own)
            xp, _ = nearest(data_after, xt)
            _, xn = nearest(data_before, xt)
        else:
            op, on = nearest(data, own)
            xp, xn = nearest(data, xt)
        for j in range(nb):
            ev = ws._events(p, j)
            e_own_p = min([d for kind, d, _ in ev if kind == "own" and d > 0], default=np.inf)
            e_own_n = min([-d for kind, d, _ in ev if kind == "own" and d < 0], default=np.inf)
            e_xt_p = min([d for kind, d, _ in ev if kind == "xt" and d > 0], default=np.inf)
            e_xt_n = min([-d for kind, d, _ in ev if kind == "xt" and d < 0], default=np.inf)
            for nm, code, true in (("own+", e_own_p, op[j]), ("own−", e_own_n, on[j]), ("xt+", e_xt_p, xp[j]),
                                   ("xt−", e_xt_n, xn[j])):
                if not math.isclose(code, true, rel_tol=1e-9) and not (np.isinf(code) and np.isinf(true)):
                    if pl == "guard_slot" and nm == "own−":
                        continue                                    # next frame's own pulse: code ignores it (≫)
                    mism.append((pl, nb, ng, k, j, nm, code / tr, true / tr))
print(f"{'placement':10s} {'N_b':>3s} {'ng':>2s} {'port':>4s} {'j':>2s} {'event':>5s} {'code/T_rep':>10s} {'true/T_rep':>10s}")
for m in mism:
    if m[3] in (0, 15, 31) and m[4] in (0, m[1] - 1):
        print(f"{m[0]:10s} {m[1]:3d} {m[2]:2d} {m[3]:4d} {m[4]:2d} {m[5]:>5s} {m[6]:10.3g} {m[7]:10.3g}")
# every finite mismatch is the code being closer (conservative): cal_frames xt = the two end ports' worst case applied
# to every port (documented as worst-port, WS-06).
assert all(m[6] <= m[7] + 1e-9 for m in mism if np.isfinite(m[6]))
# fix round (WS-06 / WS-05): guard_slot now matches the brute force exactly for every port and n_guard (was: xt− at
# N_b − j for n_guard = 2, and xt+ omitted); the only mismatches left are the documented cal_frames worst-port positions
assert not any(m[0] == "guard_slot" for m in mism)
assert all(m[0] == "cal_frames" and m[5] in ("xt+", "xt−") for m in mism)
# afterpulse undercount from omitted leaked pulses before the slot (guard_slot, every port's leaked guard pulse clicks
# when ε_dmx·Φ_pk·t_ref·QE ≳ 1): Σ_m e^{−m B T_rep/τ_ap} relative to the own pulse alone
for base, tr, B in (("F_liu20", 50e-9, 5), ("F_bor6", 1 / 6e6, 2)):
    for tau_ap in (1e-7, 1e-6):
        f = sum(math.exp(-m * B * tr / tau_ap) for m in range(32))
        print(f"  afterpulse undercount {base} τ_ap = {tau_ap:.0e}: ×{f:.2f}")

# ------------------------------------------------------------------------------------------------ part 6
print("\n== 6. Guard slot: L_k factor and the D78 Δν_max (σ = 10 mrad, R = 1, worst port)")
print(f"{'preset':8s} {'lk_factor':>9s} {'Δν_max ref':>10s} {'without':>8s} {'blind (TES 10 µs)':>17s} {'blind (43 ns)':>13s}")
for base in PRESETS:
    out = []
    for rec in (1e-5, 4.3e-8):
        s = spec_of(base, "ref/aurora_guardslot", {**REF_POINT, "det_recovery_s": rec})
        b = blk(s, "frefstray1")
        out.append(b.check())
    r = out[0]
    p = s.offline["params"]
    nb = p["N_b"]
    assert math.isclose(r["lk_factor"], (nb + 1) / nb)
    assert math.isclose(r["dnu_max_hz"], 1e-4 / (2 * math.pi * 31 * (nb + 1) / p["f_rep_hz"]), rel_tol=1e-9)
    print(f"{base:8s} {r['lk_factor']:9.3f} {r['dnu_max_hz']:10.3g} {r['dnu_max_without_ref_hz']:8.3g} "
          f"{r['blind_prob_mean']:17.3g} {out[1]['blind_prob_mean']:13.3g}")
    assert r["blind_prob_mean"] > 0.999                                    # D88-3: TES → every data slot blind
# ------------------------------------------------------------------------------------------------ part 7 (fix round)
print("\n== 7. Fix round: OPO-mode overlays (WS-01) and the calibration-frame slow gate ε_gate (WS-02)")
print("-- P1 cal_frames leak, P_pk = 1 µW, ε_mod = −50 dB; OPO presets with ref/opo_mode_*, waveguide presets Gaussian")
print(f"{'preset':8s} {'mode':>6s} " + " ".join(f"{'ε_gate=' + format(g, '.0e'):>13s}" for g in (1.0, 1e-3, 1e-5, 1e-9))
      + f" {'εε·P max [W]':>12s}")
for base in PRESETS:
    post = (f"ref/opo_mode_{OPO_OVL[base]}",) if base in OPO else ()
    vals = []
    for g in (1.0, 1e-3, 1e-5, 1e-9):
        s = spec_of(base, "ref/aurora_calframe", REF_POINT, {"fref_eps_gate": g}, post=post)
        vals.append(blk(s, "frefstray1").check()["photons_per_mode_mean"])
    lk = leak[base][2]                                                 # part 1: Φ_leak·τ_eff of the right mode
    print(f"{base:8s} {'OPO' if base in OPO else 'gauss':>6s} " + " ".join(f"{v:13.3g}" for v in vals)
          + f" {leak[base][4]:12.3g}")
    assert math.isclose(vals[0], lk, rel_tol=0.05)                     # ε_gate = 1: the leak of the right mode dominates
    assert vals[-1] < 1e-3 * vals[0]                                   # 90 dB gate: the leak is gone
print("-- F_liu20 + probe/freq_80MHz + probe/power_1uW + stray/window_tail (T = 1 ms, D = 0.65), vs stray_n_mode_max;"
      " CHK on the worst slot (statistic max)")
print(f"{'mode':>6s} {'R_refl':>7s} {'on':>10s} {'mean':>10s} {'max':>10s} {'limit':>9s} {'CHK max':>7s}")
for post, lab in (((), "gauss"), (("stray/opo_mode_liu20",), "OPO")):
    for rr in (0.0, 1e-3):
        s = spec_of("F_liu20", "stray/window_tail", {**WIN_POINT, "win_refl_r": rr, "win_limit_statistic": "max"},
                    pre=("probe/freq_80MHz", "probe/power_1uW"), post=post)
        b = blk(s, "winstray1")
        r = b.check()
        lim = s.offline["params"]["stray_n_mode_max"]
        print(f"{lab:>6s} {rr:7.0e} {r['photons_per_mode_on']:10.3g} {r['photons_per_mode_mean']:10.3g} "
              f"{r['photons_per_mode_max']:10.3g} {lim:9.2g} {'warn' if b.chk_log.events else 'pass':>7s}")
        if lab == "OPO" and rr == 0.0:
            on1 = stray_photons(blk(s, "stray1").params, blk(s, "stray1")._toggles, s.offline["params"]["p_pr_avg_w"])
            assert math.isclose(r["photons_per_mode_on"], on1["photons_per_mode"], rel_tol=1e-9)   # = stray1 (WS-01)
            # the review's "≈ 2e−3, 4× the limit" is the on value; the worst slot (just after switch-off) also exceeds
            # it and the CHK fires; the window mean (×0.017) stays below
            assert r["photons_per_mode_on"] > 3 * lim and r["photons_per_mode_max"] > lim and b.chk_log.events
            assert r["photons_per_mode_mean"] < lim
        if lab == "gauss" and rr == 0.0:
            assert r["photons_per_mode_max"] < 1e-2 * lim and not b.chk_log.events             # the old false pass
print("\nall self-checks passed")
