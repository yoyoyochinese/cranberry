"""Ref D: frozen blind prediction (tag refD-blind-freeze) vs. the measured run in refs/data (PLAN §6.5, D04).

No parameter is changed: the measured PDB450C output (DSC1 input) is analysed with the frozen prediction's own
log binning and peak detector (``run_prediction.log_bin`` / ``detect_peaks``) and set against the frozen
``results/``. Measured files used: ``spectrum.csv`` (Welch, Hann, 0.2 s, 50 % overlap → 14 segments) and
``raw.npy`` scaled by the metadata ``waveform_scaling`` (= ``waveform.csv`` to 1e-11 V) for the PDB output;
``ch1.npz`` (yzero 0) for the DSC1 output (CV), identified by its coherence with the PDB output and its
transfer function (the D53 PID; to be confirmed by the user).

    uv run --extra plot python predictions/refD/comparison/compare.py
"""

from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRED = HERE.parent
ROOT = PRED.parents[1]
sys.path.insert(0, str(PRED))
import run_prediction as rp  # noqa: E402  (frozen analysis functions)

DATA = sorted(d + "/" for d in glob.glob(str(ROOT / "refs/data/*/rounds/*")))
NSEG_MEAS = 14  # (1.5 s − 0.2 s) / 0.1 s + 1
BANDS = [(1.0, 10.0), (10.0, 100.0), (100.0, 1.0e3), (1.0e3, 1.0e4), (1.0e4, 1.0e5), (1.0e5, 2.0e5), (0.5, 2.5e5)]
K_D_PRED = 2.86  # V/rad at the zero crossing in refD_nominal (PREDICTION.md §3)


def load_round(d):
    m = json.load(open(d + "metadata.json"))["round"]
    s = np.genfromtxt(d + "spectrum.csv", delimiter=",", names=True)
    ws = m["waveform_scaling"]
    w = (np.load(d + "raw.npy") - ws["yoff"]) * ws["ymult"] + ws["yzero"]
    z = np.load(d + "ch1.npz")
    cv = (z["raw"] - float(z["yoff"])) * float(z["ymult"]) + float(z["yzero"])
    st = m["observed_state"]
    return {"seq": m["sequence"], "setpoint_v": st["setpoint_v"], "gain_enabled": st["gain_enabled"],
            "pid": [st["p"], st["i"], st["d"]], "objective_v2": m["objective_v2"],
            "f": s["frequency_hz"], "p": s["psd_v2_per_hz"], "w": w, "cv": cv}


def dsc1_c(f):
    """DSC1 discrete PID (DSC1 p.19, D53 gains, D82): C = K_p + K_i h/2 + K_i/s + K_d (1 − e^{−sh})/h."""
    kp, i, d, h = 2.41e-4, 1.406e-3, 1.108e-2, 2.0e-6
    s = 2j * np.pi * f
    return kp + i / 2 + (i / h) / s + d * (1 - np.exp(-s * h))


def cv_analysis(r):
    """DSC1 output: scope numbers, PSD on the prediction's bins, and the CV / error transfer vs. −C (H1 and H2
    estimators bracket the true gain when both channels carry uncorrelated noise)."""
    from scipy.signal import csd, welch

    v, c = r["w"], r["cv"]
    f, pcc = welch(c, fs=1.0e6, nperseg=200000)
    out = {"scope": rp.scope_features(c, 1.0e6), "f": f, "p": pcc}
    if r["gain_enabled"]:
        f2, pvc = csd(v, c, fs=1.0e6, nperseg=100000)
        _, pvv = welch(v, fs=1.0e6, nperseg=100000)
        _, pc2 = welch(c, fs=1.0e6, nperseg=100000)
        m = (f2 >= 50.0) & (f2 <= 1.0e3)
        cm = np.abs(dsc1_c(f2[m]))
        out["gain_ratio_h1"] = float(np.median(np.abs(pvc[m] / pvv[m]) / cm))
        out["gain_ratio_h2"] = float(np.median(np.abs(pc2[m] / np.conj(pvc[m])) / cm))
        out["phase_deg_vs_minus_c"] = float(np.median(np.angle(pvc[m] / pvv[m] / -dsc1_c(f2[m]), deg=True)))
        out["coherence_50_1k"] = float(np.median(np.abs(pvc[m]) ** 2 / (pvv[m] * pc2[m])))
    return out


def main():
    rounds = [load_round(d) for d in DATA]
    s = json.load(open(PRED / "results/summary.json"))
    fb_pred = np.asarray(s["f_bin_hz"])
    g = s["groups"]
    nom = np.asarray(g["nominal_q1"]["spectra"]["v.pdb"]["median"])
    nom0 = np.asarray(g["nominal_q0"]["spectra"]["v.pdb"]["median"])
    lo, hi = (np.asarray(g["lhs_q1_stable"]["spectra"]["v.pdb"][k]) for k in ("p05", "p95"))
    ulo, uhi = (np.asarray(g["lhs_q1_unstable"]["spectra"]["v.pdb"][k]) for k in ("p05", "p95"))

    free = next(r for r in rounds if not r["gain_enabled"])
    v0 = 0.5 * (free["w"].max() - free["w"].min())  # fringe amplitude from the free-running sweep
    out = {"fringe_amplitude_v": v0, "k_d_pred_v_per_rad": K_D_PRED, "rounds": []}
    # sensitivity |S|² = in-loop phase PSD / free-running phase PSD (free: fringe-averaged slope² = V0²/2)
    fbf, pbf, _ = rp.log_bin(free["f"], free["p"])
    free_phase = np.interp(fb_pred, fbf[np.isfinite(pbf)], pbf[np.isfinite(pbf)]) / (v0**2 / 2.0)
    unl = np.asarray(g["unlocked"]["spectra"]["lock.lockD.err"]["median"])
    s2_model = nom / K_D_PRED**2 / unl
    out["f_bin_hz"] = fb_pred.tolist()
    out["s2_model_nominal"] = s2_model.tolist()
    for r in rounds:
        fb, pb, cnt = rp.log_bin(r["f"], r["p"])
        peaks = rp.detect_peaks(fb, pb, cnt, NSEG_MEAS)
        locked = r["gain_enabled"] and abs(r["setpoint_v"]) <= 0.9 * v0  # S = 2.5 V ≈ V0 sits on the fringe top
        kd = math.sqrt(max(v0**2 - r["setpoint_v"] ** 2, 0.0)) if locked else math.nan
        rec = {k: r[k] for k in ("seq", "setpoint_v", "gain_enabled", "pid", "objective_v2")}
        rec["lockable_setpoint"] = locked
        rec["k_d_v_per_rad"] = kd
        rec["band_rms_v"] = {f"{a:g}-{b:g}": rp.band_rms(r["f"], r["p"], a, b) for a, b in BANDS}
        rec["peaks"] = [{k: p[k] for k in ("f_hz", "excess_db", "fwhm_hz", "psd")} | {"rms_v": math.sqrt(p["power"])}
                        for p in peaks]
        rec["scope"] = rp.scope_features(r["w"], 1.0e6)
        cva = cv_analysis(r)
        fbc, pbc, _ = rp.log_bin(cva.pop("f"), cva.pop("p"))
        cva["psd_on_pred_bins"] = np.interp(fb_pred, fbc[np.isfinite(pbc)], pbc[np.isfinite(pbc)]).tolist()
        cva["band_rms_v"] = {f"{a:g}-{b:g}": rp.band_rms(fb_pred, np.asarray(cva["psd_on_pred_bins"]), a, b)
                             for a, b in BANDS[1:6]}
        rec["dsc1_out"] = cva
        # in-loop phase-equivalent PSD (v.pdb / K_d) vs. the prediction's v.pdb / K_D_PRED, on the prediction's bins
        if locked:
            pm = np.interp(fb_pred, fb[np.isfinite(pb)], pb[np.isfinite(pb)]) / kd**2
            m = (fb_pred >= 10.0) & (fb_pred <= 2.0e5)
            rec["phase_eq_ratio_to_nominal_median"] = float(np.median(pm[m] / (nom[m] / K_D_PRED**2)))
            rec["frac_bins_in_lhs_stable_band"] = float(np.mean((pm[m] >= lo[m] / K_D_PRED**2) & (pm[m] <= hi[m] / K_D_PRED**2)))
            rec["phase_eq_rms_rad_10_200k"] = rp.band_rms(fb_pred[m], pm[m], 10.0, 2.0e5)
            rec["s2_measured"] = (pm / free_phase).tolist()
        out["rounds"].append(rec)
    (HERE / "comparison.json").write_text(json.dumps(out, indent=1, default=float))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(3, 1, figsize=(9, 13), sharex=True)
    ax[0].fill_between(fb_pred, lo, hi, color="C0", alpha=0.25, label="pred. LHS stable 5–95 %")
    ax[0].fill_between(fb_pred, ulo, uhi, color="C3", alpha=0.10, label="pred. LHS unstable 5–95 %")
    ax[0].loglog(fb_pred, nom, "C0", lw=1.5, label="pred. nominal (DAC quant. on)")
    ax[0].loglog(fb_pred, nom0, "C0--", lw=1, label="pred. nominal (quant. off)")
    for r in rounds:
        if r["setpoint_v"] == 0.0:
            fb, pb, _ = rp.log_bin(r["f"], r["p"])
            ax[0].loglog(fb, pb, "k" if r["gain_enabled"] else "0.6", lw=1,
                         label=f"meas. round {r['seq']} ({'locked, S = 0' if r['gain_enabled'] else 'free-running'})")
    ax[0].set_ylabel("v.pdb PSD [V²/Hz]")
    ax[0].legend(fontsize=8)
    ax[1].fill_between(fb_pred, lo / K_D_PRED**2, hi / K_D_PRED**2, color="C0", alpha=0.25)
    ax[1].loglog(fb_pred, nom / K_D_PRED**2, "C0", lw=1.5, label="pred. nominal")
    for i, r in enumerate(x for x in rounds if x["gain_enabled"] and abs(x["setpoint_v"]) <= 0.9 * v0):
        fb, pb, _ = rp.log_bin(r["f"], r["p"])
        kd2 = v0**2 - r["setpoint_v"] ** 2
        ax[1].loglog(fb, pb / kd2, color=plt.cm.viridis(i / 7), lw=0.8, label=f"round {r['seq']}, S = {r['setpoint_v']:+g} V")
    ax[1].set_ylabel("in-loop phase-equivalent PSD [rad²/Hz]")
    ax[1].legend(fontsize=7, ncol=2)
    ax[2].loglog(fb_pred, s2_model, "C0", lw=1.5, label="pred. nominal (locked / free-running)")
    for i, rec in enumerate(x for x in out["rounds"] if "s2_measured" in x):
        ax[2].loglog(fb_pred, rec["s2_measured"], color=plt.cm.viridis(i / 7), lw=0.8,
                     label=f"meas. round {rec['seq']}, S = {rec['setpoint_v']:+g} V")
    ax[2].axhline(1.0, color="0.5", lw=0.5)
    ax[2].set_ylabel("sensitivity |S|² (phase)")
    ax[2].set_ylim(3e-3, 50)
    ax[2].set_xlabel("f [Hz]")
    ax[2].set_xlim(1, 2.5e5)
    ax[2].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "comparison.png", dpi=130)


if __name__ == "__main__":
    main()
