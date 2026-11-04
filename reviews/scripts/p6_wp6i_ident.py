"""Phase 6 WP6-I (D85-7): validation of the in-lock identification (excitation + loop_identification).

Part 1  refD_nominal + ident/wire_refD + ident/in_lock_stepped (25 steps, 300 Hz – 100 kHz, one locked run):
        L̂ against (a) the loop_stability_check model (D82 dsc1) and (b) the exact multirate z-domain loop gain at the
        PID instants (tests/integration/refD_support.loop_parts: L_d = −g C (B A + B′A′); the loop is LTI at 500 kHz
        there), with the noise sources on (0.1 V) and off (0.02 V). Stats over 1–30 kHz at coherence ≥ 0.9.
Part 2  amplitude trade-off on refD_nominal (12 steps 1–30 kHz, corner 2 kHz): noise off isolates the fringe-sine bias
        (describing function 2J₁(A)/A ≈ 1 − A²/8, A = phase amplitude g_k·a·|T|/k_d), noise on adds the random error
        (falls as 1/a); 3 seeds per amplitude.
Part 3  multisine on refD_nominal (same power as the stepped sine).
Part 4  F_wg6_nominal: ident/wire_F (full stack at the fast controller's output) and ident/p7_literal.
Part 5  refinement sweep 14–19.5 kHz on refD_nominal: f_180 and the gain margin against the exact loop.
Usage: uv run python reviews/scripts/p6_wp6i_ident.py
"""

import dataclasses
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml
from scipy.special import j1

from cranberry import config
from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refD_support", ROOT / "tests/integration/refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)
NOISE_OFF = {**rs.NOISE_OFF, "pdb450c.per_crosscoupling": False, "pdb450c.visibility_hf_averaging": False}
TMP = Path(tempfile.mkdtemp())


def point(values):
    p = TMP / f"pt{len(list(TMP.iterdir()))}.yaml"
    p.write_text(yaml.safe_dump({"offline": {"params": values}}))
    return str(p)


def load(preset, overlays, params=None, noise=True, seed=None):
    spec = config.load(preset, overlays=[*overlays, *([point(params)] if params else [])])
    if not noise:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **NOISE_OFF})
    if seed is not None:
        spec = dataclasses.replace(spec, run=dataclasses.replace(spec.run, seed=seed))
    return spec


def l_exact(spec, f):
    fr = rs.fringe(spec)
    k_d = -fr["sign"] * fr["v0"] * fr["vis"]
    a, b, c, g = rs.loop_parts(spec, f, k_d=k_d)
    a2, b2, _, _ = rs.loop_parts(spec, f + spec.run.f_sim_hz / 2, k_d=k_d)
    return -g * c * (b * a + b2 * a2), abs(k_d)


def stats(m, lz, lo=1.0e3, hi=3.0e4, cmin=0.9):
    f = m["f_hz"]
    sel = (f >= lo) & (f <= hi) & (m["coherence"] >= cmin)
    ez = m["l_hat"][sel] / lz[sel]
    return {"n": int(sel.sum()), "n_band": int(((f >= lo) & (f <= hi)).sum()),
            "model_mag_pct": float(100 * np.max(np.abs(m["mag_err"][sel]))),
            "model_phase_deg": float(np.max(np.abs(m["phase_err_deg"][sel]))),
            "exact_mag_pct": float(100 * np.max(np.abs(np.abs(ez) - 1))),
            "exact_mag_mean_pct": float(100 * np.mean(np.abs(ez) - 1)),
            "exact_phase_deg": float(np.max(np.abs(np.degrees(np.angle(ez))))),
            "coh_min": float(np.min(m["coherence"][(f >= lo) & (f <= hi)]))}


MARG = ("f_u_hz", "phase_margin_deg", "f_180_hz", "gain_margin_db")


def part1():
    out = {}
    for noise, amp in ((True, 0.1), (False, 0.02)):
        spec = load(rs.NOMINAL, ["ident/wire_refD", "ident/in_lock_stepped"], {"ident_amplitude": amp}, noise)
        r = run(spec)
        m = r.measurements["loop_identification"]["lockD"]
        lz, k_d = l_exact(spec, m["f_hz"])
        key = "noise_on" if noise else "noise_off"
        err = np.asarray(r.series["lock.lockD.err"])[r.t_s >= 0.02]
        d = np.angle(np.exp(1j * (err - np.mean(err))))
        out[key] = {"amplitude_v": amp, **stats(m, lz), **{k: m[k] for k in MARG},
                    **{f"model_{k}": m[f"model_{k}"] for k in MARG}, "sigma_phase_rad": float(np.std(d)),
                    "cos_factor": float(np.mean(np.cos(d)))}
        z = {"f_hz": m["f_hz"].tolist(), "l_mag": m["l_mag"].tolist(), "l_model_mag": m["l_model_mag"].tolist(),
             "mag_err_pct": (100 * m["mag_err"]).tolist(), "phase_err_deg": m["phase_err_deg"].tolist(),
             "exact_err_pct": (100 * (np.abs(m["l_hat"] / lz) - 1)).tolist(),
             "exact_phase_deg": np.degrees(np.angle(m["l_hat"] / lz)).tolist(), "coherence": m["coherence"].tolist()}
        out[key]["table"] = z
        lzm = l_exact(spec, np.geomspace(300.0, 1e5, 4000))[0]
        from cranberry.blocks.meas.loop_identification import loop_margins

        out[key]["exact_margins"] = loop_margins(np.geomspace(300.0, 1e5, 4000), lzm)
    return out


def part2():
    rows = []
    for noise in (False, True):
        for amp in (0.01, 0.03, 0.1, 0.3, 1.0):
            per_seed = []
            for seed in ((0,) if not noise else (1, 2, 3)):
                spec = load(rs.NOMINAL, ["ident/wire_refD", "ident/in_lock_stepped"],
                            {"ident_amplitude": amp, "ident_n_freqs": 12, "ident_f_lo_hz": 1.0e3,
                             "ident_f_hi_hz": 3.0e4}, noise, seed=None if not noise else 20261031 + seed)
                m = run(spec).measurements["loop_identification"]["lockD"]
                lz, k_d = l_exact(spec, m["f_hz"])
                g = np.maximum(1.0, m["f_hz"] / 2.0e3)
                a_phi = g * amp * np.abs(lz / (1 + lz)) / k_d            # phase amplitude at the fringe
                pred = 2 * j1(a_phi) / a_phi - 1                          # sine describing-function bias
                per_seed.append({**stats(m, lz, cmin=0.0), "max_phase_amp_rad": float(a_phi.max()),
                                 "pred_bias_pct_at_max": float(100 * pred[np.argmax(a_phi)]),
                                 "bias_pct_at_max": float(100 * (np.abs(m["l_hat"] / lz) - 1)[np.argmax(a_phi)])})
            row = {"noise": noise, "amplitude_v": amp}
            for k in per_seed[0]:
                row[k] = float(np.mean([p[k] for p in per_seed]))
            rows.append(row)
            print(json.dumps(row), flush=True)
    return rows


def part3():
    spec = load(rs.NOMINAL, ["ident/wire_refD", "ident/in_lock_multisine"], {"ident_amplitude": 0.3})
    m = run(spec).measurements["loop_identification"]["lockD"]
    lz, _ = l_exact(spec, m["f_hz"])
    return {"amplitude_v": 0.3, "duration_s": spec.run.duration_s, **stats(m, lz), **{k: m[k] for k in MARG}}


def part4():
    out = {}
    for wire in ("wire_F", "p7_literal"):
        for noise in (True,):
            spec = load("presets/future/F_wg6_nominal.yaml", [f"ident/{wire}", "ident/in_lock_stepped"])
            m = run(spec).measurements["loop_identification"]["lock0"]
            f = m["f_hz"]
            sel = (f >= 1.0e3) & (f <= 2.0e4) & (m["coherence"] >= 0.9)
            out[wire] = {"n": int(sel.sum()), "model_mag_pct": float(100 * np.max(np.abs(m["mag_err"][sel]))),
                         "model_phase_deg": float(np.max(np.abs(m["phase_err_deg"][sel]))),
                         **{k: m[k] for k in MARG}, **{f"model_{k}": m[f"model_{k}"] for k in MARG}}
    return out


def part5():
    """Refinement sweep around the −180° crossing (the 25-point log sweep interpolates across the 915B mode)."""
    from cranberry.blocks.meas.loop_identification import loop_margins

    out = {}
    for noise, amp in ((True, 0.1), (False, 0.02)):
        spec = load(rs.NOMINAL, ["ident/wire_refD", "ident/in_lock_stepped"],
                    {"ident_amplitude": amp, "ident_n_freqs": 16, "ident_f_lo_hz": 1.4e4, "ident_f_hi_hz": 1.95e4,
                     "ident_amp_f_corner_hz": 2.0e3}, noise)
        m = run(spec).measurements["loop_identification"]["lockD"]
        lz, _ = l_exact(spec, m["f_hz"])
        fd = np.geomspace(1.4e4, 1.95e4, 4000)
        out["noise_on" if noise else "noise_off"] = {
            "f_180_hz": m["f_180_hz"], "gain_margin_db": m["gain_margin_db"],
            "exact": loop_margins(fd, l_exact(spec, fd)[0]), "model_f_180_hz": m["model_f_180_hz"],
            "model_gain_margin_db": m["model_gain_margin_db"], "coh_min": float(m["coherence"].min()),
            "exact_mag_pct": float(100 * np.max(np.abs(np.abs(m["l_hat"] / lz) - 1)))}
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["1", "2", "3", "4"]
    res = {}
    for p in parts:
        res[f"part{p}"] = globals()[f"part{p}"]()
        short = {k: v for k, v in res[f"part{p}"].items()} if isinstance(res[f"part{p}"], dict) else res[f"part{p}"]
        print(f"== part {p} ==")
        print(json.dumps(short, indent=1, default=float)[:6000], flush=True)
