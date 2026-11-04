"""Phase 6 axis 1 (d): in-lock identification of the 914 → 915B transfer on refD_lab (ident/wire_refD + a stepped sine,
D93), as a rehearsal of the bench measurement (DSC1 locked, excitation summed into the error inside the DSC1 / at its
input, L̂ = −resp/ref, D93-2).

For each "truth" (5.4 kHz side mode Q 15 / r 0.11 or Q 5 / r 0.34; D96-8: both match the measured |S|²) one locked run
per seed sweeps 1–50 kHz; the 915B response is then read as Ĥ_915B = L̂ / L_known, L_known = chk.lockD_stability's
open loop with the mode list removed (DSC1 law, h, delay, poles, K_d: the parts that do not depend on the 915B), and
compared with both mode models. The two truths are told apart when |Ĥ_q15 − Ĥ_q5| exceeds the seed spread of Ĥ.

Usage: uv run python reviews/scripts/p6_axis1_refD_ident.py [n_seeds] [amplitude_V]
Writes docs/figures/phase6/axis1/d_refD_ident_915b.png and d_refD_ident_summary.json.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import p6_axis1_lib as L  # noqa: E402

MODELS = {"q15": {"offline.params.mm2_q": 15.0, "offline.params.mm2_r": 0.11},
          "q5": {"offline.params.mm2_q": 5.0, "offline.params.mm2_r": 0.34}}
F_LO, F_HI, NF = 1.0e3, 5.0e4, 41


def one(args):
    import warnings

    from cranberry.runner.run import offline_block_specs, run

    truth, seed, amp = args
    # the ident keys go in a last overlay (not sweep axes): run.duration_s (ident_duration_s of n_freqs) is resolved
    # at load time and expand() keeps the resolved run: section
    ident = {"offline": {"params": {"ident_f_lo_hz": F_LO, "ident_f_hi_hz": F_HI, "ident_n_freqs": NF,
                                    "ident_amplitude": amp}}}
    s = L.point(**L.d_point(seed=seed, duration_s=None, axes=MODELS[truth],
                            overlays=["ident/wire_refD", "ident/in_lock_stepped"]), extra=ident)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run(s)
    m = r.measurements["loop_identification"]["lockD"]
    params = {}
    for name in MODELS:
        s2 = L.point(**L.d_point(seed=seed, duration_s=0.01, axes=MODELS[name]))
        (b,) = [x for x in offline_block_specs(s2) if x.name == "chk.lockD_stability"]
        params[name] = dict(b.params)
    return {"truth": truth, "seed": seed, "f": np.asarray(m["f_hz"]).tolist(),
            "l_hat": [np.real(m["l_hat"]).tolist(), np.imag(m["l_hat"]).tolist()],
            "coh": np.asarray(m["coherence"]).tolist(), "f_u": m.get("f_u_hz"), "pm": m.get("phase_margin_deg"),
            "f180": m.get("f_180_hz"), "gm": m.get("gain_margin_db"), "duration_s": s.run.duration_s,
            "valid": bool(r.valid), "params": params}


def main(n_seeds=5, amp=0.1):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from cranberry.blocks.chk.loop_stability import open_loop

    jobs = [(t, s, amp) for t in MODELS for s in range(1, n_seeds + 1)]
    with ProcessPoolExecutor(max_workers=min(L.N_PROC, len(jobs))) as ex:
        rows = list(ex.map(one, jobs))
    f = np.asarray(rows[0]["f"])
    fd = np.geomspace(F_LO, F_HI, 2000)
    p0 = rows[0]["params"]
    known = {n: open_loop(f, {**p, "modes": []}) for n, p in p0.items()}   # identical for both (modes only differ)
    model_h = {n: open_loop(fd, p) / open_loop(fd, {**p, "modes": []}) for n, p in p0.items()}
    model_h_f = {n: open_loop(f, p) / known[n] for n, p in p0.items()}
    summ = {"amplitude_V": amp, "n_seeds": n_seeds, "duration_s": rows[0]["duration_s"], "f_hz": f.tolist()}
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), layout="constrained", sharex=True)
    for j, truth in enumerate(MODELS):
        rr = [r for r in rows if r["truth"] == truth]
        lh = np.array([np.array(r["l_hat"][0]) + 1j * np.array(r["l_hat"][1]) for r in rr])
        h = lh / known[truth]
        hm, hs = h.mean(axis=0), h.std(axis=0, ddof=1)
        # distinguishability: |mean Ĥ − model| / spread for both models, in 4–7 kHz
        band = (f > 4.0e3) & (f < 7.0e3)
        dist = {n: float(np.sqrt(np.mean(np.abs(hm[band] - model_h_f[n][band]) ** 2 / (np.abs(hs[band]) ** 2 + 1e-30))))
                for n in MODELS}
        summ[truth] = {"f_u": [r["f_u"] for r in rr], "pm": [r["pm"] for r in rr], "f180": [r["f180"] for r in rr],
                       "gm": [r["gm"] for r in rr], "coh_min": float(np.min([r["coh"] for r in rr])),
                       "valid": all(r["valid"] for r in rr), "rms_z_vs_model_4_7kHz": dist,
                       "h_mag_mean": np.abs(hm).tolist(), "h_mag_std": np.abs(hs).tolist()}
        ax0, ax1 = axs[0, j], axs[1, j]
        for n, c in zip(MODELS, ("C0", "C3")):
            ax0.loglog(fd, np.abs(model_h[n]), color=c, lw=1, label=f"model {n}")
            ax1.semilogx(fd, np.degrees(np.unwrap(np.angle(model_h[n]))), color=c, lw=1)
        ax0.errorbar(f, np.abs(hm), yerr=np.abs(hs), fmt="k.", ms=3, capsize=1.5, label=f"in-lock Ĥ ({n_seeds} seeds)")
        ax1.errorbar(f, np.degrees(np.unwrap(np.angle(hm))), yerr=np.degrees(np.abs(hs) / np.abs(hm)), fmt="k.", ms=3,
                     capsize=1.5)
        ax0.set_title(f"truth {truth}: Ĥ_915B = L̂ / L_known (amp {amp} V, {rows[0]['duration_s']:.2f} s per sweep)")
        ax0.set_ylabel("|Ĥ_915B| (rel. to DC)")
        ax1.set_ylabel("phase [deg]")
        ax1.set_xlabel("frequency [Hz]")
        ax0.legend(fontsize=8)
        for a in (ax0, ax1):
            a.grid(True, which="both", alpha=0.25)
        print(f"truth {truth}: f_u {np.mean(summ[truth]['f_u']):.0f} Hz, PM {np.mean(summ[truth]['pm']):.1f}°, "
              f"f180 {np.mean(summ[truth]['f180']):.0f} Hz, GM {np.mean(summ[truth]['gm']):.2f} dB, coherence ≥ "
              f"{summ[truth]['coh_min']:.3f}; rms z (4–7 kHz) vs models {dist}", flush=True)
    L.FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(L.FIG / "d_refD_ident_915b.png", dpi=130)
    (L.FIG / "d_refD_ident_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    print("saved")


if __name__ == "__main__":
    main(*(int(sys.argv[1]),) if len(sys.argv) > 1 else (), *(float(sys.argv[2]),) if len(sys.argv) > 2 else ())
