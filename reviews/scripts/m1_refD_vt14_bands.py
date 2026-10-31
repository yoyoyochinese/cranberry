"""Phase 5 M1 review (Ref D): VT-14 band ratios — Welch window bias vs statistics.

Re-runs the VT-14 free-running spec (feedback off) for several seeds and compares, in each VT-14 band,
the Welch estimate with (a) the PSD model averaged over the bins (what VT-14 does) and (b) the exact
expected Welch value E[P(f_k)] = ∫ S(ν)|G_k(ν)|² dν / (f_s Σw²) (Hann window after the per-segment
linear detrend, i.e. the effective taper G_k = P_detrend(w·e_k)). Also prints the statistical 1σ of a
band mean, to show what tolerance VT-14 actually needs.

Run: uv run python reviews/scripts/m1_refD_vt14_bands.py
"""

import dataclasses
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy import signal

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("vt14", ROOT / "tests" / "vt" / "test_vt_14_refD_free_running.py")
vt = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(vt)
rs = vt.rs

FS, NPERSEG, BANDS = vt.FS, vt.NPERSEG, vt.BANDS


def expected_welch(spec, f_bins, pad=2**22):
    """Exact E[Welch] at bins f_bins (Hann, linear detrend, one-sided density)."""
    n = NPERSEG
    w = signal.get_window("hann", n)
    t = np.arange(n)
    X = np.vstack([np.ones(n), t]).T
    XtX = X.T @ X
    nu = np.fft.rfftfreq(pad, 1 / FS)
    s = np.zeros_like(nu)
    s[1:] = vt._s_expected(spec, nu[1:])[0]
    s[0] = s[1]
    out = []
    for fk in f_bins:
        a = w * np.exp(-2j * np.pi * fk * t / FS)
        g = a - X @ np.linalg.solve(XtX, X.T @ a)              # linear-detrend projector applied to the taper
        # |Σ g[n] e^{+2πiνn/fs}|² for ν ≥ 0 and ν < 0 (one-sided S: both signs of ν carry S(|ν|)/2)
        gp = np.abs(np.fft.ifft(g, pad) * pad) ** 2           # e^{+i...}: ν on the FFT grid
        gm = np.abs(np.fft.fft(g, pad)) ** 2                  # ν → −ν
        half = pad // 2 + 1
        integrand = 0.5 * s * (gp[:half] + gm[:half])
        out.append(2 * np.trapezoid(integrand, nu) / (FS * np.sum(w**2)))
    return np.array(out)


def main():
    spec0, _ = None, None
    rows = []
    for seed in (20261031, 1, 2, 3, 4):
        spec = rs.load(run={"duration_s": vt.DUR, "seed": seed}, measure={"record": {"slots": ["lock.lockD.err"], "decim": 1}},
                       toggles={"dsc1.pid.feedback": False, "pdb450c.shot_noise": False, "pdb450c.electronic_noise": False})
        x = run(spec).series["lock.lockD.err"]
        f, p = signal.welch(x, fs=FS, window="hann", nperseg=NPERSEG, noverlap=NPERSEG // 2, detrend="linear")
        n_avg = (x.size - NPERSEG) // (NPERSEG // 2) + 1
        if spec0 is None:
            spec0 = spec
            exp_bins = {}
            for lo, hi in BANDS:
                m = (f >= lo) & (f < hi) & (f > 0)
                # exact kernel only where the f⁻² curvature matters (lowest two bands); else the model
                ew = expected_welch(spec, f[m]) if hi <= 2e3 else vt._s_expected(spec, f[m])[0]
                exp_bins[(lo, hi)] = (m, ew)
        r = []
        for lo, hi in BANDS:
            m, ew = exp_bins[(lo, hi)]
            model = vt._s_expected(spec, f[m])[0].mean()
            r.append((p[m].mean() / model, p[m].mean() / ew.mean(), ew.mean() / model, int(m.sum())))
        rows.append(r)
        print(f"seed {seed}: N_avg {n_avg}; per band (Welch/model, Welch/E[Welch], E[Welch]/model, bins):")
        for (lo, hi), q in zip(BANDS, r):
            print(f"   {lo:8.0f}-{hi:8.0f} Hz  {q[0]:.3f}  {q[1]:.3f}  {q[2]:.3f}  {q[3]}")
    a = np.array([[q[0] for q in r] for r in rows])
    b = np.array([[q[1] for q in r] for r in rows])
    print("mean over seeds  Welch/model:", np.round(a.mean(0), 3), " Welch/E[Welch]:", np.round(b.mean(0), 3))
    print("seed spread (std) Welch/model:", np.round(a.std(0, ddof=1), 3))
    n_avg = (int(vt.DUR * FS) - NPERSEG) // (NPERSEG // 2) + 1
    print(f"VT-14 tolerance 3/sqrt(N_avg) = {3 / math.sqrt(n_avg):.3f} (per bin); band-mean 1σ ≈ 1/sqrt(N_avg·bins/1.5)")


if __name__ == "__main__":
    main()
