"""Phase 5 M1 framework review: D57-7 (linear fractional-delay interpolation in segment_phase /
beam_assembler) quantified.

Closed form (any stationary driver x with lag-1 autocorrelation ρ1 at the block rate):
    linear interpolation y = (1-μ)x[n-k] + μ x[n-k-1]  has  |H|² = 1 - 2μ(1-μ)(1 - cos ωh)
    ⇒ Var(y)/Var(x) = 1 - 2μ(1-μ)(1 - ρ1)                     (ρ1 = 0: D57-7's 2μ(1-μ))
The exact (band-limited) delay keeps the variance. Part 1 checks the closed form through the real
segment_phase → beam_assembler blocks; part 2 tabulates it versus f_sim; part 3 looks at what F6 is
for — the *difference* of two beams on a shared segment (counter-propagation, VT-37 V6) — in the
sub-sample regime of the reference presets, for linear and thiran1; part 4 confirms part 3 in the
time domain with the real blocks.

Run: uv run python reviews/scripts/f6_fracdelay.py
"""

import numpy as np

import cranberry.path.testing as pt
from cranberry.core.delay import DelayDesign
from cranberry.core.fft import psd_onesided
from cranberry.path.topology import C_M_PER_S

N_G = 1.468


def h_lin(d_samples, w):
    k = np.floor(d_samples)
    mu = d_samples - k
    z = np.exp(-1j * w)
    return z**k * ((1 - mu) + mu * z)


def h_thiran(d_samples, w):
    return DelayDesign(float(d_samples), "thiran1").freq_response(w / (2 * np.pi), 1.0)


def h_exact(d_samples, w):
    return np.exp(-1j * w * d_samples)


# ------------------------------------------------------------------------------------------- part 1
def part1(fs=8e6, L=1.0, n=2**18):
    """White and AR(1) driver on a lumped one-subsegment segment (the last one on the path, so only
    segment_phase interpolates: μ = (δ/2)·fs)."""
    tau_half = N_G * L / C_M_PER_S / 2
    mu = tau_half * fs
    segs = [pt.seg_spec("s", L, drivers={"strain": "uniform"})]
    beams = [pt.beam_spec("b", [pt.step_spec("s", 1)], "det")]
    topo = pt.topo_spec(segs, beams, ["det"], pt.coef_spec("s", "strain", "*", 1.0))
    rng = np.random.default_rng(57)
    out = []
    for rho in (0.0, 0.5, 0.9):
        w = rng.standard_normal(n)
        x = np.empty(n)
        x[0] = w[0]
        for i in range(1, n):  # AR(1), lag-1 autocorrelation rho
            x[i] = rho * x[i - 1] + np.sqrt(1 - rho**2) * w[i]
        pipe = [pt.pathtest("pathtest_table", "drv", "drv.s.strain", (), "1", n_rows=n),
                pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
                *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.b.phi"], n)]
        rec, _, _ = pt.run_and_record(pipe, fs, n, {"drv.s.strain_table": x}, compiled=True)
        y = rec["beam.b.phi"][10:]
        out.append((rho, np.var(y) / np.var(x[10:]), 1 - 2 * mu * (1 - mu) * (1 - rho)))
    print(f"part 1: L = {L} m, f_sim = {fs/1e6:g} MHz, μ = {mu:.4f}")
    for rho, meas, pred in out:
        print(f"  ρ1 = {rho:.1f}: Var ratio measured {meas:.5f}, closed form {pred:.5f}")


# ------------------------------------------------------------------------------------------- part 2
def part2():
    print("part 2: variance deficit 2μ(1-μ)(1-ρ1) of one interpolation stage [%]")
    print("  driver: white (ρ1=0) / OU with corner f_c (ρ1 = exp(-2π f_c / f_sim))")
    hdr = "  delay      f_sim  |  white   OU 1 MHz  OU 100 kHz  OU 10 kHz"
    print(hdr)
    for label, tau in (("1 m  δ/2", N_G * 1.0 / C_M_PER_S / 2), ("63 m δ/2", N_G * 63.0 / C_M_PER_S / 2)):
        for fs in (1e6, 2e6, 4e6, 8e6, 20e6, 40e6):
            mu = (tau * fs) % 1.0
            row = []
            for fc in (None, 1e6, 1e5, 1e4):
                rho = 0.0 if fc is None else np.exp(-2 * np.pi * fc / fs)
                row.append(100 * 2 * mu * (1 - mu) * (1 - rho))
            print(f"  {label}  {fs/1e6:5g} MHz |  " + "  ".join(f"{v:7.3f}" for v in row))


# ------------------------------------------------------------------------------------------- part 3
def counterprop_spectra(L, N, fs, interp, f):
    """S_D(f) and S_fwd(f) for N independent unit-variance white subsegment drivers (per_sub, weight
    1/N): forward leaves at x = L, backward at x = 0 (F6 timing: centre of subsegment j passed
    (N - j - 1/2)δ before the exit)."""
    d = N_G * L / C_M_PER_S / N * fs  # δ in samples
    w = 2 * np.pi * f / fs
    sd = np.zeros_like(f)
    sf = np.zeros_like(f)
    for m in range(N):
        df = (N - m - 0.5) * d  # forward passes sub m at j = m
        db = (m + 0.5) * d      # backward passes sub m at j = N-1-m
        hf, hb = interp(df, w), interp(db, w)
        sd += np.abs(hf - hb) ** 2
        sf += np.abs(hf) ** 2
    return sd / N**2, sf / N**2


def part3(L=63.0, N=20):
    tau = N_G * L / C_M_PER_S
    print(f"part 3: counter-propagation on a shared distributed segment, L = {L} m (τ = {tau*1e9:.0f} ns), "
          f"N_sub = {N}; ratio S_D(interp)/S_D(exact)")
    for fs in (1e6, 4e6, 20e6):
        fn = fs / 2
        fr = np.array([0.1, 0.25, 0.5, 0.75, 0.95]) * fn
        ex = counterprop_spectra(L, N, fs, h_exact, fr)[0]
        li = counterprop_spectra(L, N, fs, h_lin, fr)[0]
        th = counterprop_spectra(L, N, fs, h_thiran, fr)[0]
        ff = np.linspace(1, fn, 4000)
        exw, liw, thw = (counterprop_spectra(L, N, fs, h, ff)[0] for h in (h_exact, h_lin, h_thiran))
        rw = [np.trapezoid(a, ff) / np.trapezoid(exw, ff) for a in (liw, thw)]
        r2 = [np.trapezoid(a / ff**2, ff) / np.trapezoid(exw / ff**2, ff) for a in (liw, thw)]
        print(f"  f_sim = {fs/1e6:g} MHz (f_N = {fn/1e6:g} MHz)")
        print("    f/f_N       : " + "  ".join(f"{x:6.2f}" for x in fr / fn))
        print("    linear      : " + "  ".join(f"{x:6.3f}" for x in li / ex))
        print("    thiran1     : " + "  ".join(f"{x:6.3f}" for x in th / ex))
        print(f"    in-band variance of D, white driver: linear {rw[0]:.3f}, thiran1 {rw[1]:.3f}; "
              f"1/f² driver: linear {r2[0]:.3f}, thiran1 {r2[1]:.3f}")


# ------------------------------------------------------------------------------------------- part 4
def part4(L=63.0, N=20, fs=4e6, n=2**19):
    tau = N_G * L / C_M_PER_S
    segs = [pt.seg_spec("c", L, n_sub=N, spatial="distributed", drivers={"strain": "per_sub"})]
    beams = [pt.beam_spec("fwd", [pt.step_spec("c", 1)], "xL"), pt.beam_spec("bwd", [pt.step_spec("c", -1)], "x0")]
    topo = pt.topo_spec(segs, beams, ["xL", "x0"], pt.coef_spec("c", "strain", "*", 1.0))
    pipe = [pt.pathtest("pathtest_white", "noise", "drv.c.strain", (N,), "1", sigma=1.0),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.fwd.phi", "beam.bwd.phi"], n)]
    rec, _, _ = pt.run_and_record(pipe, fs, n, compiled=True, seed=37)
    fwd, bwd = rec["beam.fwd.phi"][100:], rec["beam.bwd.phi"][100:]
    f, s_d = psd_onesided(fwd - bwd, fs, nperseg=2**12)
    print(f"part 4: time domain, real blocks, f_sim = {fs/1e6:g} MHz, L = {L} m, N_sub = {N} (white per_sub)")
    print("    f/f_N   S_D measured / exact-delay prediction   (linear-interp prediction / exact)")
    for x in (0.25, 0.5, 0.75, 0.95):
        k = int(round(x * fs / 2 / (f[1] - f[0])))
        band = slice(k - 4, k + 5)
        ex = counterprop_spectra(L, N, fs, h_exact, f[band])[0] * 2 / fs  # one-sided PSD of unit white
        li = counterprop_spectra(L, N, fs, h_lin, f[band])[0] * 2 / fs
        print(f"    {x:4.2f}    {np.mean(s_d[band]) / np.mean(ex):.3f}        "
              f"                       ({np.mean(li) / np.mean(ex):.3f})")
    print(f"    (VT-37 V6 runs at f_sim = 1 GHz with integer subsegment delays; τ = {tau*1e9:.0f} ns)")


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
