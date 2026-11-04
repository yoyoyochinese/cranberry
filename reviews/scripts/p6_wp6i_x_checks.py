"""Phase 5 review of WP6-I (D93): checks not covered by the author's p6_wp6i_ident.py.

Part 1  Moku:Delta summing-point reduction L̂ = −H/(g + H) (hardware/moku.py): exactness for the true g, and the error of
        L̂ and of the margins when the control-matrix coefficient is off by 1 dB (the D93-3 "1 dB → 12 %" statement), or
        when the FRA output reaches the PID input with a small route latency (complex g e^{−sτ_r}). L = refD's
        chk.lockD_stability model on 300 Hz – 100 kHz (any L will do: the reduction is algebraic).
Part 2  p7_literal TEC hold loop: L_hold = ω_h/(s(1 + sτ_th)) (pid0s reads the error, K_i = 2π f_hold/(K k_TEC)),
        damping ζ = 1/(2√(ω_h τ_th)) and phase margin over F_common's f_12 × τ_th grid; hold time constant vs the
        0.65 s identification run.
Part 3  wire_F with controller/dsc1 (D93-4: "combinable with the controller overlays"): F_wg6_nominal + controller/dsc1
        + ident/wire_F + ident/in_lock_stepped, L̂ against chk.lock_stability (dsc1 law), 1–20 kHz, coherence ≥ 0.9.
Usage: uv run python reviews/scripts/p6_wp6i_x_checks.py [1 2 3]
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np

from cranberry import config
from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.blocks.meas.loop_identification import loop_margins
from cranberry.hardware.moku import reduce_summing_point
from cranberry.runner.run import offline_block_specs, run

ROOT = Path(__file__).resolve().parents[2]


def chk_params(spec, name):
    (b,) = [b for b in offline_block_specs(spec) if b.name == name]
    return dict(b.params)


def part1():
    spec = config.load(ROOT / "presets/reference/refD_nominal.yaml")
    mp = chk_params(spec, "chk.lockD_stability")
    f = np.geomspace(300.0, 1.0e5, 4000)
    lt = open_loop(f, mp)
    g = 1.0
    h = -lt * g / (1 + lt)                                     # FRA ratio e/x for e = −L(e + g x)
    assert np.allclose(reduce_summing_point(h, g)["l_hat"], lt, rtol=1e-12)
    m0 = loop_margins(f, lt)
    print("Moku reduction, refD model L (true g = 1):")
    print(f"  true margins: f_u {m0['f_u_hz']:.0f} Hz, PM {m0['phase_margin_deg']:.2f}°, f_180 {m0['f_180_hz']:.0f} Hz, "
          f"GM {m0['gain_margin_db']:.2f} dB")
    for what, g_true in (("coefficient +1 dB (g_true = 1.122, assumed 1)", 10 ** (1 / 20)),
                         ("coefficient −1 dB", 10 ** (-1 / 20)),
                         ("route latency 10 ns", np.exp(-2j * np.pi * f * 10e-9)),
                         ("route latency 100 ns", np.exp(-2j * np.pi * f * 100e-9))):
        h = -lt * g_true / (1 + lt)
        lh = reduce_summing_point(h, 1.0)["l_hat"]
        err = np.abs(lh / lt) - 1
        m = loop_margins(f, lh)
        rows = []
        for lmag in (7.0, 3.0, 1.0, 0.3, 0.1):
            i = int(np.argmin(np.abs(np.abs(lt) - lmag)))
            rows.append(f"|L| {abs(lt[i]):4.2f} ({f[i]:6.0f} Hz): {100 * err[i]:+6.1f} %")
        print(f"  {what}:")
        print("    " + "; ".join(rows))
        print(f"    margins: f_u {m['f_u_hz']:.0f} Hz, PM {m['phase_margin_deg']:.2f}°, f_180 {m['f_180_hz']:.0f} Hz, "
              f"GM {m['gain_margin_db']:.2f} dB")
    # analytic: (1 + L̂)/L̂ = (g'/g)(1 + L)/L  ⇒  |L| → ∞: L̂ → 1/(g'/g − 1)
    print(f"  high-gain limit of L̂ for +1 dB: 1/(10^(1/20) − 1) = {1 / (10 ** (1 / 20) - 1):.2f}")


def part2():
    print("p7_literal TEC hold (L_hold = ω_h / (s(1 + sτ_th)), ω_h = 2π f_12):")
    print(f"  {'f_hold [Hz]':>11s} {'τ_th [s]':>9s} {'ζ':>6s} {'f_c [Hz]':>9s} {'PM [°]':>7s} {'1/(2π f_c) [s]':>15s}")
    for fh in (3e-3, 3e-2, 0.3):
        for tau in (1.0, 10.0, 100.0):
            wh = 2 * math.pi * fh
            zeta = 1 / (2 * math.sqrt(wh * tau))
            # unity gain: wh/(w sqrt(1 + (w tau)^2)) = 1
            w = math.sqrt((-1 + math.sqrt(1 + 4 * (wh * tau) ** 2)) / (2 * tau**2))
            pm = 90 - math.degrees(math.atan(w * tau))
            print(f"  {fh:11.3g} {tau:9.0f} {zeta:6.3f} {w / 2 / math.pi:9.4f} {pm:7.1f} {1 / w:15.2f}")


def part3():
    spec = config.load(ROOT / "presets/future/F_wg6_nominal.yaml",
                       overlays=["controller/dsc1", "ident/wire_F", "ident/in_lock_stepped"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run(spec)
    m = r.measurements["loop_identification"]["lock0"]
    f = m["f_hz"]
    sel = (f >= 1.0e3) & (f <= 2.0e4) & (m["coherence"] >= 0.9)
    print("F_wg6_nominal + controller/dsc1 + ident/wire_F + ident/in_lock_stepped:")
    print(f"  points {int(sel.sum())}, max |mag err| {100 * np.max(np.abs(m['mag_err'][sel])):.2f} %, "
          f"max |phase err| {np.max(np.abs(m['phase_err_deg'][sel])):.2f}°")
    for k in ("f_u_hz", "phase_margin_deg", "f_180_hz", "gain_margin_db"):
        print(f"  {k}: identified {m[k]:.4g}, model {m['model_' + k]:.4g}")


def _ident(overlays, params=None, toggles=None):
    import dataclasses
    import tempfile

    import yaml

    ovs = list(overlays)
    if params:
        t = Path(tempfile.mkdtemp()) / "pt.yaml"
        t.write_text(yaml.safe_dump({"offline": {"params": params}}))
        ovs.append(str(t))
    spec = config.load(ROOT / "presets/future/F_wg6_nominal.yaml", overlays=ovs)
    if toggles:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **toggles})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(spec).measurements["loop_identification"]["lock0"]


def part4():
    """Part 3's mismatch: DSC1 12-bit quantisation (LSB ≈ 1 excitation amplitude) and the coarse sweep near f_180."""
    q = {f"{b}.quantization": False for b in ("dac0f", "adc0", "dac1f", "adc1")}
    d = ["controller/dsc1", "ident/wire_F", "ident/in_lock_stepped"]
    fine = {"ident_f_lo_hz": 1.6e4, "ident_f_hi_hz": 2.2e4, "ident_n_freqs": 16}
    for what, ovs, prm, tg in (("dsc1, quantisation off, 25-pt sweep", d, None, q),
                               ("dsc1, 16-pt sweep 16–22 kHz", d, fine, None),
                               ("dsc1, quantisation off, 16-pt sweep", d, fine, q),
                               ("moku_delta, 16-pt sweep", ["controller/moku_delta", *d[1:]], fine, None),
                               ("fpga_generic, 16-pt sweep", d[1:], fine, None)):
        m = _ident(ovs, prm, tg)
        f = m["f_hz"]
        sel = (f >= 1.0e3) & (f <= 2.2e4) & (m["coherence"] >= 0.9)
        print(f"  {what}: max |mag err| {100 * np.max(np.abs(m['mag_err'][sel])):.2f} %, f_180 {m['f_180_hz']:.0f} "
              f"(model {m['model_f_180_hz']:.0f}) Hz, GM {m['gain_margin_db']:.2f} (model {m['model_gain_margin_db']:.2f}) dB")


def part5():
    """loop_margins takes the FIRST −180° crossing; chk.loop_stability the largest |L| over all crossings."""
    from cranberry.blocks.chk.loop_stability import n_integrators, nyquist_stable

    for pre in ("presets/bench/refD_lab.yaml", "presets/reference/refD_nominal.yaml"):
        p = chk_params(config.load(ROOT / pre), "chk.lockD_stability")
        f = np.geomspace(300.0, 1.0e5, 20000)
        lg = open_loop(f, p)
        ph = np.unwrap(np.angle(lg))
        cr = np.nonzero(np.diff(np.floor((ph + np.pi) / (2 * np.pi))))[0]
        m = loop_margins(f, lg)
        print(f"  {pre}: loop_margins GM {m['gain_margin_db']:.2f} dB at {m['f_180_hz']:.0f} Hz; Nyquist max |L| at −180° "
              f"{nyquist_stable(f, lg, n_integrators(p))[1]:.3f}; crossings (Hz, |L|) "
              + ", ".join(f"({f[i]:.0f}, {abs(lg[i]):.3f})" for i in cr))


if __name__ == "__main__":
    for p in sys.argv[1:] or ["1", "2", "3", "4", "5"]:
        globals()[f"part{p}"]()
