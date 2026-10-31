"""Phase 5 M1 review (Ref D): f_sim × 2 convergence — systematic change vs seed scatter.

Same comparison as tests/integration/test_refD_f_sim_doubling (band variance 100 Hz–100 kHz of the true
phase, 1 MHz vs 2 MHz with the DSC1 blocks at their device rates), repeated over seeds, plus the
noise-free deterministic loop (P-only step response is not needed: the closed-loop |S| is compared
through a sine disturbance would need the injector; here the deterministic DAC limit cycle with ΔL = 0
and all noise off is used as a second, seed-free observable).

Run: uv run python reviews/scripts/m1_refD_fsim_doubling.py
"""

import importlib.util
from pathlib import Path

import numpy as np
from scipy import signal

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refD_support", ROOT / "tests" / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

REC = {"record": {"slots": ["lock.lockD.err"], "decim": 1}}
B2 = [{"name": "dsc1.adc", "rate_div": 2}, {"name": "dsc1.pid", "rate_div": 4},
      {"name": "dsc1.sup", "rate_div": 4}, {"name": "dsc1.dac", "rate_div": 4}]


def band_var(spec, f_lo=100.0, f_hi=1.0e5):
    fs = spec.run.f_sim_hz
    x = run(spec).series["lock.lockD.err"][int(0.02 * fs):]
    f, p = signal.welch(x, fs=fs, window="hann", nperseg=int(8192 * fs / 1e6), detrend="linear")
    m = (f >= f_lo) & (f <= f_hi)
    return float(np.trapezoid(p[m], f[m])), x


def main():
    r = []
    for seed in (20261031, 11, 12, 13):
        v1, _ = band_var(rs.load(run={"duration_s": 0.32, "seed": seed}, measure=REC))
        v2, _ = band_var(rs.load(run={"duration_s": 0.32, "f_sim_hz": 2e6, "seed": seed}, blocks=B2, measure=REC))
        r.append((v1, v2))
        print(f"seed {seed}: σ(1 MHz) {np.sqrt(v1) * 1e3:.3f} mrad, σ(2 MHz) {np.sqrt(v2) * 1e3:.3f} mrad, "
              f"ratio σ2/σ1 − 1 = {np.sqrt(v2 / v1) - 1:+.4f}")
    a = np.array(r)
    s1, s2 = np.sqrt(a[:, 0]), np.sqrt(a[:, 1])
    print(f"seed scatter of σ at 1 MHz: {s1.std(ddof=1) / s1.mean():.4f} (rel); mean ratio {np.mean(s2 / s1) - 1:+.4f}")
    # seed-free observable: the DAC limit cycle (ΔL = 0, all noise off, rigid / nominal resonance)
    tog = {**rs.NOISE_OFF, "dsc1.dac.quantization": True, "dsc1.adc.quantization": True}
    for res in (False, True):
        t = {**tog, "stretcher915b.resonance": res}
        out = []
        for fs, blocks in ((1e6, None), (2e6, B2)):
            spec = rs.load(params={"dL_m": 0.0}, run={"duration_s": 0.1, "f_sim_hz": fs}, blocks=blocks, measure=REC,
                           toggles=t)
            e = run(spec).series["lock.lockD.err"][int(0.05 * fs):]
            f, p = signal.welch(e, fs=fs, nperseg=int(32768 * fs / 1e6))
            out.append((np.std(e), f[np.argmax(p)]))
        print(f"DAC limit cycle {'20 kHz Q10' if res else 'rigid'}: rms {out[0][0] * 1e3:.3f} → {out[1][0] * 1e3:.3f} mrad, "
              f"peak {out[0][1] / 1e3:.2f} → {out[1][1] / 1e3:.2f} kHz")


if __name__ == "__main__":
    main()
