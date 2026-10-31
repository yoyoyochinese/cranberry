"""Phase 5 M1 review (Ref D): independent check of the DAC limit cycle (D61-3, D56-7, T15 §2.11 / test 10).

Independent of the cranberry code: a quasi-continuous model (dt = 0.1 µs) of the DSC1 loop built only
from the documents (DSC1 p.1-2, p.17-20; EO p.1-2; T15 §2.7; D53 gains; D61-2 delay accounting):

    v = −V0 sin e                      (fringe at +π/2, e = true phase error, T15 eq. 3 with S = 0)
    → LPF 120 kHz (1st order) → ADC 1 MHz, 2-point average, 16 bit ±4 V
    → PID every 2 µs: u = K_p e_D + h Σ K_i e_D + (K_d/h)(e_D[n] − e_D[n−1]),  e_D = S − v (p.18-19)
    → DAC 12 bit ±4 V, held 2 µs → LPF 100 kHz → pure delay 4.2 µs (5.7 − 1.5, D61-2)
    → 914 poles 29.5 / 112 kHz → 915B (rigid, or 20 kHz / Q 10) → e = φ_dist − k_act·g·y

All analogue sections are one bilinear state-space at 10 MHz (no half-sample ZOH artefact).
Also: a describing-function estimate (relay biased at a rounding threshold) of the limit-cycle
amplitude, and the "ideal" case (no delay, static actuator) for the T15 ±LSB/2 statement.

Run: uv run python reviews/scripts/m1_refD_dac_limit_cycle.py
"""

from __future__ import annotations

import math

import numba
import numpy as np
from scipy import signal
from scipy.optimize import brentq

DT = 1.0e-7            # 10 MHz quasi-continuous step
N_ADC = 10             # ADC sample every 1 µs
N_PID = 20             # PID / DAC every 2 µs
V0 = 3.33              # G ℜ P_in, T_opt = 1 (T15 §6.1)
K_ACT = 2 * math.pi * 16e-6 / (1.5488e-6 * 10.0)   # 6.491 rad/V (EO p.1, conv §4: δOPL, no ξ_eo)
G_DAC = 5.0e5 / (5.0e5 + 220.0)                  # scope 1 MΩ ‖ 914 1 MΩ, R_out 220 Ω (DSC1 p.2 note b)
KP, KI, KD = 2.41e-4, 1.406e-3 / 2e-6, 1.108e-2 * 2e-6   # D53 entries → p.17 units via p.20
H = 2.0e-6
Q_ADC, Q_DAC = 8.0 / 2**16, 8.0 / 2**12
LSB_PHASE = K_ACT * G_DAC * Q_DAC                # 12.67 mrad


def chain_ss(resonance: bool, poles: bool = True, lpf_out: bool = True):
    """Bilinear (A, B, C, D) at DT of LPF_out · 914 · 915B (DC gain 1)."""
    num, den = np.array([1.0]), np.array([1.0])
    ps = ([1.0e5] if lpf_out else []) + ([2.95e4, 1.12e5] if poles else [])
    for p in ps:
        w = 2 * math.pi * p
        num, den = np.polymul(num, [w]), np.polymul(den, [1.0, w])
    if resonance:
        w = 2 * math.pi * 2.0e4
        num, den = np.polymul(num, [w * w]), np.polymul(den, [1.0, w / 10.0, w * w])
    A, B, C, D = signal.tf2ss(num, den)
    Ad, Bd, Cd, Dd, _ = signal.cont2discrete((A, B, C, D), DT, method="bilinear")
    return (np.ascontiguousarray(Ad), np.ascontiguousarray(Bd.ravel()), np.ascontiguousarray(Cd.ravel()),
            float(Dd.ravel()[0]))


@numba.njit(cache=True)
def _run(n, phi0, ramp, a_in, Ad, Bd, Cd, Dd, n_delay, quant, use_in_lpf, use_avg, stride):
    nx = Bd.shape[0]
    x = np.zeros(nx)
    tmp = np.zeros(nx)
    dl = np.zeros(max(n_delay, 1))
    vf = 0.0
    adc_prev = 0.0
    adc_cur = 0.0
    integ = 0.0
    e_prev = 0.0
    dac = 0.0
    y = 0.0
    out = np.empty(n // stride)
    j = 0
    for i in range(n):
        e = phi0 + ramp * i * DT - K_ACT * G_DAC * y
        v = -V0 * math.sin(e)
        if use_in_lpf:
            vf += a_in * (v - vf)          # 1st-order LPF (exact for a held input at DT ≪ 1/ω)
        else:
            vf = v
        if i % N_ADC == 0:
            adc_prev = adc_cur
            s = vf
            if quant:
                s = Q_ADC * math.floor(s / Q_ADC + 0.5)
            adc_cur = s
        if i % N_PID == 0:
            pv = 0.5 * (adc_cur + adc_prev) if use_avg else adc_cur
            ed = 0.0 - pv                    # DSC1: e = S − v, S = 0 (p.18)
            integ += H * KI * ed
            u = KP * ed + integ + KD / H * (ed - e_prev)
            e_prev = ed
            u = min(max(u, -4.0), 4.0)
            if quant:
                u = Q_DAC * math.floor(u / Q_DAC + 0.5)
                u = min(max(u, -4.0), 4.0 - Q_DAC)
            dac = u
        # pure delay (ring buffer), then the analogue chain
        if n_delay > 0:
            k = i % n_delay
            d_out = dl[k]
            dl[k] = dac
        else:
            d_out = dac
        yy = Dd * d_out
        for p in range(nx):
            yy += Cd[p] * x[p]
        for p in range(nx):
            acc = Bd[p] * d_out
            for q in range(nx):
                acc += Ad[p, q] * x[q]
            tmp[p] = acc
        for p in range(nx):
            x[p] = tmp[p]
        y = yy
        if i % stride == 0 and j < out.size:
            out[j] = e
            j += 1
    return out


def simulate(resonance: bool, ramp_rad_s: float, dur_s: float, phi0: float = 0.3, ideal: bool = False):
    if ideal:  # no delay, static actuator, no filters, no 2-point average
        Ad, Bd, Cd, Dd = (np.zeros((1, 1)), np.zeros(1), np.zeros(1), 1.0)
        n_delay, use_lpf, use_avg = 0, False, False
    else:
        Ad, Bd, Cd, Dd = chain_ss(resonance)
        n_delay, use_lpf, use_avg = int(round(4.2e-6 / DT)), True, True
    a_in = 1.0 - math.exp(-2 * math.pi * 1.2e5 * DT)
    n = int(dur_s / DT)
    e = _run(n, phi0, ramp_rad_s, a_in, Ad, Bd, Cd, Dd, n_delay, True, use_lpf, use_avg, N_ADC)
    tail = e[int(0.05 / 1e-6):]
    return {"band_lsb": float(np.ptp(tail) / LSB_PHASE), "maxabs_lsb": float(np.max(np.abs(tail)) / LSB_PHASE),
            "rms_mrad": float(np.std(tail) * 1e3), "mean_mrad": float(np.mean(tail) * 1e3)}, tail


def loop_gain(f, resonance: bool):
    """Continuous L(iω) seen by the DAC quantiser (break at the DAC; DSC1 sign absorbed: L > 0 = negative fb)."""
    s = 2j * np.pi * np.asarray(f, float)
    c = KP + KI / s + KD * s
    p = np.exp(-s * 5.7e-6) / ((1 + s / (2 * np.pi * 1.2e5)) * (1 + s / (2 * np.pi * 1.0e5)))
    p = p / ((1 + s / (2 * np.pi * 2.95e4)) * (1 + s / (2 * np.pi * 1.12e5)))
    if resonance:
        w = 2 * np.pi * 2.0e4
        p = p * w * w / (s * s + s * w / 10 + w * w)
    return V0 * K_ACT * G_DAC * c * p, p


def describing_function(resonance: bool):
    """Biased relay (DAC input sitting on a rounding threshold): output ±q/2 square wave, DF 2q/(πA).

    Harmonic balance N(A)·L(iω₁₈₀) = −1 → A = 2q|L(iω₁₈₀)|/π; the phase error's fundamental is
    (4/π)(q/2)·k·g·|P(iω₁₈₀)| (P = actuator path from DAC output to phase, DC gain 1)."""
    f = np.geomspace(1e3, 2e5, 200001)
    lg, p = loop_gain(f, resonance)
    ph = np.unwrap(np.angle(lg))
    i = np.where(ph < -np.pi)[0][0]
    f180 = float(f[i])
    a_over_q = 2 * abs(lg[i]) / math.pi
    fund = (4 / math.pi) * 0.5 * abs(p[i])          # in LSB of phase (k g q = 1 LSB)
    return {"f180_khz": f180 / 1e3, "L_at_f180": float(abs(lg[i])), "gain_margin_db": float(-20 * np.log10(abs(lg[i]))),
            "relay_amp_over_q": a_over_q, "phase_fund_amp_lsb": fund, "phase_ptp_fund_lsb": 2 * fund}


def main():
    print(f"LSB (phase) = {LSB_PHASE * 1e3:.3f} mrad; K_p = {KP}, K_i = {KI:.1f} /s, K_d = {KD:.3e} s")
    for res in (False, True):
        print(f"DF, {'915B 20 kHz Q10' if res else 'rigid 915B'}:", describing_function(res))
    cases = [("ideal (no delay/filters), fixed", dict(resonance=False, ramp_rad_s=0.0, dur_s=0.2, ideal=True)),
             ("ideal (no delay/filters), ramp 0.4 rad/s", dict(resonance=False, ramp_rad_s=0.4, dur_s=1.2, ideal=True)),
             ("rigid, fixed", dict(resonance=False, ramp_rad_s=0.0, dur_s=0.2)),
             ("rigid, ramp 0.4 rad/s", dict(resonance=False, ramp_rad_s=0.4, dur_s=1.2)),
             ("20 kHz Q10, fixed", dict(resonance=True, ramp_rad_s=0.0, dur_s=0.2)),
             ("20 kHz Q10, ramp 0.4 rad/s", dict(resonance=True, ramp_rad_s=0.4, dur_s=1.2))]
    for name, kw in cases:
        r, tail = simulate(**kw)
        f, pxx = signal.welch(tail, fs=1e6, nperseg=65536)
        fpk = f[np.argmax(pxx[f > 500]) + np.argmax(f > 500)]
        print(f"{name:42s} band {r['band_lsb']:.3f} LSB, max|e| {r['maxabs_lsb']:.3f} LSB, "
              f"rms {r['rms_mrad']:.2f} mrad, mean {r['mean_mrad']:+.3f} mrad, PSD peak {fpk / 1e3:.1f} kHz")
    # fixed points: the band depends on where u* sits between two codes
    bands = []
    for phi0 in np.linspace(0.3, 0.3 + LSB_PHASE, 9):
        r, _ = simulate(False, 0.0, 0.12, phi0=float(phi0))
        bands.append((r["band_lsb"], r["maxabs_lsb"], r["rms_mrad"]))
    b = np.array(bands)
    print("rigid, 9 fixed points across one LSB: band", np.round(b[:, 0], 3), "max|e|", np.round(b[:, 1], 3),
          "rms mrad", np.round(b[:, 2], 2))


if __name__ == "__main__":
    main()
