"""VT-13 (Ref D part, WP-I2): closed-loop transfer function of the refD preset = analytic (prompt §7-7;
T15 §7 test 4; notes/00_summary.md §4 VT-13).

Linear test case: ``refD_nominal`` (D53 DSC1 gains through the F5 mapping, 914 two poles, 915B
resonance 20 kHz / Q 10, τ_L = 5.7 µs, 120 kHz / 100 kHz poles, PID / DAC every 2 samples) with every
noise source, drift, quantiser, PER term and γ_coh switched off. A sinusoidal phase disturbance d is
injected on arm 2 (strain driver of segment a2, in place of the mechanical-noise block) and the true
φ_D(f)/d(f) is demodulated from the run.

(a) It equals the exact multirate z-domain transfer function built from the same discretisations
    (``refD_support.sensitivity_z``: C(z₂) of the dsc1 PID, 2-point average, decimation by 2 with its
    image f + 250 kHz... folded back, ZOH hold, delay-buffer samples, zoh-discretised 914, Tustin 915B
    mode, bilinear LPFs, segment-transit interpolation) to relative 1e−3 and 0.1° at f = 10 Hz, 100 Hz,
    1 kHz, f_c/2, f_c, 2 f_c and 20 kHz (T15: 1e−3, 0.1°).
(b) It equals the continuous 1/(1 + L(i2πf)) of T15 eq. (10) (C(s) = K_p + K_i/s + K_d s, first-order
    120 kHz / 100 kHz poles, e^{−sτ_L} with τ_L = 5.7 µs, 914 poles, 915B mode) within 5 % and 3° for
    f ≤ f_c (T15).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parents[1] / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

FS = 1.0e6
A_PHI = 1.0e-3  # rad, disturbance amplitude on φ_D (fringe nonlinearity ~ A²/6 = 2e−7)
T_SETTLE = 0.02  # s, lock acquisition + transients (integrator pole e^{−2π f_c t}, resonance Q/(π f_r))
N_REC = 100_000  # 0.1 s: integer periods of f and of the folded image f + 500 kHz for f ∈ 10 Hz · ℤ
TOG = {**rs.NOISE_OFF, "pdb450c.per_crosscoupling": False, "pdb450c.visibility_hf_averaging": False}


def _base():
    return rs.load(run={"duration_s": T_SETTLE + N_REC / FS},
                   measure={"record": {"slots": ["lock.lockD.err", "drv.a2.strain"], "decim": 1}}, toggles=TOG)


def _continuous_s(spec, f, k_d):
    p = rs.block_params(spec, "dsc1.pid")
    st = rs.block_params(spec, "stretcher915b")
    pr = spec.offline["params"]
    s = 2j * np.pi * f
    c = p["kp"] + p["ki_per_s"] / s + p["kd_s"] * s
    k_act = 2 * np.pi * st["dopl_full_m"] / (st["lambda0_m"] * st["v_full_v"]) * st["polarity"]
    rl = 1.0 / (1.0 / pr["r_scope_ohm"] + 1.0 / pr["r_914_ohm"])
    g = rl / (rl + 220.0)
    amp = 1.0
    for fp in (2.95e4, 1.12e5):
        amp = amp / (1 + s / (2 * np.pi * fp))
    wr, q = 2 * np.pi * pr["f_res_hz"], pr["q_res"]
    h915 = wr**2 / (s**2 + s * wr / q + wr**2)
    lpf = 1 / ((1 + s / (2 * np.pi * 1.2e5)) * (1 + s / (2 * np.pi * 1.0e5)))
    loop = g * c * k_d * lpf * np.exp(-s * pr["latency_s"]) * amp * (-k_act) * h915
    return 1.0 / (1.0 - loop), -loop


@pytest.fixture(scope="module")
def setup():
    spec = _base()
    fr = rs.fringe(spec)
    k_d = -fr["sign"] * fr["v0"] * fr["vis"]  # ∂v/∂φ_D at +π/2 (T15 eq. 5)
    fgrid = np.geomspace(100.0, 2e4, 4000)
    _, loop = _continuous_s(spec, fgrid, k_d)
    f_c = float(fgrid[np.argmin(np.abs(np.abs(loop) - 1.0))])
    return spec, k_d, f_c


def _measure(spec, f_hz):
    a2 = rs.segment(spec, "a2")
    k_eps = spec.offline["params"]["k_eps_rad_per_m"]
    amp_strain = A_PHI / (k_eps * a2["L_m"])
    inj = [{"block": "sine_injection", "name": "inj", "params": {"in": None, "out": "drv.a2.strain", "unit": "1",
                                                                 "freq_hz": f_hz, "amplitude": amp_strain,
                                                                 "toggles": {"enabled": True}}},
           {"block": "sine_injection", "name": "zero_a1", "params": {"in": None, "out": "drv.a1.strain", "unit": "1",
                                                                     "freq_hz": 1.0, "amplitude": 0.0,
                                                                     "toggles": {"enabled": True}}}]
    sp = rs.replace_block(spec, "env.mech", inj)
    res = run(sp)
    n = np.arange(res.t_s.size)[-N_REC:]
    err = res.series["lock.lockD.err"][-N_REC:]
    strain = res.series["drv.a2.strain"][-N_REC:]
    from cranberry.path.segments import frac_delay

    k2, mu2 = frac_delay(0.5 * a2["n_g"] * a2["L_m"] / rs.C_M_PER_S, 1.0 / FS)
    zi = np.exp(-2j * np.pi * f_hz / FS)
    d = k_eps * a2["L_m"] * rs.demod(strain, n, f_hz, FS) * zi**k2 * ((1 - mu2) + mu2 * zi)
    return rs.demod(err, n, f_hz, FS) / d


def _freqs(f_c):
    r10 = lambda x: float(round(x / 10.0) * 10.0)  # noqa: E731
    return [10.0, 100.0, 1000.0, r10(f_c / 2), r10(f_c), r10(2 * f_c), 20000.0]


def test_vt13_refD_matches_z_domain_and_continuous(setup):
    spec, k_d, f_c = setup
    assert 1.5e3 < f_c < 3.0e3  # K_0 K_i / 2π with D53 K_i = 703 s⁻¹ (T15 §2.9 scale)
    rows = []
    for f in _freqs(f_c):
        sm = _measure(spec, f)
        sz = complex(rs.sensitivity_z(spec, np.array([f]), k_d=k_d)[0])
        sc = complex(_continuous_s(spec, np.array([f]), k_d)[0][0])
        rows.append((f, sm, sz, sc))
        assert abs(sm / sz - 1) < 1e-3, (f, sm, sz)
        assert abs(math.degrees(np.angle(sm / sz))) < 0.1, (f, sm, sz)
        if f <= f_c:
            assert abs(abs(sm) / abs(sc) - 1) < 0.05, (f, sm, sc)
            assert abs(math.degrees(np.angle(sm / sc))) < 3.0, (f, sm, sc)
    # strong rejection inside the loop band, |S| ≈ f/f_c for an integrator-dominated loop
    assert abs(rows[0][1]) < 2 * 10.0 / f_c
