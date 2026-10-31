"""VT-34 (WP-F6 part; the laser-source parts are WP-B2's ``test_vt_34_laser.py``): emission-time
evaluation of the source phase (conv §7, T05 §3.5, V3–V5).

* A pure laser phase step reaches the observer exactly τ_g later: integer delay → a step at n0 + k;
  fractional delay D = k + μ → value (1 − μ) at n0 + k, full at n0 + k + 1 (linear interpolation,
  ``core.delay``), for a multi-segment path (τ_b = Σ segment group delays).
* T05 V3: two beams with equal total group delay along different segment splits see the same laser
  phase: |Δφ| < 1e−12 (laser common path).
* T05 V5 (Ref D, ΔL = 10 m, τ = 48.97 ns, f_sim = 100 MHz): a +1 MHz frequency step (φ_L = −2π∫δν,
  conv §4) changes φ_long − φ_short by +2π·δν·τ = +0.3077 rad.
* T05 V4 (sub-sample delay, dt = 1 µs, τ = 49 ns): with white FM noise the in-band PSD of
  φ_L(t) − φ_L(t − τ) through the beam assembler equals (2πτ)²·S_ν within ±10 % (f < f_Nyq/2).
* Frequency-offset beam: static phase 2π·f_off·τ_b (beam in its own rotating frame, conv §4).
"""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.core.fft import psd_onesided
from cranberry.path.topology import C_M_PER_S

N_G = 1.468


def length_for(samples, fs):
    return samples * C_M_PER_S / (N_G * fs)


def laser_run(topo, phi_l, fs, beams, **toggles):
    n = len(phi_l)
    pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=n),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}, **toggles), pt.record([f"beam.{b}.phi" for b in beams], n)]
    rec, _, _ = pt.run_and_record(pipe, fs, n, {"laser.las.phi_table": phi_l})
    return rec


@pytest.mark.parametrize("delay_samples", [7.0, 7.3])
def test_laser_phase_step_arrives_after_tau_g(delay_samples):
    fs = 1e8
    L = length_for(delay_samples, fs)
    segs = [pt.seg_spec("a", 0.4 * L, n_sub=3), pt.seg_spec("b", 0.6 * L, n_sub=2)]
    topo = pt.topo_spec(segs, [pt.beam_spec("d", [pt.step_spec("a"), pt.step_spec("b")], "det")], ["det"], [])
    n, n0 = 60, 20
    phi = np.zeros(n)
    phi[n0:] = 1.0
    out = laser_run(topo, phi, fs, ["d"])["beam.d.phi"]
    k = int(math.floor(delay_samples))
    mu = delay_samples - k
    assert np.all(out[: n0 + k] == 0.0)
    assert out[n0 + k] == pytest.approx(1.0 - mu, abs=1e-12)
    np.testing.assert_allclose(out[n0 + k + (1 if mu else 0):], 1.0, atol=1e-12)


def test_common_path_laser_cancels():
    fs = 1e8
    L = length_for(5.37, fs)
    segs = [pt.seg_spec("x", 0.3 * L), pt.seg_spec("y", 0.7 * L, n_sub=4), pt.seg_spec("z", L, n_sub=2)]
    beams = [pt.beam_spec("a", [pt.step_spec("x"), pt.step_spec("y")], "det"),
             pt.beam_spec("b", [pt.step_spec("z")], "det")]
    topo = pt.topo_spec(segs, beams, ["det"], [])
    phi = np.cumsum(np.random.default_rng(1).normal(size=400))
    rec = laser_run(topo, phi, fs, ["a", "b"])
    assert np.std(rec["beam.a.phi"]) > 1.0
    assert np.max(np.abs(rec["beam.a.phi"] - rec["beam.b.phi"])) < 1e-12 * np.max(np.abs(phi))


def test_refD_frequency_step():
    fs, dL, L2 = 1e8, 10.0, 0.5
    segs = [pt.seg_spec("A1", L2 + dL, n_sub=4), pt.seg_spec("A2", L2, n_sub=1)]
    beams = [pt.beam_spec("pr", [pt.step_spec("A1")], "pd"), pt.beam_spec("lo", [pt.step_spec("A2")], "pd")]
    topo = pt.topo_spec(segs, beams, ["pd"], [])
    n, n0, dnu = 200, 50, 1e6
    t = np.arange(n) / fs
    phi = np.where(np.arange(n) >= n0, -2 * np.pi * dnu * (t - t[n0]), 0.0)  # φ_L = −2π∫δν (conv §4)
    rec = laser_run(topo, phi, fs, ["pr", "lo"])
    d = rec["beam.pr.phi"] - rec["beam.lo.phi"]
    tau = N_G * dL / C_M_PER_S
    assert tau == pytest.approx(48.97e-9, rel=1e-3)
    assert d[-1] - d[0] == pytest.approx(2 * np.pi * dnu * tau, rel=1e-9)
    assert d[-1] - d[0] == pytest.approx(0.3077, abs=2e-4)


def test_sub_sample_delay_white_fm():
    fs, h0, n = 1e6, 8e3, 2**17
    L = 49e-9 * C_M_PER_S / N_G
    segs = [pt.seg_spec("long", L), pt.seg_spec("zero", 0.0)]
    beams = [pt.beam_spec("a", [pt.step_spec("long")], "det"), pt.beam_spec("b", [pt.step_spec("zero")], "det")]
    topo = pt.topo_spec(segs, beams, ["det"], [])
    rng = np.random.default_rng(34)
    dnu = rng.normal(size=n) * math.sqrt(h0 * fs / 2)  # white S_ν = h0 (one-sided) up to f_Nyq
    phi = -2 * np.pi * np.cumsum(dnu) / fs
    rec = laser_run(topo, phi, fs, ["a", "b"])
    f, s = psd_onesided(rec["beam.b.phi"] - rec["beam.a.phi"], fs, nperseg=4096)
    band = (f > 1e3) & (f < fs / 4)
    assert np.mean(s[band]) == pytest.approx((2 * np.pi * 49e-9) ** 2 * h0, rel=0.10)


def test_frequency_offset_static_phase():
    fs = 1e8
    L = length_for(3.0, fs)
    topo = pt.topo_spec([pt.seg_spec("c", L)], [pt.beam_spec("p", [pt.step_spec("c")], "det", freq_offset_hz=80e6)],
                        ["det"], [])
    rec = laser_run(topo, np.zeros(10), fs, ["p"])
    tau = N_G * L / C_M_PER_S
    assert rec["beam.p.phi"][-1] == pytest.approx(2 * np.pi * 80e6 * tau, rel=1e-12)
    rec_off = laser_run(topo, np.zeros(10), fs, ["p"], emission_time=False)
    assert rec_off["beam.p.phi"][-1] == 0.0
