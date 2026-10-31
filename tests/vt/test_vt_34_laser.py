"""VT-34 (WP-B2 part; emission-time / path parts are WP-F6's ``test_vt_34_f6_emission_time.py``): laser
source, delay transfer, ledger, filter cavity, ASE, RIN through CMRR (T05 §10).

* V1 delay transfer (white FM Δν = 25 kHz, dt = 1 ns, τ = 49 dt): Welch PSD of Δφ = φ(t) − φ(t − τ)
  vs 4 sin²(πfτ) h0/f² within ±10 % for f ≤ 0.1/dt, away from the nulls; up to 0.4/dt against the
  same formula with the aliases of the sampled process included (the generator samples an exact
  Wiener phase, so the aliased sum Σ_k S(f + k/dt) is the PSD of the samples; the note's 0.4/dt for
  the bare formula would need a band-limited φ). Variance 2πΔντ = 7.70e−3 rad² ±5 %.
* V2 coherence: τ = 10 µs, ⟨cos Δφ⟩ = exp(−πΔντ) = 0.456 ± 0.01 (dt = 0.1 µs: the white-FM phase is an
  exact random walk at any step, so the statistic does not depend on dt).
* V3 common path: equal delays ⇒ laser part of Δφ exactly 0 (< 1e−12), through segment_phase /
  beam_assembler.
* V6 Ref D floor (+ V4 sub-sample delay): ΔL = 10 m (τ = 49.0 ns), f_sim = 1 MHz, only laser noise:
  PSD of φ_pr − φ_lo = (2πτ)² h0 = 7.5e−10 rad²/Hz ±10 %; in-band variance + the ledger's out-of-band
  variance (the block's ``oob_pairs`` row, 4 sin²(πfτ) folded in) = 2πΔντ ±2 %; ``laser_oob_variance``
  agrees with the ledger row.
* V8 port 31 (no lock): Δν = 1 kHz, T_rep = 10 ns, k = 31: Var = 2πΔν·31·T_rep = 1.95e−3 rad²
  (σ = 44.1 mrad) ±5 %.
* V11 filter cavity: δ = 0.1 Δν_FC → arg T = +0.1974 rad, |T|²/T0² = 0.9615; |H_FC(Δν_FC/2)|² = 0.5;
  the block writes the same ψ / transmission.
* V12 ASE: S_ASE = 0 changes nothing (exact); ε_ASE = 1e−3 with Ref A R± (−10.58 / +20.56 dB) → −8.41 dB
  (±0.01 dB); a cavity of finesse F passes π/(2F) of broadband ASE (B_sq ≫ FSR, ±2 %).
* V13 RIN through CMRR: P_LO = 15.4 mW, η = 0.99, RIN −145 dB/Hz, ε = 0.01 → excess/shot = 1.88e−2.
"""

import math

import numpy as np
import pytest

import cranberry.blocks.noise.testing as nt
import cranberry.path.testing as pt
from cranberry.blocks.err.testing import record, rig
from cranberry.blocks.noise.laser import laser_oob_variance, rin_cmrr_excess, s_nu_model
from cranberry.blocks.spf.filter_cavity import (ase_epsilon, ase_mixed_variance, ase_suppression, phase_transfer,
                                                transmission)
from cranberry.core.compile import run_reference
from cranberry.core.fft import psd_onesided
from cranberry.core.pipeline import iter_leaves
from cranberry.offline.band_ledger import BandLedger
from cranberry.path.topology import C_M_PER_S

pytestmark = pytest.mark.numba
DNU = 25e3
H0 = DNU / math.pi


def white(h0=H0, **over):
    return nt.laser(on=("fm_white",), h0_hz2_per_hz=h0, **over)


def phi_series(fs, n, h0=H0, seed=11):
    return record([white(h0)], fs, n, ["laser.dx1.phi"], seed=seed, external={}, chunk=1 << 20)["laser.dx1.phi"]


def test_v1_delay_transfer_psd_and_variance():
    dt, k = 1e-9, 49
    tau = k * dt
    phi = phi_series(1 / dt, 4_000_000)
    d = phi[k:] - phi[:-k]
    assert d.var() == pytest.approx(2 * math.pi * DNU * tau, rel=0.05)
    assert 2 * math.pi * DNU * tau == pytest.approx(7.70e-3, rel=1e-3)
    f, s = psd_onesided(d, 1 / dt, nperseg=1024)
    s4 = 4 * np.sin(np.pi * f * tau) ** 2
    ok = (f > 0) & (s4 > 1.0)
    bare = s4 * H0 / np.where(f > 0, f, 1.0) ** 2
    aliased = s4 * H0 * np.pi**2 * dt**2 / np.where(f > 0, np.sin(np.pi * f * dt), 1.0) ** 2
    lo = ok & (f <= 0.1 / dt)
    hi = ok & (f <= 0.4 / dt)
    np.testing.assert_allclose(s[lo] / bare[lo], 1.0, atol=0.10)
    np.testing.assert_allclose(s[hi] / aliased[hi], 1.0, atol=0.10)


def test_v2_coherence_at_10us():
    dt, k = 1e-7, 100
    phi = phi_series(1 / dt, 4_000_000, seed=12)
    c = np.mean(np.cos(phi[k:] - phi[:-k]))
    assert c == pytest.approx(math.exp(-math.pi * DNU * 1e-5), abs=0.01)
    assert math.exp(-math.pi * DNU * 1e-5) == pytest.approx(0.456, abs=5e-4)


def _two_arm(L_lo, L_pr, n, fs, laser, seed=5):
    segs = [pt.seg_spec("s_lo", L_lo), pt.seg_spec("s_pr", L_pr)]
    beams = [pt.beam_spec("lo", [pt.step_spec("s_lo", 1)], "det", source="dx1"),
             pt.beam_spec("pr", [pt.step_spec("s_pr", 1)], "det", source="dx1")]
    topo = pt.topo_spec(segs, beams, ["det"], [])
    pipe = [laser, *pt.path_pipeline(topo, {"dx1": "laser.dx1.phi"}), pt.record(["beam.lo.phi", "beam.pr.phi"], n)]
    rec, st, _ = pt.run_and_record(pipe, fs, n, compiled=True, seed=seed)
    return rec["beam.pr.phi"] - rec["beam.lo.phi"], st


def test_v3_common_path_cancels_exactly():
    d, _ = _two_arm(3.0, 3.0, 20000, 1e6, nt.laser(on=("fm_white", "fm_flicker", "drift", "dither")))
    assert np.max(np.abs(d)) < 1e-12


def test_v6_refD_floor_and_v4_band_split():
    fs, n, dL = 1e6, 1 << 20, 10.0
    tau = 1.468 * dL / C_M_PER_S
    las = white(oob_pairs=[{"target": "laser_pr_lo", "tau_s": tau}])
    d, _ = _two_arm(0.5, 0.5 + dL, n, fs, las)
    f, s = psd_onesided(d, fs, nperseg=1024)
    band = (f > 2 * fs / 1024) & (f < fs / 4)  # first Welch bins carry the Hann/detrend bias
    floor = (2 * math.pi * tau) ** 2 * H0
    assert floor == pytest.approx(7.5e-10, rel=0.02)
    assert np.mean(s[band]) == pytest.approx(floor, rel=0.10)
    np.testing.assert_allclose(s[band] / floor, 1.0, atol=0.15)
    # V4: in-band variance + ledger out-of-band variance = 2πΔντ
    root, _, _ = rig([las], fs, external={})
    led = BandLedger(root.bands(), beams={"pr": ("laser_pr_lo",), "lo": ()})
    led.define_metric("refD", f_max_hz=1e3 / tau, weight=None)
    oob = led.hf_variance("refD", {"pr": 1.0, "lo": -1.0})
    assert d.var() + oob == pytest.approx(2 * math.pi * DNU * tau, rel=0.02)
    helper = laser_oob_variance(s_nu_model(las["params"]), tau, fs / 2, 1e3 / tau)
    assert oob == pytest.approx(helper, rel=0.01)


def test_v8_port31_decorrelation():
    fs, k = 1e8, 31
    phi = phi_series(fs, 3_000_000, h0=1e3 / math.pi, seed=13)
    v = np.var(phi[k:] - phi[:-k])
    assert v == pytest.approx(2 * math.pi * 1e3 * k / fs, rel=0.05)
    assert math.sqrt(2 * math.pi * 1e3 * k / fs) == pytest.approx(44.1e-3, rel=2e-3)


def test_v11_filter_cavity():
    t = transmission(0.1, 1.0, 0.9)
    assert np.angle(t) == pytest.approx(0.1974, abs=1e-4) and np.angle(t) > 0
    assert abs(t) ** 2 / 0.81 == pytest.approx(0.9615, abs=1e-4)
    assert abs(phase_transfer(0.5e6, 1e6)) ** 2 == pytest.approx(0.5, rel=1e-12)
    blk = {"block": "filter_cavity", "name": "fc", "params": {
        "segment": "fc", "beam": "lo", "dnu_fc_hz": 1e6, "t0": 0.9, "source": "V11 test", "in_detuning": "v.det",
        "toggles": {"lock_residual": True, "laser_lowpass": False}}}
    root, st, pa = rig([blk], 1e7, external={"v.det": "Hz"})
    st["v.det"] = 1e5
    run_reference(root, st, pa, 1)
    assert float(st["psi.fc.lo"]) == pytest.approx(np.angle(transmission(1e5, 1e6)), rel=1e-14)
    assert float(st["psi.fc.lo.trans"]) == pytest.approx(abs(transmission(1e5, 1e6, 0.9)) ** 2, rel=1e-14)


def test_v12_ase():
    vm, vp = 10 ** (-1.058), 10 ** (2.056)
    assert ase_mixed_variance(vm, vp, vm, ase_epsilon(0.0, 15.4e-3)) == vm
    assert 10 * math.log10(ase_mixed_variance(vm, vp, vm, 1e-3)) == pytest.approx(-8.41, abs=0.01)
    for finesse in (100, 1000, 1e4):
        fsr = 1e9
        r = ase_suppression(fsr, fsr / finesse, 2.5e12 + 0.37 * fsr)
        assert r == pytest.approx(math.pi / (2 * finesse), rel=0.02)


def test_v13_rin_through_cmrr():
    x = rin_cmrr_excess(15.4e-3, 0.99, 10 ** (-14.5), 0.01, 1548.8e-9)
    assert x == pytest.approx(1.88e-2, rel=5e-3)


def test_laser_slots_wired_to_path():
    _, st = _two_arm(1.0, 2.0, 10, 1e6, white())
    assert "laser.dx1.phi" in st.layout and any(leaf.kind == "laser_source" for leaf in iter_leaves(
        rig([white()], 1e6, external={})[0]))
