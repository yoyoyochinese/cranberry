"""VT-39: probe-cancellation breakdown (T10 §7 tests 2–13; 00_summary §4 VT-39; WP-B8 with F6, B7, B3, B10).

Already covered elsewhere (not repeated here):
* T10 tests 2, 3, 4, 9 (point / distributed counter-propagating noise, PM axis difference 102.5 rad with
  TEC → 0 and 90° splice → 0, reciprocity) — ``test_vt_39_f6_counterprop_lock.py`` (WP-F6);
* T10 test 6 (stray total 302.8 photons / 10 ns, CW 1 mW, 63 m) — ``test_vt_36_fiber_quantum_extra.py`` (B10).

Here:
* test 5 (frequency offset): Δν_pr = 1 GHz, dτ_g/dT = 39 ps/(km·K), 63 m, 1 K, ideal static lock through the
  F6 path (the offset enters the probe's dT coefficient, segment_phase docstring): residual
  −2πΔν·L·(dτ_g/dT)·ΔT = −0.0154 rad (0.1 %).
* test 7 (gating, ``probe_pulser`` on B10's table): 100 MHz, T_f = 320 ns, τ_pr = 10 ns, T_gate = 1 ns,
  ε_AOM = 0: port 15 (30.6 m) has an offset with exactly 0 distributed stray photons in the gate; port 16
  has none (scan of 640 offsets, minimum > 0); the offset average equals the CW count (1e−3).
* test 8 (SMF orthogonal probe): Haar-random M, ⟨|⟨p|MᵀM|s⟩|²⟩ = 1/3 ± 1 % (10⁵ samples); diagonal (PM) M → 0.
* test 10 (laser phase noise, two ports, ideal locks with the actuator at the upstream end, through the F6
  path with emission-time laser phases): unlocked Θ_k = eq. (14), 4 sin²(πf(τ_k − τ_0)) S_φ (1e−9); locked
  Θ_k = 4 sin²(πf(τ_k − τ_0))·|e^{−iωT_a} + e^{−iωT_ref} − e^{−iωT_b}(e^{−iωτ_k} + e^{−iωτ_0})|² S_φ (1e−9),
  derived from T10 §2.3 keeping the reference phase and the τ_k transit of the correction (see
  ``_locked_laser_tf``). It reduces to T10 eq. (13), 16 sin²(πf(τ_k − τ_0)) sin²(πf(T_a − T_b′)) S_φ, with
  T_b′ = T_b + τ_k (probe delay to the detector **including** the common segment) only when T_ref = T_b + τ_0;
  the O(f²) in-band cancellation (|H|² ∝ f⁴) holds for any delays. **Finding (Phase 5):** T10 eq. (12)/(13)
  with T_b excluding the common segment and the claim "T_a = T_b ⇒ 0" do not follow from the path model.
* test 11 (probe shot noise): ``liu_probe_interferometer`` homodyne, P_det = 1 nW, η = 1 (ℜ = e/hν),
  strong reference: S_φ = hν0/(2P_det) = 6.43e−11 rad²/Hz (2 %). An ideal lock writes this sensor noise
  1:1 into the actuator, hence into the squeezed phase (in band).
* tests 12, 13: IL 0.3 dB on −10 dB → −7.96 dB; V = 0.05 + n_mode 0.03 → −7.96 dB.
* Liu specifics (T10 §3.2, D40): the actuator acts on the probe after the transit n_g x_a/c
  (``liu_loop_delay_s``; upstream 0, downstream ≈ n_g L/c), computed by the F6 subsegment model; detection
  after the demux (probe through the ``demux_switch`` segment) cancels the switch phase exactly in the
  squeezed light, detection before the demux leaves it whole in Θ_k (component 4).
"""

import math

import numpy as np
import pytest

import cranberry.blocks.dmx.testing as dt
import cranberry.path.testing as pt
from cranberry.blocks.dmx.liu import liu_loop_delay_s
from cranberry.blocks.err.testing import liu as liu_spec
from cranberry.blocks.err.testing import record as err_record
from cranberry.blocks.fib.scattering import H_PLANCK, stray_photons
from cranberry.blocks.mod.probe_pulser import gated_stray_photons
from cranberry.offline.gaussian import add_noise, cov_from_noise_levels, pure_loss
from cranberry.path.topology import C_M_PER_S

pytestmark = pytest.mark.numba
N_G = dt.N_G


def seg_len(samples, fs):
    return samples * C_M_PER_S / (N_G * fs)


def dtft(h, f, dt_s):
    n = np.arange(len(h))
    return np.exp(-2j * np.pi * np.outer(np.atleast_1d(f), n) * dt_s) @ h


# ----------------------------------------------------------------------------- test 5: frequency offset


def test_5_frequency_offset_probe():
    L, C_T, dnu, k_tau = 63.0, 44.0, 1e9, 39e-15  # dτ_g/dT = 39 ps/(km·K) [AE, T08]
    seg = pt.seg_spec("c", L, n_sub=4, drivers={"dT": "uniform"})
    beams = [pt.beam_spec("sq", [pt.step_spec("c", 1)], "o"),
             pt.beam_spec("pr", [pt.step_spec("c", -1)], "o", freq_offset_hz=dnu)]
    coef = pt.coef_spec("c", "dT", "sq", C_T * L) + pt.coef_spec("c", "dT", "pr", C_T * L + 2 * np.pi * dnu * k_tau * L)
    topo = pt.topo_spec([seg], beams, ["o"], coef)
    n = 300
    tab = np.zeros(n)
    tab[50:] = 1.0
    pipe = [pt.pathtest("pathtest_table", "t", "drv.c.dT", (), "K", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "beam.pr.phi"], n)]
    rec, _, _ = pt.run_and_record(pipe, 1e8, n, {"drv.c.dT_table": tab})
    d = rec["beam.sq.phi"] - rec["beam.pr.phi"]
    res = d[-1] - d[0]  # ideal static lock holds the probe: the squeezed light keeps Δψ_sq − Δψ_pr
    assert res == pytest.approx(-2 * np.pi * dnu * L * k_tau, rel=1e-9)
    assert res == pytest.approx(-0.0154, abs=5e-5)  # T10 quotes 0.0154 (rounded)


# ---------------------------------------------------------------------------------------- test 7: gating

TG = {"stray_rayleigh": True, "stray_brillouin": False, "stray_raman": False, "stray_reflections": False}


def port_stray(k, offsets):
    L = k * C_M_PER_S / (100e6 * N_G)
    p = dt.stray_params(L)
    kw = dict(p_avg_w=1e-3 * 10 / 320, duty=10 / 320, period_s=320e-9, eps_aom=0.0, gating=True, n_g=N_G,
              reflection_z_m=[])
    return p, np.array([gated_stray_photons(p, TG, t_offset_s=o, **kw)["rayleigh_gate"] for o in offsets])


def test_7_gating_k15_zero_k16_never():
    offs = np.arange(640) * 0.5e-9
    p15, n15 = port_stray(15, offs)
    assert p15["L_m"] == pytest.approx(30.6, abs=0.05)
    assert n15.min() == 0.0 and n15[offs == 5e-9][0] == 0.0  # 2·15·T_rep + τ_pr + T_gate = 311 ns < 320 ns
    p16, n16 = port_stray(16, offs)
    assert n16.min() > 0.0  # 2·16·T_rep alone fills the frame
    cw16 = stray_photons(p16, TG, 1e-3 * 10 / 320)["rayleigh_gate"]
    assert n16.mean() == pytest.approx(cw16, rel=1e-3)  # distributed: average = CW of the same P̄ (T10 §3.4.4)


def test_7_probe_pulser_block_gating():
    pp = dt.make(dt.pulser(), 1e6)[0]
    p = dt.stray_params(15 * C_M_PER_S / (100e6 * N_G))
    assert pp.stray(p, TG, n_g=N_G, reflection_z_m=[])["clicks_per_gate"] == 0.0
    cw = dt.make(dt.pulser(toggles={"gating": False}), 1e6)[0]
    assert cw.stray(p, TG, n_g=N_G, reflection_z_m=[])["clicks_per_gate"] > 0.0


# ------------------------------------------------------------------------------- test 8: SMF polarization


def test_8_smf_orthogonal_probe_one_third():
    rng = np.random.default_rng(8)
    z = rng.normal(size=(100_000, 2, 2)) + 1j * rng.normal(size=(100_000, 2, 2))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    m = q * (d / np.abs(d))[:, None, :]  # Haar U(2)
    mtm = np.swapaxes(m, 1, 2) @ m  # round trip MᵀM (reciprocity, T09)
    val = np.abs(mtm[:, 0, 1]) ** 2  # |⟨p|MᵀM|s⟩|², p = x, s = y
    assert val.mean() == pytest.approx(1 / 3, rel=0.01)
    pm = np.diag(np.exp(1j * np.array([0.3, -1.2])))
    assert abs((pm.T @ pm)[0, 1]) ** 2 == 0.0


# --------------------------------------------------------------------------------- test 10: laser noise


def two_port_topology(Ta, Tb, Tref, t0, tk, fs):
    segs = [pt.seg_spec("up", seg_len(Ta, fs)), pt.seg_spec("pb", seg_len(Tb, fs)),
            pt.seg_spec("rf", seg_len(Tref, fs)), pt.seg_spec("c0", seg_len(t0, fs)), pt.seg_spec("ck", seg_len(tk, fs)),
            pt.seg_spec("a0", 0.0, drivers={"stretcher": "uniform"}),
            pt.seg_spec("ak", 0.0, drivers={"stretcher": "uniform"})]
    s = pt.step_spec
    beams = [pt.beam_spec("sq0", [s("up"), s("a0"), s("c0")], "g0"),
             pt.beam_spec("sqk", [s("up"), s("ak"), s("ck")], "gk"),
             pt.beam_spec("pr0", [s("pb"), s("c0", -1), s("a0", -1)], "d0"),
             pt.beam_spec("prk", [s("pb"), s("ck", -1), s("ak", -1)], "dk"),
             pt.beam_spec("ref", [s("rf")], "d0")]
    coef = pt.coef_spec("a0", "stretcher", "*", 1.0) + pt.coef_spec("ak", "stretcher", "*", 1.0)
    return pt.topo_spec(segs, beams, ["g0", "gk", "d0", "dk"], coef)


def laser_lock_responses(Ta, Tb, Tref, t0, tk, fs=1e9, n=400, n0=5):
    """Impulse φ_L(n0) = 1; returns (Θ unlocked, Θ locked) impulse responses. Ideal locks (actuator at the
    upstream end x_a = 0, so the probe sees a(t) at detection): a_j = φ_ref − φ_pr,j,free (replay)."""
    topo = two_port_topology(Ta, Tb, Tref, t0, tk, fs)
    beams = ["sq0", "sqk", "pr0", "prk", "ref"]
    pipe = [pt.pathtest("pathtest_table", "act0", "drv.a0.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "actk", "drv.ak.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=n),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record([f"beam.{b}.phi" for b in beams], n)]
    las = np.zeros(n)
    las[n0] = 1.0
    free, _, _ = pt.run_and_record(pipe, fs, n, {"laser.las.phi_table": las}, compiled=True)
    tabs = {"laser.las.phi_table": las,
            "drv.a0.stretcher_table": free["beam.ref.phi"] - free["beam.pr0.phi"],
            "drv.ak.stretcher_table": free["beam.ref.phi"] - free["beam.prk.phi"]}
    locked, _, _ = pt.run_and_record(pipe, fs, n, tabs, compiled=True)
    assert np.allclose(locked["beam.prk.phi"], locked["beam.ref.phi"], atol=1e-12)  # the locks hold
    return free["beam.sqk.phi"] - free["beam.sq0.phi"], locked["beam.sqk.phi"] - locked["beam.sq0.phi"]


def _locked_laser_tf(f, Ta, Tb, Tref, t0, tk, dt_s):
    """Θ̃_k/φ̃_L for ideal locks with the actuator upstream: θ_k(t) = φ_L(t − T_a − τ_k) + a_k(t − τ_k),
    a_k(t) = φ_L(t − T_ref) − φ_L(t − T_b − τ_k) (T10 §2.3 with the reference term kept)."""
    e = lambda d: np.exp(-2j * np.pi * f * d * dt_s)  # noqa: E731
    return (e(tk) - e(t0)) * (e(Ta) + e(Tref) - e(Tb) * (e(tk) + e(t0)))


@pytest.mark.parametrize("delays", [(30, 20, 47, 5, 36), (30, 20, 25, 5, 36), (12, 12, 40, 3, 31)])
def test_10_laser_phase_noise_through_locks(delays):
    Ta, Tb, Tref, t0, tk = delays
    fs = 1e9
    free, locked = laser_lock_responses(*delays, fs=fs)
    f = np.array([1e5, 1e6, 3e6, 1e7, 2.3e7, 4.1e7])
    hf, hl = np.abs(dtft(free, f, 1 / fs)) ** 2, np.abs(dtft(locked, f, 1 / fs)) ** 2
    np.testing.assert_allclose(hf, 4 * np.sin(np.pi * f * (tk - t0) / fs) ** 2, rtol=1e-9)  # eq. (14)
    want = np.abs(_locked_laser_tf(f, *delays, 1 / fs)) ** 2
    np.testing.assert_allclose(hl, want, rtol=1e-9, atol=1e-24)
    # in band the lock removes the O(f²) term: |H|² ∝ f⁴ at low f
    lo = np.abs(dtft(locked, np.array([1e4, 2e4]), 1 / fs)) ** 2
    assert lo[1] / lo[0] == pytest.approx(16.0, rel=1e-3)
    if Tref == Tb + t0:  # T10 eq. (13) with T_b′ = T_b + τ_k
        t10 = 16 * np.sin(np.pi * f * (tk - t0) / fs) ** 2 * np.sin(np.pi * f * (Ta - Tb - tk) / fs) ** 2
        np.testing.assert_allclose(hl, t10, rtol=1e-9, atol=1e-24)


# ---------------------------------------------------------------------------- test 11: probe shot noise


def test_11_probe_shot_noise():
    lam, p_det, fs, n = dt.LAMBDA0, 1e-9, 1e6, 200_000
    hnu = H_PLANCK * C_M_PER_S / lam
    resp = 1.602176634e-19 / hnu  # η = 1
    blk = liu_spec(noise=True, gating=False, duty=1.0, p_pr_avg_w=p_det, p_ref_w=1e-3, responsivity_a_per_w=resp,
                   i_n_a_rthz=0.0, toggles={"electronic_noise": False})
    v = err_record([blk], fs, n, ["v.liu"], init={"beam.pr.phi": math.pi / 2})["v.liu"]
    k = 1e4 * 2 * resp * math.sqrt(p_det * 1e-3)  # fringe slope at X = π/2 (G = 1e4 V/A, visibility 1)
    s_phi = v.var() / (fs / 2) / k**2  # white, one-sided PSD = var/(f_s/2)
    assert hnu / (2 * p_det) == pytest.approx(6.43e-11, rel=1e-3)
    assert s_phi == pytest.approx(hnu / (2 * p_det) * (1 + p_det / 1e-3), rel=0.02)


# --------------------------------------------------------------------------------------- tests 12, 13


def test_12_13_loss_and_stray_noise():
    db = lambda v: 10 * math.log10(v / 0.5)  # noqa: E731
    v = cov_from_noise_levels(0.05, 20.0, 0.0)
    assert db(pure_loss(v, 10 ** (-0.03))[0, 0]) == pytest.approx(-7.96, abs=0.01)
    assert db(add_noise(v, 0.03)[0, 0]) == pytest.approx(-7.96, abs=0.01)


# ------------------------------------------------------------------------------------- Liu specifics


@pytest.mark.parametrize("sub", [0, 9])
def test_liu_actuator_position_sets_probe_delay(sub):
    N, L = 10, 63.0
    fs = 2 * N * C_M_PER_S / (N_G * L)  # subsegment delay δ = 2 samples
    seg = pt.seg_spec("c", L, n_sub=N, spatial="mixed", drivers={"stretcher": "actuator"}, actuator_sub=sub)
    topo = pt.topo_spec([seg], [pt.beam_spec("pr", [pt.step_spec("c", -1)], "det")], ["det"],
                        pt.coef_spec("c", "stretcher", "*", 1.0))
    n, n0 = 120, 10
    imp = np.zeros(n)
    imp[n0] = 1.0
    pipe = [pt.pathtest("pathtest_table", "a", "drv.c.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.pr.phi"], n)]
    rec, _, _ = pt.run_and_record(pipe, fs, n, {"drv.c.stretcher_table": imp})
    h = rec["beam.pr.phi"]
    delay = (np.arange(n) * h).sum() / h.sum() - n0
    x_a = (sub + 0.5) * L / N  # subsegment centre (F6 placement)
    assert delay / fs == pytest.approx(liu_loop_delay_s(0.0, x_a, L, N_G), rel=1e-9)


def demux_lock_residual(after_demux):
    """Switch phase (bias random walk) on a lumped demux segment upstream of port k's fiber; ideal lock on
    the probe (actuator upstream). Returns the squeezed phase at the GBS (laser off)."""
    fs, n, tk = 1e8, 600, 20
    s = pt.step_spec
    segs = [pt.seg_spec("dmx", 0.0, drivers={"eom": "per_sub"}), pt.seg_spec("act", 0.0, drivers={"stretcher": "uniform"}),
            pt.seg_spec("c", seg_len(tk, fs))]
    pr = [s("c", -1), s("act", -1)] + ([s("dmx", -1)] if after_demux else [])
    beams = [pt.beam_spec("sq", [s("dmx"), s("act"), s("c")], "gbs"), pt.beam_spec("pr", pr, "det")]
    coef = pt.coef_spec("dmx", "eom", "*", 1.0) + pt.coef_spec("act", "stretcher", "*", 1.0)
    topo = pt.topo_spec(segs, beams, ["gbs", "det"], coef)
    sw = dt.switch(n_stages=3, toggles={"drive_phase_noise": False, "bias_drift": True}, bias_rw_rad_per_rts=3.0)
    pipe = [sw, pt.pathtest("pathtest_table", "a", "drv.act.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "beam.pr.phi", "act.sw.dpsi"], n)]
    free, _, _ = pt.run_and_record(pipe, fs, n, {}, compiled=True)
    locked, _, _ = pt.run_and_record(pipe, fs, n, {"drv.act.stretcher_table": -free["beam.pr.phi"]}, compiled=True)
    np.testing.assert_array_equal(locked["act.sw.dpsi"], free["act.sw.dpsi"])  # same switch noise (same seed)
    return free, locked, tk


def test_liu_detection_after_demux_cancels_switch_phase():
    free, locked, tk = demux_lock_residual(after_demux=True)
    assert np.std(free["beam.sq.phi"]) > 1e-4
    assert np.max(np.abs(locked["beam.sq.phi"][tk + 2:])) < 1e-12


def test_liu_detection_before_demux_keeps_switch_phase():
    free, locked, tk = demux_lock_residual(after_demux=False)
    np.testing.assert_allclose(locked["beam.sq.phi"][tk + 2:], free["act.sw.dpsi"][2:-tk], rtol=0, atol=1e-12)
