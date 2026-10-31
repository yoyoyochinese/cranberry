"""VT-23 (WP-B6 part; WP-I4 adds the Ref B preset, WP-B11 the pulsed OPO): no cancellation on a shared
segment with unequal segment × beam coefficients — the OPO cavity (prompt §7-11 added by D51-4;
00_summary §4 VT-23; T02 V-B7; D47).

Topology: squeezed beam ``sq`` and co-propagating probe ``pr`` both traverse the OPO segment (driver
``detuning``, coefficients from ``opo_cw_spectrum.coef_entries``: sq −1/(κ(1 + x²)), pr −1/(2κ)) and a
shared fibre segment with equal coefficients. Lock 2 holds the LO on the probe, so the lock combination
φ_pr − φ_sq is the residual φ_LO − θ_sq.

* Ref B, x = 0.745, constant Δ/κ = 1e−3: residual = (Δ/κ)(1 − x²)/(2(1 + x²)) (relative 1e−12 against the
  coefficients; relative 1e−3 against the Langevin ellipse angle minus the probe singular-vector angle).
* With opo_detuning_noise → opo_length_lock driving Δ(t) and white noise on the shared fibre: the probe
  does not fully follow — φ_pr − φ_sq = (c_pr − c_sq)/c_pr · ψ_opo,pr at every sample, read at the
  emission time through the shared fibre's group delay (the shared fibre and laser cancel to 1e−12), i.e. a fraction (1 − x²)/(1 + x²) = 0.286 of the probe's OPO phase remains.
* Control: equal OPO coefficients (probe given the squeezed coefficient) → residual 0 (1e−12), as VT-22.
"""

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.source.opo_testing import REF_B, opo
from cranberry.offline.source_output import wrap_half_pi

X = 0.745
FS = 1.0e6
N = 3000


def topology(block, *, equal=False):
    segs = [pt.seg_spec("opo", 0.5, kind="opo", n_g=1.0, drivers={"detuning": "uniform"}),
            pt.seg_spec("s1", 3.0, n_sub=3, drivers={"dT": "uniform"})]
    path = [pt.step_spec("opo"), pt.step_spec("s1")]
    beams = [pt.beam_spec("sq", path, "bhd"), pt.beam_spec("pr", path, "bhd")]
    if equal:  # control: the probe gets the squeezed-beam coefficient
        coef = block.coef_entries(["sq", "pr"], [])
    else:
        coef = block.coef_entries(["sq"], ["pr"])
    coef = coef + pt.coef_spec("s1", "dT", "*", 300.0)
    locks = [{"id": "lock2", "combo": {"pr": 1.0, "sq": -1.0}, "setpoint_rad": 0.0, "observer": "bhd"}]
    return pt.topo_spec(segs, beams, ["bhd"], coef, locks)


def refb():
    return opo(REF_B, loss_ic=0.0038, pump_mode="x", x=X)


def run_constant(block, delta, **kw):
    topo = topology(block, **kw)
    pipe = [pt.pathtest("pathtest_table", "det", "drv.opo.detuning", (), "rad/s", n_rows=N),
            pt.pathtest("pathtest_white", "fib", "drv.s1.dT", (), "K", sigma=0.0),
            pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=0.0),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True),
            pt.record(["lock.lock2.err", "psi.opo.sq", "psi.opo.pr"], N)]
    rec, _, _ = pt.run_and_record(pipe, FS, N, {"drv.opo.detuning_table": np.full(N, delta)})
    return rec


def test_constant_detuning_residual_matches_formula_and_langevin():
    b = refb()
    k = b.kappa(X)
    d = 1e-3
    rec = run_constant(b, d * k)
    err = rec["lock.lock2.err"][10:]
    want = d * (1 - X**2) / (2 * (1 + X**2))
    np.testing.assert_allclose(err, want, rtol=1e-12)
    # independent: Langevin ellipse angle (θ_off of the detuned spectrum) and probe singular vector
    th = float(opo(REF_B, loss_ic=0.0038, pump_mode="x", x=X, f_s_hz=[0.0], detuning_rad_s=d * k,
                   toggles={"detuning": True}).source_output.project(None)[2][0])
    u, _, _ = np.linalg.svd(np.linalg.inv(np.array([[k - X * k, -d * k], [d * k, k + X * k]])))
    phi_pr = float(wrap_half_pi(np.arctan2(u[1, 0], u[0, 0])))
    assert err[-1] == pytest.approx(phi_pr - th, rel=1e-3)
    # equal coefficients: cancels
    assert np.max(np.abs(run_constant(b, d * k, equal=True)["lock.lock2.err"])) < 1e-12 * d


@pytest.mark.parametrize("compiled", [False, True])
def test_noisy_detuning_probe_does_not_follow(compiled):
    b = refb()
    topo = topology(b)
    pipe = [{"block": "opo_detuning_noise", "name": "dn", "params": {
                "out": "opo.dn.delta_free", "segment": "opo", "nu0_hz": 3.49e14, "round_trip_m": 0.5,
                "h_m2_hz": 1e-24, "f_lo_hz": 10.0, "f_c_hz": 2e4,
                "peaks": [{"f_hz": 3e3, "q": 20.0, "area_m2": 1e-20}],
                "toggles": {"length_noise": True, "peaks": True, "thermal": False}}},
            {"block": "opo_length_lock", "name": "ll", "params": {
                "in": "opo.dn.delta_free", "out": "drv.opo.detuning", "f_ugf_hz": 2e3, "lock_light_ratio": 1.0,
                "offset_rad_s": 0.0, "sensor_psd": 1e4, "toggles": {"lock": True, "sensor_noise": True}}},
            pt.pathtest("pathtest_white", "fib", "drv.s1.dT", (), "K", sigma=1e-3),
            pt.pathtest("pathtest_white", "laser", "laser.las.phi", (), "rad", sigma=0.3),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}, locks=True),
            pt.record(["lock.lock2.err", "psi.opo.sq", "psi.opo.pr", "beam.pr.phi", "drv.opo.detuning"], N)]
    rec, _, _ = pt.run_and_record(pipe, FS, N, seed=23, compiled=compiled)
    err, pr, sq = rec["lock.lock2.err"], rec["psi.opo.pr"], rec["psi.opo.sq"]
    k = b.kappa(X)
    assert np.std(rec["drv.opo.detuning"]) > 1e-4 * k  # Δ(t) is really there
    assert np.std(rec["beam.pr.phi"]) > 10 * np.std(pr)  # the shared fibre dominates the beam phase ...
    scale = np.max(np.abs(rec["beam.pr.phi"]))
    # the OPO phase reaches the observer after half the OPO segment (μ₁, to its exit) plus the s1 group
    # delay (μ₂; emission time, conv §7), through ONE linear interpolation at μ₁ + μ₂ < 1 sample (Phase 5
    # F6-2: segment_phase's psi.<seg>.<b>.obs; path.segments.frac_delay)
    c = 299_792_458.0
    mu1, mu2 = 0.5 * 0.5 * 1.0 / c * FS, 3.0 * 1.468 / c * FS
    lin = lambda y, mu: (1 - mu) * y[1:] + mu * y[:-1]  # noqa: E731
    drv = rec["drv.opo.detuning"]
    l1 = lin(drv, mu1)  # the segment-exit reading of the (lumped) OPO subsegment
    c_pr = float(np.dot(pr[1:], l1) / np.dot(l1, l1))
    np.testing.assert_allclose(pr[1:], c_pr * l1, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(sq, pr * 2 / (1 + X**2), rtol=1e-12, atol=1e-15)
    c_res = c_pr * (1 - 2 / (1 + X**2))  # = c_pr − c_sq: the probe does not fully follow the ellipse
    np.testing.assert_allclose(err[1:], lin(c_res * drv, mu1 + mu2), rtol=1e-9, atol=1e-12 * scale)  # fibre cancels
    resid = pr * (1 - 2 / (1 + X**2))
    assert np.std(resid) / np.std(pr) == pytest.approx((1 - X**2) / (1 + X**2), rel=1e-9)
