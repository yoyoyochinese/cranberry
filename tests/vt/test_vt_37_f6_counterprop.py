"""VT-37 (WP-F6 part; the noise-source parts are WP-B2's ``test_vt_37_environmental_noise.py``): the
path model's handling of environmental drivers (T08 §9 V2, V6, V11, V12; conv §9).

* V6 counter-propagation (T08 eq. 2.11): with N_sub independent white subsegment drivers
  (distributed), the difference D(t) = ψ_fwd(t) − ψ_bwd(t) of the forward beam leaving at x = L and the
  backward beam leaving at x = 0 at the same time has S_D = 2S[1 − sinc(2πfτ)] (S: PSD of ψ_fwd).
  Checked at 2πfτ = 0.3, 1, π, 2π: exactly from impulse responses (tolerance 2 %, tighter than the
  note's ±10 %) and statistically with a Welch estimate of a seeded white-noise run (±10 %). A lumped
  (uniform) driver gives S_D = 0 (exact to rounding). This observable is separate from the T10 ideal-
  lock residual 2[1 − sinc(4πf n_g L/c)] (VT-39 F6 file).
* V2: L = 63 m, ΔT = 1 K, κ_τ = 39 ps/(km·K) → Δτ_g = 2.457 ps.
* V11: segment_phase toggles off → ψ ≡ 0 (exact), τ_g constant (nominal).
* V12: same ΔT on the slow and fast axis → ψ_slow − ψ_fast = K_T·L·ΔT exactly (conv §9 sign, K_T < 0);
  a driver not in ``axis_drivers`` (strain) gives only the difference of its coefficients (0 here).
"""

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.core.fft import psd_onesided
from cranberry.path.topology import C_M_PER_S

pytestmark = pytest.mark.numba
N_G = 1.468


def dtft(h, f, dt):
    n = np.arange(len(h))
    return np.exp(-2j * np.pi * np.outer(np.atleast_1d(f), n) * dt) @ h


def two_way(L, N, placement, spatial):
    segs = [pt.seg_spec("c", L, n_sub=N, spatial=spatial, drivers={"strain": placement})]
    beams = [pt.beam_spec("fwd", [pt.step_spec("c", 1)], "xL"), pt.beam_spec("bwd", [pt.step_spec("c", -1)], "x0")]
    return pt.topo_spec(segs, beams, ["xL", "x0"], pt.coef_spec("c", "strain", "*", 1.0))


def pipe_for(topo, driver_block, n):
    return [driver_block, pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.fwd.phi", "beam.bwd.phi"], n)]


def expected(x):
    return 2 * (1 - np.sin(x) / x)


X = np.array([0.3, 1.0, np.pi, 2 * np.pi])  # 2πfτ


def test_v6_counterprop_exact_impulse_responses():
    L, N = 63.0, 100
    tau = N_G * L / C_M_PER_S
    fs = 2 * N / tau  # subsegment centres on integer sample delays
    G, n0 = 2 * N + 6, 3
    n = n0 + N * G
    tab = np.zeros((n, N))
    for m in range(N):
        tab[n0 + m * G, m] = 1.0
    topo = two_way(L, N, "per_sub", "distributed")
    drv = pt.pathtest("pathtest_table", "noise", "drv.c.strain", (N,), "1", n_rows=n)
    rec, _, _ = pt.run_and_record(pipe_for(topo, drv, n), fs, n, {"drv.c.strain_table": tab}, compiled=True)
    f = X / (2 * np.pi * tau)
    num = den = 0.0
    for m in range(N):
        w = slice(n0 + m * G, n0 + (m + 1) * G)
        hf, hb = rec["beam.fwd.phi"][w], rec["beam.bwd.phi"][w]
        num = num + np.abs(dtft(hf - hb, f, 1 / fs)) ** 2
        den = den + np.abs(dtft(hf, f, 1 / fs)) ** 2
    np.testing.assert_allclose(num / den, expected(X), rtol=0.02)


def test_v6_counterprop_welch_white_noise():
    N, n = 20, 2**20
    fs = 1.0e9
    L = 2 * N * C_M_PER_S / (N_G * fs)  # τ = 2N samples
    tau = N_G * L / C_M_PER_S
    topo = two_way(L, N, "per_sub", "distributed")
    drv = pt.pathtest("pathtest_white", "noise", "drv.c.strain", (N,), "1", sigma=1.0)
    rec, _, _ = pt.run_and_record(pipe_for(topo, drv, n), fs, n, compiled=True, seed=37)
    fwd, bwd = rec["beam.fwd.phi"][1000:], rec["beam.bwd.phi"][1000:]
    nper = 2**13
    f, s_d = psd_onesided(fwd - bwd, fs, nperseg=nper)
    _, s_f = psd_onesided(fwd, fs, nperseg=nper)
    for x in X:
        k = int(round(x / (2 * np.pi * tau) / (f[1] - f[0])))
        band = slice(k - 3, k + 4)
        ratio = np.mean(s_d[band]) / np.mean(s_f[band])
        assert ratio == pytest.approx(float(np.mean(expected(2 * np.pi * f[band] * tau))), rel=0.10)


def test_v6_lumped_cancels_exactly():
    N, n, fs = 25, 3000, 1.0e9
    L = 3 * N * C_M_PER_S / (N_G * fs) * 1.137  # fractional subsegment delays
    topo = two_way(L, N, "uniform", "lumped")
    drv = pt.pathtest("pathtest_white", "noise", "drv.c.strain", (), "1", sigma=1.0)
    rec, _, _ = pt.run_and_record(pipe_for(topo, drv, n), fs, n, compiled=True, seed=5)
    assert np.std(rec["beam.fwd.phi"]) > 0.05
    assert np.max(np.abs(rec["beam.fwd.phi"] - rec["beam.bwd.phi"])) < 1e-13


def axis_topology(K_T):
    pm = {"K_T_rad_per_m_K": K_T, "K_tau_s_per_m_K": -1.78e-15, "b_g_s_per_m": 1.36e-12, "axis_drivers": ["dT"]}
    segs = [pt.seg_spec("c", 63.0, n_sub=2, drivers={"dT": "uniform", "strain": "uniform"}, pol_model="pm_diagonal",
                        pol_params=pm)]
    beams = [pt.beam_spec("s", [pt.step_spec("c", 1, axis="slow")], "o"),
             pt.beam_spec("f", [pt.step_spec("c", 1, axis="fast")], "o")]
    coef = pt.coef_spec("c", "dT", "*", 44.0 * 63.0, 39e-15 * 63.0) + pt.coef_spec("c", "strain", "*", 7.0e6)
    return pt.topo_spec(segs, beams, ["o"], coef)


def run_axis(topo, dT, strain, n=60, **toggles):
    pipe = [pt.pathtest("pathtest_table", "t", "drv.c.dT", (), "K", n_rows=n),
            pt.pathtest("pathtest_table", "e", "drv.c.strain", (), "1", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}, **toggles),
            pt.record(["psi.c.s", "psi.c.f", "tau_g.c.s", "tau_g.c.f", "beam.s.phi"], n)]
    return pt.run_and_record(pipe, 1e8, n, {"drv.c.dT_table": np.full(n, dT), "drv.c.strain_table": np.full(n, strain)})[0]


def test_v12_axis_coefficients_and_v2_group_delay():
    K_T = -2.05
    r = run_axis(axis_topology(K_T), 1.0, 0.0)
    assert r["psi.c.s"][-1] - r["psi.c.f"][-1] == pytest.approx(K_T * 63.0, rel=1e-12)
    assert r["psi.c.s"][-1] == pytest.approx(44.0 * 63.0, rel=1e-12)
    # V2: Δτ_g = κ_τ L ΔT on the slow (reference) axis
    tau0 = 1.468 * 63.0 / C_M_PER_S
    assert r["tau_g.c.s"][-1] - tau0 == pytest.approx(2.457e-12, rel=1e-3)
    # fast axis: nominal DGD −b_g L and its drift −K_τ L ΔT
    assert r["tau_g.c.f"][-1] - r["tau_g.c.s"][-1] == pytest.approx(-1.36e-12 * 63.0 + 1.78e-15 * 63.0, rel=1e-9)
    rs = run_axis(axis_topology(K_T), 0.0, 1e-6)
    assert rs["psi.c.s"][-1] == rs["psi.c.f"][-1] == pytest.approx(7.0, rel=1e-12)


def test_v11_toggles_off():
    r = run_axis(axis_topology(-2.05), 1.0, 1e-6, phase=False, group_delay=False, pm_axis_drift=False)
    assert np.all(r["psi.c.s"] == 0.0) and np.all(r["psi.c.f"] == 0.0) and np.all(r["beam.s.phi"] == 0.0)
    tau0 = 1.468 * 63.0 / C_M_PER_S
    assert np.all(r["tau_g.c.s"] == r["tau_g.c.s"][0]) and r["tau_g.c.f"][0] == pytest.approx(tau0, rel=1e-15)
