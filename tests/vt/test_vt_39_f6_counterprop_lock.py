"""VT-39 (WP-F6 groundwork; the full VT-39 is WP-B8's ``test_vt_39_probe_cancellation.py``):
probe-cancellation breakdown computed by the path model (T10 §3.1, §3.2, §7 tests 2–4 and 9; T09 §5
tests 1–2; conv §9).

* T10 test 2 (point noise, ideal lock): residual PSD / free PSD = 4 sin²(2πf n_g|x − x_a|/c), tolerance
  2 % (achieved: exact to rounding, integer delays); x = x_a → 0.
* T10 test 3 (distributed noise, actuator at the end, L = 63 m, n_g = 1.468): ratio = 2[1 − sinc(4πf
  n_g L/c)]; ratio 1 at f₁ = 0.489 MHz ± 5 %; −53 dB ± 1 dB at 1 kHz.
* T10 test 4 (axis difference, PM, orthogonal probe, ideal PZT lock, dB/dT = −4×10⁻⁷ /K, L = 63 m,
  λ0 = 1545.32 nm, ΔT = 1 K): squeezed phase change K_T·L·ΔT = −102.5 rad ± 0.1 % (negative: conv §9);
  TEC actuator → 0; 90° splice at the midpoint → 0; first half heated → −51.2 rad.
* T09 test 1: k_T = −1.37×10⁻³ /K, L_b = 4.2 mm, 1 mK → −0.1291 rad (1e−6); same-polarization probe 0
  (1e−12).
* T09 test 2 / T10 test 9 (reciprocity): 100 random SU(2) segments; same-polarization backward probe
  phase = forward phase (1e−12), |t_s| = |t_p|; orthogonal probe error = 2·arg(a), a = xᵀUx (1e−12);
  PM (diagonal U) → 2·arg(a) = Δψ_b.

Ideal lock (bandwidth ∞, delay 0) is realised by replay: the actuator drive is a(t) = −φ_pr,free(t +
τ_{x_a}) (the probe passes the actuator τ_{x_a} before it leaves the segment), fed back through the
actuator's subsegment placement in a second run. Noise is a unit impulse per subsegment (staggered in
time) so the transfer functions are computed exactly from the block outputs (DTFT), instead of with a
noisy Welch estimate.
"""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.path.topology import C_M_PER_S

pytestmark = pytest.mark.numba
N_G = 1.468


def dtft(h, f, dt):
    n = np.arange(len(h))
    return np.exp(-2j * np.pi * np.outer(np.atleast_1d(f), n) * dt) @ h


def counterprop_topology(L, N, noise_placement, act_seg_len=None, act_sub=None):
    """Common segment c (x = 0 upstream … L at the GBS end); squeezed forward, probe backward. The
    actuator is either a zero-length segment `act` at x = 0 (`act_seg_len=0`) or a point on c."""
    drivers = {"strain": noise_placement}
    segs = []
    if act_seg_len is not None:
        segs.append(pt.seg_spec("act", act_seg_len, drivers={"stretcher": "uniform"}))
        coef_act = pt.coef_spec("act", "stretcher", "*", 1.0)
    else:
        drivers["stretcher"] = "actuator"
        coef_act = pt.coef_spec("c", "stretcher", "*", 1.0)
    segs.insert(0, pt.seg_spec("c", L, n_sub=N, spatial="mixed", drivers=drivers, actuator_sub=act_sub))
    act = [pt.step_spec("act", 1)] if act_seg_len is not None else []
    act_b = [pt.step_spec("act", -1)] if act_seg_len is not None else []
    beams = [pt.beam_spec("sq", [*act, pt.step_spec("c", 1)], "gbs"),
             pt.beam_spec("pr", [pt.step_spec("c", -1), *act_b], "det")]
    return pt.topo_spec(segs, beams, ["gbs", "det"], pt.coef_spec("c", "strain", "*", 1.0) + coef_act)


def pipeline(topo, N, n_rows, strain_shape, act_slot):
    return [pt.pathtest("pathtest_table", "noise", "drv.c.strain", strain_shape, "1", n_rows=n_rows),
            pt.pathtest("pathtest_table", "act", act_slot, (), "V", n_rows=n_rows),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}),
            pt.record(["beam.sq.phi", "beam.pr.phi"], n_rows)]


def ideal_lock(pipe, fs, n, noise_table, act_slot, k_a):
    free, _, _ = pt.run_and_record(pipe, fs, n, {"drv.c.strain_table": noise_table}, compiled=True)
    a = np.zeros(n)
    a[: n - k_a] = -free["beam.pr.phi"][k_a:]
    locked, _, _ = pt.run_and_record(pipe, fs, n, {"drv.c.strain_table": noise_table, f"{act_slot}_table": a},
                                     compiled=True)
    return free, locked


def test_t10_point_noise_ideal_lock():
    N, fs = 20, 1.0e9
    L = N * 2 * C_M_PER_S / (N_G * fs)  # δ = 2 samples
    tau_sub = N_G * L / C_M_PER_S / N
    for x, xa in [(14, 3), (5, 11), (7, 7)]:
        topo = counterprop_topology(L, N, {"point": x}, act_sub=xa)
        n = 300
        pipe = pipeline(topo, N, n, (), "drv.c.stretcher")
        noise = np.zeros(n)
        noise[80] = 1.0
        k_a = round((xa + 0.5) * tau_sub * fs)
        free, locked = ideal_lock(pipe, fs, n, noise, "drv.c.stretcher", k_a)
        assert np.max(np.abs(locked["beam.pr.phi"][: n - k_a - 1])) < 1e-12  # the probe is held
        dt = 1 / fs
        sep = abs(x - xa) * tau_sub  # n_g|x − x_a|/c
        if x == xa:
            assert np.max(np.abs(locked["beam.sq.phi"])) < 1e-12
            continue
        f = np.linspace(0.01, 0.3, 30) / sep
        ratio = np.abs(dtft(locked["beam.sq.phi"], f, dt)) ** 2 / np.abs(dtft(free["beam.sq.phi"], f, dt)) ** 2
        np.testing.assert_allclose(ratio, 4 * np.sin(2 * np.pi * f * sep) ** 2, rtol=0.02)


def staggered_impulses(N, G, n0):
    n = n0 + N * G
    tab = np.zeros((n, N))
    for m in range(N):
        tab[n0 + m * G, m] = 1.0
    return tab


def windows(x, N, G, n0):
    return [x[n0 + m * G: n0 + (m + 1) * G] for m in range(N)]


def test_t10_distributed_noise_ideal_lock():
    L, N = 63.0, 50
    tau_L = N_G * L / C_M_PER_S
    fs = 2 * N / tau_L  # δ = 2 samples: subsegment centres on integer delays
    G, n0 = 4 * N + 10, 5
    topo = counterprop_topology(L, N, "per_sub", act_seg_len=0.0)
    noise = staggered_impulses(N, G, n0)
    n = noise.shape[0]
    pipe = pipeline(topo, N, n, (N,), "drv.act.stretcher")
    free, locked = ideal_lock(pipe, fs, n, noise, "drv.act.stretcher", 0)
    dt = 1 / fs
    h_free = windows(free["beam.sq.phi"], N, G, n0)
    h_res = windows(locked["beam.sq.phi"], N, G, n0)

    def ratio(f):
        num = sum(np.abs(dtft(h, f, dt)) ** 2 for h in h_res)
        den = sum(np.abs(dtft(h, f, dt)) ** 2 for h in h_free)
        return num / den

    def expected(f):
        y = 4 * np.pi * f * tau_L
        return 2 * (1 - np.sin(y) / y)

    f = np.array([1e3, 1e4, 1e5, 3e5, 6e5, 1e6])
    np.testing.assert_allclose(ratio(f), expected(f), rtol=0.01)
    assert 10 * np.log10(ratio(1e3)[0]) == pytest.approx(-53.0, abs=1.0)
    lo, hi = 0.2e6, 0.8e6  # ratio crosses 1 at f₁ = 0.151 c/(n_g L)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ratio(mid)[0] < 1 else (lo, mid)
    assert lo == pytest.approx(0.489e6, rel=0.05)
    assert lo == pytest.approx(0.151 * C_M_PER_S / (N_G * L), rel=0.01)


# ------------------------------------------------------------------------------- axis difference


def pm_params(K_T, axis_drivers=("dT", "tec")):
    return {"K_T_rad_per_m_K": K_T, "K_tau_s_per_m_K": 0.0, "b_g_s_per_m": 0.0, "axis_drivers": list(axis_drivers)}


def static_axis_residual(segments, beams, coef, steps_dT, fs=1e8, extra=()):
    """Static ideal PZT lock with ρ_act = 1: residual = Δφ_sq − Δφ_pr after the drivers settle."""
    topo = pt.topo_spec(segments, beams, ["o"], coef)
    n = 400
    pipe = [pt.pathtest("pathtest_table", f"t_{slot}", slot, (), "K", n_rows=n) for slot in steps_dT]
    pipe += [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
             *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(["beam.sq.phi", "beam.pr.phi"], n)]
    tables = {}
    for slot, val in steps_dT.items():
        tab = np.zeros(n)
        tab[50:] = val
        tables[f"{slot}_table"] = tab
    rec, _, _ = pt.run_and_record(pipe, fs, n, tables)
    d = rec["beam.sq.phi"] - rec["beam.pr.phi"]
    return d[-1] - d[0], rec


def test_t10_axis_difference_pm_orthogonal_probe():
    lam0, L, C_T = 1545.32e-9, 63.0, 44.0
    K_T = 2 * np.pi / lam0 * (-4e-7)  # (2π/λ0)·dB/dT, dB/dT < 0 for stress birefringence (conv §9)
    seg = pt.seg_spec("c", L, n_sub=4, drivers={"dT": "uniform", "tec": "uniform"}, pol_model="pm_diagonal",
                      pol_params=pm_params(K_T))
    beams = [pt.beam_spec("sq", [pt.step_spec("c", 1, axis="slow")], "o"),
             pt.beam_spec("pr", [pt.step_spec("c", -1, axis="fast")], "o")]
    coef = pt.coef_spec("c", "dT", "*", C_T * L) + pt.coef_spec("c", "tec", "*", C_T * L)
    res, _ = static_axis_residual([seg], beams, coef, {"drv.c.dT": 1.0, "drv.c.tec": 0.0})
    assert res == pytest.approx(K_T * L, rel=1e-9) and res < 0
    assert abs(res) == pytest.approx(102.5, rel=1e-3)
    # same polarization (probe on the slow axis too): 0
    same = [beams[0], pt.beam_spec("pr", [pt.step_spec("c", -1, axis="slow")], "o")]
    res_same, _ = static_axis_residual([seg], same, coef, {"drv.c.dT": 1.0, "drv.c.tec": 0.0})
    assert abs(res_same) < 1e-12 * C_T * L
    # TEC as actuator: it writes temperature, so the ideal TEC lock (ΔT_tec = −ΔT) leaves 0 on both axes
    res_tec, rec = static_axis_residual([seg], beams, coef, {"drv.c.dT": 1.0, "drv.c.tec": -1.0})
    assert abs(res_tec) < 1e-12 * C_T * L and abs(rec["beam.pr.phi"][-1] - rec["beam.pr.phi"][0]) < 1e-9
    # 90° splice at the midpoint: slow/fast swap between the halves
    halves = [pt.seg_spec(h, L / 2, n_sub=2, drivers={"dT": "uniform"}, pol_model="pm_diagonal",
                          pol_params=pm_params(K_T)) for h in ("h1", "h2")]
    sw_beams = [pt.beam_spec("sq", [pt.step_spec("h1", 1, axis="slow"), pt.step_spec("h2", 1, axis="fast")], "o"),
                pt.beam_spec("pr", [pt.step_spec("h2", -1, axis="slow"), pt.step_spec("h1", -1, axis="fast")], "o")]
    sw_coef = pt.coef_spec("h1", "dT", "*", C_T * L / 2) + pt.coef_spec("h2", "dT", "*", C_T * L / 2)
    res_sw, _ = static_axis_residual(halves, sw_beams, sw_coef, {"drv.h1.dT": 1.0, "drv.h2.dT": 1.0})
    assert abs(res_sw) < 1e-12 * C_T * L
    res_half, _ = static_axis_residual(halves, sw_beams, sw_coef, {"drv.h1.dT": 1.0, "drv.h2.dT": 0.0})
    assert res_half == pytest.approx(K_T * L / 2, rel=1e-9)
    assert abs(res_half) == pytest.approx(51.2, rel=2e-3)


def test_t09_pm_differential_drift():
    L, L_b, k_T = 63.0, 4.2e-3, -1.37e-3  # T09 test 1 (λ0 = 1545.32 nm enters through L_b)
    K_T = k_T * 2 * np.pi / L_b
    seg = pt.seg_spec("c", L, n_sub=3, drivers={"dT": "uniform"}, pol_model="pm_diagonal", pol_params=pm_params(K_T, ["dT"]))
    beams = [pt.beam_spec("sq", [pt.step_spec("c", 1, axis="slow")], "o"),
             pt.beam_spec("pr", [pt.step_spec("c", -1, axis="fast")], "o")]
    coef = pt.coef_spec("c", "dT", "*", 44.0 * L)
    res, _ = static_axis_residual([seg], beams, coef, {"drv.c.dT": 1e-3})
    assert res == pytest.approx(-0.1291, abs=1e-4)
    assert res == pytest.approx(K_T * L * 1e-3, rel=1e-6)
    same = [beams[0], pt.beam_spec("pr", [pt.step_spec("c", -1, axis="slow")], "o")]
    res_same, _ = static_axis_residual([seg], same, coef, {"drv.c.dT": 1e-3})
    assert abs(res_same) < 1e-12


# ------------------------------------------------------------------------------------ reciprocity


def haar_su2(rng):
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    q = q * (np.diag(r) / np.abs(np.diag(r)))
    return q / np.sqrt(np.linalg.det(q))


def reciprocity_phases(mats):
    names = [f"s{i}" for i in range(len(mats))]
    segs = [pt.seg_spec(nm, 0.63, pol_model="jones", pol_params={"J": m.tolist()}, kind="smf") for nm, m in zip(names, mats)]
    fwd = [pt.step_spec(nm, 1) for nm in names]
    bwd = [pt.step_spec(nm, -1) for nm in reversed(names)]
    x, y = [1, 0], [0, 1]
    beams = [pt.beam_spec("sq", fwd, "gbs", jones_in=x, jones_out=x),
             pt.beam_spec("pr_same", bwd, "src", jones_in=x, jones_out=x),
             pt.beam_spec("pr_orth", bwd, "src", jones_in=y, jones_out=y)]
    topo = pt.topo_spec(segs, beams, ["gbs", "src"], [])
    slots = [f"beam.{b}.{f}" for b in ("sq", "pr_same", "pr_orth") for f in ("phi", "power")]
    pipe = [pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1),
            *pt.path_pipeline(topo, {"las": "laser.las.phi"}), pt.record(slots, 3)]
    rec, _, _ = pt.run_and_record(pipe, 1e6, 3)
    return {k: v[-1] for k, v in rec.items()}


def wrap(x):
    return math.remainder(x, 2 * math.pi)


def test_t09_reciprocity_random_su2():
    rng = np.random.default_rng(20261031)
    mats = [haar_su2(rng) for _ in range(100)]
    r = reciprocity_phases(mats)
    U = np.eye(2)
    for m in mats:
        U = m @ U
    a = U[0, 0]
    assert abs(wrap(r["beam.sq.phi"] - r["beam.pr_same.phi"])) < 1e-12
    assert r["beam.sq.power"] == pytest.approx(r["beam.pr_same.power"], rel=1e-12)
    assert abs(wrap(r["beam.sq.phi"] - r["beam.pr_orth.phi"] - 2 * np.angle(a))) < 1e-12
    assert r["beam.sq.power"] == pytest.approx(1e-3 * abs(a) ** 2, rel=1e-12)
    # PM (diagonal U, axes aligned): 2·arg(a) = Δψ_b
    dpsi = 0.7
    r_pm = reciprocity_phases([np.diag([np.exp(0.5j * dpsi), np.exp(-0.5j * dpsi)])])
    assert abs(wrap(r_pm["beam.sq.phi"] - r_pm["beam.pr_orth.phi"] - dpsi)) < 1e-12
    assert abs(wrap(r_pm["beam.sq.phi"] - r_pm["beam.pr_same.phi"])) < 1e-12
