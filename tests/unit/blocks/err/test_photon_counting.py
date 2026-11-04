"""``photon_counting_fringe`` (Phase 6 WP6-C, D85-5): reference vs njit (VT-44, 1e-12), Poisson / click statistics,
saturation, dark counts, gate mask, dead time, windows and hold, the estimator's bias and variance against the
Cramér-Rao bound, and the photon-budget numbers of notes/lit-pittaluga2021.md §3.3."""

import math

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from cranberry.blocks.err.photon_counting import fisher_per_gate, sqvac_background
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

FS = 1.0e6
ON = {"counting_noise": True, "saturation": True, "dark": True, "background": True, "dead_time": False,
      "squeezed_vacuum": True}


def phase(name, x):
    return {"block": "ramp", "name": name, "params": {"out": f"v.{name}", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                     "period_s": 1.0, "offset_v": x, "toggles": {"enabled": True}}}


def pcf(toggles=None, mode="gates", win=1000, f_gate=FS, cal=1.0, in_b="v.b", gate_in=None, **kw):
    p = {"in_a": "v.a", "in_b": in_b, "bg_in": None, "gate_in": gate_in, "out": "v.pcf", "f_gate_hz": f_gate,
         "cal_fraction": cal, "window_mode": mode, "mu_per_gate": 3.2, "visibility": 0.95, "qe": 0.93,
         "dark_per_gate": 1e-3, "bg_per_gate": 0.0, "phi_offset_rad": 0.0, "setpoint_rad": math.pi / 2,
         "dead_time_s": 4.3e-8, "recovery_s": 1e-8, "sqvac_per_gate": 0.0, "toggles": {**ON, **(toggles or {})}}
    p["window_t_s" if mode == "time" else "window_n_gates"] = win
    p.update(kw)
    return {"block": "photon_counting_fringe", "name": "pcf", "params": p}


def make(objs, seed=7):
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(objs), lay, seed)
    return blk, State(lay), Params(lay)


def windows(objs, n_windows, steps_per_window, compiled=True):
    """(n0, n1, n_gates, x_hat) of every completed window."""
    blk, st, pr = make(objs)
    step = run_compiled if compiled else run_reference
    out = []
    for _ in range(n_windows):
        step(blk, st, pr, steps_per_window)
        out.append([float(st[f"cnt.pcf.{k}"]) for k in ("n0", "n1", "n_gates", "x_hat")])
    return np.array(out)


def pipe(x=math.pi / 2, **kw):
    return [phase("a", x), phase("b", 0.0), pcf(**kw)]


@pytest.mark.numba
@pytest.mark.parametrize("kw", [dict(), dict(toggles={"saturation": False}), dict(toggles={"dead_time": True}),
                                dict(mode="time", win=3.7e-4, f_gate=2.3e6, cal=0.37),
                                dict(toggles={"counting_noise": False, "dead_time": True}, in_b=None)])
def test_reference_matches_njit(kw):
    objs = pipe(x=1.3, **{"win": 37, **kw})
    a, sa, pa = make(objs)
    b, sb, pb = make(objs)
    for n in (1, 500, 4000):
        run_reference(a, sa, pa, n)
        run_compiled(b, sb, pb, n)
    for name in sa.layout.slot_names():
        np.testing.assert_allclose(sb[name], sa[name], rtol=1e-12, atol=1e-300, err_msg=name)
    assert int(sa["cnt.pcf.window_index"]) > 10


@pytest.mark.numba
def test_poisson_statistics_pnr():
    g = 200
    w = windows(pipe(x=1.0, toggles={"saturation": False}, win=g), 3000, g)
    m = 0.5 * 0.93 * 3.2
    for col, s in ((1, 1.0), (0, -1.0)):
        lam = g * (m * (1 + s * 0.95 * math.cos(1.0)) + 1e-3)
        assert w[:, col].mean() == pytest.approx(lam, rel=0.01)
        assert w[:, col].var() == pytest.approx(lam, rel=0.07)   # Poisson: variance = mean
    assert np.all(w[:, 2] == g)


@pytest.mark.numba
def test_click_saturation_and_dark():
    g = 1000
    for mu in (0.5, 30.0):
        w = windows(pipe(x=math.pi / 2, win=g, mu_per_gate=mu), 400, g)
        lam = 0.5 * 0.93 * mu + 1e-3
        p = -math.expm1(-lam)
        assert w[:, 1].mean() / g == pytest.approx(p, abs=4 * math.sqrt(p * (1 - p) / (g * 400)) + 1e-12)
        assert w[:, 1].max() <= g                          # at most one click per gate
    # dark counts alone (μ = 0): rate 1 − e^{−d}; dark toggle off → no clicks at all
    w = windows(pipe(win=g, mu_per_gate=0.0, dark_per_gate=0.02), 300, g)
    assert w[:, :2].mean() / g == pytest.approx(-math.expm1(-0.02), rel=0.05)
    w = windows(pipe(win=g, mu_per_gate=0.0, toggles={"dark": False}), 20, g)
    assert np.all(w[:, :2] == 0.0)
    # background adds qe·b per output; background toggle off removes it
    w = windows(pipe(win=g, mu_per_gate=0.0, dark_per_gate=0.0, bg_per_gate=0.05), 300, g)
    assert w[:, :2].mean() / g == pytest.approx(-math.expm1(-0.93 * 0.05), rel=0.05)
    w = windows(pipe(win=g, mu_per_gate=0.0, dark_per_gate=0.0, bg_per_gate=0.05, toggles={"background": False}), 5, g)
    assert np.all(w[:, :2] == 0.0)


def test_expected_counts_and_exact_estimate():
    """counting_noise off: expected counts G·p and X̂ = X exactly (MLE inverse of the click model)."""
    for x in (0.4, 1.2, math.pi / 2, 2.6):
        w = windows(pipe(x=x, win=50, toggles={"counting_noise": False}), 2, 50, compiled=False)
        lam1 = 0.5 * 0.93 * 3.2 * (1 + 0.95 * math.cos(x)) + 1e-3
        assert w[-1, 1] == pytest.approx(50 * -math.expm1(-lam1), rel=1e-12)
        assert w[-1, 3] == pytest.approx(x, abs=1e-9)


def test_window_hold_mask_and_offset():
    objs = [phase("a", 0.2), phase("m", 0.0), pcf(win=1e-4, mode="time", f_gate=2e5, cal=0.5, in_b=None,
                                                  gate_in="v.m", phi_offset_rad=1.0,
                                                  toggles={"counting_noise": False})]
    blk, st, pr = make(objs)
    run_reference(blk, st, pr, 250)                   # mask 0: windows close with no gate, out held at 0
    assert int(st["cnt.pcf.window_index"]) == 2 and float(st["cnt.pcf.n_gates"]) == 0.0
    assert float(st["v.pcf"]) == 0.0
    pr["m.offset_v"] = 1.0
    run_reference(blk, st, pr, 150)
    assert float(st["cnt.pcf.n_gates"]) == pytest.approx(10.0, abs=1.0)   # 1e5 gates/s × 100 µs
    assert float(st["cnt.pcf.x_hat"]) == pytest.approx(1.2, abs=1e-9)   # X = φ_a + offset (one input)
    assert float(st["v.pcf"]) == pytest.approx(1.2 - math.pi / 2, abs=1e-9)
    held = float(st["v.pcf"])
    run_reference(blk, st, pr, 99)                    # inside the next window: held
    assert float(st["v.pcf"]) == held


@pytest.mark.numba
def test_dead_time():
    """Gates 1 µs apart: a dead time > 1 µs blocks the gate after every click; a short one changes nothing."""
    g = 2000
    long = windows(pipe(win=g, mu_per_gate=30.0, dead_time_s=1.5e-6, toggles={"dead_time": True}), 50, g)
    p = -math.expm1(-(0.5 * 0.93 * 30.0 + 1e-3))
    assert long[:, 1].mean() / g == pytest.approx(p / (1 + p), rel=0.02)   # alternate click / dead gate
    short = windows(pipe(win=g, mu_per_gate=30.0, dead_time_s=1e-9, recovery_s=1e-9, toggles={"dead_time": True}), 5, g)
    assert short[:, 1].mean() / g == pytest.approx(p, rel=1e-3)


@pytest.mark.numba
@pytest.mark.parametrize("mu,x,g", [(3.2, math.pi / 2, 2000), (0.3, math.pi / 2, 2000), (3.2, math.pi / 2 + 0.3, 2000),
                                    (12.0, math.pi / 2, 20000)])
def test_estimator_bias_and_variance_vs_cramer_rao(mu, x, g):
    """Asymptotic regime (≳ 50 empty gates per output and window): the MLE is unbiased and reaches the CRB. Deep in
    saturation (μ = 12: e^{−5.6} empty gates) G = 2000 leaves ≈ 7 empties per window and the variance is ≈ 1.23 × the
    CRB (small-count regime), hence G = 20000 there. The bound is the μ-unknown one (the normalised estimator does not
    know m = qe μ/2): 2×2 Fisher matrix in (m, X); at π/2 it equals the μ-known 1/F_XX, off π/2 it is larger (+7.5 % at
    +0.3 rad, review C-5), which the 10 % tolerance would otherwise hide."""
    w = windows(pipe(x=x, win=g, mu_per_gate=mu), 1500, g)
    crb = cramer_rao_mu_unknown(0.5 * 0.93 * mu, 0.95, x, 1e-3) / g
    if abs(x - math.pi / 2) > 0.1:
        assert crb > 1.05 / (g * fisher_per_gate(0.5 * 0.93 * mu, 0.95, x, 1e-3))
    err = w[:, 3] - x
    sem = math.sqrt(crb / len(err))
    assert abs(err.mean()) < 4 * sem + 0.05 * math.sqrt(crb)        # unbiased within statistics
    assert err.var() == pytest.approx(crb, rel=0.1)                # efficient: reaches the CRB


def cramer_rao_mu_unknown(m, vis, x, dark):
    """Var(X̂) bound per gate with m unknown: [F⁻¹]_XX of the click model λ± = m(1 ± V cos X) + d."""
    f = np.zeros((2, 2))
    for s in (1.0, -1.0):
        lam = m * (1.0 + s * vis * math.cos(x)) + dark
        g = np.array([1.0 + s * vis * math.cos(x), -s * m * vis * math.sin(x)])
        f += np.outer(g, g) / math.expm1(lam)
    return float(np.linalg.inv(f)[1, 1])


def test_squeezed_vacuum_background():
    """b' of the review's table (C-2: r = 1.8, η·qe 0.59 / 0.68 → 1.678 / 1.773; 0 without squeezing or overlap); in the
    block it adds b' clicks per output and, being known, is subtracted: X̂ = X exactly with expected counts."""
    assert sqvac_background(1.8, 0.59) == pytest.approx(1.678, abs=1e-3)
    assert sqvac_background(1.8, 0.68) == pytest.approx(1.773, abs=1e-3)
    assert sqvac_background(0.0, 0.7) == pytest.approx(0.0, abs=1e-15)
    assert sqvac_background(1.8, 0.7, 0.0) == pytest.approx(0.0, abs=1e-15)
    bp = sqvac_background(1.8, 0.63)
    g = 400
    w = windows(pipe(x=1.2, win=g, mu_per_gate=0.0, dark_per_gate=0.0, sqvac_per_gate=bp,
                     toggles={"counting_noise": False}), 2, g, compiled=False)
    assert w[-1, 1] / g == pytest.approx(-math.expm1(-bp), rel=1e-12)
    w = windows(pipe(x=1.2, win=g, sqvac_per_gate=bp, toggles={"counting_noise": False}), 2, g, compiled=False)
    assert w[-1, 3] == pytest.approx(1.2, abs=1e-9)
    w = windows(pipe(x=1.2, win=g, mu_per_gate=0.0, dark_per_gate=0.0, sqvac_per_gate=bp,
                     toggles={"counting_noise": False, "squeezed_vacuum": False}), 2, g, compiled=False)
    assert np.all(w[-1, :2] == 0.0)


def test_unequal_seed_intensities():
    """ia_in / ib_in: λ± = qe μ(I_a + I_b ± 2V√(I_a I_b) cos X)/2 + dark; balanced (½, ½) equals the unwired block; the
    estimator, not told the contrast 2√(I_a I_b)/(I_a + I_b), reads cos X̂ scaled by it (no bias at π/2)."""
    def run(ia, ib, x):
        objs = [phase("a", x), phase("b", 0.0), phase("ia", ia), phase("ib", ib),
                pcf(win=50, ia_in="v.ia", ib_in="v.ib", toggles={"counting_noise": False})]
        blk, st, pr = make(objs)
        run_reference(blk, st, pr, 100)
        return float(st["cnt.pcf.n1"]), float(st["cnt.pcf.x_hat"])
    n_bal, x_bal = run(0.5, 0.5, 1.0)
    w = windows(pipe(x=1.0, win=50, toggles={"counting_noise": False}), 2, 50, compiled=False)
    assert n_bal == pytest.approx(w[-1, 1], rel=1e-12) and x_bal == pytest.approx(1.0, abs=1e-9)
    ia, ib = 0.9, 0.1
    n1, xh = run(ia, ib, 1.0)
    lam = 0.93 * 3.2 * (0.5 * (ia + ib) + math.sqrt(ia * ib) * 0.95 * math.cos(1.0)) + 1e-3
    assert n1 == pytest.approx(50 * -math.expm1(-lam), rel=1e-12)
    c = 2 * math.sqrt(ia * ib) / (ia + ib)
    assert math.cos(xh) == pytest.approx(c * math.cos(1.0), rel=1e-9)
    assert run(ia, ib, math.pi / 2)[1] == pytest.approx(math.pi / 2, abs=1e-9)


def test_photon_budget_lit_pittaluga():
    """lit-pittaluga2021 §3.3: information per gate peaks at μ_tot ≈ 3.2 detected photons (0.648 V²); 10 mrad at
    V = 0.95 needs ≈ 8.5e3 calibration gates → 45–140 ms when every port is calibrated in every frame."""
    res = minimize_scalar(lambda m: -fisher_per_gate(m, 1.0, math.pi / 2), bounds=(0.1, 10.0), method="bounded")
    assert 2 * res.x == pytest.approx(3.19, abs=0.01)
    assert fisher_per_gate(res.x, 1.0, math.pi / 2) / 2 == pytest.approx(0.648, abs=1e-3)
    # the GBS SNSPD (qe 0.93) moves the optimum to μ = 3.43 photons at the detector, the maximum is unchanged
    n_g = 1.0 / (fisher_per_gate(res.x, 0.95, math.pi / 2) * 0.01**2)
    assert n_g == pytest.approx(8.5e3, rel=0.01)
    t_frame = {"F_liu20": 32 * 4 / 20e6, "F_wg6": 32 * 1 / 6e6, "F_bor6": 32 * 1 / 6e6, "F_wg100": 32 * 50 / 100e6}
    t_int = {k: n_g * v for k, v in t_frame.items()}
    assert t_int["F_wg6"] == pytest.approx(0.0456, rel=0.01)
    assert t_int["F_liu20"] == pytest.approx(0.0547, rel=0.01)
    assert t_int["F_wg100"] == pytest.approx(0.137, rel=0.01)
    # μ_tot = 0.2 (m = 0.1): 5.8e4 gates (×6.8)
    assert 1.0 / (fisher_per_gate(0.1, 0.95, math.pi / 2) * 0.01**2) == pytest.approx(5.83e4, rel=0.01)
