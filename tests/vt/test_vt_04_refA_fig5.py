"""VT-04 (prompt §7-2, WP-B1 part): Ref A Fig. 5, four curves × 8 pump powers and the optima.

Sources: S §4 VT-04; T01 V-T01-2, §7.2 table; T04 V2; D48 (§7-2 reading). Single-mode reference
model, α = 9.06 /W, lumped loss, small-angle average (the §7-1 form the figure uses).
Tolerances: ±0.01 dB to the T01 §7.2 table (values given to 0.01 dB); optima P_opt = 612 mW at
−10.152 dB (θ 9 mrad, L 8 %) and 503 mW at −8.397 dB (14 mrad, 12 %) to ±1 mW / ±1e−3 dB.
Excluded (D48): the right end (1200 mW) of the "previous work" dashed curve — the figure shows
≈ −4.8 dB against the computed −6.07 dB, so it is not a shape target.

Runner part (WP-I1): the four curves through the full time-loop simulation, presets/reference/refA.yaml
with the overlays presets/axes/refA_curve/{previous_work, blue, green}.yaml, pump power set on the
offline OPA (lock-2 gain ∝ A+(P) follows, as with fixed electronics). Tolerance ±0.10 dB (D48, as VT-02)
for P ≤ 1000 mW; at 1200 mW ±0.20 dB: the simulated θ is not constant in P — lock 2's loop gain grows
with A+ and a delayed loop raises the variance of the flat laser/filter-cavity noise it acts on (Bode
sensitivity integral), +1 % (9 mrad) / +3.4 % (14 mrad) at 1.2 W, which the constant-θ curves of the
figure do not contain (finding for Phase 5). The optimum of the simulated 9 mrad curve lies between
400 and 800 mW (P_opt = 612 mW analytically).

Digitised data points (WP-I1): Fig. 5 is an embedded 891 × 484 px raster (200 ppi) on PDF page 7,
extracted with ``pdfimages``. Axes calibrated on the gridlines (5 dB = 43.0 px, 200 mW = 86.4 px;
residuals < 0.5 px); a point is the bounding-box centre of a marker's blue outline (b > 100,
b − r > 60, b − g > 60, g < 110). Only markers whose outline is complete (≥ 11 × 11 px) are kept:
7 of 12 squeezed and 9 of 12 anti-squeezed (the rest are cut by the theory lines and would be biased
by up to 3 px). Uncertainty ±1 px = ±0.12 dB, small against the ±0.35 dB check (S §4 VT-04), which
the simulated 9 mrad / 8 % curve meets at every kept point (max |Δ| 0.26 dB, anti-squeezed at 551 mW).
"""

import copy

import numpy as np
import pytest

from cranberry.blocks.source.single_mode import make_single_mode_reference, optimal_pump_w

P_MW = [40, 100, 200, 400, 640, 800, 1000, 1200]
CURVES = {  # (θ [rad], L): [(R′− dB, R′+ dB)] per P, T01 §7.2
    "this_work": (0.009, 0.08, [(-4.48, 4.98), (-6.62, 7.96), (-8.43, 11.36), (-9.84, 16.18), (-10.15, 20.56),
                                (-9.98, 23.02), (-9.47, 25.78), (-8.67, 28.28)]),
    "previous_work": (0.014, 0.12, [(-4.15, 4.85), (-5.98, 7.80), (-7.40, 11.18), (-8.32, 15.99), (-8.29, 20.36),
                                    (-7.91, 22.83), (-7.14, 25.59), None]),  # right end excluded (D48)
    "blue": (0.0, 0.08, [(-4.49, 4.98), (-6.63, 7.96), (-8.47, 11.36), (-9.98, 16.18), (-10.58, 20.56),
                         (-10.75, 23.02), (-10.85, 25.78), (-10.90, 28.28)]),
    "green": (0.014, 0.08, [(-4.48, 4.98), (-6.61, 7.96), (-8.39, 11.35), (-9.64, 16.18), (-9.60, 20.56),
                            (-9.08, 23.02), (-8.06, 25.78), (-6.71, 28.28)]),
}


def curve_db(theta, loss, p_w):
    b = make_single_mode_reference(dict(loss=loss, alpha_per_w=9.06, pump_w=p_w, theta_rms_rad=theta,
                                        phi_opa0_rad=0.0, average="small_angle"))
    return b.r_prime_db()


@pytest.mark.parametrize("name", list(CURVES))
def test_fig5_curves(name):
    theta, loss, pts = CURVES[name]
    for p_mw, target in zip(P_MW, pts, strict=True):
        if target is None:
            continue
        rm, rp = curve_db(theta, loss, p_mw * 1e-3)
        assert rm == pytest.approx(target[0], abs=0.01), (name, p_mw)
        assert rp == pytest.approx(target[1], abs=0.01), (name, p_mw)


@pytest.mark.parametrize("theta, loss, p_opt_mw, rmin_db", [(0.009, 0.08, 612, -10.152), (0.014, 0.12, 503, -8.397)])
def test_optima(theta, loss, p_opt_mw, rmin_db):
    p = optimal_pump_w(9.06, theta)
    assert p * 1e3 == pytest.approx(p_opt_mw, abs=1.0)
    assert curve_db(theta, loss, p)[0] == pytest.approx(rmin_db, abs=1e-3)
    assert 10 * np.log10(loss + (1 - loss) * np.sin(2 * theta)) == pytest.approx(rmin_db, abs=1e-3)


# --------------------------------------------------------------------------- runner part (WP-I1)

OVERLAY = {"this_work": None, "previous_work": "refA_curve/previous_work", "blue": "refA_curve/blue",
           "green": "refA_curve/green"}
#: digitised Fig. 5 points (module docstring): (P [mW], dB)
DATA_SQZ = [(40, -4.65), (74, -5.69), (198, -8.31), (276, -9.18), (379, -9.76), (853, -9.88), (987, -9.41)]
DATA_ANTI = [(40, 4.89), (124, 8.90), (198, 11.29), (277, 13.44), (379, 15.65), (453, 17.28), (551, 19.31),
             (853, 23.67), (986, 25.59)]


def sim_db(p_mw, overlay=None, duration_s=0.03):
    """(R′− dB, R′+ dB, θ) of the refA time-loop simulation at pump power P (offline OPA, output reference)."""
    from cranberry import config
    from cranberry.config.yaml_io import load_yaml
    from cranberry.runner.run import run

    t = load_yaml("presets/reference/refA.yaml")
    t["run"]["duration_s"] = duration_s
    t["measure"].pop("esa")  # short record: no zero-span traces
    for b in t["offline"]["blocks"]:
        if b["name"] == "sqz_opa":
            b["params"]["pump_w"] = p_mw * 1e-3
    r = run(config.load(copy.deepcopy(t), overlays=[overlay] if overlay else []), keep_state=False)
    nl = r.metrics.noise_levels["i.port0"]
    return nl["r_minus_db"], nl["r_plus_db"], nl["sigma_total_rad"]


@pytest.mark.numba
@pytest.mark.integration
@pytest.mark.parametrize("name", list(CURVES))
def test_fig5_curves_through_runner(name):
    theta, _, pts = CURVES[name]
    sims = [sim_db(p, OVERLAY[name]) for p in P_MW]
    for p_mw, target, (rm, rp, th) in zip(P_MW, pts, sims, strict=True):
        if target is None:
            continue  # D48: right end of the dashed curve excluded
        tol = 0.10 if p_mw <= 1000 else 0.20
        assert rm == pytest.approx(target[0], abs=tol), (name, p_mw)
        assert rp == pytest.approx(target[1], abs=tol), (name, p_mw)
        if theta > 0 and p_mw >= 400:
            assert th == pytest.approx(theta, rel=0.05), (name, p_mw)
    if name == "this_work":
        rms = [s[0] for s in sims]
        assert min(rms) == rms[P_MW.index(640)] and rms[P_MW.index(400)] > rms[P_MW.index(640)]


@pytest.mark.numba
@pytest.mark.integration
def test_fig5_digitised_points_through_runner():
    for data, k in ((DATA_SQZ, 0), (DATA_ANTI, 1)):
        for p_mw, db in data:
            assert sim_db(p_mw)[k] == pytest.approx(db, abs=0.35), (k, p_mw)
