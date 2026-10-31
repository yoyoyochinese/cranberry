"""VT-18 (WP-B6 unit part; WP-I4 adds the preset run): Ref B definition discrimination and back-solve
(prompt §7-9 per D48; 00_summary §4 VT-18; T02 V-B4, V-B5, V-B6).

* V-B4: x = 0.999, θ̃ = 1.5°: R′+ = +21.4 dB with Ω_B = 2πf/γ_B; the paper's literal Ω = f/γ (2π missing,
  emulated by evaluating at f/2π) gives ≈ +37 dB — the discriminator of the Ω definition.
* V-B5: L_ic = 0.0038 (measured at 100 mW), θ̃ = 1.5°: R′∓ = −9.44 / +14.95 dB (±0.02); with the
  circuit noise −21.7 dB (pure loss η_el = 1/(1 + C)) −9.21 / +14.92 dB, within 0.25 dB of the observed
  −9.01 / +15.12 dB. Back-solve from the reported generated R∓ = −10.12 / +15.15 dB: x = 0.753 ± 0.002,
  total loss 0.0708 (paper 0.0709).
* V-B6: lossless, Ω = 0 (toggle frequency_dependence off): min_x R′− = 10 log10(sin 2θ̃) = −8.67 / −12.81 dB
  for 3.9° / 1.5° (paper −8.66 / −12.81, ±0.02).
"""

import numpy as np
import pytest
from scipy.optimize import brentq, minimize_scalar

from cranberry.blocks.source.opo_cw import closed_form
from cranberry.blocks.source.opo_testing import REF_B, REF_B_ETA_D, db, measured_db, opo

TH15 = np.radians(1.5)


def test_v_b4_two_pi_in_omega():
    good = measured_db(opo(REF_B, pump_mode="x", x=0.999), REF_B_ETA_D, TH15)[1][0]
    bad = measured_db(opo(REF_B, pump_mode="x", x=0.999, f_s_hz=[1e6 / (2 * np.pi)]), REF_B_ETA_D, TH15)[1][0]
    assert good == pytest.approx(21.4, abs=0.05)
    assert bad == pytest.approx(37.4, abs=0.2)


def test_v_b5_100mw_point_and_circuit_noise():
    rm, rp = measured_db(opo(REF_B, loss_ic=0.0038), REF_B_ETA_D, TH15)
    assert (float(rm[0]), float(rp[0])) == pytest.approx((-9.44, 14.95), abs=0.02)
    c = 10 ** (-21.7 / 10)  # circuit noise / shot noise
    with_c = [db((10 ** (v / 10) + c) / (1 + c)) for v in (float(rm[0]), float(rp[0]))]
    assert with_c == pytest.approx([-9.21, 14.92], abs=0.02)
    assert with_c == pytest.approx([-9.01, 15.12], abs=0.25)  # observed [Ref B §4]


def test_v_b5_back_solve():
    b = opo(REF_B, loss_ic=0.0038)
    w = 2 * np.pi * 1e6 / b.kappa(b.x)  # κ independent of x for a constant L_ic
    rm, rp = 10 ** (-10.12 / 10), 10 ** (15.15 / 10)

    def ratio(x):  # (R+ − 1)/(1 − R−) does not depend on η_tot
        m, p = closed_form(x, 1.0, w)
        return (p - 1) / (1 - m)

    x = brentq(lambda x: ratio(x) - (rp - 1) / (1 - rm), 0.3, 0.95, xtol=1e-12)
    eta = (rp - 1) / (closed_form(x, 1.0, w)[1] - 1)
    assert x == pytest.approx(0.753, abs=0.002)
    assert 1 - eta == pytest.approx(0.0708, abs=0.0005)


@pytest.mark.parametrize("th_deg, target", [(3.9, -8.67), (1.5, -12.81)])
def test_v_b6_ideal_limit_sin_2theta(th_deg, target):
    th = np.radians(th_deg)

    def rmin(x):
        b = opo(REF_B, loss_ic=0.0, pump_mode="x", x=x, toggles={"frequency_dependence": False})
        return float(measured_db(b, 1.0, th)[0][0])

    res = minimize_scalar(rmin, bounds=(0.0, 0.9999), method="bounded", options={"xatol": 1e-9})
    assert res.fun == pytest.approx(db(np.sin(2 * th)), abs=1e-4)
    assert res.fun == pytest.approx(target, abs=0.02)


# ------------------------------------------------------------------------------------------------------------------
# WP-I4: V-B4 / V-B5 through the refB preset and the runner (prompt §7-9 per D48).
#
# * V-B4 with the preset's opo block and the run's LossRegistry η: x = 0.999 → R′+ = +21.4 dB; the literal
#   Ω = f/γ (evaluated at f/2π) → +37.4 dB;
# * V-B5 with presets/axes/refB_loss/measured_100mw.yaml (L = 0.0038 [Ref B §3]): the run's jitter-free V∓ mapped
#   by Eq. (1) at 1.5° give −9.44 / +14.95 dB; with the HD circuit noise kept (toggle hd_det.subtract_dark off,
#   −21.7 dB as the pure loss η_el = C/(1 + C), T13) −9.21 / +14.92 dB (±0.02); the closed loop (θ̃ ≈ 1.5°) and its
#   ESA zero span are within 0.25 dB of the observed −9.01 / +15.12 dB [Ref B §4].
# ------------------------------------------------------------------------------------------------------------------

import importlib.util  # noqa: E402
from pathlib import Path  # noqa: E402

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parents[1] / "integration" / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)


@pytest.mark.numba
@pytest.mark.integration
def test_refB_preset_v_b4():
    from cranberry.runner.run import run

    res = run(rs.load(run={"duration_s": 0.012}))
    eta = res.offline.loss_registry.eta("sqz", "bhd")
    params = rs.offline_block_params(res.spec, "opo")
    p999 = 0.999**2 * rs.params_of(res.spec)["p_th_w"]
    good = measured_db(opo(params, pump_w=p999), eta, TH15)[1][rs.I_1MHZ]
    bad = measured_db(opo(params, pump_w=p999, f_s_hz=[1e6 / (2 * np.pi)]), eta, TH15)[1][0]
    assert good == pytest.approx(21.4, abs=0.05)
    assert bad == pytest.approx(37.4, abs=0.2)


@pytest.mark.numba
@pytest.mark.integration
def test_refB_preset_v_b5_measured_loss_and_circuit_noise():
    from cranberry.runner.run import run

    sub = run(rs.load(run={"duration_s": 0.3}, overlays=(rs.MEASURED_LOSS,)))
    raw = run(rs.load(run={"duration_s": 0.3}, overlays=(rs.MEASURED_LOSS,), toggles={"hd_det.subtract_dark": False},
                      esa=True))
    prov = sub.offline.source_outputs["opo"].provenance
    assert prov["loss_ic"] == 0.0038
    # total loss budget 1 − ρζξ²η = 0.0645 [Ref B §4] (VT-41 Ref B part: escape efficiency in the source, the rest
    # in the LossRegistry)
    assert 1 - prov["escape_efficiency"] * sub.offline.loss_registry.eta("sqz", "bhd") == pytest.approx(0.0645, abs=5e-5)
    nl_s, nl_r = (r.metrics.noise_levels["i.port0"] for r in (sub, raw))
    assert rs.eq1_db(nl_s["v_minus_db"][rs.I_1MHZ], nl_s["v_plus_db"][rs.I_1MHZ], TH15) == pytest.approx((-9.44, 14.95), abs=0.02)
    assert rs.eq1_db(nl_r["v_minus_db"][rs.I_1MHZ], nl_r["v_plus_db"][rs.I_1MHZ], TH15) == pytest.approx((-9.21, 14.92), abs=0.02)
    observed = (-9.01, 15.12)   # [Ref B §4]
    assert (nl_r["r_minus_db"][rs.I_1MHZ], nl_r["r_plus_db"][rs.I_1MHZ]) == pytest.approx(observed, abs=0.25)
    esa = raw.measurements["esa"]
    assert (esa["squeezed"]["mean_db"], esa["anti_squeezed"]["mean_db"]) == pytest.approx(observed, abs=0.25)
