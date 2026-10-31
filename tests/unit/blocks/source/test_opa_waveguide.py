"""opa_waveguide block (WP-B1): models, closed forms, toggles-off identity, contract (T01 §2, §9)."""

import numpy as np
import pytest

from cranberry.blocks.source.opa_waveguide import db_per_cm_to_per_m, make_opa_waveguide
from cranberry.core.block import MissingParameterError
from cranberry.offline import gaussian as g
from cranberry.offline.band_ledger import BandLedger
from cranberry.offline.loss_registry import LossRegistry

OFF = dict(pump_loss=False, phase_mismatch=False, pump_depletion=False, pump_rin=False)
BASE = dict(alpha_per_w=9.06, pump_w=0.64, pump_ref="output", length_m=0.045, phi_opa0_rad=0.0,
            segment="xtal.sqz", psi_p_xtal_rad=0.0, gain_table_pump_w=None, toggles=OFF)
LUMPED = dict(BASE, model="lumped", waveguide_loss=0.02, loss_source="[Ref A §4.2]")
DIST = dict(BASE, model="distributed", loss_signal_db_per_cm=0.1)
CW = dict(BASE, model="coupled_wave", loss_signal_db_per_cm=0.1, seed_power_w=0.0)


def opa(params, **kw):
    toggles = dict(params["toggles"], **kw.pop("toggles", {}))
    return make_opa_waveguide(dict(params, **kw, toggles=toggles))


def r_pm(block):
    vm, vp, _ = block.source_output.project(None)
    return 2 * float(vm), 2 * float(vp)


def closed_form(a_s_db, pump_w, ell=0.045, alpha=9.06):
    """T01 §2.3: R± = [±2g e^{(±2g − a_s)ℓ} − a_s] / (±2g − a_s)."""
    a, gg = db_per_cm_to_per_m(a_s_db), np.sqrt(alpha * pump_w) / ell
    return tuple((s * 2 * gg * np.exp((s * 2 * gg - a) * ell) - a) / (s * 2 * gg - a) for s in (-1, 1))


@pytest.mark.parametrize("a_s_db, pump_w, l_wg", [(0.1, 0.1, 0.0353), (0.1, 0.3, 0.0268), (0.1, 0.64, 0.0203),
                                                  (0.1, 1.0, 0.0167), (0.1, 1.2, 0.0153), (0.07, 0.64, 0.0143)])
def test_distributed_closed_form_and_l_wg_table(a_s_db, pump_w, l_wg):
    """V-T01-3: ODE = closed form (1e−6), effective L_wg of T01 §2.3 table (±0.01 %)."""
    b = opa(DIST, loss_signal_db_per_cm=a_s_db, pump_w=pump_w)
    np.testing.assert_allclose(r_pm(b), closed_form(a_s_db, pump_w), rtol=1e-6)
    assert b.source_output.provenance["L_wg"] == pytest.approx(l_wg, abs=1e-4)
    assert b.losses() == []  # distributed loss is inside cov


def test_limits_and_vacuum():
    """a_s → 0 gives e^{±2r}; P = 0 gives vacuum (VT-06 OPA part, prompt §7-4)."""
    r = np.sqrt(9.06 * 0.64)
    np.testing.assert_allclose(r_pm(opa(DIST, loss_signal_db_per_cm=0.0)), [np.exp(-2 * r), np.exp(2 * r)], rtol=1e-9)
    for p in (DIST, LUMPED, CW):
        np.testing.assert_allclose(opa(p, pump_w=0.0).source_output.cov, g.vacuum(1), atol=1e-14)


def test_lumped_registers_loss_and_callable():
    b = opa(LUMPED)
    r = np.sqrt(9.06 * 0.64)
    np.testing.assert_allclose(b.source_output.cov, g.squeezed_vacuum(r, np.pi / 2), rtol=1e-9, atol=1e-12)
    (el,) = b.losses()
    assert el.eta == pytest.approx(0.98) and el.segment == "xtal.sqz" and el.kind == "component"
    reg = LossRegistry(b.losses(), paths={"sqz": {"segments": ["xtal.sqz"]}})
    assert reg.eta("sqz", "gbs") == pytest.approx(0.98)
    b2 = opa(LUMPED, waveguide_loss=lambda P_w: 0.01 + 0.01 * P_w)  # pump-induced loss L_wg(P) (T01 §2.8)
    assert b2.losses()[0].eta == pytest.approx(1 - 0.0164)
    with pytest.raises(ValueError):
        opa(LUMPED, waveguide_loss=1.5).losses()


def test_coupled_wave_equals_distributed_without_depletion():
    np.testing.assert_allclose(opa(CW, seed_power_w=1e-3).source_output.cov, opa(DIST).source_output.cov, rtol=1e-10)


def test_toggles_off_are_identity():
    """Toggle off ⇒ same result as the effect's parameter absent / zero (PLAN §7.3)."""
    ref = opa(LUMPED).source_output.cov
    for extra in (dict(delta_k_per_m=50.0), dict(loss_pump_db_per_cm=0.3)):
        np.testing.assert_array_equal(opa(LUMPED, **extra).source_output.cov, ref)
    assert opa(LUMPED, pump_rin=dict(f_lo_hz=1.0, f_hi_hz=1e5, variance_outside=1e-6)).bands() == []
    a = opa(CW, seed_power_w=0.1, gain_table_pump_w=[0.1, 0.64]).source_output.gain_table.values
    b = opa(CW, seed_power_w=0.0, gain_table_pump_w=[0.1, 0.64]).source_output.gain_table.values
    np.testing.assert_allclose(a, b, rtol=1e-12)  # no depletion: seed power irrelevant
    # and on ⇒ the effect acts
    assert not np.allclose(opa(LUMPED, delta_k_per_m=50.0, toggles=dict(phase_mismatch=True)).source_output.cov, ref)
    c = opa(CW, seed_power_w=0.1, gain_table_pump_w=[0.1, 0.64], toggles=dict(pump_depletion=True))
    assert c.source_output.gain_table.values[1, 0] < 0.9 * b[1, 0]


def test_pump_loss_reference_faces():
    """T01 §2.4: r = √(αP_in)·(2/(a_pℓ))(1 − e^{−a_pℓ/2}); output face P_in = P_out e^{a_pℓ}."""
    a_p, ell = db_per_cm_to_per_m(0.2), 0.045
    on = dict(toggles=dict(pump_loss=True), loss_pump_db_per_cm=0.2)
    r_in = opa(LUMPED, pump_ref="input", **on).source_output.provenance["r_eff"]
    assert r_in == pytest.approx(np.sqrt(9.06 * 0.64) * 2 / (a_p * ell) * (1 - np.exp(-a_p * ell / 2)), rel=1e-9)
    prov = opa(LUMPED, pump_ref="output", **on).source_output.provenance
    assert prov["pump_out_w"] == pytest.approx(0.64, rel=1e-9)
    assert prov["r_eff"] == pytest.approx(np.sqrt(9.06 * 0.64) * 2 / (a_p * ell) * (np.exp(a_p * ell / 2) - 1), rel=1e-9)


def test_gain_table_lumped():
    grid = [0.0, 0.1, 0.64, 1.2]
    t = opa(LUMPED, gain_table_pump_w=grid).source_output.gain_table
    r = np.sqrt(9.06 * np.array(grid))
    np.testing.assert_allclose(t.values, np.stack([np.exp(r), np.exp(-r)], axis=1), rtol=1e-9)
    assert t.interp(pump_w=0.64)[0] == pytest.approx(np.exp(r[2]), rel=1e-9)


def test_pump_rin_band_declared_not_added():
    b = opa(LUMPED, toggles=dict(pump_rin=True), pump_rin=dict(f_lo_hz=1.0, f_hi_hz=1e5, psd_outside=1e-14))
    (d,) = b.bands()
    assert d.quantity == "rin" and d.targets == ("xtal.sqz",) and d.source == b.name
    BandLedger().declare(d)
    assert b.source_output.provenance["dr_d_rin"] == pytest.approx(np.sqrt(9.06 * 0.64) / 2)
    with pytest.raises(MissingParameterError):
        opa(LUMPED, toggles=dict(pump_rin=True))


def test_crystal_is_pump_path_segment():
    """ψ_p,xtal is reported, not baked into cov (enters θ_sq as ½ψ via the time loop, conv §4)."""
    a, b = opa(LUMPED, psi_p_xtal_rad=0.0).source_output, opa(LUMPED, psi_p_xtal_rad=1.3).source_output
    np.testing.assert_array_equal(a.cov, b.cov)
    assert b.pump_path_phase_rad == {"xtal.sqz": 1.3}


def test_validation():
    with pytest.raises(ValueError):
        opa(LUMPED, model="magic")
    with pytest.raises(ValueError):
        opa(DIST, toggles=dict(pump_depletion=True))
    with pytest.raises(MissingParameterError):
        opa(DIST, toggles=dict(phase_mismatch=True))
    with pytest.raises(MissingParameterError):
        make_opa_waveguide({k: v for k, v in DIST.items() if k != "loss_signal_db_per_cm"})
    with pytest.raises(MissingParameterError):
        make_opa_waveguide({k: v for k, v in LUMPED.items() if k != "toggles"})
    with pytest.raises(ValueError):
        opa(LUMPED, gain_table_pump_w=[0.5, 0.1])


@pytest.mark.parametrize("params", [LUMPED, DIST, CW])
def test_physical_states(params):
    for kw in (dict(), dict(delta_k_per_m=80.0, toggles=dict(phase_mismatch=True))):
        assert g.is_physical(opa(params, **kw).source_output.cov)
