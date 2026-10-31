"""core.units: dB conversions (conv §2, §7) and the unit-suffix convention (G16-i)."""

import numpy as np
import pytest

from cranberry.core import units as u


def test_power_and_amplitude_db_round_trip():
    x = np.array([1e-3, 0.5, 1.0, 2.0, 1e4])
    np.testing.assert_allclose(u.db_to_lin(u.lin_to_db(x)), x, rtol=1e-14)
    np.testing.assert_allclose(u.db_to_amp(u.amp_to_db(x)), x, rtol=1e-14)
    assert u.db_to_lin(10.0) == pytest.approx(10.0)
    assert u.db_to_amp(20.0) == pytest.approx(10.0)
    assert u.lin_to_db(0.0) == -np.inf
    with pytest.raises(ValueError):
        u.lin_to_db(-1.0)


def test_insertion_loss():
    assert u.il_db_to_eta(0.0) == pytest.approx(1.0)
    assert u.il_db_to_eta(10.0) == pytest.approx(0.1)
    assert u.il_db_to_eta(3.0) == pytest.approx(0.501187, rel=1e-6)
    assert u.eta_to_il_db(0.1) == pytest.approx(10.0)
    with pytest.raises(ValueError):
        u.eta_to_il_db(0.0)


def test_noise_level_relative_to_shot_noise():
    """S_dB = 10 log10(V / V_SN), V_SN = V_vac = 1/2 (conv §2)."""
    assert u.V_VAC == 0.5 and u.V_SN == 0.5
    assert u.noise_db(0.5) == pytest.approx(0.0)
    assert u.noise_db(0.05) == pytest.approx(-10.0)
    assert u.noise_db_to_variance(-10.0) == pytest.approx(0.05)
    assert u.noise_db_to_variance(u.noise_db(0.37)) == pytest.approx(0.37)


def test_phase_noise_and_rin():
    """ℒ = S_φ/2 (one-sided S_φ, conv §7)."""
    assert u.dbc_per_hz_to_s_phi(-100.0) == pytest.approx(2e-10)
    assert u.s_phi_to_dbc_per_hz(2e-10) == pytest.approx(-100.0)
    assert u.db_per_hz_to_lin(-140.0) == pytest.approx(1e-14)
    assert u.lin_to_db_per_hz(1e-15) == pytest.approx(-150.0)


@pytest.mark.parametrize(
    "name, unit",
    [
        ("f_sim_hz", "Hz"),
        ("duration_s", "s"),
        ("tau_s", "s"),
        ("L_m", "m"),
        ("il_db", "dB"),
        ("gain_v_per_a", "V/A"),
        ("k_s_rad_per_v", "rad/V"),
        ("ki_per_s", "1/s"),
        ("kt_rad_per_m_k", "rad/(m·K)"),
        ("rin_db_per_hz", "dB/Hz"),
        ("l_dbc_per_hz", "dBc/Hz"),
        ("h0_hz2_per_hz", "Hz²/Hz"),
        ("power_w", "W"),
        ("dL", None),
        ("s", None),
    ],
)
def test_unit_of(name, unit):
    assert u.unit_of(name) == unit


def test_is_db_name():
    assert u.is_db_name("il_db") and u.is_db_name("rin_db_per_hz") and u.is_db_name("l_dbc_per_hz")
    assert not u.is_db_name("gain_v_per_a") and not u.is_db_name("_db")
