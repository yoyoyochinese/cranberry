"""Unit tests of the WP-B9 fiber blocks: contract (D34), toggles (off ⇒ identity / term removed),
tangent-linear step (B10 interface), step-halving control, offline tables and their cache."""

import numpy as np
import pytest

from cranberry.blocks.fib.gnlse import TOGGLES, GnlseModel, propagate, pulse_field, time_axis
from cranberry.blocks.fib.tables import LO_COLUMNS
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.registry import make_block
from cranberry.offline.tables import TableCache

B2, GAMMA, T0 = -21.7e-27, 1.3e-3, 100e-15
LD = T0**2 / abs(B2)
P1 = abs(B2) / (GAMMA * T0**2)
FIBER = dict(L_m=2 * LD, alpha_per_m=0.05, beta2_s2_per_m=B2, beta3_s3_per_m=0.126e-39, beta_higher=[1e-55],
             gamma_per_w_m=GAMMA, lambda0_m=1550e-9, raman_model="blow_wood", f_R=0.18, xpm_b=2.0, xpm_power_w=50.0)
ON = dict.fromkeys(TOGGLES, True)


def pulsed(toggles, **over):
    p = dict(FIBER, toggles=toggles, laser_mode="pulsed", pulse={"shape": "sech", "peak_power_w": 2 * P1, "t0_s": T0},
             grid={"n_t": 1024, "t_window_s": 40 * T0}, steps={"n_steps": 200, "n_save": 4, "control": "fixed"})
    p.update(over)
    return make_block(BlockSpec("fiber_gnlse", "fib", p))


# --- contract ---------------------------------------------------------------------------------------

def test_missing_toggle_and_params_raise():
    with pytest.raises(MissingParameterError):
        pulsed({k: v for k, v in ON.items() if k != "raman"})
    p = dict(FIBER, toggles=ON, laser_mode="pulsed")
    with pytest.raises(MissingParameterError):  # pulse/grid/steps missing
        make_block(BlockSpec("fiber_gnlse", "fib", p))
    q = {k: v for k, v in FIBER.items() if k != "f_R"}
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("fiber_gnlse", "fib", dict(q, toggles=ON, laser_mode="cw", power_w=1.0)))
    with pytest.raises(ValueError, match="raman_model"):
        pulsed(ON, raman_model="hollenbeck")
    # CW with the shortcut needs no grid
    b = make_block(BlockSpec("fiber_gnlse", "fib", dict(FIBER, toggles=ON, laser_mode="cw", power_w=1e-3)))
    assert b.reads == () and b.writes == () and b.toggles == ON


# --- toggles ----------------------------------------------------------------------------------------

def test_all_toggles_off_is_identity():
    r = pulsed(dict.fromkeys(TOGGLES, False)).result
    np.testing.assert_allclose(r.a_out, r.a_in, rtol=0, atol=1e-13 * np.abs(r.a_in).max())
    assert np.all(r.phi_nl_peak_rad == 0.0)


def test_toggle_equals_parameter_zero():
    def out(tg, **over):
        return pulsed(tg, **over).result.a_out

    np.testing.assert_allclose(out(dict(ON, loss=False)), out(ON, alpha_per_m=0.0), rtol=0, atol=1e-12)
    np.testing.assert_allclose(out(dict(ON, raman=False)), out(ON, f_R=0.0), rtol=0, atol=1e-12)  # D49
    np.testing.assert_allclose(out(dict(ON, xpm_external=False)), out(ON, xpm_power_w=0.0), rtol=0, atol=1e-12)
    np.testing.assert_allclose(out(dict(ON, beta_higher=False)), out(ON, beta_higher=[0.0]), rtol=0, atol=1e-12)
    # kerr and raman off: purely linear propagation e^{L D}
    blk = pulsed(dict(ON, kerr=False, raman=False, xpm_external=False))
    r = blk.result
    m = r.model
    lin = np.roll(m.to_t(m.linear(m.to_w(np.roll(r.a_in, -512)), FIBER["L_m"])), 512)
    np.testing.assert_allclose(r.a_out, lin, rtol=0, atol=1e-12 * np.abs(lin).max())
    # every term changes the answer when switched on (each toggle is live)
    base = out(ON)
    for k in ("loss", "beta2", "beta3", "kerr", "raman", "self_steepening", "xpm_external"):
        assert np.max(np.abs(out(dict(ON, **{k: False})) - base)) > 1e-6 * np.abs(base).max(), k


# --- tangent-linear step (WP-B10 interface) --------------------------------------------------------

def test_tangent_matches_finite_difference():
    t, dt = time_axis(512, 40 * T0)
    m = GnlseModel(512, dt, lambda0_m=1550e-9, alpha_per_m=0.05, betas=(B2, 1e-40), gamma_per_w_m=GAMMA,
                   c_kerr=0.82, c_raman=0.18, raman_model="lin_agrawal", self_steepening=True, xpm_rate_per_m=0.1,
                   length_m=LD)
    at0 = m.to_w(np.roll(pulse_field("sech", 2 * P1, T0, t), -256).astype(complex))
    rng = np.random.default_rng(1)
    band = np.abs(at0) > 1e-6 * np.abs(at0).max()
    da = (rng.normal(size=(2, 512)) + 1j * rng.normal(size=(2, 512))) * band * 1e-2 * np.abs(at0).max()
    pr = propagate(m, at0, LD, n_steps=40, n_save=2, dat0=da)
    assert pr.dat_z.shape == (3, 2, 512)
    eps = 1e-4
    for j in range(2):
        fp = propagate(m, at0 + eps * da[j], LD, n_steps=40, n_save=1).at_z[-1]
        fm = propagate(m, at0 - eps * da[j], LD, n_steps=40, n_save=1).at_z[-1]
        fd = (fp - fm) / (2 * eps)
        assert np.linalg.norm(fd - pr.dat_z[-1, j]) <= 1e-6 * np.linalg.norm(fd)


# --- step-halving control --------------------------------------------------------------------------

def test_halving_control_matches_fine_fixed():
    t, dt = time_axis(1024, 40 * T0)
    m = GnlseModel(1024, dt, lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1.0,
                   c_raman=0.0, raman_model=None, self_steepening=False)
    at0 = m.to_w(np.roll(pulse_field("sech", 4 * P1, T0, t), -512).astype(complex))
    z0 = np.pi / 2 * LD
    ref = propagate(m, at0, z0, n_steps=3200, n_save=1).at_z[-1]
    ad = propagate(m, at0, z0, n_steps=10, n_save=3, control="halving", rtol=1e-7)
    assert np.linalg.norm(ad.at_z[-1] - ref) / np.linalg.norm(ref) < 1e-5
    assert ad.n_rejected >= 1 and ad.n_steps < 3200 and ad.at_z.shape == (4, 1024)
    with pytest.raises(ValueError):
        propagate(m, at0, z0, n_steps=10, n_save=3)  # fixed: n_steps must be a multiple of n_save


# --- tables -----------------------------------------------------------------------------------------

def test_fiber_linear_table_and_cache(tmp_path):
    p = dict(toggles={"loss": True, "beta2": True, "beta3": False}, lengths_m=[0.0, 10.0, 1000.0],
             alpha_per_m=1e-4, n_g=1.4682, beta2_s2_per_m=B2, beta3_s3_per_m=1e-40)
    cache = TableCache(tmp_path)
    t = make_block(BlockSpec("fiber_linear_table", "lin", p)).tables(cache)["fiber_linear"]
    np.testing.assert_allclose(t.values[:, 0], np.exp(-1e-4 * np.array([0, 10, 1000.0])), rtol=1e-15)
    assert t.values[2, 1] == pytest.approx(1.4682 * 1000 / 299792458.0)
    assert t.values[2, 2] == pytest.approx(B2 * 1000) and np.all(t.values[:, 3] == 0.0)
    assert len(list(tmp_path.glob("*.npz"))) == 1
    make_block(BlockSpec("fiber_linear_table", "lin", p)).tables(cache)  # cache hit: no new entry
    assert len(list(tmp_path.glob("*.npz"))) == 1
    make_block(BlockSpec("fiber_linear_table", "lin", dict(p, alpha_per_m=2e-4))).tables(cache)
    assert len(list(tmp_path.glob("*.npz"))) == 2


def test_lo_nonlinear_table(tmp_path):
    f_rep, tau = 100e6, 1e-12
    p = dict(FIBER, L_m=20.0, alpha_per_m=4e-5, xpm_power_w=0.0, toggles={k: v for k, v in ON.items() if k != "cw_shortcut"},
             pulse={"shape": "sech", "t_fwhm_s": tau}, f_rep_hz=f_rep, p_avg_w=[0.99e-3, 1e-3, 1.01e-3],
             grid={"n_t": 1024, "t_window_s": 100e-12}, steps={"n_steps": 100, "n_save": 4, "control": "fixed"})
    tabs = make_block(BlockSpec("lo_nonlinear_table", "lo", p)).tables(TableCache(tmp_path))
    lo, sp = tabs["lo_nonlinear"], tabs["gnlse_output_spectrum"]
    assert lo.axes == ("p_avg_w", "z_m") and lo.values.shape == (3, 5, len(LO_COLUMNS))
    assert sp.values.shape == (3, 5, 1024) and np.all(np.diff(sp.grid["omega_rad_s"]) > 0)
    v = lo.values
    np.testing.assert_allclose(v[:, 0, :2], [[0.0, 1.0]] * 3, atol=1e-12)  # z = 0: ov = 1
    assert np.all(np.diff(v[1, :, 3]) > 0)  # φ_NL,peak accumulates
    # T06 §8.4 order of magnitude (100 MHz, 1 ps, 1 mW): arg ov ~ 0.1 rad at tens of m, |ov|² near 1
    assert 0.01 < v[1, -1, 0] < 0.3 and v[1, -1, 1] > 0.99
    # tangent-linear ∂ arg ov/∂ ln P vs finite difference over ±1 % power
    fd = (v[2, -1, 0] - v[0, -1, 0]) / (np.log(1.01e-3) - np.log(0.99e-3))
    assert v[1, -1, 2] == pytest.approx(fd, rel=1e-3)
    assert lo.provenance["columns"] == list(LO_COLUMNS) and "spec_hash" in lo.provenance
