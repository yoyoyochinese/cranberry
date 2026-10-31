"""Unit tests of WP-B10 ``fiber_quantum_linearized`` / ``fiber_quantum_tw``: contract (D34), toggles
(all off ⇒ identity; each noise toggle switches only its own contribution), outputs and input embedding."""

import numpy as np
import pytest

from cranberry.blocks.fib.quantum import NOISE_SOURCES, TW_TOGGLES, fiber_quantum_linear, project_cov
from cranberry.blocks.fib.quantum_testing import linear_params, toggles
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.registry import make_block
from cranberry.offline import gaussian
from cranberry.offline.source_output import SourceOutput

SQ = dict(kind="squeezed", r=1.15, theta_sq_rad=np.pi / 2, mode=dict(shape="sech", t0_s=1e-12))


def _blk(tg, **over):
    return make_block(BlockSpec("fiber_quantum_linearized", "fq", dict(linear_params(**over), toggles=tg)))


def test_contract_missing_and_invalid():
    tg = toggles()
    with pytest.raises(MissingParameterError):
        _blk({k: v for k, v in tg.items() if k != "gawbs_noise"})
    p = linear_params()
    del p["gawbs"]
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("fiber_quantum_linearized", "fq", dict(p, toggles=tg)))
    with pytest.raises(ValueError, match="fixed"):
        _blk(tg, steps=dict(n_steps=10, n_save=3, control="fixed"))
    with pytest.raises(ValueError, match="polarization"):
        _blk(tg, probe=dict(linear_params()["probe"], polarization="circular"))
    with pytest.raises(MissingParameterError):  # cw shortcut needs the sideband grid
        _blk(toggles(cw_shortcut=True), laser_mode="cw", power_w=1e-3)


def test_all_toggles_off_is_identity_temporal_and_cw():
    b = _blk(toggles(), input_state=SQ)
    r = b.result
    t0 = make_block(BlockSpec("fiber_quantum_linearized", "fq", dict(linear_params(input_state=SQ),
                                                                       toggles=toggles(),
                                                                       steps=dict(n_steps=4, n_save=1,
                                                                                  control="fixed")))).result
    np.testing.assert_allclose(r.transfer, np.eye(128), atol=1e-13)
    np.testing.assert_allclose(r.cov, t0.cov, atol=1e-13)
    assert all(not np.any(v) for v in r.noise_cov.values())
    u = np.cosh(r.t_s / 1e-12) ** -1
    vm, vp, th = b.source_output.project(u)
    assert 10 * np.log10(vm / 0.5) == pytest.approx(-20 * 1.15 / np.log(10), abs=1e-9)
    assert th == pytest.approx(0.0, abs=1e-9)
    c = _blk(toggles(cw_shortcut=True), laser_mode="cw", power_w=1e-3, f_s_hz=[1e6, 1e9],
             input_state=SQ).result
    np.testing.assert_allclose(c.cov, np.broadcast_to(gaussian.squeezed_vacuum(1.15, np.pi / 2), (2, 2, 2)),
                               atol=1e-14)


@pytest.mark.parametrize("cw", [False, True])
def test_each_noise_toggle_only_its_own_source(cw):
    full = toggles(True, cw_shortcut=cw, xpm_external=False)
    over = dict(alpha_per_m=0.01)
    if cw:
        over.update(laser_mode="cw", power_w=0.05, f_s_hz=[1e6, 3.2e8, 1e10],
                    steps=dict(n_steps=20, n_save=1, control="fixed"))
    names = dict(loss="loss_noise", raman="raman_noise", gawbs="gawbs_noise", brillouin="brillouin_backscatter",
                 rayleigh="rayleigh_backscatter", raman_bs="raman_backscatter")
    r_all = _blk(full, **over).result
    for src in NOISE_SOURCES:
        assert np.max(np.abs(r_all.noise_cov[src])) > 0, src
        off = _blk(dict(full, **{names[src]: False}), **over).result
        assert not np.any(off.noise_cov[src]), src
        np.testing.assert_allclose(off.cov + r_all.noise_cov[src], r_all.cov, atol=1e-12, rtol=1e-9)


def test_kerr_fluct_off_gives_linear_map_and_source_output_projection():
    r = _blk(toggles(beta2=True, kerr=True))
    m = r.result.transfer
    np.testing.assert_allclose(m @ m.T, np.eye(128), atol=1e-12)  # dispersion only: orthogonal
    on = _blk(toggles(beta2=True, kerr=True, kerr_fluct=True)).result.transfer
    assert np.max(np.abs(on @ on.T - np.eye(128))) > 1e-2  # Kerr squeezes the fluctuations
    b = _blk(toggles(beta2=True, kerr=True, kerr_fluct=True))
    u = np.exp(-0.5 * (b.result.t_s / 1e-12) ** 2)
    so = b.source_output
    assert so.basis == "temporal" and so.cov.shape == (128, 128)
    np.testing.assert_allclose(so.lo_mode_cov(u), b.result.lo_variances(u), atol=1e-13)


def test_source_output_input_embedding():
    b = _blk(toggles(), input_state=SQ)
    so = b.source_output
    r = fiber_quantum_linear(linear_params(), toggles(beta2=True), v0=so)
    r_ref = _blk(toggles(beta2=True), input_state=SQ).result
    np.testing.assert_allclose(r.cov, r_ref.cov, atol=1e-12)
    single = SourceOutput.single_mode_reference(loss=0.0, alpha_per_w=1.0, pump_w=1.0, phi_opa0_rad=0.0)
    with pytest.raises(ValueError, match="temporal"):
        fiber_quantum_linear(linear_params(), toggles(), v0=single)
    assert project_cov(r.noise_cov["loss"], np.ones(64), r.dt).shape == (2, 2)


def test_tw_block_small_run():
    p = linear_params(laser_mode="cw", power_w=1.0, L_m=10.0, grid=dict(n_t=4, t_window_s=1e-12),
                      steps=dict(n_steps=10, n_save=1, control="fixed"), n_traj=200, seed=3)
    tg = dict.fromkeys(TW_TOGGLES, False)
    b = make_block(BlockSpec("fiber_quantum_tw", "tw", dict(p, toggles=tg)))
    assert b.cov.shape == (8, 8)
    np.testing.assert_allclose(np.diag(b.cov), 0.5, rtol=0.35)  # vacuum, 200 samples
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("fiber_quantum_tw", "tw",
                             dict(p, toggles={k: v for k, v in tg.items() if k != "loss_noise"})))
