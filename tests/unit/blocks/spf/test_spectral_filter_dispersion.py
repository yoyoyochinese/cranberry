"""spectral_filter and dispersion_comp (WP-B11): filter action on SourceOutput, registered losses, toggles."""

import numpy as np
import pytest

from cranberry.blocks.fib.mode_overlap import gaussian_source
from cranberry.blocks.source.pulsed_testing import BOREALIS, pulsed_opo
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline import gaussian
from cranberry.offline.loss_registry import LossRegistry


def filt(shape="gauss", bw=10e9, on=True, center=0.0):
    return make_block(BlockSpec("spectral_filter", "sf", dict(
        shape=shape, bandwidth_hz=bw, center_hz=center, insertion_eta=0.98, segment="gbs.in", applies_to=["gbs"],
        source_tag="[T03 §5.4]", pad=2, toggles={"enabled": on})))


def test_wide_filter_leaves_ns_mode_alone_and_narrow_filter_loses():
    b = pulsed_opo(BOREALIS, x_ss=2.4, dt_kappa=1 / 40)
    so, u = b.source_output, b.optimal_mode[0]
    v0 = float(so.project(u)[0])
    wide = filt("gauss", 10e9).apply(so)  # 10 GHz ≫ 66 MHz linewidth
    assert float(wide.project(u)[0]) == pytest.approx(v0, rel=1e-3)
    narrow = filt("lorentz", 66e6).apply(so)
    assert float(narrow.project(u)[0]) > v0 + 0.01
    # the discretised source obeys the uncertainty relation to O((κδt)²) only (V-P2); the passive filter
    # must not make that worse
    nu = lambda c: gaussian.symplectic_eigenvalues(np.asarray(c)).min()  # noqa: E731
    assert nu(narrow.cov) >= nu(so.cov) - 1e-6


def test_single_mode_filter_is_pure_loss_of_its_overlap():
    so = gaussian_source(1e-12, 1.0, 512, 40e-12, 0.0)
    f = filt("gauss", 0.5e12)
    out = f.apply(so)
    # one mode in, one mode out: A = ⟨u|h∗u⟩, so the filtered state is a pure loss |A|² (conv §2)
    g, u = so.grid, so.modes[0]
    from cranberry.core import fft
    dt = g[1] - g[0]
    a = np.zeros(3 * u.size, complex)
    a[u.size : 2 * u.size] = u
    t0 = g[0] - u.size * dt
    hu = fft.ifft_agrawal(fft.fft_agrawal(a, dt, t0=t0) * f.t_of(fft.omega_grid(a.size, dt)), dt, t0=t0)[u.size : 2 * u.size]
    eta = abs(np.sum(u * hu) * dt) ** 2
    np.testing.assert_allclose(out.cov, gaussian.pure_loss(so.cov, eta), atol=1e-12)


def test_cw_sideband_and_toggle_off():
    cw = make_block(BlockSpec("opo_cw_spectrum", "cw", dict(
        cavity_model="linewidth", fwhm_hz=66e6, escape_efficiency=0.95, pump_mode="x", x=0.5, detuning_rad_s=0.0,
        f_s_hz=[0.0, 50e6, 200e6], phi_opa0_rad=0.0, segment="opo", gain_table_grid=None, coef_table_x=None,
        toggles={"frequency_dependence": True, "pump_dependent_loss": False, "detuning": False}))).source_output
    f = filt("rect", 200e6)
    out = f.apply(cw)
    np.testing.assert_array_equal(out.cov[:2], cw.cov[:2])
    np.testing.assert_allclose(out.cov[2], 0.5 * np.eye(2), atol=1e-15)  # 200 MHz outside the ±100 MHz band
    assert filt(on=False).apply(cw) is cw
    with pytest.raises(ValueError):
        filt(center=1e9).apply(cw)


def test_spectral_filter_loss_registration():
    reg = LossRegistry(filt().losses(), paths={"sqz": {"segments": ["gbs.in"]}})
    assert reg.eta("sqz", "gbs") == pytest.approx(0.98) and reg.eta("sqz", "bhd") == 1.0
    reg = LossRegistry(filt(on=False).losses(), paths={"sqz": {"segments": ["gbs.in"]}})
    assert reg.eta("sqz", "gbs") == 1.0


def comp(on=True):
    return make_block(BlockSpec("dispersion_comp", "dc", dict(
        ports={"p1": {"gdd_s2": 1.37e-24, "tod_s3": -8e-39, "delay_s": -2e-13, "segment": "port1.dcf", "insertion_eta": 0.95}},
        applies_to=["gbs", "bhd"], source_tag="[D44]", toggles={"enabled": on})))


def test_dispersion_comp_residual_and_toggle():
    phi2 = -21.7e-27 * 63.0
    r = comp().residual("p1", phi2, 8e-39, 2e-13)
    assert r == pytest.approx((phi2 + 1.37e-24, 0.0, 0.0), abs=1e-40)
    assert comp(False).residual("p1", phi2, 8e-39, 2e-13) == (phi2, 8e-39, 2e-13)
    reg = LossRegistry(comp().losses(), paths={"sqz1": {"segments": ["port1.dcf"]}})
    assert reg.eta("sqz1", "gbs") == pytest.approx(0.95)
    assert LossRegistry(comp(False).losses(), paths={"sqz1": {"segments": ["port1.dcf"]}}).eta("sqz1", "gbs") == 1.0
