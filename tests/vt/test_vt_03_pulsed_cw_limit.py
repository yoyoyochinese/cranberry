"""VT-03 pulsed part (WP-B11; prompt §7-1 CW limit, T03 V5/§6, T17 §8.1 ``opo.pulsed``, PLAN G07).

The GNLSE part lives in ``test_vt_03_cw_limit.py`` (WP-B9).

* ``opa_supermodes`` with ``multimode`` off is the single-mode model: after the §7-1 loss and jitter the
  numbers (−10.148 / +20.556 dB) are **bit-identical** to ``single_mode_reference``.
* ``pulsed_opo_source`` with ``pulsed`` off returns the CW ``opo_cw_spectrum`` covariance bit for bit
  (linewidth model, x = x_ss, Δ = Δ_s; also with the detuning on).
* Matched pulsed LO and zero timing / length error leave the source untouched: η_t = 1, θ_off = 0
  (``lo_temporal_shape: optimal``, ``mode_overlap_table`` at (0, 0)).
"""

import numpy as np
import pytest

from cranberry.blocks.source.single_mode import r_prime
from cranberry.blocks.source.pulsed_testing import BOREALIS, TOGGLES_ON, pulsed_opo
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline import gaussian

REF_A = dict(loss=0.08, alpha_per_w=9.06, pump_w=0.64, theta_rms_rad=0.009, phi_opa0_rad=0.0)


def test_opa_supermodes_cw_is_single_mode_bit_identical():
    ref = make_block(BlockSpec("single_mode_reference", "ref", dict(REF_A, average="exact")))
    opa = make_block(BlockSpec("opa_supermodes", "opa", dict(
        length_m=0.045, alpha_per_w=9.06, peak_pump_w=0.64, phi_opa0_rad=0.0, segment="xtal", psi_p_xtal_rad=0.0,
        toggles={"multimode": False, "gvm": True, "gvd": True})))
    so = opa.source_output
    assert so.basis == "single_mode"
    vm, vp, th = so.project(None)
    vm, vp, _ = gaussian.noise_levels_from_cov(gaussian.pure_loss(so.cov, 1 - REF_A["loss"]))  # LossRegistry action
    a, b = r_prime(2 * vm, 2 * vp, REF_A["theta_rms_rad"], "exact")
    assert (a, b) == ref.r_prime()  # bit for bit
    assert 10 * np.log10(a) == pytest.approx(-10.148, abs=1e-3) and 10 * np.log10(b) == pytest.approx(20.556, abs=1e-3)
    assert float(th) == 0.0 and opa.supermode_r[0] == np.sqrt(9.06 * 0.64)


@pytest.mark.parametrize("detuning", [False, True])
def test_pulsed_opo_off_is_cw_opo_bit_identical(detuning):
    fs = [0.0, 1e6, 30e6, 100e6]
    kc = np.pi * BOREALIS["fwhm_hz"]
    b = pulsed_opo(BOREALIS, x_ss=0.6, delta_s_rad_s=0.02 * kc, f_s_hz=fs, toggles={"pulsed": False, "detuning": detuning})
    cw = make_block(BlockSpec("opo_cw_spectrum", "cw", dict(
        cavity_model="linewidth", fwhm_hz=BOREALIS["fwhm_hz"], escape_efficiency=0.95, pump_mode="x", x=0.6,
        detuning_rad_s=0.02 * kc, f_s_hz=fs, phi_opa0_rad=0.0, segment="opo", gain_table_grid=None, coef_table_x=None,
        toggles={"frequency_dependence": True, "pump_dependent_loss": False, "detuning": detuning})))
    so = b.source_output
    assert so.basis == "cw_sideband"
    assert np.array_equal(so.cov, cw.source_output.cov)


def test_matched_pulsed_lo_is_identity():
    src = {"block": "pulsed_opo_source", "name": "opo", "params": dict(BOREALIS, toggles=TOGGLES_ON)}
    lo = make_block(BlockSpec("lo_temporal_shape", "lo", dict(source=src, shape="optimal", fit=False, observer="bhd0",
                                                              source_tag="[T17 §6.3]", toggles={"shape_mismatch": True})))
    assert lo.eta_t == pytest.approx(1.0, abs=1e-12)
    assert lo.projection[2] == pytest.approx(0.0, abs=1e-12)
    ov = make_block(BlockSpec("mode_overlap_table", "ov", dict(source=src, dt_grid_s=[0.0], dL_grid_m=[0.0],
                                                               beta2_s2_per_m=-21.7e-27, beta3_s3_per_m=0.0, pad=2,
                                                               toggles={"dispersion": True, "jitter_fast": False})))
    vm, vp, th, eta = ov.table.values[0, 0]
    assert eta == pytest.approx(1.0, abs=1e-12) and th == pytest.approx(0.0, abs=1e-12)
    assert 2 * vm == pytest.approx(pulsed_opo(BOREALIS).optimal_mode[1], rel=1e-9)
