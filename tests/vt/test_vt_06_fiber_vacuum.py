"""VT-06, linearized-fiber part (PLAN §6.2; prompt §7-4; T07 §8 V1, V1b).

* V1: no nonlinearity, no loss (dispersion β2, β3 only): V(0) = 𝟙/2 → V(L) = 𝟙/2, max|ΔV| < 1e−12.
* V1b: loss with its vacuum noise (toggle ``loss_noise``): same tolerance.
* Squeezed-vacuum path: Kerr on but no bright field (A = 0) — the fluctuation map stays linear.
* CW sideband path (``cw_shortcut``): dispersion + loss, vacuum at every f_s.
Other VT-06 parts (OPA g → 0, OPO x = 0, pulsed OPO) belong to B1 / B6 / B11.
"""

import numpy as np
import pytest

from cranberry.blocks.fib.quantum import fiber_quantum_linear
from cranberry.blocks.fib.quantum_testing import linear_params, toggles

VAC = 0.5 * np.eye(128)
DARK = dict(shape="sech", peak_power_w=0.0, t0_s=1e-12)


@pytest.mark.parametrize("loss", [False, True])
def test_v1_vacuum_stays_vacuum(loss):
    tg = toggles(beta2=True, beta3=True, loss=loss, loss_noise=loss)
    r = fiber_quantum_linear(linear_params(alpha_per_m=0.02, beta3_s3_per_m=1.2e-40), tg)
    assert np.max(np.abs(r.cov - VAC)) < 1e-12


def test_v1_squeezed_vacuum_path_with_kerr_and_no_bright_field():
    tg = toggles(True, cw_shortcut=False, xpm_external=False)  # every term and noise source on
    r = fiber_quantum_linear(linear_params(alpha_per_m=0.02, pulse=DARK, probe=dict(linear_params()["probe"],
                                                                                    p_peak_w=0.0)), tg)
    assert np.max(np.abs(r.cov - VAC)) < 1e-12


def test_v1_cw_sideband_path():
    p = linear_params(laser_mode="cw", power_w=0.0, f_s_hz=[1e6, 1e9, 1e11], alpha_per_m=0.01,
                      steps=dict(n_steps=7, n_save=1, control="fixed"))
    r = fiber_quantum_linear(p, toggles(beta2=True, loss=True, loss_noise=True, kerr=True, kerr_fluct=True,
                                        cw_shortcut=True))
    assert np.max(np.abs(r.cov - 0.5 * np.eye(2))) < 1e-12
