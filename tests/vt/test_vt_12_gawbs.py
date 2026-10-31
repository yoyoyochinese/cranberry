"""VT-12 ``gawbs_table`` (PLAN §6.2; prompt §7-6 "GAWBS peak frequencies consistent with cladding diameter";
T07 §8 V7, V7b, §4.1).

* V7: R0,m peaks for d = 125 µm, V_L = 5996, V_S = 3740 m/s: 30.5, 82.1, 130.7, 179.0, 227.2, 275.3,
  323.3, 371.3, 419.4, 467.4, 515.4, 563.4 MHz (±0.1 MHz); the 7th within 1 % of the measured
  321.75 MHz (London et al. 2018); d = 80 µm gives exactly ×125/80 = ×1.5625 (1e−12).
* V7b: σ_G² ∝ L (1e−12); ∝ T up to the quantum correction (n_th + ½ = k_BT/ħΩ·(1 + O((ħΩ/k_BT)²)),
  ~2e−10 at 295 K, so 1e−9 here, not T07's 1e−12); T → 0: the thermal part → 0 while the ½ remains
  at the fraction ħΩ/(2k_BT) ≈ 5e−5 of the 295 K value.
* Magnitude: R0,7 (2 /(W km), 4.3 MHz) gives σ² ≈ 1.33e−10 rad²/m and S_φ ≈ 1.96e−17 rad²/(Hz m) at
  the peak, 295 K (T07 §4.1), 2 %.
"""

import numpy as np
import pytest

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.quantum_testing import GAWBS_SMF, LAMBDA0
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline.tables import TableCache

PEAKS_125 = [30.5, 82.1, 130.7, 179.0, 227.2, 275.3, 323.3, 371.3, 419.4, 467.4, 515.4, 563.4]
KW = dict(lambda0_m=LAMBDA0)


def test_v7_peak_frequencies_and_cladding_scaling(tmp_path):
    t = make_block(BlockSpec("gawbs_table", "g", dict(gawbs=dict(GAWBS_SMF), temperature_k=295.0, L_m=1.0,
                                                     f_hz=[1e6], **KW))).tables(TableCache(tmp_path))
    f = t["gawbs_peaks"].values[:, 0] / 1e6
    np.testing.assert_allclose(f, PEAKS_125, atol=0.1)
    assert f[6] == pytest.approx(321.75, rel=0.01)
    f80 = sc.gawbs_peak_frequencies(80e-6, 5996.0, 3740.0, 12) / 1e6
    np.testing.assert_allclose(f80 / f, 125 / 80, rtol=1e-12)


def test_v7b_variance_scaling_and_zero_temperature():
    var = lambda L, T, th=False: sc.gawbs_variances(GAWBS_SMF, temperature_k=T, length_m=L, thermal_only=th, **KW)[6]  # noqa: E731
    assert var(126.0, 295.0) == pytest.approx(2 * var(63.0, 295.0), rel=1e-12)
    assert var(63.0, 590.0) == pytest.approx(2 * var(63.0, 295.0), rel=1e-9)
    assert var(63.0, 0.0, th=True) == 0.0 and var(63.0, 1e-3, th=True) < 1e-6 * var(63.0, 295.0, th=True)
    half = var(63.0, 0.0) / var(63.0, 295.0)
    f7 = sc.gawbs_peak_frequencies(125e-6, 5996.0, 3740.0, 7)[6]
    assert half == pytest.approx(sc.H_PLANCK * f7 / (2 * sc.K_B * 295.0), rel=1e-3)


def test_r07_magnitude_and_psd_consistency():
    s2 = sc.gawbs_variances(GAWBS_SMF, temperature_k=295.0, length_m=1.0, **KW)
    assert s2[6] == pytest.approx(1.33e-10, rel=0.02) and np.count_nonzero(s2) == 1
    f7 = sc.gawbs_peak_frequencies(125e-6, 5996.0, 3740.0, 7)[6]
    psd = lambda f: sc.gawbs_phase_psd(f, GAWBS_SMF, temperature_k=295.0, length_m=1.0, **KW)  # noqa: E731
    assert psd(f7) == pytest.approx(1.96e-17, rel=0.02)
    f = np.linspace(f7 - 43e6, f7 + 43e6, 200001)  # ±10 linewidths: Lorentz area fraction (2/π)·atan(20)
    assert np.trapezoid(psd(f), f) == pytest.approx(s2[6] * 2 / np.pi * np.arctan(20), rel=1e-3)
    assert sc.gawbs_correlation(0.0, GAWBS_SMF, temperature_k=295.0, **KW) == pytest.approx(s2[6], rel=1e-12)
