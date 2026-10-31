"""VT-11 ``raman_noise_table`` (PLAN §6.2; prompt §7-6 read by D48: only the thermal part → 0 as T → 0;
T07 §8 V6, V6b, V6c, §3.4). Blow–Wood, γ = 1.3 /(W km), f_R = 0.18.

* V6: Ω/2π = 100 GHz, PL = 1 W·m. Thermal part g_R n_th at 1 mK < 1e−6 × its 295 K value; at 1 mK
  anti-Stokes = 0 and Stokes = g_R(Ω) = 2.39e−6 /(W m) (295 K Stokes 1.48e−4), 1 %.
  1 ps Gaussian mode, 63 m, 1 mW, same-frequency probe: N_R(1 mK)/N_R(295 K) ≈ 0.012 (T07 §3.4:
  9.3e−6 → 1.1e−7), 5 %.
* V6b detailed balance: anti-Stokes/Stokes = exp(−ħΩ/k_BT) to 1e−6 relative (0.117 at 13.2 THz, 295 K).
* V6c g_R peak 0.68 /(W km) at 13.1 THz (±2 %).
* Linear layer: the Raman noise of a bright pulse drops with T but does not vanish (the ½ term).
"""

import numpy as np
import pytest

from cranberry.blocks.fib import scattering as sc
from cranberry.blocks.fib.quantum import fiber_quantum_linear
from cranberry.blocks.fib.quantum_testing import GAMMA, LAMBDA0, linear_params, toggles
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline.tables import TableCache

KW = dict(gamma_per_w_m=GAMMA, f_R=0.18, raman_model="blow_wood")
W100 = 2 * np.pi * 100e9


def test_v6_thermal_part_vanishes_spontaneous_stokes_remains():
    th = lambda T, w=W100: float(sc.raman_scatter_density(w, temperature_k=T, thermal_only=True, **KW))  # noqa: E731
    assert th(1e-3) < 1e-6 * th(295.0)
    assert float(sc.raman_scatter_density(-W100, temperature_k=1e-3, **KW)) == 0.0
    assert float(sc.raman_scatter_density(W100, temperature_k=1e-3, **KW)) == pytest.approx(2.39e-6, rel=0.01)
    assert float(sc.raman_scatter_density(W100, temperature_k=295.0, **KW)) == pytest.approx(1.48e-4, rel=0.01)
    assert float(sc.raman_scatter_density(W100, temperature_k=0.0, **KW)) == pytest.approx(
        float(sc.raman_gain(W100, **KW)), rel=1e-15)


def test_v6_mode_ratio_1ps():
    mode = dict(shape="gaussian", t_fwhm_s=1e-12)

    def n_mode(T):  # scattered offset f from the common carrier: Ω = −2πf
        return sc.mode_average(lambda f: 1e-3 * 63 * sc.raman_scatter_density(-2 * np.pi * f, temperature_k=T, **KW),
                               mode)

    hot, cold = n_mode(295.0), n_mode(1e-3)
    assert hot == pytest.approx(9.3e-6, rel=0.05) and cold == pytest.approx(1.1e-7, rel=0.05)
    assert cold / hot == pytest.approx(0.012, rel=0.05)


def test_v6b_detailed_balance():
    w = 2 * np.pi * 13.2e12
    r = sc.raman_scatter_density(-w, temperature_k=295.0, **KW) / sc.raman_scatter_density(w, temperature_k=295.0, **KW)
    assert r == pytest.approx(np.exp(-sc.HBAR * w / (sc.K_B * 295.0)), rel=1e-6)
    assert r == pytest.approx(0.117, abs=5e-4)


def test_v6c_gain_peak_and_table(tmp_path):
    f = np.linspace(1e12, 30e12, 29001)
    g = sc.raman_gain(2 * np.pi * f, **KW)
    i = int(np.argmax(g))
    assert g[i] * 1e3 == pytest.approx(0.68, rel=0.02) and f[i] == pytest.approx(13.1e12, rel=0.02)
    t = make_block(BlockSpec("raman_noise_table", "rt", dict(KW, lambda0_m=LAMBDA0, temperature_k=1e-3,
                                                            f_hz=[-100e9, 100e9]))).tables(TableCache(tmp_path))
    v = t["raman_noise"].values
    assert v[0, 1] == 0.0 and v[1, 1] == pytest.approx(2.39e-6, rel=0.01) and v[1, 2] < 1e-300


def test_linear_layer_raman_noise_keeps_zero_point_part():
    tg = toggles(beta2=True, kerr=True, raman=True, kerr_fluct=True, raman_noise=True)

    def tr(T):
        return np.trace(fiber_quantum_linear(linear_params(temperature_k=T), tg).noise_cov["raman"])

    hot, cold = tr(295.0), tr(1e-3)
    assert 0 < cold < hot
