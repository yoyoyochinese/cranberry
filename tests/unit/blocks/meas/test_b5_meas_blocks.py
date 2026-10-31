"""WP-B5 MEAS blocks: esa_zero_span / esa_swept_spectrum, pulse_sample, mode_overlap (contracts + formulas)."""

import math

import numpy as np
import pytest
from scipy.special import spence

import cranberry.blocks.det.testing as T
from cranberry.blocks.meas import esa, mode_overlap as mo
from cranberry.blocks.meas.pulse_sample import expected_r_meas
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.fft import fft_agrawal, omega_grid
from cranberry.core.registry import make_block
from cranberry.core.rng import stream
from cranberry.offline import DoubleRegistrationError
from cranberry.offline.loss_registry import LossRegistry


def mk(obj):
    return make_block(BlockSpec.from_obj(obj))


# --------------------------------------------------------------------------------------------- ESA


def test_log_var_ratio_constants():
    """∫Li₂(|ρ|²)/∫|ρ|² for rect (|ρ| = sinc) and Gaussian (|ρ|² = e^{−2πx²}) RBW shapes."""
    li2 = lambda x: spence(1 - x)  # noqa: E731
    x = np.linspace(-3000, 3000, 6_000_001)
    s2 = np.sinc(x) ** 2
    assert li2(s2).sum() / s2.sum() == pytest.approx(esa.LOG_VAR_RATIO["rect"], rel=2e-4)
    g = np.exp(-2 * np.pi * x[np.abs(x) < 10] ** 2)
    assert li2(g).sum() / g.sum() == pytest.approx(esa.LOG_VAR_RATIO["gaussian"], rel=1e-6)
    assert esa.LOG_BIAS_DB == pytest.approx(-2.507, abs=1e-3) and esa.VOLTAGE_BIAS_DB == pytest.approx(-1.049, abs=1e-3)


def test_zero_span_noise_off_is_expectation():
    for avg in esa.AVERAGE_TYPES:
        r = mk(T.zero_span(noise=False, average=avg)).evaluate(0.05)  # −10 dB
        np.testing.assert_allclose(r["trace_db"], -10.0, atol=1e-9)
    r = mk(T.zero_span(noise=False, average="log")).evaluate(0.5)
    np.testing.assert_allclose(r["raw_db"], esa.LOG_BIAS_DB, atol=1e-12)


def test_zero_span_contract():
    obj = T.zero_span()
    del obj["params"]["vbw_hz"]
    with pytest.raises(MissingParameterError):
        mk(obj)
    with pytest.raises(ValueError):
        mk(T.zero_span(average="peak"))
    with pytest.raises(ValueError):  # radiometer noise needs the block's rng stream
        mk(T.zero_span()).evaluate(0.5)


def test_rbw_filter_shape_and_normalization():
    blk = mk(T.zero_span(rbw_shape="gaussian", oversample=6, sweep_time_s=0.02))
    z2 = blk._z2(stream(1, "esa"))
    assert z2.mean() == pytest.approx(1.0, rel=0.02)
    h2 = blk._rbw_power(np.linspace(-5e6, 5e6, 200001))
    assert h2.sum() * 50.0 == pytest.approx(blk.enbw_hz, rel=1e-6)


def test_swept_noise_off_and_subtraction():
    sa = mk(T.swept())
    c, r = 10 ** 2.8, 0.1
    out = sa.evaluate(lambda f: r + 1 / c, lambda f: 1 + 1 / c, lambda f: 1 / c + 0 * f)
    np.testing.assert_allclose(10 ** (out["normalized_db"] / 10), (c * r + 1) / (c + 1), rtol=1e-12)
    np.testing.assert_allclose(10 ** (out["corrected_db"] / 10), r, rtol=1e-12)
    assert mk(T.swept(subtract=False)).evaluate(lambda f: 1.0, lambda f: 1.0, lambda f: 0.0 * f)["corrected_db"] is None


# ------------------------------------------------------------------------------------- pulse_sample


def test_pulse_sample_expectations():
    r = np.ones(32)
    r[0] = 0.1
    c = np.array([1.0, 0.3, 0.05])
    ps = mk(T.pulse_sample(on=False))
    assert ps.evaluate(r, c=c, clearance_lin=100.0)["expected"][0] == pytest.approx(0.1)  # toggles off: identity
    assert expected_r_meas(r, "integrating", [0], c, 0.0)[0] == pytest.approx((0.1 + 31) / 32)
    # crosstalk mixes vacuum neighbours in (pulse-resolved), gating removes them
    exp = (0.1 + 0.09 + 0.0025) / (1 + 0.09 + 0.0025)
    assert expected_r_meas(r, "pulse_resolved", [0], c, 0.0)[0] == pytest.approx(exp)
    assert expected_r_meas(r, "gated", [0], c, 0.0)[0] == pytest.approx(0.1)
    # electronic noise = loss C/(1+C) for a single slot without crosstalk
    assert expected_r_meas(r, "gated", [0], [1.0], 1 / 99)[0] == pytest.approx((99 * 0.1 + 1) / 100)
    with pytest.raises(ValueError):
        mk(T.pulse_sample(target_slots=[32]))


# ------------------------------------------------------------------------------------- mode_overlap


def test_mode_overlap_registration_once_per_path_pair():
    a = mk(T.overlap(name="mm1"))
    assert a.eta == pytest.approx(0.99097, abs=1e-5)
    (e,) = a.losses()
    assert e.kind == "mode_overlap_spatial" and e.observer == "bhd"
    reg = LossRegistry(paths={"sqz": {"segments": [], "observers": ["bhd"]}})
    reg.add_all(a.losses())
    with pytest.raises(DoubleRegistrationError):  # second spatial overlap on (sqz, bhd), G06 / VT-F3
        reg.add_all(mk(T.overlap(name="mm2", model="eta", eta=0.98)).losses())
    reg.add_all(mk(T.overlap(name="lo", overlap="lo_shape", model="eta", eta=0.951)).losses())  # other kind: ok
    assert reg.eta("sqz", "bhd") == pytest.approx(a.eta * 0.951)
    off = mk(T.overlap(name="mm3", on=False))
    assert not off.losses()[0].enabled


def test_mode_overlap_contract():
    with pytest.raises(ValueError):
        mk(T.overlap(model="gaussian_pulse"))  # temporal model on a spatial overlap
    with pytest.raises(ValueError):
        mk(T.overlap(segment="s1"))  # both segment and observer
    obj = T.overlap()
    del obj["params"]["d_m"]
    with pytest.raises(MissingParameterError):
        mk(obj)


def test_gaussian_pulse_kappa_numeric_t03_v1():
    """T03 V1: analytic κ vs numerical ∫u*u_s dT (squeezed delayed by Δt with GDD φ₂, conv §8)."""
    t0, n, dt = 1.0, 2**14, 0.01
    t = (np.arange(n) - n // 2) * dt
    u = np.exp(-t**2 / (2 * t0**2))
    om = omega_grid(n, dt)
    from cranberry.core.fft import ifft_agrawal

    for d, chi in [(0.0, 0.0), (1.3, 0.0), (0.0, 2.0), (-2.1, -4.5), (3.0, 5.0)]:
        us = ifft_agrawal(fft_agrawal(u + 0j, dt, t0=t[0]) * np.exp(1j * om * d + 1j * chi * t0**2 * om**2),
                          dt, t0=t[0])
        k_num = mo.temporal_overlap(u, us, dt)
        k_ana = mo.gaussian_pulse_kappa(d, 2 * chi * t0**2, t0)
        assert abs(k_num - k_ana) < 1e-6
    blk = mk(T.overlap(overlap="temporal", model="gaussian_pulse", dt_s=1e-12, phi2_s2=0.0, t0_s=0.6e-12))
    assert blk.eta == pytest.approx(math.exp(-1 / (2 * 0.36)), rel=1e-12) and blk.losses()[0].kind == "mode_overlap_temporal"


def test_lo_shape_overlap_analytic():
    """Rectangular LO vs one-sided exponential mode √(2κ)e^{−κt}: max 2(1 − e^{−x})²/x = 0.814529 at x = 1.2564."""
    t = np.linspace(-5, 20, 250001)
    u = np.where(t >= 0, np.sqrt(2) * np.exp(-np.clip(t, 0, None)), 0.0)
    best = mo.best_shape_overlap(u, t, "rect", (0.1, 1.0))
    assert best["eta"] == pytest.approx(0.814529, abs=2e-5) and best["b"] == pytest.approx(1.2564, abs=2e-3)
    assert mo.best_shape_overlap(u, t, "one_sided_exp", (0.2, 0.7))["eta"] == pytest.approx(1.0, abs=1e-6)
