"""opo_detuning_noise / opo_length_lock (WP-B6): toggles-off identity, compiled = reference, PSD and
ledger, thermal sign, loop suppression (T02 §2.4, §6.1)."""

import numpy as np
import pytest
from scipy import signal

import cranberry.path.testing as pt
from cranberry.blocks.source.opo_detuning import make_opo_detuning_noise, make_opo_length_lock
from cranberry.core.block import MissingParameterError
from cranberry.core.layout import SlotLayout
from cranberry.core.rng import stream

FS = 1.0e6
NU0, L = 3.49e14, 0.5
GAIN = -2 * np.pi * NU0 / L


def noise_params(**tg):
    return {"out": "opo.dn.delta_free", "segment": "opo", "nu0_hz": NU0, "round_trip_m": L,
            "h_m2_hz": 1e-24, "f_lo_hz": 20.0, "f_c_hz": 5e4, "peaks": [{"f_hz": 4e3, "q": 10.0, "area_m2": 4e-20}],
            "temperature_slot": "drv.opo.dT", "ddelta_dT_rad_s_per_k": -1.2e10,
            "toggles": {"length_noise": False, "peaks": False, "thermal": False, **tg}}


def lock_params(**kw):
    tg = kw.pop("toggles", {})
    return {"in": "opo.dn.delta_free", "out": "drv.opo.detuning", "f_ugf_hz": 1e3, "lock_light_ratio": 1.0,
            "offset_rad_s": 0.0, "sensor_psd": 0.0, **kw, "toggles": {"lock": True, "sensor_noise": False, **tg}}


def run(pipe, n, tables=None, compiled=False, seed=1):
    return pt.run_and_record(pipe, FS, n, tables, compiled=compiled, seed=seed)[0]


def noise_pipe(n, **tg):
    return [pt.pathtest("pathtest_table", "temp", "drv.opo.dT", (), "K", n_rows=n),
            {"block": "opo_detuning_noise", "name": "dn", "params": noise_params(**tg)},
            pt.record(["opo.dn.delta_free"], n)]


def test_all_off_is_zero_and_thermal_sign():
    n = 50
    temp = np.linspace(0.0, 1e-3, n)
    assert np.all(run(noise_pipe(n), n, {"drv.opo.dT_table": temp})["opo.dn.delta_free"] == 0.0)
    out = run(noise_pipe(n, thermal=True), n, {"drv.opo.dT_table": temp})["opo.dn.delta_free"]
    np.testing.assert_allclose(out, -1.2e10 * temp, rtol=1e-15)
    with pytest.raises(MissingParameterError):
        make_opo_detuning_noise({k: v for k, v in noise_params(thermal=True).items() if k != "temperature_slot"})


@pytest.mark.numba
def test_noise_compiled_equals_reference_and_psd():
    n = 1 << 16
    tab = {"drv.opo.dT_table": np.zeros(n)}
    pipe = noise_pipe(n, length_noise=True, peaks=True)
    ref = run(pipe, 4000, tab)["opo.dn.delta_free"]
    cmp_ = run(pipe, 4000, tab, compiled=True)["opo.dn.delta_free"]
    np.testing.assert_allclose(cmp_, ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(ref)))
    y = run(pipe, n, tab, compiled=True, seed=7)["opo.dn.delta_free"]
    f, s = signal.welch(y, fs=FS, nperseg=4096)  # one-sided (conv §7)
    blk = make_opo_detuning_noise(noise_params(length_noise=True, peaks=True))
    band = (f > 2e3) & (f < 2e4)  # ~150 bins × 31 segments: 3σ ≪ 10 %
    assert np.mean(s[band]) / np.mean(blk.psd_delta(f[band])) == pytest.approx(1.0, rel=0.1)


def test_bands_declared():
    blk = make_opo_detuning_noise(noise_params(length_noise=True))
    blk.init(SlotLayout(FS), stream(0, "dn"))
    (b,) = blk.bands()
    assert (b.quantity, b.targets, b.unit, b.f_hi_hz) == ("detuning", ("opo",), "rad/s", FS / 2)
    assert b.psd_outside(np.array([1e6]))[0] == pytest.approx(GAIN**2 * 1e-24 / 20**2 / ((1 + 5e4**2) * (1 + 20**2)), rel=1e-9)
    off = make_opo_detuning_noise(noise_params())
    off.init(SlotLayout(FS), stream(0, "dn"))
    assert off.bands()[0].variance_outside == 0.0


def lock_pipe(n, **kw):
    return [pt.pathtest("pathtest_table", "free", "opo.dn.delta_free", (), "rad/s", n_rows=n),
            {"block": "opo_length_lock", "name": "ll", "params": lock_params(**kw)},
            pt.record(["drv.opo.detuning", "opo.ll.err"], n)]


def test_lock_off_identity_and_steady_state():
    n = 6000  # (1 − 2π f_ugf h)^n ≈ e^{−38}: converged
    free = np.full(n, 2.0e5)
    out = run(lock_pipe(n, toggles={"lock": False}), n, {"opo.dn.delta_free_table": free})["drv.opo.detuning"]
    np.testing.assert_array_equal(out, free)
    rec = run(lock_pipe(n, lock_light_ratio=0.9, offset_rad_s=1e3), n, {"opo.dn.delta_free_table": free})
    assert rec["drv.opo.detuning"][0] == 2.0e5  # first sample before the integrator acts
    assert rec["drv.opo.detuning"][-1] == pytest.approx((1 - 0.9) * 2.0e5 + 1e3, rel=1e-9)  # (1 − r)Δ + Δ_off
    assert abs(rec["opo.ll.err"][-1]) < 1e-6
    with pytest.raises(ValueError, match="f_ugf_hz too high"):
        run(lock_pipe(10, f_ugf_hz=2e5), 10, {"opo.dn.delta_free_table": np.zeros(10)})


@pytest.mark.parametrize("f", [20.0, 100.0])
def test_lock_suppression_below_ugf(f):
    n = int(5 * FS / f)
    t = np.arange(n) / FS
    free = 1e5 * np.sin(2 * np.pi * f * t)
    out = run(lock_pipe(n), n, {"opo.dn.delta_free_table": free})["drv.opo.detuning"]
    tail = out[n // 2:]
    assert np.max(np.abs(tail)) / 1e5 == pytest.approx(f / 1e3, rel=0.02)


@pytest.mark.numba
def test_lock_sensor_noise_compiled_equals_reference():
    n = 2000
    tab = {"opo.dn.delta_free_table": np.zeros(n)}
    pipe = lock_pipe(n, sensor_psd=1e4, toggles={"sensor_noise": True})
    a = run(pipe, n, tab)["drv.opo.detuning"]
    b = run(pipe, n, tab, compiled=True)["drv.opo.detuning"]
    np.testing.assert_allclose(b, a, rtol=1e-12, atol=1e-12)
    assert np.std(a) > 0
    assert make_opo_length_lock(lock_params()).bands() == []
