"""VT-F2 (WP-F1 part): chunk-length invariance (bit-identical results for any ``chunk_samples``, G11,
PLAN §6.4) and f_sim doubling convergence (§6.4: f_sim ×2, ×4; change shrinks with the step), run
through ``runner.run``. Also: compiled = reference through the runner (VT-44, 1e−12), toggle
overrides, the seed record. Noise blocks of WP-B2 add their own VT-F2 cases.
"""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.pipeline import spec_from_obj
from cranberry.core.rates import with_f_sim
from cranberry.core.rng import stream_key
from cranberry.core.testing import toy
from cranberry.runner.run import run

pytestmark = pytest.mark.numba
FS = 1.0e6
REC = ["v.n1", "v.n2", "v.y", "ctrl.pid.u", "psi.plant.x", "act.u", "n", "t"]


def noisy_loop_obj(kp=0.4, ki=3e3):
    return [
        toy("toy_noise_ar1", "n1", out="v.n1", unit="V", a=0.999, b=0.05),
        toy("toy_noise_ar1", "n2", 3, out="v.n2", unit="V", a=0.9, b=0.2),
        toy("toy_sum", "sens", ins=["psi.plant.x", "v.n1", "v.n2"], gains=[1.0, 1.0, 1.0], out="v.y", unit="V"),
        toy("toy_pid", "pid", 2, in_="v.y", out="ctrl.pid.u", unit="V", kp=kp, ki=ki, kd=0.0, setpoint=0.0),
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out",
                                                           "tau_s": 3.37e-6, "method": "thiran1"}},
        toy("toy_sat", "act", in_="delay.dl.out", out="act.u", unit="V", gain=1.0, limit=0.8),
        toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="V", k=2e4),
    ]


def spec(obj, chunk, n=9000, fs=FS, decim=7, slots=REC, toggles=None, seed=20261031):
    return RunSpec(run=RunSettings(duration_s=n / fs, f_sim_hz=fs, chunk_samples=chunk, seed=seed),
                   pipeline=spec_from_obj(obj), toggles=toggles or {},
                   measure={"record": {"slots": list(slots), "decim": decim}})


def test_chunk_length_invariance_bitwise():
    base = run(spec(noisy_loop_obj(), 65536))
    assert base.mode == "compiled"
    assert base.series["v.n1"].shape == (len(base.t_s),) and len(base.t_s) == -(-9000 // 7)
    np.testing.assert_array_equal(base.series["n"], np.arange(0, 9000, 7))
    np.testing.assert_array_equal(base.series["t"], base.t_s)
    assert np.std(base.series["v.n1"]) > 0
    for chunk in (1, 7, 1000, 4095, 4096, 4097, 8999):
        res = run(spec(noisy_loop_obj(), chunk))
        for k in REC:
            np.testing.assert_array_equal(res.series[k], base.series[k], err_msg=f"chunk={chunk} slot={k}")
        # every model slot (draw / record buffers are scratch space whose size depends on the chunk)
        for name in base.final_state.layout.slot_names():
            if not name.startswith(("rng.", "rec.")):
                np.testing.assert_array_equal(res.final_state[name], base.final_state[name], err_msg=name)


def test_reference_runner_matches_compiled():
    a = run(spec(noisy_loop_obj(), 1000, n=3000), mode="compiled")
    b = run(spec(noisy_loop_obj(), 999, n=3000), mode="reference")
    assert b.mode == "reference"
    for k in REC:
        np.testing.assert_allclose(b.series[k], a.series[k], rtol=1e-12, atol=1e-14)


def test_seeds_and_noise_source_independence():
    """Switching one noise source off leaves the other's realisation unchanged (G11, G13)."""
    on = run(spec(noisy_loop_obj(), 2048, n=4000))
    off = run(spec(noisy_loop_obj(), 2048, n=4000, toggles={"n2.enabled": False}))
    np.testing.assert_array_equal(off.series["v.n1"], on.series["v.n1"])
    assert np.all(off.series["v.n2"] == 0.0) and np.any(on.series["v.n2"] != 0.0)
    assert on.seeds["n1"] == list(stream_key("n1")) and on.seeds["root"] == 20261031
    other = run(spec(noisy_loop_obj(), 2048, n=4000, seed=1))
    assert not np.array_equal(other.series["v.n1"], on.series["v.n1"])


def deterministic_loop_obj():
    """Sine disturbance into a PI loop with a delay: deterministic, so f_sim runs are comparable."""
    obj = noisy_loop_obj(kp=0.3, ki=2e3)
    obj = [toy("toy_sine", "dist", out="v.d", unit="V", amp=0.5, freq_hz=3e3)] + obj[2:]
    obj[1] = toy("toy_sum", "sens", ins=["psi.plant.x", "v.d"], gains=[1.0, 1.0], out="v.y", unit="V")
    return obj


def test_f_sim_doubling_converges():
    """f_sim ×1, ×2, ×4 (PID update rate kept fixed, a device property): the recorded loop response
    at common times converges with the step size."""
    n, decim = 4000, 10
    y, t0 = [], None
    for f in (1, 2, 4):
        # duration fixed; decim scales with f_sim so the record times coincide
        s = with_f_sim(spec(deterministic_loop_obj(), 4096, n=n, decim=decim * f, slots=["v.y"]), f,
                       keep_rate=["pid"])
        r = run(s)
        assert r.spec.run.f_sim_hz == FS * f and r.spec.run.n_samples == n * f
        t0 = r.t_s if t0 is None else t0
        np.testing.assert_allclose(r.t_s, t0, rtol=1e-12)
        y.append(r.series["v.y"])
    amp = np.max(np.abs(y[2]))
    e1, e2 = np.max(np.abs(y[0] - y[2])), np.max(np.abs(y[1] - y[2]))
    assert amp > 0.05
    assert e1 < 0.1 * amp  # already close at the preset f_sim
    # vs the ×4 run: first order gives e2/e1 = (1/4)/(3/4) = 1/3, second order 0.2
    assert 0.15 * e1 < e2 < 0.45 * e1


def test_with_f_sim_scales_rate_divs():
    s = spec(noisy_loop_obj(), 100)
    s2 = with_f_sim(s, 4, keep_rate=["pid", "n2"])
    divs = {c.name: c.rate_div for c in s2.pipeline.children}
    assert divs["pid"] == 8 and divs["n2"] == 12 and divs["n1"] == 1 and s2.run.f_sim_hz == 4 * FS
    with pytest.raises(KeyError):
        with_f_sim(s, 2, keep_rate=["nope"])
