"""ACT blocks: ``pzt_fiber_stretcher``, ``pzt_mirror``, ``hv_amp_914``, ``thermal_actuator`` — gain modes
(conv §4), driver slot (G08, F6 names), hysteresis, saturation, noise, losses, toggles off = identity."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

FS = 1.0e6
LAM = 1548.8e-9
OFF = {"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": False}


def make(obj, seed=2, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def tri(amp=6.0, period=1e-3):
    return {"block": "ramp", "name": "u", "params": {"out": "v.u", "waveform": "triangle", "amplitude_v": amp,
                                                     "period_s": period, "offset_v": 0.0, "toggles": {"enabled": True}}}


def stretcher(toggles=None, kind="pzt_fiber_stretcher", **kw):
    p = {"in": "v.u", "gain_mode": "direct", "k_rad_per_v": 1.0, "polarity": 1.0, "modes": [],
         "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -5.0, "drive_max_v": 5.0,
         "hysteresis_frac": 0.1, "segment": "arm1", "driver": "stretcher", "n_sub": 4, "sub": [1, 3],
         "ring_len": 0, "insertion_loss_db": None, "toggles": dict(toggles or OFF)}
    p.update(kw)
    return {"block": kind, "name": "a", "params": p}


def trace(obj, n, slots):
    blk, st, pr = make(obj)
    out = {s: [] for s in slots}
    for _ in range(n):
        run_reference(blk, st, pr, 1)
        for s in slots:
            out[s].append(np.array(st[s], dtype=float))
    return blk, st, pr, {s: np.array(v) for s, v in out.items()}


def test_toggles_off_identity_and_driver_split():
    _, st, _, tr = trace([tri(), stretcher()], 300, ["v.u", "act.a.dpsi", "drv.arm1.stretcher"])
    np.testing.assert_array_equal(tr["act.a.dpsi"], tr["v.u"])
    d = tr["drv.arm1.stretcher"]
    np.testing.assert_array_equal(d[:, [0, 3]], 0.0)
    np.testing.assert_allclose(d[:, 1] + d[:, 2], tr["v.u"], rtol=1e-15)  # split over sub [1, 3)


@pytest.mark.parametrize("mode, keys, k", [
    ("dopl", dict(dopl_full_m=16e-6, v_full_v=10.0, lambda0_m=LAM), 2 * math.pi * 16e-6 / (LAM * 10)),
    ("elongation", dict(dl_full_m=16e-6, v_full_v=10.0, lambda0_m=LAM, n_eff=1.468, xi_eo=0.784),
     2 * math.pi * 1.468 * 0.784 * 16e-6 / (LAM * 10)),
])
def test_stretcher_gain_modes(mode, keys, k):
    """δψ = K u with K from δOPL (no ξ_eo) or from δL (ξ_eo β δL), conv §4; polarity flips the sign."""
    obj = stretcher(gain_mode=mode, polarity=-1.0, **keys)
    del obj["params"]["k_rad_per_v"]
    blk, st, pr, tr = trace([tri(), obj], 50, ["v.u", "act.a.dpsi"])
    np.testing.assert_allclose(tr["act.a.dpsi"], -k * tr["v.u"], rtol=1e-12)
    assert blk.children[1].gain_rad_per_v(pr) == pytest.approx(k, rel=1e-15)
    if mode == "dopl":
        assert k == pytest.approx(6.49, abs=0.005)
    else:
        assert k == pytest.approx(7.47, abs=0.01)  # T16 §3.1.3: ≈ 7.5 rad/V for the elongation reading


def test_mirror_gain():
    obj = stretcher(kind="pzt_mirror", gain_mode="mirror", d_x_m_per_v=1e-8, theta_i_rad=0.3, lambda0_m=1064e-9)
    del obj["params"]["k_rad_per_v"]
    _, _, _, tr = trace([tri(), obj], 20, ["v.u", "act.a.dpsi"])
    k = 2 * (2 * math.pi / 1064e-9) * math.cos(0.3) * 1e-8
    np.testing.assert_allclose(tr["act.a.dpsi"], k * tr["v.u"], rtol=1e-12)


def test_hysteresis_play_operator_and_saturation():
    """Play operator of width w = 0.1·10 V: the output lags a triangle by w/2 on each branch; the drive
    is clipped at ±5 V with the saturation flag."""
    tog = {**OFF, "hysteresis": True, "saturation": True}
    _, _, _, tr = trace([tri(amp=6.0), stretcher(tog)], 1000, ["v.u", "act.a.dpsi", "act.a.sat"])
    u, y, sat = tr["v.u"], tr["act.a.dpsi"], tr["act.a.sat"]
    uc = np.clip(u, -5, 5)
    rising = np.diff(u, prepend=u[0]) > 0
    sel = rising & (np.arange(u.size) > 600) & (np.abs(uc) < 4.4)
    np.testing.assert_allclose(y[sel], uc[sel] - 0.5, atol=1e-12)
    assert np.max(y) == pytest.approx(4.5) and np.min(y) == pytest.approx(-4.5)
    np.testing.assert_array_equal(sat, np.where(u > 5, 1, np.where(u < -5, -1, 0)))


def test_driver_ring_history():
    """ring_len > 0: each step pushes a row (head = newest), older rows keep the history (F6 reads it)."""
    blk, st, pr = make([tri(), stretcher(ring_len=3)])
    vals = []
    for _ in range(5):
        run_reference(blk, st, pr, 1)
        vals.append(float(st["act.a.dpsi"]))
    ring = st["drv.arm1.stretcher"]
    head = st.ring_head("drv.arm1.stretcher")
    for d in range(3):
        assert ring[(head - d) % 3, 1] * 2 == pytest.approx(vals[-1 - d], rel=1e-15)
    with pytest.raises(ValueError, match="rate_div"):
        make([tri(), {**stretcher(ring_len=3), "rate_div": 2}])


def test_losses_registered_with_toggle():
    obj = stretcher(insertion_loss_db=0.5, loss_applies_to=["gbs", "bhd"], loss_source="UNVERIFIED (915B)")
    blk, _, _ = make([tri(), obj])
    (el,) = blk.losses()
    assert el.segment == "arm1" and el.eta == pytest.approx(10 ** -0.05) and not el.enabled and el.kind == "component"
    obj["params"]["toggles"] = {**OFF, "insertion_loss": True}
    blk, _, _ = make([tri(), obj])
    assert blk.losses()[0].enabled
    with pytest.raises(MissingParameterError):
        bad = stretcher()
        del bad["params"]["polarity"]
        make([tri(), bad])
    with pytest.raises(ValueError, match="gain_mode"):
        make([tri(), stretcher(gain_mode="mirror")])


def amp(toggles, **kw):
    p = {"in": "v.u", "out": "v.out", "poles_hz": [29.5e3, 112e3], "discretisation": "zoh", "prewarp_hz": None,
         "input_min_v": -5.0, "input_max_v": 5.0, "noise_v_per_rthz": 0.0, "slew_v_per_s": 1e5, "toggles": toggles}
    p.update(kw)
    return {"block": "hv_amp_914", "name": "amp", "params": p}


def test_914_toggles_and_noise():
    off = {"poles": False, "slew": False, "saturation": False, "noise": False}
    _, _, _, tr = trace([tri(), amp(off)], 200, ["v.u", "v.out"])
    np.testing.assert_array_equal(tr["v.out"], tr["v.u"])
    _, _, _, tr = trace([tri(), amp({**off, "slew": True})], 200, ["v.u", "v.out"])
    assert np.max(np.abs(np.diff(tr["v.out"]))) <= 1e5 / FS + 1e-15
    _, _, _, tr = trace([tri(), amp({**off, "saturation": True})], 1000, ["v.u", "v.out"])
    np.testing.assert_array_equal(tr["v.out"], np.clip(tr["v.u"], -5, 5))
    const = {"block": "ramp", "name": "u", "params": {"out": "v.u", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                      "period_s": 1.0, "offset_v": 0.0, "toggles": {"enabled": True}}}
    blk, st, pr = make([const, amp({**off, "noise": True}, noise_v_per_rthz=1e-6)])
    xs = []
    for _ in range(200):
        run_compiled(blk, st, pr, 100)
        xs.append(float(st["v.out"]))
    # white density e_n up to f_sim/2: per-sample std = e_n·√(f_sim/2)
    assert np.std(xs) == pytest.approx(1e-6 * math.sqrt(FS / 2), rel=0.25)


def test_914_dc_gain_and_response():
    blk, st, pr = make([tri(), amp({"poles": True, "slew": False, "saturation": False, "noise": False})])
    h = blk.children[1].freq_response([1.0, 29.5e3])
    assert abs(h[0]) == pytest.approx(1.0, rel=1e-6)
    assert abs(h[1]) == pytest.approx(1 / math.sqrt(2) / abs(1 + 1j * 29.5 / 112), rel=0.01)


def tec(**kw):
    p = {"in": "v.u", "segment": "spool", "driver": "tec", "n_sub": 2, "sub": [0, 2], "ring_len": 0,
         "k_t_k_per_v": 2.0, "tau_th_s": 1e-4, "dt_min_k": -0.5, "dt_max_k": 0.5, "toggles": {"saturation": True}}
    p.update(kw)
    return {"block": "thermal_actuator", "name": "tec", "params": p}


def test_tec_saturation_copy_and_decimation():
    _, _, _, tr = trace([tri(amp=1.0, period=4e-3), tec()], 4000, ["v.u", "act.tec.dT", "act.tec.sat",
                                                                   "drv.spool.tec"])
    assert np.max(tr["act.tec.dT"]) == 0.5 and np.min(tr["act.tec.dT"]) == -0.5 and np.any(tr["act.tec.sat"] == 1)
    np.testing.assert_array_equal(tr["drv.spool.tec"][:, 0], tr["act.tec.dT"])  # temperature copied, not split
    np.testing.assert_array_equal(tr["drv.spool.tec"][:, 1], tr["act.tec.dT"])
    blk, st, pr = make([tri(amp=0.1), {**tec(tau_th_s=1e-3), "rate_div": 10}])
    run_reference(blk, st, pr, 1)
    x1 = float(st["act.tec.dT"])
    a = math.exp(-10 / FS / 1e-3)
    assert x1 == pytest.approx((1 - a) * 2.0 * -0.1, rel=1e-12)  # first step at n = 0 with h = 10/f_sim
