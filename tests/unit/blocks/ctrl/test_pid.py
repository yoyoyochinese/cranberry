"""``pid``: DSC1 law and sign (D53), discretisations vs C(s), anti-windup, bumpless, toggles."""

import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy
from cranberry.hardware.mapping import PidCoefficients, map_to_canonical

FS = 1.0e6
ZERO = dict(kp=0.0, ki_per_s=0.0, kd_s=0.0, omega_is_rad_per_s=0.0, omega_ds_rad_per_s=0.0, ki2_per_s2=0.0,
            setpoint_v=0.0, offset_v=0.0, u_min_v=-4.0, u_max_v=4.0)
ON = {"feedback": True, "p_term": True, "i_term": True, "d_term": True}


def pid_obj(disc="tustin", aw="none", toggles=None, device=None, rate_div=1, inp="v.x", **kw):
    p = {"in": inp, "out": "ctrl.p.u", "device": device, "anti_windup": aw, "toggles": dict(toggles or ON),
         **{**ZERO, **kw}}
    if disc is not None:
        p["discretisation"] = disc
    o = {"block": "pid", "name": "p", "params": p}
    if rate_div != 1:
        o["rate_div"] = rate_div
    return o


def src(amp=0.3, a=0.99):
    return toy("toy_noise_ar1", "x", out="v.x", unit="V", a=a, b=amp)


def make(obj, seed=1, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def trace(obj, n, slots, fs=FS):
    blk, st, pr = make(obj, fs=fs)
    out = {s: np.empty(n) for s in slots}
    for i in range(n):
        run_reference(blk, st, pr, 1)
        for s in slots:
            out[s][i] = float(st[s])
    return out


D53 = {"p": 0.000241, "i": 0.001406, "d": 0.01108, "p_enabled": True, "i_enabled": True, "d_enabled": True}


def test_dsc1_d53_settings_map_to_canonical():
    c, _ = map_to_canonical(D53, "dsc1")
    assert c.kp == pytest.approx(-2.41e-4, rel=1e-12)
    assert c.ki_per_s == pytest.approx(-703.0, rel=1e-12)
    assert c.kd_s == pytest.approx(-2.216e-8, rel=1e-12)


def test_dsc1_preset_reproduces_native_dsc1_law_and_sign():
    """Canonical gains (e = v − S) in the block = DSC1's own law with e_DSC1 = S − v (p.18–19), entries in
    p.20 units (D53), PID every 2 samples at 1 MHz (h = 2 µs, device), CV limits ±4 V."""
    c, _ = map_to_canonical(D53, "dsc1")
    sp = 0.05
    obj = [src(amp=0.4, a=0.98), pid_obj(disc=None, device="dsc1", rate_div=2, setpoint_v=sp, **c.as_params())]
    tr = trace(obj, 4000, ["v.x", "ctrl.p.u"])
    acc = ed1 = 0.0
    u_ref = np.empty(2000)
    for k, v in enumerate(tr["v.x"][::2]):
        ed = sp - v  # DSC1 error
        acc += D53["i"] * ed
        u = D53["p"] * ed + acc + D53["d"] * (ed - ed1)
        ed1 = ed
        u_ref[k] = min(max(u, -4.0), 4.0)
    np.testing.assert_allclose(tr["ctrl.p.u"][::2], u_ref, rtol=1e-12, atol=1e-15)
    np.testing.assert_array_equal(tr["ctrl.p.u"][1::2], tr["ctrl.p.u"][::2])  # held between PID updates
    # sign: input above the setpoint (e_DSC1 < 0) with positive DSC1 gains → negative P action
    blk, st, pr = make([toy("toy_gain", "x", in_="t", out="v.x", unit="V", gain=0.0),
                        pid_obj(disc=None, device="dsc1", setpoint_v=-1.0, **c.as_params())])
    run_reference(blk, st, pr, 1)
    assert float(st["ctrl.p.u"]) < 0


@pytest.mark.parametrize("disc", ["dsc1", "tustin", "forward_euler"])
@pytest.mark.parametrize("terms", [dict(ki_per_s=2e3, omega_is_rad_per_s=50.0), dict(ki2_per_s2=1e5),
                                   dict(kd_s=2e-6), dict(kd_s=2e-6, omega_ds_rad_per_s=2e5)])
def test_discretisations_match_continuous_pid_below_fs_over_20(disc, terms):
    """V8, term by term: |C(e^{iωh})| vs |C(s)| within 1 % for f ≤ f_s/20; Tustin also in phase
    (complex 1 %); the Euler forms differ in phase by at most half a sample (ωh/2) plus 0.2°."""
    if disc != "tustin" and "omega_ds_rad_per_s" in terms:
        pytest.skip("Euler forms of a filtered derivative near its pole (32 kHz) are not 1 % accurate at f_s/20")
    blk, st, pr = make([src(), pid_obj(disc=disc, **terms)])
    leaf = blk.children[1]
    g = {**ZERO, **terms}
    c = PidCoefficients(0.0, g["ki_per_s"], g["kd_s"], g["omega_is_rad_per_s"], g["omega_ds_rad_per_s"] or None,
                        True, g["ki2_per_s2"])
    f = np.logspace(1, np.log10(FS / 20), 30)
    cd, cc = leaf.freq_response(f, pr), c.transfer(f)
    tol = 0.02 if "ki2_per_s2" in terms else 0.01  # two cascaded integrators: (x/tan x)² = 0.984 at f_s/20
    np.testing.assert_allclose(np.abs(cd), np.abs(cc), rtol=tol)
    if disc == "tustin":
        np.testing.assert_allclose(cd, cc, rtol=tol)
    else:
        n_int = 2 if "ki2_per_s2" in terms else 1
        assert np.all(np.abs(np.angle(cd / cc)) <= n_int * np.pi * f / FS + np.radians(0.2))


def test_freq_response_matches_simulation():
    """freq_response (used for z-domain loop analytics) equals the simulated sine response."""
    g = dict(kp=0.2, ki_per_s=5e3, kd_s=1e-6, omega_ds_rad_per_s=3e5, ki2_per_s2=2e5)
    f = 7.3e3
    obj = [{"block": "sine_injection", "name": "x", "params": {"in": None, "out": "v.x", "unit": "V",
                                                               "freq_hz": f, "amplitude": 0.01,
                                                               "toggles": {"enabled": True}}},
           pid_obj(disc="tustin", u_min_v=-1e9, u_max_v=1e9, **g)]
    blk, st, pr = make(obj)
    run_compiled(blk, st, pr, 200_000)
    n0 = st.n
    n = int(round(50 * FS / f))
    x, u = np.empty(n), np.empty(n)
    for i in range(n):
        run_compiled(blk, st, pr, 1)
        x[i], u[i] = float(st["v.x"]), float(st["ctrl.p.u"])
    t = (n0 + np.arange(n)) / FS
    w = np.exp(-2j * np.pi * f * t)
    # remove the integrators' ramp/offset drift by fitting the sine component only
    a = np.column_stack([np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t), np.ones(n), t, t**2])
    cu = np.linalg.lstsq(a, u, rcond=None)[0]
    h_sim = (cu[0] - 1j * cu[1]) / (np.sum(x * w) * 2 / n)
    np.testing.assert_allclose(h_sim, blk.children[1].freq_response(f, pr)[0], rtol=2e-3)


def test_anti_windup_modes():
    """Constant error +1 V against u_max = 1 V for 10 ms: ``none`` winds up, ``clamp`` holds the
    integrator at the limit, ``back_calculation`` settles it near the limit; all outputs clip."""
    res = {}
    for aw in ("none", "clamp", "back_calculation"):
        kw = {"aw_tracking_s": 1e-4} if aw == "back_calculation" else {}
        obj = [{"block": "ramp", "name": "x", "params": {"out": "v.x", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                         "period_s": 1.0, "offset_v": 1.0, "toggles": {"enabled": True}}},
               pid_obj(disc="dsc1", aw=aw, kp=0.1, ki_per_s=1e3, u_min_v=-1.0, u_max_v=1.0, **kw)]
        blk, st, pr = make(obj)
        run_compiled(blk, st, pr, 10_000)
        res[aw] = float(st["ctrl.p.x"][1])
        assert 0.99 <= float(st["ctrl.p.u"]) <= 1.0
    assert res["none"] == pytest.approx(10.0, rel=1e-3)  # h Σ K_i e = 1e3 · 10 ms
    assert 0.89 < res["clamp"] <= 0.9  # integration stops when it would push past the limit
    assert 0.85 < res["back_calculation"] < 1.0


def test_bumpless_integrator_switch_on():
    """V9: K_i switched from 0 to 703 s⁻¹: the next output changes by h·K_i·e only (< 1 DSC1 LSB)."""
    obj = [{"block": "ramp", "name": "x", "params": {"out": "v.x", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                     "period_s": 1.0, "offset_v": 0.5, "toggles": {"enabled": True}}},
           pid_obj(disc="dsc1", kp=-2.41e-4, rate_div=2)]
    blk, st, pr = make(obj)
    run_compiled(blk, st, pr, 1000)
    u0 = float(st["ctrl.p.u"])
    pr["p.ki_per_s"] = -703.0
    run_compiled(blk, st, pr, 2)
    du = float(st["ctrl.p.u"]) - u0
    assert du == pytest.approx(2e-6 * -703.0 * 0.5, rel=1e-12)
    assert abs(du) < 8.0 / 2**12


def test_toggles():
    g = dict(kp=0.5, ki_per_s=1e3, kd_s=1e-6, offset_v=0.25)
    off = {"feedback": False, "p_term": True, "i_term": True, "d_term": True}
    tr = trace([src(), pid_obj(toggles=off, **g)], 50, ["ctrl.p.u"])
    np.testing.assert_array_equal(tr["ctrl.p.u"], 0.25)  # loop open: u = clip(offset)
    only_p = {"feedback": True, "p_term": True, "i_term": False, "d_term": False}
    tr = trace([src(), pid_obj(toggles=only_p, **g)], 50, ["v.x", "ctrl.p.u"])
    np.testing.assert_allclose(tr["ctrl.p.u"], np.clip(0.5 * tr["v.x"] + 0.25, -4, 4), rtol=1e-15)


def test_contract_and_validation():
    blk, st, pr = make([src(), pid_obj()])
    leaf = blk.children[1]
    from cranberry.hardware.mapping import CANONICAL_KEYS

    assert set(CANONICAL_KEYS) <= set(leaf.runtime_params) and "feedback" in leaf.toggles
    bad = pid_obj()
    del bad["params"]["anti_windup"]
    with pytest.raises(MissingParameterError):
        make([src(), bad])
    with pytest.raises(MissingParameterError):
        make([src(), pid_obj(disc=None)])  # no device → discretisation required (D34)
    with pytest.raises(ValueError, match="tustin"):
        make([src(), pid_obj(disc="moku")])
    assert math.isfinite(leaf.h) and leaf.h == 1 / FS
