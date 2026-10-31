"""VT-44 (block part, WP-B4): actuator / controller kernels = pure-Python reference (relative 1e−12,
PLAN §6.3), ξ_eo, mechanical-resonance discretisation, 914 response, stacked = parallel with the real
blocks, sign check, re-acquisition, AOM, TEC, f_sim convergence.

Sources (notes/00_summary.md §4 VT-44 row → T16 §8): V1 static gain 6.491 rad/V and range 51.93 rad,
V2 ξ_eo = 0.7837, V3 |H(f_m)| = Q ± 1 % and −90° ± 0.5° (f_m = 20 kHz, Q = 10, f_sim = 1 MHz),
V4 914 six points ± 5 % (2-pole) and slew model three points ± 10 %, V13 stacked = parallel (1e−9) and
offload u_1 → 0, V14 sign flip diverges, V16 fringe hop by n·2π (< 1e−6 rad), V17 AOM, V18 TEC,
V20 compiled = reference (1e−12) and f_sim doubling (5 %). The F1 framework part is
test_vt_44_delay_stack_compile.py.
"""

import math

import numpy as np
import pytest
from scipy import signal

import cranberry.core.testing  # noqa: F401  (toy blocks)
from cranberry.blocks.act.pzt import xi_eo
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj, spec_to_obj
from cranberry.core.rates import PIECE_SAMPLES
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba
FS = 1.0e6
LAM = 1548.8e-9
PID_GAINS = dict(kp=0.0, ki_per_s=0.0, kd_s=0.0, omega_is_rad_per_s=0.0, omega_ds_rad_per_s=0.0, ki2_per_s2=0.0,
                 setpoint_v=0.0, offset_v=0.0, u_min_v=-4.0, u_max_v=4.0)
ON = {"feedback": True, "p_term": True, "i_term": True, "d_term": True}


def make(obj, seed=7, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, seed)
    return blk, State(lay), Params(lay)


def pid(name, inp, out=None, rate_div=1, disc="tustin", aw="none", toggles=None, **gains):
    p = {"in": inp, "out": out or f"ctrl.{name}.u", "device": None, "discretisation": disc, "anti_windup": aw,
         "toggles": dict(toggles or ON), **{**PID_GAINS, **gains}}
    o = {"block": "pid", "name": name, "params": p}
    if rate_div != 1:
        o["rate_div"] = rate_div
    return o


def stretcher(name, inp, segment="arm1", modes=None, disc="zoh", prewarp=None, toggles=None, ring=0, **kw):
    p = {"in": inp, "gain_mode": "dopl", "dopl_full_m": 16e-6, "v_full_v": 10.0, "lambda0_m": LAM, "polarity": 1.0,
         "modes": modes if modes is not None else [{"f_hz": 20e3, "q": 10.0, "r": 1.0}],
         "discretisation": disc, "prewarp_hz": prewarp, "drive_min_v": -5.0, "drive_max_v": 5.0,
         "hysteresis_frac": 0.1, "segment": segment, "driver": "stretcher", "n_sub": 4, "sub": [1, 3],
         "ring_len": ring, "insertion_loss_db": None,
         "toggles": toggles or {"resonance": True, "hysteresis": False, "saturation": True, "insertion_loss": True}}
    p.update(kw)
    return {"block": "pzt_fiber_stretcher", "name": name, "params": p}


def gain_act(name, inp, k):
    """Static actuator δψ = k·u (pzt_fiber_stretcher in direct mode, no dynamics / limits)."""
    return stretcher(name, inp, segment=None, modes=[], gain_mode="direct", k_rad_per_v=k,
                     toggles={"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": False})


def const(name, out, value, unit="V"):
    """Constant source: a ramp with zero amplitude."""
    return {"block": "ramp", "name": name, "params": {"out": out, "waveform": "sawtooth", "amplitude_v": 0.0,
                                                      "period_s": 1.0, "offset_v": value, "toggles": {"enabled": True}}}


def dl(name, inp, tau=0.0):
    return {"block": "delay_buffer", "name": name, "params": {"in": inp, "out": f"delay.{name}.out", "tau_s": tau,
                                                              "method": "linear"}}


def all_blocks_obj():
    """Every WP-B4 kind with every effect switched on, in one closed loop."""
    sos = signal.butter(2, 30e3, fs=FS, output="sos").tolist()
    return [
        toy("toy_noise_ar1", "n1", out="v.n1", unit="V", a=0.9995, b=0.02),
        {"block": "sine_injection", "name": "inj", "params": {"in": "v.n1", "out": "v.inj", "unit": "V",
                                                              "freq_hz": 3.1e3, "amplitude": 0.05,
                                                              "toggles": {"enabled": True}}},
        {"block": "adder", "name": "sum", "params": {"ins": ["v.inj", "psi.plant.x"], "gains": [1.0, 0.3],
                                                     "out": "v.sum", "unit": "V"}},
        {"block": "adc", "name": "adc", "params": {"in": "v.sum", "out": "v.adc", "device": "dsc1",
                                                   "noise_vrms": 2e-4, "toggles": {"quantization": True,
                                                                                   "noise": True, "saturation": True}}},
        {"block": "lpf", "name": "lpf", "params": {"in": "v.adc", "out": "v.lp", "unit": "V", "order": 1,
                                                   "f_c_hz": 120e3, "method": "bilinear", "toggles": {"enabled": True}}},
        {"block": "band_split", "name": "bs", "params": {"in": "v.lp", "out_low": "v.lo", "out_high": "v.hi",
                                                         "unit": "V", "order": 2, "f_c_hz": 5e3, "method": "zoh",
                                                         "toggles": {"enabled": True}}},
        {"block": "iir", "name": "iir", "params": {"in": "v.hi", "out": "v.iir", "unit": "V", "sos": sos,
                                                   "toggles": {"enabled": True}}},
        {"block": "hpf", "name": "hpf", "params": {"in": "v.lo", "out": "v.hp", "unit": "V", "order": 1,
                                                   "f_c_hz": 10.0, "method": "bilinear", "toggles": {"enabled": True}}},
        {"block": "adder", "name": "err", "params": {"ins": ["v.hp", "v.iir"], "gains": [1.0, 1.0], "out": "v.err",
                                                     "unit": "V"}},
        pid("pid", "v.err", rate_div=2, disc="dsc1", aw="back_calculation", aw_tracking_s=2e-4, kp=-0.02,
            ki_per_s=-800.0, kd_s=-1e-8, omega_ds_rad_per_s=2e5, ki2_per_s2=-1e4),
        {"block": "lock_supervisor", "name": "sup", "params": {
            "u_in": "ctrl.pid.u", "e_in": "v.err", "true_in": "psi.plant.x", "out": "ctrl.sup.u",
            "relock": "fringe_hop", "u_min_v": -0.3, "u_max_v": 0.3, "t_sat_s": 0.0, "e_th": 0.5, "t_e_s": 1e-4,
            "t_hop_s": 2e-5, "u_2pi_v": 0.2, "toggles": {"relock": True}}},
        {"block": "dac", "name": "dac", "params": {"in": "ctrl.sup.u", "out": "v.dac", "device": "dsc1",
                                                   "noise_vrms": 1e-4, "r_load_ohm": 1e5,
                                                   "toggles": {"quantization": True, "noise": True,
                                                               "saturation": True}}, "rate_div": 2},
        {"block": "delay_buffer", "name": "dl", "params": {"in": "v.dac", "out": "delay.dl.out", "tau_s": 3.3e-6,
                                                           "method": "thiran1"}},
        {"block": "hv_amp_914", "name": "amp", "params": {
            "in": "delay.dl.out", "out": "v.amp", "poles_hz": [29.5e3, 112e3], "discretisation": "zoh",
            "prewarp_hz": None, "input_min_v": -5.0, "input_max_v": 5.0, "noise_v_per_rthz": 1e-7,
            "slew_v_per_s": 4.3e5, "toggles": {"poles": True, "slew": True, "saturation": True, "noise": True}}},
        {"parallel": [
            stretcher("str", "v.amp", modes=[{"f_hz": 20e3, "q": 10.0, "r": 1.2}, {"f_hz": 45e3, "q": 30.0, "r": -0.2}],
                      ring=6, hysteresis_frac=1e-4, toggles={"resonance": True, "hysteresis": True,
                                                             "saturation": True, "insertion_loss": True}),
            {"block": "pzt_mirror", "name": "mir", "params": {
                "in": "v.amp", "gain_mode": "mirror", "d_x_m_per_v": 2e-8, "theta_i_rad": 0.1, "lambda0_m": LAM,
                "polarity": -1.0, "modes": [{"f_hz": 14e3, "q": 8.0, "r": 1.0}], "discretisation": "tustin",
                "prewarp_hz": 14e3, "drive_min_v": -5.0, "drive_max_v": 5.0, "hysteresis_frac": 0.05,
                "segment": "fs1", "driver": "stretcher", "n_sub": 1, "sub": [0, 1], "ring_len": 0,
                "insertion_loss_db": None,
                "toggles": {"resonance": True, "hysteresis": True, "saturation": True, "insertion_loss": True}}},
            {"block": "eom_phase", "name": "eom", "params": {
                "in": "delay.dl.out", "gain_mode": "vpi", "v_pi_v": 4.0, "polarity": 1.0, "poles_hz": [200e3],
                "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -1.0, "drive_max_v": 1.0,
                "segment": "lo1", "driver": "eom", "n_sub": 2, "sub": [0, 2], "ring_len": 0,
                "insertion_loss_db": 3.0, "loss_applies_to": ["gbs", "bhd"], "loss_source": "UNVERIFIED test",
                "toggles": {"bandwidth": True, "saturation": True, "insertion_loss": True}}},
            {"block": "thermal_actuator", "name": "tec", "rate_div": 4, "params": {
                "in": "delay.dl.out", "segment": "arm2", "driver": "tec", "n_sub": 3, "sub": [0, 2], "ring_len": 0,
                "k_t_k_per_v": 0.5, "tau_th_s": 2e-4, "dt_min_k": -0.01, "dt_max_k": 0.01,
                "toggles": {"saturation": True}}},
            {"block": "aom_shift", "name": "aom", "params": {
                "in": "delay.dl.out", "order": 1, "input": "frequency", "k_hz_per_v": 1e3, "polarity": 1.0,
                "segment": "pr1", "driver": "eom", "n_sub": 1, "sub": [0, 1], "ring_len": 4,
                "insertion_loss_db": 1.0, "loss_applies_to": ["probe"], "loss_source": "UNVERIFIED test",
                "toggles": {"enabled": True, "insertion_loss": True}}},
        ]},
        {"block": "dither_modulator", "name": "dith", "params": {
            "dither_mode": "timeseries", "m_rad": 0.0785, "f_m_hz": 1.1e5, "phase_rad": 0.3, "segment": "arm1",
            "driver": "dither", "n_sub": 4, "sub": [0, 4], "ring_len": 0, "toggles": {"enabled": True}}},
        {"block": "ramp", "name": "ramp", "params": {"out": "v.ramp", "waveform": "triangle", "amplitude_v": 0.2,
                                                     "period_s": 1.3e-3, "offset_v": 0.01, "toggles": {"enabled": True}}},
        {"block": "adder", "name": "plant", "params": {
            "ins": ["act.str.dpsi", "act.mir.dpsi", "act.eom.dpsi", "act.tec.dT", "act.aom.dpsi", "act.dith.dpsi",
                    "v.ramp"], "gains": [1.0, 1.0, 0.5, 30.0, 0.01, 1.0, 0.5], "out": "psi.plant.x", "unit": "rad"}},
    ]


def assert_states_close(a, b, rtol=1e-12):
    np.testing.assert_array_equal(a.I, b.I)
    np.testing.assert_allclose(a.S, b.S, rtol=rtol, atol=rtol * max(1.0, float(np.max(np.abs(b.S)))))


def test_v20_compiled_equals_reference_all_blocks():
    n = PIECE_SAMPLES + 1500
    b1, s1, p1 = make(all_blocks_obj())
    b2, s2, p2 = make(all_blocks_obj())
    run_reference(b1, s1, p1, n)
    run_compiled(b2, s2, p2, n)
    for slot in ("act.str.dpsi", "act.mir.dpsi", "act.eom.dpsi", "act.tec.dT", "act.aom.dpsi", "v.dac", "v.amp"):
        assert abs(float(s1[slot])) > 0, slot  # every block is alive
    assert_states_close(s1, s2)
    p1.set({"pid.kp": -0.05, "str.hysteresis_frac": 0.2})
    p2.set({"pid.kp": -0.05, "str.hysteresis_frac": 0.2})
    run_reference(b1, s1, p1, 700)
    run_compiled(b2, s2, p2, 700)
    assert_states_close(s1, s2)
    assert [e.id for e in b2.losses()] == ["eom.il", "aom.il"]
    assert spec_from_obj(spec_to_obj(b1.spec)) == spec_from_obj(all_blocks_obj())  # exact spec round trip


def test_v1_static_gain_and_range():
    """u = +1 V → δψ = +6.491 rad (1e−9 vs the formula); DSC1 ±4 V clip → 51.93 rad span (T16 V1)."""
    k = 2 * math.pi * 16e-6 / (LAM * 10.0)
    assert k == pytest.approx(6.491, abs=5e-4)
    out = []
    for u in (1.0, 10.0, -10.0):
        obj = [const("src", "v.u", u),
               {"block": "dac", "name": "dac", "params": {"in": "v.u", "out": "v.dac", "device": "dsc1",
                                                          "noise_vrms": 0.0, "r_load_ohm": float("inf"),
                                                          "toggles": {"quantization": False, "noise": False,
                                                                      "saturation": True}}},
               stretcher("s", "v.dac", modes=[], toggles={"resonance": False, "hysteresis": False,
                                                          "saturation": True, "insertion_loss": False})]
        blk, st, pr = make(obj)
        run_reference(blk, st, pr, 1)
        out.append(float(st["act.s.dpsi"]))
    assert out[0] == pytest.approx(k, rel=1e-9)
    assert out[1] - out[2] == pytest.approx(8.0 * k, rel=1e-9) and 8.0 * k == pytest.approx(51.93, abs=0.01)


def test_v2_xi_eo():
    assert xi_eo(0.121, 0.270, 0.17, 1.458) == pytest.approx(0.7837, abs=5e-5)


def lockin(x, t, f):
    return 2.0 / x.size * np.sum(x * np.exp(-2j * np.pi * f * t))


def drive_sine(obj_act, f, n_settle, n_meas, fs=FS, amp=1e-3):
    obj = [{"block": "sine_injection", "name": "src", "params": {"in": None, "out": "v.u", "unit": "V",
                                                                 "freq_hz": f, "amplitude": amp,
                                                                 "toggles": {"enabled": True}}}, obj_act]
    blk, st, pr = make(obj, fs=fs)
    run_compiled(blk, st, pr, n_settle)
    t0 = st.n
    xs, us = np.empty(n_meas), np.empty(n_meas)
    out_slot = "act.a.dpsi" if obj_act["block"] != "hv_amp_914" else "v.out"
    for i in range(n_meas):
        run_compiled(blk, st, pr, 1)
        xs[i], us[i] = float(st[out_slot]), float(st["v.u"])
    t = (t0 + np.arange(n_meas)) / fs
    return lockin(xs, t, f) / lockin(us, t, f)


@pytest.mark.parametrize("disc, prewarp", [("tustin", 20e3), ("zoh", None)])
def test_v3_resonance_discretisation(disc, prewarp):
    """1 mode 20 kHz, Q = 10, f_sim 1 MHz: |H| = Q ± 1 %, arg H = −90° ± 0.5° (Tustin prewarped at f_m);
    ZOH adds the held-input half-sample lag −180°·f·h on top (same tolerances)."""
    obj = stretcher("a", "v.u", segment=None, modes=[{"f_hz": 20e3, "q": 10.0, "r": 1.0}], disc=disc,
                    prewarp=prewarp, gain_mode="direct", k_rad_per_v=1.0,
                    toggles={"resonance": True, "hysteresis": False, "saturation": False, "insertion_loss": False})
    h = drive_sine(obj, 20e3, 4000, 2000)  # 2000 samples = 40 periods
    lag = 0.0 if disc == "tustin" else -180.0 * 20e3 / FS
    assert abs(h) == pytest.approx(10.0, rel=0.01)
    assert np.degrees(np.angle(h)) == pytest.approx(-90.0 + lag, abs=0.5)


def test_v4_914_two_pole_and_slew():
    """Fig. 1 of the 914 (T16 §3.1.4): 2-pole 29.5/112 kHz within 5 % at six points; pure slew model
    (SR = 0.432 V/µs input-referred) within 10 % at 30, 50, 100 kHz."""
    table = {7.06e3: 51.7, 10e3: 50.1, 19.7e3: 42.3, 29.9e3: 36.2, 50e3: 23.7, 100e3: 11.0}

    def amp(poles, slew_on):
        return {"block": "hv_amp_914", "name": "a", "params": {
            "in": "v.u", "out": "v.out", "poles_hz": poles, "discretisation": "zoh", "prewarp_hz": None,
            "input_min_v": -5.0, "input_max_v": 5.0, "noise_v_per_rthz": 0.0, "slew_v_per_s": 4.32e5,
            "toggles": {"poles": True, "slew": slew_on, "saturation": True, "noise": False}}}

    for f, v in table.items():
        n_per = FS / f
        h = drive_sine(amp([29.5e3, 112e3], False), f, int(20 * n_per) + 200, int(round(10 * n_per)), amp=5.0)
        assert abs(h) * 52.1 == pytest.approx(v, rel=0.05), f
    for f, v in ((29.9e3, 36.2), (50e3, 23.7), (100e3, 11.0)):
        blk, st, pr = make([{"block": "sine_injection", "name": "src", "params": {
            "in": None, "out": "v.u", "unit": "V", "freq_hz": f, "amplitude": 5.0, "toggles": {"enabled": True}}},
            amp([], True)])
        run_compiled(blk, st, pr, 2000)
        ys = []
        for _ in range(int(4 * FS / f)):
            run_compiled(blk, st, pr, 1)
            ys.append(float(st["v.out"]))
        peak = 0.5 * (max(ys) - min(ys))
        assert peak * 52.1 / 5.0 == pytest.approx(v, rel=0.10), f


def stack_objs():
    """Disturbance → error → stacked c1 → c2 → c3 driving three static actuators, and the parallel
    form F_k = C_k⋯C_1 (copies), both with the real pid / actuator blocks (T16 §4.4, V13)."""
    gains = [dict(kp=0.5, ki_per_s=2e3), dict(kp=0.3, ki_per_s=5e2), dict(kp=0.2, ki_per_s=1e2)]
    src = toy("toy_noise_ar1", "e", out="v.e", unit="V", a=0.95, b=1.0)
    stages = [pid("c1", "v.e", **gains[0]), pid("c2", "ctrl.c1.u", **gains[1]), pid("c3", "ctrl.c2.u", **gains[2])]
    stages = [{**s, "params": {**s["params"], "u_min_v": -1e9, "u_max_v": 1e9}} for s in stages]
    acts = [gain_act(f"a{k}", f"ctrl.c{k}.u", float(k)) for k in (1, 2, 3)]
    stacked = [src, {"stack": {"stages": stages, "actuators": acts}}]
    branches = []
    for k in (1, 2, 3):
        chain, inp = [], "v.e"
        for j in range(1, k + 1):
            nm = f"b{k}c{j}"
            s = pid(nm, inp, **gains[j - 1])
            s["params"].update(u_min_v=-1e9, u_max_v=1e9)
            chain.append(s)
            inp = f"ctrl.{nm}.u"
        chain.append(gain_act(f"b{k}a", inp, float(k)))
        branches.append({"seq": chain, "name": f"br{k}"})
    return stacked, [src, {"parallel": branches}]


def test_v13_stacked_equals_parallel_real_blocks():
    stacked, parallel = stack_objs()
    b1, s1, p1 = make(stacked)
    b2, s2, p2 = make(parallel)
    run_compiled(b1, s1, p1, 3000)
    run_compiled(b2, s2, p2, 3000)
    tot1 = sum(float(s1[f"act.a{k}.dpsi"]) for k in (1, 2, 3))
    tot2 = sum(float(s2[f"act.b{k}a.dpsi"]) for k in (1, 2, 3))
    assert tot1 == pytest.approx(tot2, rel=1e-9)
    for k in (1, 2, 3):
        assert float(s1[f"act.a{k}.dpsi"]) == pytest.approx(float(s2[f"act.b{k}a.dpsi"]), rel=1e-9)


def offload_loop(sign_c2):
    """Two-stage stack on a constant disturbance d = 0.7 rad: e = d − φ_1 − φ_2, φ_k = u_k; fast c1
    (K_i = 2e4 s⁻¹) on e, integrating c2 (K_i = ±2e3 s⁻¹) on u_1."""
    obj = [const("d", "v.d", 0.7),
           {"block": "adder", "name": "err", "params": {"ins": ["v.d", "act.a1.dpsi", "act.a2.dpsi"],
                                                        "gains": [1.0, -1.0, -1.0], "out": "v.e", "unit": "rad"}},
           {"stack": {"stages": [pid("c1", "v.e", ki_per_s=2e4, u_min_v=-1e6, u_max_v=1e6),
                                 pid("c2", "ctrl.c1.u", ki_per_s=sign_c2 * 2e3, u_min_v=-1e6, u_max_v=1e6)],
                      "actuators": [{"seq": [dl("d1", "ctrl.c1.u"), gain_act("a1", "delay.d1.out", 1.0)]},
                                    {"seq": [dl("d2", "ctrl.c2.u"), gain_act("a2", "delay.d2.out", 1.0)]}]}}]
    blk, st, pr = make(obj)
    run_compiled(blk, st, pr, 20000)
    return st


def test_v13_offload_and_v14_sign():
    """Constant disturbance: with an integrating second stage of the right sign the fast actuator is
    unloaded (u_1 → 0) and the slow one carries d (T16 V13); with C_2's sign flipped it diverges (V14)."""
    st = offload_loop(+1.0)
    assert abs(float(st["ctrl.c1.u"])) < 1e-6
    assert float(st["act.a2.dpsi"]) == pytest.approx(0.7, rel=1e-6)
    st = offload_loop(-1.0)
    assert abs(float(st["ctrl.c1.u"])) > 1e3


def test_v16_fringe_hop_by_two_pi():
    """Fringe hop: the output jumps by exactly −n·u_2π; with a static actuator K the phase moves by
    −n·2π (error < 1e−6 rad) and the hop is counted."""
    k = 6.49
    u2pi = 2 * math.pi / k
    obj = [toy("toy_sine", "ramp", out="v.u", unit="V", amp=3.5, freq_hz=50.0),
           {"block": "adder", "name": "e", "params": {"ins": ["v.u"], "gains": [0.0], "out": "v.e", "unit": "V"}},
           {"block": "lock_supervisor", "name": "sup", "params": {
               "u_in": "v.u", "e_in": "v.e", "true_in": None, "out": "v.o", "relock": "fringe_hop",
               "u_min_v": -2.5, "u_max_v": 2.5, "t_sat_s": 0.0, "e_th": float("inf"), "t_e_s": 0.0,
               "t_hop_s": 1e-4, "u_2pi_v": u2pi, "toggles": {"relock": True}}},
           gain_act("a", "v.o", k)]
    blk, st, pr = make(obj)
    jumps = []
    prev_in, prev = None, None
    for _ in range(20000):
        run_compiled(blk, st, pr, 1)
        u_in, y = float(st["v.u"]), float(st["act.a.dpsi"])
        if prev is not None:
            dj = (y - prev) - k * (u_in - prev_in)
            if abs(dj) > 1.0:
                jumps.append(dj)
        prev_in, prev = u_in, y
    assert int(st["lock.sup.n_relock"]) == len(jumps) >= 2
    for dj in jumps:
        n = round(dj / (2 * math.pi))
        assert n != 0 and abs(dj - 2 * math.pi * n) < 1e-6
    assert abs(float(st["v.o"])) <= 2.5 + u2pi


def test_v17_aom():
    """δf_RF = 1 Hz (m = +1) → ψ = −2π t; a φ_RF step jumps ψ by −m φ_RF (conv §4, T16 V17)."""
    def aom(mode, order):
        key = "k_hz_per_v" if mode == "frequency" else "k_rad_per_v"
        return {"block": "aom_shift", "name": "a", "params": {
            "in": "v.c", "order": order, "input": mode, key: 1.0, "polarity": 1.0, "segment": None,
            "insertion_loss_db": None, "toggles": {"enabled": True, "insertion_loss": False}}}

    blk, st, pr = make([const("c", "v.c", 1.0), aom("frequency", 1)], fs=1e3)
    run_compiled(blk, st, pr, 1000)
    assert float(st["act.a.dpsi"]) == pytest.approx(-2 * math.pi * 1.0, rel=1e-12)  # t = 1000 · 1 ms
    blk, st, pr = make([const("c", "v.c", 0.4), aom("phase", -1)], fs=1e3)
    run_compiled(blk, st, pr, 1)
    assert float(st["act.a.dpsi"]) == pytest.approx(+0.4, rel=1e-12)  # −m φ_RF with m = −1


def test_v18_tec_step_and_offload():
    """TEC step: ΔT(t) = K_T u (1 − e^{−t/τ}) (exact ZOH), DC gain K_T; driver copied to the subsegments.
    PZT + TEC stack under a ramp drift: the PZT command returns to the centre (T16 V18)."""
    fs, tau, kt = 1e3, 0.05, 0.8
    obj = [const("c", "v.c", 1.0),
           {"block": "thermal_actuator", "name": "tec", "params": {
               "in": "v.c", "segment": "spool", "driver": "tec", "n_sub": 5, "sub": [1, 4], "ring_len": 0,
               "k_t_k_per_v": kt, "tau_th_s": tau, "dt_min_k": -10.0, "dt_max_k": 10.0,
               "toggles": {"saturation": True}}}]
    blk, st, pr = make(obj, fs=fs)
    for n in (1, 50, 400):
        run_compiled(blk, st, pr, n - st.n)
        assert float(st["act.tec.dT"]) == pytest.approx(kt * (1 - math.exp(-n / fs / tau)), rel=1e-12)
    np.testing.assert_allclose(st["drv.spool.tec"], [0.0, *[float(st["act.tec.dT"])] * 3, 0.0], rtol=0)
    # PZT (fast, range-limited) + TEC (slow) stack on a ramp drift r [rad/s]
    r, c_t_l = 0.5, 30.0  # c_T·L_h [rad/K] (coefficient side, F6)
    obj = [{"block": "ramp", "name": "drift", "params": {"out": "v.d", "waveform": "sawtooth", "amplitude_v": r * 50,
                                                         "period_s": 100.0, "offset_v": r * 50,
                                                         "toggles": {"enabled": True}}},
           {"block": "adder", "name": "err", "params": {"ins": ["v.d", "act.pzt.dpsi", "act.tec.dT"],
                                                        "gains": [1.0, -1.0, -c_t_l], "out": "v.e", "unit": "rad"}},
           {"stack": {"stages": [pid("c1", "v.e", ki_per_s=-50.0, kp=0.0), pid("c2", "ctrl.c1.u", ki_per_s=-0.5, kp=0.0)],
                      "actuators": [{"seq": [dl("d1", "ctrl.c1.u"), gain_act("pzt", "delay.d1.out", -1.0)]},
                                    {"seq": [dl("d2", "ctrl.c2.u"), {"block": "thermal_actuator", "name": "tec", "params": {
                                        "in": "delay.d2.out", "segment": "spool", "driver": "tec", "n_sub": 1,
                                        "sub": [0, 1], "ring_len": 0, "k_t_k_per_v": 1.0, "tau_th_s": 1.0,
                                        "dt_min_k": -100.0, "dt_max_k": 100.0, "toggles": {"saturation": True}}}]}]}}]
    for s in obj[2]["stack"]["stages"]:
        s["params"].update(u_min_v=-1e6, u_max_v=1e6)
    blk, st, pr = make(obj, fs=fs)
    run_compiled(blk, st, pr, int(60 * fs))
    # stable for K_i2 = 0.5 (Routh: s³ + 51 s² + 50 s + 1500·K_i2); the PZT ends near r/(c_T L·K_i2) = 0.033 V
    assert abs(float(st["ctrl.c1.u"])) < 0.05 * r * 60  # PZT carries < 5 % of the accumulated drift
    assert float(st["act.tec.dT"]) * c_t_l == pytest.approx(r * 60, rel=0.05)


def test_f_sim_convergence_closed_loop():
    """Closed loop (PI + 914 + resonant stretcher), sine disturbance at 500 Hz: the residual amplitude at
    f_sim = 1 MHz and 2 MHz agrees within 5 % (PLAN §6.4, T16 V20)."""
    res = []
    for fs in (1e6, 2e6):
        obj = [{"block": "sine_injection", "name": "d", "params": {"in": None, "out": "v.d", "unit": "rad",
                                                                   "freq_hz": 500.0, "amplitude": 0.1,
                                                                   "toggles": {"enabled": True}}},
               {"block": "adder", "name": "err", "params": {"ins": ["v.d", "act.s.dpsi"], "gains": [1.0, 1.0],
                                                            "out": "v.e", "unit": "rad"}},
               pid("pid", "v.e", kp=-0.002, ki_per_s=-2000.0),
               dl("dl", "ctrl.pid.u", 5.7e-6),
               {"block": "hv_amp_914", "name": "amp", "params": {
                   "in": "delay.dl.out", "out": "v.amp", "poles_hz": [29.5e3, 112e3], "discretisation": "zoh",
                   "prewarp_hz": None, "input_min_v": -5.0, "input_max_v": 5.0, "noise_v_per_rthz": 0.0,
                   "slew_v_per_s": 1e9, "toggles": {"poles": True, "slew": False, "saturation": True, "noise": False}}},
               stretcher("s", "v.amp", segment=None, modes=[{"f_hz": 40e3, "q": 10.0, "r": 1.0}])]
        blk, st, pr = make(obj, fs=fs)
        run_compiled(blk, st, pr, int(5e-3 * fs))
        n = int(10e-3 * fs)
        e = np.empty(n)
        for i in range(n):
            run_compiled(blk, st, pr, 1)
            e[i] = float(st["v.e"])
        res.append(abs(lockin(e, np.arange(n) / fs, 500.0)))
    assert res[1] == pytest.approx(res[0], rel=0.05)
    assert res[0] < 0.05  # the loop (f_c ≈ 2 kHz) does suppress the disturbance
