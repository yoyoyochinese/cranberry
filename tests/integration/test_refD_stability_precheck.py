"""Phase 5 M1D-02: per-sweep-point linear-stability pre-check of refD (``chk.lockD_stability``,
``loop_stability_check``).

* the check sits in refD's ``offline.blocks``, so every sweep point evaluates it from its own
  ``offline.params``; unstable points get a ``loop_unstable`` CHK warning and are NOT pruned (D04: pruning
  is an open user decision);
* over the 915B-resonance grid of refD (f_res × Q, 16 points, other axes at refD_nominal) the unstable
  set is the M1 review's (``reviews/scripts/m1_refD_stability_grid.py``): f_res = 5 kHz at every Q and
  11.45 kHz at Q ≥ 31 (V_0 = 2.90 V); with T_opt = 1 (V_0 = 3.33 V) also 11.45 kHz at Q 9.65;
* the analytic verdict agrees with the time loop: an unstable point does not hold the lock, the
  nominal point does;
* the check's Nyquist logic on textbook loops (integrator + delay: K_u = π/(2τ));
* (Ref D prediction follow-up, reviews/phase5_refD_stability_mismatch.md) the check models the DSC1's
  discrete PID law, and agrees with the time loop at two points where the continuous law did not;
  the lock supervisor counts fringe slips about the lock point the loop holds (−π/2 here).
"""

import dataclasses
import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.chk.loop_stability import make_loop_stability_check
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import ensure_loaded, make_block
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parent / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

REFD = "presets/reference/refD.yaml"
F_RES = [round(x, 2) for x in np.geomspace(5.0e3, 6.0e4, 4)]
Q_RES = [round(x, 2) for x in np.geomspace(3.0, 100.0, 4)]


def _resonance_grid(**pinned_over):
    """(f_res, Q) → (stable, max |L| at −180°, CHK codes) over refD's resonance grid (other axes nominal)."""
    nominal = config.load(rs.NOMINAL).offline["params"]
    refd = config.load(REFD)
    axes = {}
    for path in config.sweep_markers(refd):
        if not path.startswith("offline.params."):
            continue
        name = path.split(".")[-1]
        lin = name[: -len("_db")] + "_lin" if name.endswith("_db") else name
        v = pinned_over.get(name, nominal[name] if name in nominal else 10 * math.log10(nominal[lin]))
        axes[path] = [v]
    axes["offline.params.f_res_hz"] = list(np.geomspace(5.0e3, 6.0e4, 4))  # the preset's own axis values
    axes["offline.params.q_res"] = list(np.geomspace(3.0, 100.0, 4))
    ensure_loaded()
    out = {}
    for pt in expand(refd, axes):
        prm = pt.spec.offline["params"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            (chk,) = [make_block(b) for b in offline_block_specs(pt.spec) if b.kind == "loop_stability_check"]
        codes = [e.code for e in chk.chk_log.events]
        assert sum(issubclass(x.category, ChkWarning) for x in w) == len(codes)
        out[(round(prm["f_res_hz"], 2), round(prm["q_res"], 2))] = (chk.result["stable"],
                                                                    chk.result["max_l_at_180"], codes)
    return out


@pytest.mark.parametrize("v0_case", ["nominal", "t_opt_1"])
def test_refD_resonance_grid_unstable_points(v0_case, record_property):
    grid = _resonance_grid(**({"coupler_excess_db": 0.0} if v0_case == "t_opt_1" else {}))
    assert len(grid) == 16  # every point kept (no pruning)
    unstable = sorted(k for k, (ok, _, _) in grid.items() if not ok)
    record_property("unstable_points_f_res_hz_q", unstable)
    want = [(F_RES[0], q) for q in Q_RES] + [(F_RES[1], q) for q in Q_RES[2:]]
    if v0_case == "t_opt_1":
        want.append((F_RES[1], Q_RES[1]))
    assert unstable == sorted(want)
    for k, (ok, mx, codes) in grid.items():
        assert codes == ([] if ok else ["loop_unstable"]), k
        assert (mx > 1.0) == (not ok), (k, mx)  # one dominant −180° crossing: |L| > 1 ⇔ unstable
    # refD_nominal's own point (20 kHz, Q 10) is stable: GM 12.8 dB with the DSC1's discrete PID law
    # (14.0 dB, max |L| 0.199, with the continuous K_p + K_i/s + K_d s of the M1 review)
    ensure_loaded()
    (chk,) = [make_block(b) for b in offline_block_specs(config.load(rs.NOMINAL))]
    assert chk.result["stable"] and chk.result["max_l_at_180"] == pytest.approx(0.229, abs=0.002)
    assert chk.chk_log.events == []


@pytest.mark.numba
def test_refD_precheck_agrees_with_time_loop():
    """Noise off: the nominal point holds +π/2 to < 1 µrad; an unstable grid point (5 kHz, Q 9.65,
    max |L| at −180° = 3.7) does not hold the lock (oscillation of tenths of a rad or fringe slips)."""
    tails = {}
    for label, prm in {"nominal": {}, "unstable": {"f_res_hz": 5.0e3, "q_res": 9.65}}.items():
        spec = rs.load(params=prm, run={"duration_s": 0.03}, toggles=rs.NOISE_OFF,
                       measure={"record": {"slots": ["lock.lockD.err"], "decim": 1}})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            e = run(spec).series["lock.lockD.err"]
        flagged = any("loop_unstable" in str(x.message) for x in w if issubclass(x.category, ChkWarning))
        assert flagged == (label == "unstable")
        tails[label] = e[-10000:]
    assert np.ptp(tails["nominal"]) < 1e-6
    assert np.ptp(tails["unstable"]) > 0.1


def _textbook(k, tau, **kw):
    p = {"k0": 1.0, "kp": 0.0, "ki_per_s": k, "kd_s": 0.0, "delay_s": tau, "poles_hz": [math.inf], "modes": [],
         "f_grid_hz": [1.0e-2, 1.0e8, 200001], "toggles": {"enabled": True}, **kw}
    return make_loop_stability_check(p, name="t")


def test_nyquist_logic_textbook_loops():
    tau = 1.0e-5
    ku = math.pi / (2 * tau)  # L = K e^{−sτ}/s: |L| = 1 at the −180° frequency 1/(4τ) when K = π/(2τ)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ChkWarning)
        assert _textbook(0.95 * ku, tau).result["stable"]
        bad = _textbook(1.05 * ku, tau)
        assert not bad.result["stable"] and bad.chk_log.events[0].code == "loop_unstable"
        assert bad.result["f_180_hz"] == pytest.approx(1 / (4 * tau), rel=1e-3)
        assert bad.result["max_l_at_180"] == pytest.approx(1.05, rel=1e-3)
        # P only (no integrator, n_int = 0): K e^{−sτ}/(1 + s/ω_p) unstable for K > √(1 + (ω_180/ω_p)²)
        p_only = dict(ki_per_s=0.0, poles_hz=[1.0e4], delay_s=tau)
        assert _textbook(0.0, tau, kp=0.9, **p_only).result["stable"]
        assert not _textbook(0.0, tau, kp=5.0, **p_only).result["stable"]
        # grid that does not reach |L| ≪ 1 (derivative term without poles): inconclusive, warned
        inc = _textbook(1.0, 0.0, kd_s=1.0, poles_hz=[], f_grid_hz=[1.0, 1.0e3, 1001])
        assert inc.chk_log.events[0].code == "loop_check_grid"
    assert _textbook(1.0, tau, toggles={"enabled": False}).result is None
    with pytest.raises(KeyError):
        _textbook(1.0, tau, discretisation="dsc1")  # h_s required
    with pytest.raises(ValueError):
        _textbook(1.0, tau, discretisation="zoh")


# Ref D blind prediction points (predictions/refD/results/runs.json, LHS seed 20261101) where the continuous
# check disagreed with the time loop (reviews/phase5_refD_stability_mismatch.md)
LHS_002 = {  # 30.8 kHz, Q 66.5 near the PID's I/D notch: continuous check "stable" (0.538), time loop oscillates
    "dL_m": 13.083376963983474, "kappa1": 0.5087424772903861, "kappa2": 0.5096350222052648,
    "coupler_excess_db": 0.047060158254620295, "il_915b_db": 0.5326794061792505, "h0_hz2_per_hz": 26056.762395588285,
    "h_m1_hz2": 415864319.7121703, "k_T_rad_per_m_K": 44.755189267388914, "r_T_k_per_s": 0.0003259494080580608,
    "tau_th_s": 1.5490118533138657, "h_ref_rad2_hz": 0.0010729462965244955, "mech_gamma": 1.0,
    "mech_f_c_hz": 167.4048261741738, "mech_rho_0": 0.3875538140414559, "D_pol_rad2_per_s": 0.005869462104284686,
    "polarizer_er_db": 26.47340813075201, "K_T_rad_per_m_K": -0.6790631330567094, "dpsi_b0_rad": 2.3897185276183115,
    "theta_coupler_rad": 0.028908451953870467, "responsivity_a_per_w": 1.015166000143706,
    "eps_r": -0.00879056815085185, "i_n_a_rthz": 4.6968587806065973e-11, "v_off0_v": 0.007821598844350939,
    "adc_noise_vrms": 0.00035610346319716793, "latency_s": 5.7e-06, "in_out_poles": False,
    "r_914_ohm": 190943.64569792442, "anti_windup": "clamp", "f_res_hz": 30826.634451786762,
    "q_res": 66.52630362089371, "loop_sign": -1.0, "coupler_sign": -1.0}
LHS_044 = {  # 53.9 kHz, Q 81.6: continuous check "unstable" (1.132), time loop locks (at −π/2)
    "dL_m": 1.9315176745679858, "kappa1": 0.49107831434524213, "kappa2": 0.5071392648927708,
    "coupler_excess_db": 0.012238854796321523, "il_915b_db": 0.4081586743331438, "h0_hz2_per_hz": 14453.418190799126,
    "h_m1_hz2": 352463203.2346414, "k_T_rad_per_m_K": 46.06042837720213, "r_T_k_per_s": 0.00027640482970346685,
    "tau_th_s": 99.30864511715767, "h_ref_rad2_hz": 0.008477982966639697, "mech_gamma": 2.0,
    "mech_f_c_hz": 126.04694718435022, "mech_rho_0": 0.3629566352520164, "D_pol_rad2_per_s": 1.371493993036598e-07,
    "polarizer_er_db": 26.190052237961112, "K_T_rad_per_m_K": -1.9526473778397215, "dpsi_b0_rad": 1.4721000095525423,
    "theta_coupler_rad": 0.02573682749615287, "responsivity_a_per_w": 1.0174381172287847,
    "eps_r": -0.021516280580504837, "i_n_a_rthz": 3.124515516046214e-11, "v_off0_v": -0.00545987761219649,
    "adc_noise_vrms": 0.0003823419576098157, "latency_s": 5.7e-06, "in_out_poles": True,
    "r_914_ohm": 181157.77555648226, "anti_windup": "clamp", "f_res_hz": 53931.99359972174,
    "q_res": 81.61485388934297, "loop_sign": 1.0, "coupler_sign": -1.0}


def _point(values, **tree):
    (pt,) = expand(config.load({"extends": REFD, **tree}), {f"offline.params.{k}": [v] for k, v in values.items()})
    return pt.spec


def _chk(spec, **over):
    ensure_loaded()
    (b,) = [b for b in offline_block_specs(spec) if b.kind == "loop_stability_check"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ChkWarning)
        return make_loop_stability_check({**b.params, **over}, name="t").result


def test_check_controller_is_the_dsc1_discrete_law():
    """The check's ``discretisation: dsc1`` controller equals the time loop's PID (pid.freq_response, DSC1
    p.19) up to the Tustin warp of the integrator, (ωh)²/12 of the I term (< 3 % of |C| to 60 kHz, largest at
    the shallow notch of the discrete law near 28 kHz); the
    continuous K_p + K_i/s + K_d s is off by up to ×6 near its I/D notch √(K_i/K_d)/2π = 28 kHz."""
    from cranberry.blocks.chk.loop_stability import _ctrl

    spec = rs.load()
    leaves, params = rs.built_leaves(spec)
    f = np.geomspace(1.0e2, 6.0e4, 200)
    c_pid = -leaves["dsc1.pid"].freq_response(f, params)  # canonical gains carry the DSC1 sign (e = S − v)
    ensure_loaded()
    (b,) = [b for b in offline_block_specs(spec) if b.kind == "loop_stability_check"]
    s = 2j * np.pi * f
    assert b.params["discretisation"] == "dsc1"
    assert np.max(np.abs(_ctrl(s, b.params) / c_pid - 1)) < 0.03
    cont = _ctrl(s, {**b.params, "discretisation": "continuous"})
    assert np.max(np.abs(c_pid / cont)) > 5.0  # the notch the device does not have
    f_notch = math.sqrt(b.params["ki_per_s"] / b.params["kd_s"]) / (2 * math.pi)
    assert f_notch == pytest.approx(28.3e3, rel=0.01)


@pytest.mark.numba
@pytest.mark.parametrize("label, values, stable", [("lhs_002", LHS_002, False), ("lhs_044", LHS_044, True)])
def test_check_agrees_with_time_loop_at_prediction_mismatch_points(label, values, stable):
    """Regression (Ref D prediction §6.3): with the DSC1's discrete PID law the check's verdict is the
    time loop's; with the continuous law (the frozen prediction's) it was the opposite. Noise and both
    quantisers off: lhs_002 grows into a ≈ 1.3 rad rms oscillation at f_res, lhs_044 settles to < 1 µrad."""
    spec = _point(values, run={"duration_s": 0.03},
                  measure={"record": {"slots": ["lock.lockD.err"], "decim": 1}})
    assert _chk(spec)["stable"] is stable
    assert _chk(spec, discretisation="continuous")["stable"] is (not stable)
    spec = dataclasses.replace(spec, toggles={**spec.toggles, **rs.NOISE_OFF})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tail = run(spec).series["lock.lockD.err"][-10000:]
    if stable:
        assert np.ptp(tail) < 1e-6
    else:
        assert np.std(tail) > 1.0


@pytest.mark.numba
def test_refD_slips_counted_about_the_held_lock_point():
    """Regression (Ref D prediction §10-2): with loop_sign × coupler_sign = −1 the loop holds −π/2
    (lock.lockD.err ≈ ±π); noise makes err cross ±π, which the old setpoint-centred count took for slips.
    The supervisor's slip_center_rad (π here) puts the boundaries on the unstable points: 0 slips."""
    spec = rs.load(params={"coupler_sign": -1.0}, run={"duration_s": 0.05},
                   measure={"record": {"slots": ["lock.lockD.err"], "decim": 1}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run(spec, keep_state=True)
    e = np.asarray(res.series["lock.lockD.err"])
    e = e[len(e) // 2:]
    assert abs(abs(np.mean(e)) - math.pi) < 0.05 and np.std(e) < 0.1  # held at −π/2
    old = np.floor(e / (2 * math.pi) + 0.5)
    assert np.count_nonzero(np.diff(old)) > 10  # what the setpoint-centred count saw
    assert int(res.final_state["lock.dsc1.sup.n_slip"]) == 0
