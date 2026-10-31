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
* the check's Nyquist logic on textbook loops (integrator + delay: K_u = π/(2τ)).
"""

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
    # refD_nominal's own point (20 kHz, Q 10) is stable with GM 14 dB (M1 review)
    ensure_loaded()
    (chk,) = [make_block(b) for b in offline_block_specs(config.load(rs.NOMINAL))]
    assert chk.result["stable"] and chk.result["max_l_at_180"] == pytest.approx(0.199, abs=0.002)
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
