"""Phase 5 M2 Ref C review fixes (reviews/phase5_M2_refC.md, D69): M2C-01 (zero span vs record anti-aliasing),
M2C-02 (stability of the three-stage stack), M2C-03 (actuator clamp fraction as validity).

* M2C-01: the zero span reads cos 2δ / sin 2δ formed at the full rate and then decimated (``Recorder``
  ``cos_sin2``). With the time-series dither (1 MHz, m/√2 = 0.0555 rad, all in δ) the 30 MHz squeezed level is
  the same for ``antialias`` none / boxcar / fir within 0.02 dB and equals the full-rate moments (swept-spectrum
  bin); before the fix boxcar gave −5.626 dB against −5.164 dB (0.46 dB too good).
* M2C-02: ``chk.lockC_stability`` (``loop_stability_check`` with ``stages``) over the six stack axes of refC
  (f_u1 × f_12 × f_23 × mirror_f_res × mirror_q × latency = 972 points, others at refC_nominal) flags exactly
  the review's 84 unstable points (closed-loop roots with a (6,6) Padé delay, ``m2_refC_stack_stability.py``):
  all at f_12 = 3 kHz with f_u1 ≥ 9.65 kHz, every Q ≥ 9.5 point at f_u1 = 30 kHz, none with Q = 3 below 30 kHz. Every
  point is kept (no pruning) and the nominal point is stable. In the time loop an unstable point oscillates at
  tenths of a rad while the nominal point holds; its CHK warning reaches ``RunResult.events``.
* M2C-03: ``measure.validity`` reports the clamp fraction of the three controller outputs (``RunResult.validity``,
  ``summary()['validity.<slot>.frac']``, ``valid``): 0 at the nominal point; with the T08 placeholder fibre
  noise (pre-M2C-04 nominal: ρ = 0, f_lo = 0.1 Hz, f_u1 = 10 kHz) and h_ref_conduit = 40 the EOM sits at its
  clamp ≈ 10.6 % of the time (review) and the non-dither residual is 82 mrad; a ±5 V EOM drive range
  (``eom_range_v`` axis) removes the clamp and the residual returns to the linear ≈ 14 mrad.
"""

import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.core.checks import ChkWarning
from cranberry.core.registry import ensure_loaded, make_block
from cranberry.runner.run import offline_block_specs, run
from cranberry.runner.sweep import expand

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parent / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

M = math.pi * 0.1 / 4.0


def _db(x):
    return 10 * math.log10(float(x))


@pytest.mark.numba
@pytest.mark.integration
def test_m2c01_zero_span_independent_of_record_antialias():
    out = {}
    for aa in ("none", "boxcar", "fir"):
        spec = rs.load(params={"x": rs.X_T04}, toggles=rs.NOISE_OFF, overlays=(rs.TIMESERIES,), esa=True,
                       measure={"record": {"slots": ["lock.lockC.err"], "decim": 25, "antialias": aa}})
        res = run(spec)
        sp = res.measurements["esa"]["spectrum_sqz"]
        i30 = int(np.argmin(np.abs(np.asarray(sp["f_hz"]) - 3e7)))
        out[aa] = (_db(res.measurements["esa"]["sqz_30m"]["r_lin"]), _db(np.atleast_1d(sp["r_lin"])[i30]))
    ref = out["none"][1]                                  # full-rate moments (swept-spectrum bin at 30 MHz)
    assert ref == pytest.approx(-5.164, abs=0.01)         # the dither is in the zero span (−5.63 dB without it)
    for aa, (zs, swept) in out.items():
        assert swept == ref, aa
        assert zs == pytest.approx(ref, abs=0.02), aa     # review: boxcar −5.626 dB before the fix


STACK = ("f_u1_hz", "f_12_hz", "f_23_hz", "mirror_f_res_hz", "mirror_q", "latency_s")


def _stack_grid():
    """{(f_u1, f_12, f_23, f_res, Q, τ_L): (stable, CHK codes)} over refC's stack axes, other axes nominal."""
    nominal = config.load(rs.NOMINAL).offline["params"]
    refc = config.load(rs.PRESET)
    axes = {path: [nominal[path.split(".")[-1]]] for path in config.sweep_markers(refc)
            if path.startswith("offline.params.") and path.split(".")[-1] not in STACK}
    ensure_loaded()
    out = {}
    for pt in expand(refc, axes):
        prm = pt.spec.offline["params"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            (chk,) = [make_block(b) for b in offline_block_specs(pt.spec) if b.kind == "loop_stability_check"]
        codes = [e.code for e in chk.chk_log.events]
        assert sum(issubclass(x.category, ChkWarning) for x in w) == len(codes)
        out[tuple(prm[k] for k in STACK)] = (chk.result["stable"], codes)
    return out


@pytest.mark.integration
def test_m2c02_stack_stability_grid(record_property):
    grid = _stack_grid()
    assert len(grid) == 4 * 3 * 3 * 3 * 3 * 3            # every sweep point kept (no pruning, D34)
    unstable = sorted(k for k, (ok, _) in grid.items() if not ok)
    record_property("n_unstable", len(unstable))
    assert len(unstable) == 84                            # the review's closed-loop roots (m2_refC_stack_stability.py)
    for k, (ok, codes) in grid.items():
        assert codes == ([] if ok else ["loop_unstable"]), k
    f_u1 = sorted({k[0] for k in grid})
    q = sorted({k[4] for k in grid})
    assert all(k[1] == pytest.approx(3.0e3) and k[0] >= f_u1[2] for k in unstable)   # f_12 = 3 kHz, f_u1 ≥ 9.65 kHz
    assert {k for k in grid if k[1] == pytest.approx(3.0e3) and k[0] == f_u1[-1] and k[4] >= q[1]} <= set(unstable)
    assert not any(k[4] == q[0] and k[0] < f_u1[-1] for k in unstable)
    # refC_nominal's own point
    ensure_loaded()
    (chk,) = [make_block(b) for b in offline_block_specs(config.load(rs.NOMINAL)) if b.kind == "loop_stability_check"]
    assert chk.result["stable"] and chk.chk_log.events == []


@pytest.mark.numba
@pytest.mark.integration
def test_m2c02_unstable_point_oscillates_in_time_loop():
    """Review point f_u1 30 kHz, f_12 3 kHz, f_res 10 kHz, Q 30 (time loop 0.5 rad there): warned and not locked."""
    stds = {}
    for label, prm in {"nominal": {}, "unstable": {"f_u1_hz": 3.0e4, "f_12_hz": 3.0e3, "mirror_f_res_hz": 1.0e4,
                                                   "mirror_q": 30.0}}.items():
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = run(rs.load(params=prm, run={"duration_s": 0.06}))
        flagged = any("loop_unstable" in str(x.message) for x in w if issubclass(x.category, ChkWarning))
        assert flagged == (label == "unstable")
        assert [e["code"] for e in res.events] == (["loop_unstable"] if label == "unstable" else [])
        stds[label] = rs.at(res.metrics.noise_levels["i.port0"], "std_loop_rad", 3e7)
    assert stds["nominal"] < 0.03
    assert stds["unstable"] > 0.1


@pytest.mark.numba
@pytest.mark.integration
def test_m2c03_clamp_fraction_is_reported():
    old = {"rho_conduit": 0.0, "mech_f_lo_hz": 0.1, "f_u1_hz": 1.0e4, "h_ref_conduit_rad2_hz": 40.0}
    out = {}
    for label, prm in {"nominal": {}, "clamped": old, "wide_eom": dict(old, eom_range_v=5.0)}.items():
        res = run(rs.load(params=prm, run={"duration_s": 0.3}))
        tot = rs.at(res.metrics.noise_levels["i.port0"], "sigma_total_rad", 3e7)
        out[label] = (res, math.sqrt(max(tot**2 - M**2 / 2, 0.0)))
    slots = ["ctrl.pid1.u", "ctrl.pid2.u", "ctrl.pid3.u"]
    for label, (res, _) in out.items():
        assert sorted(res.validity) == slots
        assert all(res.summary()[f"validity.{s}.frac"] == res.validity[s] for s in slots)
    nom, clamped, wide = (out[k][0] for k in ("nominal", "clamped", "wide_eom"))
    assert nom.valid and all(v == 0.0 for v in nom.validity.values())
    assert not clamped.valid
    assert clamped.validity["ctrl.pid1.u"] == pytest.approx(0.106, abs=0.02)        # review: 10.6 %
    assert clamped.validity["ctrl.pid2.u"] == 0.0 and clamped.validity["ctrl.pid3.u"] == 0.0
    assert wide.valid
    assert out["clamped"][1] > 0.07 and out["wide_eom"][1] < 0.02                   # 82 vs ≈ 14 mrad (linear)

