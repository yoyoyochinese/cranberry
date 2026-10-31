"""WP-I4: numerical convergence of the refB preset (PLAN §6.4, prompt §8; ``runner.convergence``).

* chunk length: bit-identical recorded series and final physical state for two chunk lengths (G11);
* f_sim × 2 (2 → 4 MHz; the FPGA-stand-in ADC / PID / DAC keep their 2 MHz device rate through
  ``measure.convergence.keep_rate``, D64-3, and τ_L − 1/f_sim is re-evaluated, D55-1):
  (a) a seed-free loop observable — a 1 kHz sine on the probe phase (PZT2 chain), every noise source off — the
      lock-1 error and Δ amplitudes change by ≤ 1 % (tolerance 3 %, criterion tol/3);
  (b) the closed-loop θ̃ at the nominal point changes by ≤ 1/3 of VT-17's 0.1° (0.5 s records);
* seed ensemble (10 seeds): spread of θ̃ ≪ the 0.1° tolerance of VT-17 (a conclusion must exceed it, prompt §8).
"""

import dataclasses
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.convergence import convergence
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refB_support", Path(__file__).parent / "refB_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

REC = {"record": {"slots": ["lock.lock1.err", "lock.lock2.err", "metric.mi.port0"], "decim": 1}}


def test_refB_chunk_invariance():
    spec = rs.load(run={"duration_s": 0.05, "chunk_samples": 30011}, measure=REC)
    a = run(spec)
    b = run(dataclasses.replace(spec, run=dataclasses.replace(spec.run, chunk_samples=7777)))
    for k in a.series:
        assert a.series[k].tobytes() == b.series[k].tobytes(), k
    names = [n for n in a.final_state.layout.slot_names() if not n.startswith(("rec.", "rng."))]
    assert len(names) > 50
    for n in names:
        assert np.asarray(a.final_state[n]).tobytes() == np.asarray(b.final_state[n]).tobytes(), n


def test_refB_f_sim_doubling_loop_observable():
    amp = 0.02   # rad on the probe phase (PZT2 chain, pure gain)
    pipe = [{"name": "dither", "params": {"freq_hz": 1.0e3, "amplitude": amp / 30.0}}]
    spec = rs.load(run={"duration_s": 0.03}, toggles={**rs.NOISE_OFF, "dither.enabled": True}, pipeline=pipe,
                   measure={"record": {"slots": ["lock.lock1.err", "metric.mi.port0"], "decim": 1,
                                       "antialias": "none"}})
    q = {"lock.lock1.err.std": 0.03 * 0.02, "i.port0.std_rad": 0.03 * 1e-3}
    rep = convergence(spec, ("f_sim",), quantities=q, f_sim_factors=(2,))
    vals = rep.checks["f_sim"]["values"]
    assert rep.checks["f_sim"]["passed"], vals
    ref = run(spec).summary()
    assert ref["lock.lock1.err.std"] > 1e-3 and ref["i.port0.std_rad"] > 1e-5   # the sine is really in both


def test_refB_f_sim_doubling_theta():
    spec = rs.load(run={"duration_s": 0.5})
    rep = convergence(spec, ("f_sim",), quantities={"i.port0.rms_rad": math.radians(0.1)}, f_sim_factors=(2,))
    assert rep.checks["f_sim"]["passed"], rep.checks["f_sim"]["values"]


def test_refB_seed_ensemble(record_property):
    spec = rs.load(run={"duration_s": 0.15})
    rep = convergence(spec, ("seed",), quantities={"i.port0.rms_rad": math.radians(0.1)}, seeds=list(range(1, 11)))
    seed = rep.checks["seed"]
    assert seed["passed"]
    v = seed["values"]["i.port0.rms_rad"]
    record_property("theta_deg_mean", math.degrees(v["mean"]))
    record_property("theta_deg_spread", math.degrees(v["std"]))
    assert math.degrees(v["mean"]) == pytest.approx(1.5, abs=0.1)
    assert math.degrees(v["std"]) < 0.1 / 1.5   # ≪ VT-17's 0.1°
