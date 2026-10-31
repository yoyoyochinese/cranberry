"""WP-I2: numerical convergence of the refD preset (PLAN §6.4, prompt §8; ``runner.convergence``).

* chunk length: bit-identical series and final state for two chunk lengths (G11);
* f_sim × 2 (1 → 2 MHz, the DSC1 blocks keep their device rates: ADC 1 MHz, PID / DAC 500 kHz, and the
  loop delay τ_L − 1/f_sim − explicit delays is re-evaluated): the residual phase variance in the band
  100 Hz–100 kHz and the lock point change by ≤ 1/3 of a 10 % tolerance (the in-loop band above the
  old Nyquist frequency is the BandLedger's business, D57-6: compared in a common band);
* seed ensemble (10 seeds): mean and spread of the residual (a conclusion must exceed the spread).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from cranberry.runner.convergence import convergence
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parent / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

REC = {"record": {"slots": ["lock.lockD.err", "v.dsc1.cv"], "decim": 1}}


def test_refD_chunk_invariance():
    """Bit-identical series and physics state. ``runner.convergence('chunk')`` also compares the
    DrawBuffer internals ``rng.<block>.buf/.idx``, whose fill level depends on where the chunk ends
    (the random stream itself is chunk-invariant): reported to WP-F4; compared here without them."""
    import dataclasses

    spec = rs.load(run={"duration_s": 0.03, "chunk_samples": 7000}, measure=REC)
    a = run(spec)
    b = run(dataclasses.replace(spec, run=dataclasses.replace(spec.run, chunk_samples=2334)))
    for k in a.series:
        assert a.series[k].tobytes() == b.series[k].tobytes(), k
    names = [n for n in a.final_state.layout.slot_names() if not n.startswith(("rec.", "rng."))]
    assert len(names) > 50
    for n in names:
        assert np.asarray(a.final_state[n]).tobytes() == np.asarray(b.final_state[n]).tobytes(), n


def test_refD_seed_ensemble(record_property):
    """10 seeds (``run.seed``; G11 streams per block): ``runner.convergence('seed')`` runs the ensemble;
    the residual is taken after the lock acquisition (last 20 ms), so it is computed here per seed."""
    import dataclasses

    spec = rs.load(run={"duration_s": 0.04}, measure=REC)
    rep = convergence(spec, ["seed"], seeds=list(range(1, 11)))
    assert rep.checks["seed"]["passed"]
    tails = []
    for seed in range(1, 11):
        e = run(dataclasses.replace(spec, run=dataclasses.replace(spec.run, seed=seed))).series["lock.lockD.err"]
        tails.append(float(np.std(e[-20000:])))
    m, s = float(np.mean(tails)), float(np.std(tails, ddof=1))
    record_property("residual_std_mean_rad", m)
    record_property("residual_std_spread_rad", s)
    assert len(set(tails)) == 10  # independent realisations
    assert 0.0 < s < 0.1 * m  # 20 ms of laser-dominated residual: defined to a few %


def _band_variance(spec, f_lo=100.0, f_hi=1.0e5):
    res = run(spec)
    fs = spec.run.f_sim_hz
    x = res.series["lock.lockD.err"][int(0.02 * fs):]
    f, p = signal.welch(x, fs=fs, window="hann", nperseg=int(8192 * fs / 1e6), detrend="linear")
    m = (f >= f_lo) & (f <= f_hi)
    return float(np.trapezoid(p[m], f[m])), float(np.mean(x))


def test_refD_f_sim_doubling():
    dur = 0.32
    v1, m1 = _band_variance(rs.load(run={"duration_s": dur}, measure=REC))
    blocks2 = [{"name": "dsc1.adc", "rate_div": 2}, {"name": "dsc1.pid", "rate_div": 4},
               {"name": "dsc1.sup", "rate_div": 4}, {"name": "dsc1.dac", "rate_div": 4}]
    spec2 = rs.load(run={"duration_s": dur, "f_sim_hz": 2.0e6}, blocks=blocks2, measure=REC)
    assert rs.block_params(spec2, "dsc1.latency")["tau_s"] == pytest.approx(5.7e-6 - 0.5e-6 - 1.5e-6, abs=1e-15)
    v2, m2 = _band_variance(spec2)
    tol = 0.10
    assert abs(np.sqrt(v2 / v1) - 1) <= tol / 3, (v1, v2)
    assert abs(m2 - m1) <= tol / 3 * np.sqrt(v1)
