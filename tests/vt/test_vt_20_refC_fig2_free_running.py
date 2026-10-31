"""VT-20 (d) (Phase 5 M2C-04, D69-1): the free-running fibre phase of refC against [Ref C Fig. 2 right].

Fig. 2 right is the only measurement of the fibre phase noise in Ref C: a 30 MHz zero span through the 1 km
fibres with the lock OFF over 1 s. The trace stays near the squeezing level for the whole second, with ms-long
excursions to +2.5…+7 dB (review §1 M2C-04: |δ| ≲ 0.1–0.2 rad most of the time, bursts of 0.4–0.6 rad).
Qualitative-but-numeric criteria on the differential phase δ(t) (all three PIDs open, the last 1 s of 1.2 s):

1. time near the squeezing angle: fraction of the second with |δ − median δ| < 0.2 rad ≥ 0.9;
2. rms of δ in 1–300 Hz (what a VBW-300 Hz trace shows within 1 s) ≤ 0.2 rad;
3. 1 s random-walk excursion: peak-to-peak of δ smoothed over 10 ms (the ms bursts removed) ≤ 0.4 rad
   (|δ| ≲ 0.2 rad on either side).

``refC_nominal`` (rho_conduit = 0.99 below f_rho, mech_f_lo_hz = 100 Hz, chosen by these criteria and
documented there) passes. The T08 placeholders it replaced (uncorrelated fibres, f⁻² from 0.1 Hz) fail all
three by far (review: 5.5 rad random walk, 6 % of the time near squeezing), and so do the neighbouring choices
that justify the nominal point: mech_f_lo_hz = 10 Hz fails (3) and rho_conduit = 0 fails (1). The criteria
constrain only the low-frequency noise; the kHz band (mech_f_c_hz, h_ref) that sets the locked residual is
not fitted to anything — in particular not to the paper's 0.1 rad (D69-1).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parents[1] / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

OPEN = {"pid1.feedback": False, "pid2.feedback": False, "pid3.feedback": False}
DECIM = 40
F_SIM = 1.0e6                        # the free-running phase is slow; 1 MHz suffices (f_N ≫ f_c of the fibre noise)


def fig2_stats(params):
    """(time fraction near squeezing, rms 1–300 Hz [rad], 10 ms-smoothed 1 s peak-to-peak [rad]) of δ."""
    spec = rs.load(params=params, run={"duration_s": 1.2, "f_sim_hz": F_SIM}, toggles=OPEN,
                   measure={"record": {"slots": ["lock.lockC.err"], "decim": DECIM, "antialias": "boxcar"}})
    err = np.asarray(run(spec).series["lock.lockC.err"], dtype=float)
    fr = F_SIM / DECIM
    d = (np.unwrap(2 * err) / 2)[int(0.2 * fr):]          # quadrature angle is mod π; last 1 s
    w = np.hanning(d.size)
    psd = 2 * np.abs(np.fft.rfft((d - d.mean()) * w)) ** 2 / (fr * np.sum(w**2))
    f = np.fft.rfftfreq(d.size, 1 / fr)
    band = (f >= 1.0) & (f <= 300.0)
    rms_band = math.sqrt(np.sum(psd[band]) * (f[1] - f[0]))
    near = float(np.mean(np.abs(d - np.median(d)) < 0.2))
    k = int(0.01 * fr)
    slow = np.convolve(d, np.ones(k) / k, mode="valid")
    return near, rms_band, float(np.ptp(slow))


def test_vt20d_nominal_free_running_phase_matches_fig2():
    near, rms_band, excursion = fig2_stats({})
    assert near >= 0.9
    assert rms_band <= 0.2
    assert excursion <= 0.4


def test_vt20d_t08_placeholders_violate_fig2():
    near, rms_band, excursion = fig2_stats({"rho_conduit": 0.0, "mech_f_lo_hz": 0.1})
    assert near < 0.2 and rms_band > 0.5 and excursion > 2.0


@pytest.mark.parametrize("params, failing", [({"mech_f_lo_hz": 10.0}, "excursion"), ({"rho_conduit": 0.0}, "near")])
def test_vt20d_neighbouring_points_fail(params, failing):
    near, rms_band, excursion = fig2_stats(params)
    passed = {"near": near >= 0.9, "rms": rms_band <= 0.2, "excursion": excursion <= 0.4}
    assert not passed[failing], (near, rms_band, excursion)
