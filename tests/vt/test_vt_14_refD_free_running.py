"""VT-14 (WP-I2): feedback off ⇒ the Ref D phase returns to its free-running noise (prompt §7-7; T15 §7
test 5; notes/00_summary.md §4 VT-14).

``refD_nominal`` with the DSC1 PID's ``feedback`` toggle off (u = clip(O) = 0: the 915B is static).
Sources on: DX1 white FM (h0 = 2.546e4 Hz²/Hz), RIN, mechanical noise on both arms (independent,
T08 eq. 2.7), thermal drift, SMF polarization drift, PDB offset. Checks:

(a) The Welch PSD of the true φ_D (``lock.lockD.err``, ≥ 50 segments, Hann, linear detrend for the
    thermal ramp) equals the sum of the injected disturbances
        S_laser + S_env = 4 sin²(πfτ)·h0/f² + Σ_arms (L/L_ref)^γ h_ref /f_lo² /((1 + (f/f_lo)²)(1 + (f/f_c)²)),
    τ = n_g ΔL / c (T15 eq. 8, T08 eq. 2.7; S_act = 0 with the loop open), within ±3/√N_avg
    (relative) in each of five bands from 200 Hz to 200 kHz (T15: per bin; here per band average,
    which is the stricter statistic).
(b) The detector output v (the DSC1 input, whose in-loop error is e_DSC1 = S − v) equals the nonlinear
    T15 eq. (3) evaluated on the recorded beam phases, (x, Δχ), δp and V_off, to 1e−12 V (sensor noise
    off for this check).
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parents[1] / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

FS = 1.0e6
DUR = 1.2  # s
NPERSEG = 16384
BANDS = [(200.0, 500.0), (500.0, 2e3), (2e3, 1e4), (1e4, 5e4), (5e4, 2e5)]


@pytest.fixture(scope="module")
def free_run():
    slots = ["lock.lockD.err", "beam.lo.phi", "beam.pr.phi", "pol.per.x", "pol.per.dchi", "v.power.dp",
             "v.pdb.voff", "v.pdb"]
    spec = rs.load(run={"duration_s": DUR}, measure={"record": {"slots": slots, "decim": 1}},
                   toggles={"dsc1.pid.feedback": False, "pdb450c.shot_noise": False,
                            "pdb450c.electronic_noise": False})
    return spec, run(spec)


def _s_expected(spec, f):
    pr = spec.offline["params"]
    a1, a2 = rs.segment(spec, "a1"), rs.segment(spec, "a2")
    tau = pr["n_g"] * (a1["L_m"] - a2["L_m"]) / rs.C_M_PER_S
    s_laser = 4 * np.sin(np.pi * f * tau) ** 2 * pr["h0_hz2_per_hz"] / f**2
    mech = rs.block_params(spec, "env.mech")
    two_pole = mech["h_ref_rad2_hz"] / mech["f_lo_hz"] ** 2 / ((1 + (f / mech["f_lo_hz"]) ** 2) * (1 + (f / mech["f_c_hz"]) ** 2))
    s_env = sum((seg["L_m"] / mech["L_ref_m"]) ** mech["gamma"] for seg in (a1, a2)) * two_pole
    return s_laser + s_env, s_laser, s_env


def test_vt14_free_running_psd_equals_injected_sum(free_run):
    spec, res = free_run
    assert int(res.final_state["lock.dsc1.sup.n_loss"]) == 0  # nothing reported as "lock" logic with the loop open
    x = res.series["lock.lockD.err"]
    f, p = signal.welch(x, fs=FS, window="hann", nperseg=NPERSEG, noverlap=NPERSEG // 2, detrend="linear")
    n_avg = (x.size - NPERSEG) // (NPERSEG // 2) + 1
    assert n_avg >= 50
    tol = 3.0 / math.sqrt(n_avg)
    s_tot, s_las, s_env = _s_expected(spec, f[1:])
    p = p[1:]
    fr = f[1:]
    for lo, hi in BANDS:
        m = (fr >= lo) & (fr < hi)
        ratio = p[m].mean() / s_tot[m].mean()
        assert abs(ratio - 1) < tol, (lo, hi, ratio, tol)
    # both contributions are exercised: mechanical noise dominates at 200 Hz, laser × ΔL above 2 kHz
    m_lo, m_hi = (fr >= 200) & (fr < 500), (fr >= 1e4) & (fr < 5e4)
    assert s_env[m_lo].mean() > 0.5 * s_las[m_lo].mean() and s_las[m_hi].mean() > 100 * s_env[m_hi].mean()


def test_vt14_detector_output_follows_eq3(free_run):
    spec, res = free_run
    s = res.series
    p = rs.block_params(spec, "pdb450c")
    fr = rs.fringe(spec)
    phid = s["beam.lo.phi"] - s["beam.pr.phi"]
    x, dchi = s["pol.per.x"], s["pol.per.dchi"]
    mu = p["gamma_coh"] * ((1 - x) + x * np.exp(1j * dchi))
    p1 = p["kappa1"] * p["t_arm1"] * p["p_in_w"]
    p2 = (1 - p["kappa1"]) * p["t_arm2"] * p["p_in_w"]
    eps = fr["sign"] * (2 * p["kappa2"] - 1) * (p1 - p2) / fr["pdet"] + p["eps_r"]
    v = fr["v0"] * (1 + s["v.power.dp"]) * (eps + fr["sign"] * fr["vis"] * np.real(mu * np.exp(1j * phid))) \
        + p["v_off_v"] + s["v.pdb.voff"]
    assert np.max(np.abs(np.clip(v, -p["v_sat_v"], p["v_sat_v"]) - s["v.pdb"])) < 1e-12
    # φ_D wanders by tenths of a rad around the fringe top (u = 0): the cos is far from linear there
    assert np.ptp(phid) > 0.2 and np.max(np.abs(phid)) < 1.0
