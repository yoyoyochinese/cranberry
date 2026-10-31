"""VT-19: dither contribution (S §4 VT-19 ← T04 V5, T14 V8, T16 V15; prompt §7-10).

WP-B3 part (this file so far): with the Ref C LO dither m = 0.0785 rad (±4.5°) as a time series at the
loop rate, (a) the measured-quadrature angle Δ(t) has rms m/√2 = 0.0555 rad, and (b) the square-law
quantum-noise detector of ``quantum_noise_power`` averages to P̄ − ΔP·χ·cos 2Δ0 with χ = ⟨cos 2δ⟩ =
J0(2m) = 0.99384 (conv §6), i.e. the block's noise power responds to the dither exactly as the phase →
squeezing map assumes. Tolerances: 1e−6 (deterministic, integer number of dither periods).

Still to add by their owners: the Ref C source mapping −5.63 → −5.18 dB and the time-domain PSD check
(±0.05 dB) via ``phase_to_quadrature`` (WP-F4) and the Ref C preset (WP-I3).
"""

import math

import numpy as np
import pytest
from scipy.special import j0

import cranberry.core.testing  # noqa: F401
from cranberry.blocks.err.testing import REF_C, qnl, qnl_flat, record
from cranberry.core.testing import toy

pytestmark = pytest.mark.numba


def test_vt19_dither_rms_and_chi_from_the_qnl_detector():
    m, f_m, fs = REF_C["m"], 1e6, 16e6
    n = 16 * 5000  # 5000 dither periods
    objs = [toy("toy_sine", "dith", out="beam.lo.phi", unit="rad", amp=m, freq_hz=f_m), qnl(rf=True)]
    rec = record(objs, fs, n, ["beam.lo.phi", "v.qnl"], external={"beam.p.phi": "rad"}, init={"beam.p.phi": -math.pi})
    delta = rec["beam.lo.phi"]  # Δ = φ_LO (φ_p = −π, φ_OPA = 0)
    assert math.sqrt(np.mean(delta**2)) == pytest.approx(m / math.sqrt(2), rel=1e-6)
    assert m / math.sqrt(2) == pytest.approx(0.0555, abs=5e-5)
    q = qnl_flat()
    chi = (q["p_mean"] - np.mean(rec["v.qnl"])) / q["p_delta"]  # Δ0 = 0: ⟨P⟩ = P̄ − ΔP χ
    assert chi == pytest.approx(j0(2 * m), rel=1e-9)
    assert j0(2 * m) == pytest.approx(0.99384, abs=1e-5)
