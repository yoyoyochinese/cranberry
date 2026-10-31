"""VT-10 (prompt §7-5, T06 V4): photon number conserved with Raman + self-steepening.

N = 2 soliton, α = 0, 3 z0, 800 fixed steps, Blow–Wood (f_R = 0.18) and Lin–Agrawal (f_R = 0.245).
Tolerance: |ΔN_ph/N_ph| ≤ 1e−5 (prototype 8e−7 / 7e−7); sign tests ΔE < 0 (energy to phonons) and
spectral-centroid shift < 0 (Stokes side, conv §8 sign; prototype −6.9 / −8.8 THz).
Also checks that with self-steepening off the energy is conserved instead (T06 §2.4).
"""

import numpy as np
import pytest

from cranberry.blocks.fib.gnlse import GnlseModel, pulse_field, solve_pulsed, time_axis

B2, GAMMA, T0 = -21.7e-27, 1.3e-3, 100e-15
LD = T0**2 / abs(B2)
P1 = abs(B2) / (GAMMA * T0**2)
Z0 = np.pi / 2 * LD


def _run(model, shock):
    t, dt = time_axis(4096, 80 * T0)
    fr = {"blow_wood": 0.18, "lin_agrawal": 0.245}[model]
    m = GnlseModel(4096, dt, lambda0_m=1550e-9, alpha_per_m=0.0, betas=(B2,), gamma_per_w_m=GAMMA, c_kerr=1 - fr,
                   c_raman=fr, raman_model=model, self_steepening=shock)
    return m, solve_pulsed(m, pulse_field("sech", 4 * P1, T0, t), t, 3 * Z0, n_steps=800, n_save=1,
                           control="fixed", rtol=None)


@pytest.mark.parametrize("model", ["blow_wood", "lin_agrawal"])
def test_photon_number_with_raman_and_shock(model):
    m, r = _run(model, True)
    n0, n1 = r.photon_number
    e0, e1 = r.energy_j
    assert abs(n1 - n0) / n0 <= 1e-5
    assert e1 < e0 and (e0 - e1) / e0 > 1e-2
    shift = m.centroid(r.spectrum(-1)) - m.centroid(r.spectrum(0))
    assert shift < -2 * np.pi * 1e12  # several THz to the red


def test_energy_conserved_without_shock():
    _, r = _run("blow_wood", False)
    e0, e1 = r.energy_j
    assert abs(e1 - e0) / e0 <= 1e-5
