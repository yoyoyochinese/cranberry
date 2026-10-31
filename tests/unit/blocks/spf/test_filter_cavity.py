"""Unit tests of ``filter_cavity`` (WP-B2): contract (D34, toggles, losses), static detuning phase and
transmission, laser-phase low pass (delay τ_eff), compiled = reference. T05 V11/V12 numbers: VT-34."""

import math

import numpy as np
import pytest

from cranberry.blocks.err.testing import record, rig
from cranberry.blocks.spf.filter_cavity import ase_suppression, phase_transfer, tau_eff_s, transmission
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference

pytestmark = pytest.mark.numba
EXT = {"v.det": "Hz", "laser.dx1.phi": "rad"}


def fc(lock=True, lp=True, **over):
    p = {"segment": "fc", "beam": "lo", "dnu_fc_hz": 1e5, "t0": 0.95, "source": "test value",
         "in_detuning": "v.det", "in_phase": "laser.dx1.phi", "toggles": {"lock_residual": lock, "laser_lowpass": lp}}
    p.update(over)
    return {"block": "filter_cavity", "name": "fc1", "params": p}


def test_missing_input_raises():
    obj = fc()
    del obj["params"]["in_detuning"]
    with pytest.raises(MissingParameterError):
        rig([obj], 1e6, external=EXT)
    rig([fc(lock=False)], 1e6, external=EXT)  # not needed when its toggle is off


def test_static_detuning_and_toggles_off():
    root, st, pa = rig([fc(lp=False)], 1e6, external=EXT)
    st["v.det"] = 0.1 * 1e5
    run_reference(root, st, pa, 1)
    assert float(st["psi.fc.lo"]) == pytest.approx(math.atan(0.2), rel=1e-14)
    assert float(st["psi.fc.lo.trans"]) == pytest.approx(0.95**2 / 1.04, rel=1e-14)
    root, st, pa = rig([fc(lock=False, lp=False)], 1e6, external=EXT)
    st["v.det"], st["laser.dx1.phi"] = 1e4, 3.0
    run_reference(root, st, pa, 5)
    assert float(st["psi.fc.lo"]) == 0.0 and float(st["psi.fc.lo.trans"]) == 0.95**2


def test_losses_register_t0():
    root, _, _ = rig([fc()], 1e6, external=EXT)
    (el,) = root.losses()
    assert el.eta == pytest.approx(0.95**2) and el.segment == "fc" and el.applies_to == {"lo"}


def test_lowpass_delay_is_tau_eff():
    """A laser-phase ramp φ_L = s·t: the LO carries H_FC φ_L, i.e. ψ → −s·τ_eff in steady state."""
    fs, n, slope = 1e7, 4000, 1e3
    rec = record([fc(lock=False)], fs, n, ["psi.fc.lo"], external=EXT, init={})
    assert np.all(rec["psi.fc.lo"] == 0.0)  # φ_L ≡ 0 → ψ ≡ 0
    root, st, pa = rig([fc(lock=False)], fs, external=EXT)
    for k in range(n):
        st["laser.dx1.phi"] = slope * k / fs
        run_reference(root, st, pa, 1)
    assert float(st["psi.fc.lo"]) == pytest.approx(-slope * tau_eff_s(1e5), rel=1e-4)


def test_compiled_equals_reference_and_helpers():
    a = rig([fc()], 1e6, external=EXT)
    b = rig([fc()], 1e6, external=EXT)
    for (root, st, pa), runner in ((a, run_reference), (b, run_compiled)):
        for k in range(50):
            st["v.det"], st["laser.dx1.phi"] = 1e3 * math.sin(k), 0.01 * k * k
            runner(root, st, pa, 1)
    np.testing.assert_allclose(b[1].S, a[1].S, rtol=1e-12, atol=1e-300)
    assert abs(phase_transfer(5e4, 1e5)) ** 2 == pytest.approx(0.5)
    assert np.angle(transmission(1e4, 1e5)) == pytest.approx(math.atan(0.2))
    assert ase_suppression(1e9, 1e9 / 100, 0.5e9) < 1.0  # partial period path
