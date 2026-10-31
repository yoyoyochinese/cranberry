"""VT-05 (basic parts, WP-F3 / M0): prompt §7-3 pure loss and §7-4 vacuum stays vacuum.

Sources: prompt §7-3, §7-4; T12 §9 V1–V3, V8, V13; T07 §8 V1, V1b, V5, V11; PLAN §6.2 VT-05, §6.3
(analytic identities: relative 1e−12). The OPO / pulsed-OPO parts of VT-05 belong to B6 / B11; the
canonical file location is tests/vt/test_vt_05_pure_loss.py (PLAN §6.1) — kept here because WP-F3 owns
tests/unit/offline only; proposed move in the WP-F3 report.
"""

import numpy as np
import pytest

from cranberry.core.rng import stream
from cranberry.core.units import db_to_lin, il_db_to_eta
from cranberry.offline import gaussian as g
from cranberry.offline.loss_registry import LossElement, LossRegistry

pytestmark = pytest.mark.vt

REL = 1e-12


@pytest.mark.parametrize("eta", [0.0, 0.3, 0.92, 1.0])
def test_pure_loss_exact_prompt_7_3(eta):
    """V_out = η V_in + (1 − η) V_vac exactly, for squeezed, anti-squeezed and thermal inputs (T12 V1)."""
    r, nbar = 1.7, 3.2
    states = [
        g.squeezed_vacuum(r, 0.0),
        g.squeezed_vacuum(r, 0.9),
        (2 * nbar + 1) / 2 * np.eye(2),
    ]
    theta = np.linspace(0, np.pi, 25)
    for cov in states:
        v_in = g.homodyne_variance(cov, theta)
        v_out = g.homodyne_variance(g.pure_loss(cov, eta), theta)
        np.testing.assert_allclose(v_out, eta * v_in + (1 - eta) * 0.5, rtol=REL, atol=0)
    # V− = e^{−2r}/2 and V+ = e^{2r}/2 explicitly
    vm, vp, _ = g.noise_levels_from_cov(g.pure_loss(g.squeezed_vacuum(r, 0.0), eta))
    assert vm == pytest.approx(eta * np.exp(-2 * r) / 2 + (1 - eta) / 2, rel=REL)
    assert vp == pytest.approx(eta * np.exp(2 * r) / 2 + (1 - eta) / 2, rel=REL)


def test_series_losses_multiply_and_commute_with_rotation():
    """T12 V2: 10 random losses in series (rotations in between) = one loss Π η_i."""
    rng = stream(5, "vt05.series")
    etas = rng.uniform(0.5, 1.0, 10)
    rots = rng.uniform(-1, 1, 10)
    cov = g.squeezed_vacuum(2.0, 0.3)
    seq = cov
    for e, psi in zip(etas, rots, strict=True):
        seq = g.rotate(g.pure_loss(seq, e), psi)
    one = g.rotate(g.pure_loss(cov, float(np.prod(etas))), float(np.sum(rots)))
    np.testing.assert_allclose(seq, one, rtol=REL, atol=1e-15)


def test_distributed_fiber_loss_equals_lumped():
    """T12 V3: 0.5 dB/km × 633 m split into M = 1, 10, 1000 pieces → η = 10^(−0.3165/10)."""
    eta_total = 10 ** (-0.5 * 633 / 1e4)
    assert eta_total == pytest.approx(0.9297, abs=5e-5)
    cov = g.squeezed_vacuum(1.5, 0.0)
    lumped = g.pure_loss(cov, eta_total)
    for m in (1, 10, 1000):
        piece = float(il_db_to_eta(0.5 * 633 / 1e3 / m))
        out = cov
        for _ in range(m):
            out = g.pure_loss(out, piece)
        np.testing.assert_allclose(out, lumped, rtol=1e-12)
        # the same through the registry: M fiber elements on one path
        reg = LossRegistry(paths={"sqz": {"segments": ["fib"], "observers": []}})
        for k in range(m):
            reg.add(
                LossElement.from_il_db(
                    f"fib.part{k}", 0.5 * 633 / 1e3 / m, kind="fiber", applies_to={"gbs"}, source="[T12 V3]",
                    enabled=True, segment="fib",
                )
            )
        assert reg.eta("sqz", "gbs") == pytest.approx(eta_total, rel=1e-12)


def test_vacuum_stays_vacuum_prompt_7_4():
    """No nonlinearity, no loss (passive gates only) — and loss alone — keep vacuum = I/2 (T07 V1, V1b)."""
    rng = stream(6, "vt05.vacuum")
    for n in (1, 2, 4):
        cov = g.vacuum(n)
        for _ in range(20):
            if n > 1 and rng.uniform() < 0.5:
                i, j = rng.choice(n, 2, replace=False)
                cov = g.beam_splitter(cov, rng.uniform(), int(i), int(j))
            else:
                cov = g.rotate(cov, rng.uniform(-np.pi, np.pi), mode=int(rng.integers(n)))
            cov = g.pure_loss(cov, rng.uniform(0, 1, n))
        assert np.max(np.abs(cov - g.vacuum(n))) < 1e-12
        assert g.is_physical(cov)


def test_lossy_squeezed_state_stays_physical():
    """T07 V3: physicality V + iΩ/2 ≥ 0 after loss and added noise."""
    for eta in (0.0, 0.01, 0.5, 0.999, 1.0):
        cov = g.add_noise(g.pure_loss(g.squeezed_vacuum(3.0, 0.4), eta), 0.0)
        assert g.is_physical(cov)


def test_circuit_noise_is_pure_loss():
    """T12 V8: unsubtracted white circuit noise at clearance C = 28 dB equals a pure loss
    η_el = 1/(1 + 10^(−2.8)); L_el = 0.158 %. R_meas = (R + E)/(1 + E)."""
    e_el = 10 ** (-2.8)
    eta_el = 1 / (1 + e_el)
    assert 1 - eta_el == pytest.approx(0.00158, abs=5e-6)
    cov = g.squeezed_vacuum(2.4, 0.0)
    theta = np.linspace(0, np.pi, 11)
    r_in = g.homodyne_variance(cov, theta) / 0.5
    r_meas = (r_in + e_el) / (1 + e_el)
    r_loss = g.homodyne_variance(g.pure_loss(cov, eta_el), theta) / 0.5
    np.testing.assert_allclose(r_meas, r_loss, rtol=1e-9)
    # conv §5 form: η_el = C/(1 + C)
    c_lin = float(db_to_lin(28.0))
    assert c_lin / (1 + c_lin) == pytest.approx(eta_el, rel=1e-12)


def test_added_noise_phase_insensitive():
    """T07 V11: added N photons/mode raises Var(x_θ) by N at every θ (ΔV/V_SN = 2N)."""
    cov = g.pure_loss(g.squeezed_vacuum(1.0, 0.2), 0.9)
    theta = np.linspace(0, np.pi, 19)
    out = g.add_noise(cov, 1e-3)
    np.testing.assert_allclose(g.homodyne_variance(out, theta) - g.homodyne_variance(cov, theta), 1e-3, rtol=1e-9)


def test_registry_eta_applied_once_equals_elementwise():
    """The registry's total η applied once by pure_loss = each element applied in turn (G06 + conv §2)."""
    etas = [0.98, 0.99, 0.99, 0.998, 0.96]
    reg = LossRegistry(paths={"sqz": {"segments": ["s"], "observers": ["bhd1"]}})
    for k, e in enumerate(etas):
        reg.add(
            LossElement(
                id=f"e{k}", eta=e, kind="component", applies_to={"bhd"}, source="[Ref A]", enabled=True, segment="s"
            )
        )
    cov = g.squeezed_vacuum(2.408, 0.0)
    seq = cov
    for e in etas:
        seq = g.pure_loss(seq, e)
    np.testing.assert_allclose(g.pure_loss(cov, reg.eta("sqz", "bhd")), seq, rtol=1e-12)
