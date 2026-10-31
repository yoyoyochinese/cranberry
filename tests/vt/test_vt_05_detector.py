"""VT-05 detector parts (WP-B5): circuit noise = pure loss, QE = pure loss, GBS input vs BHD separation.

Sources: T13 §6 tests 1–2, T12 §9 V8 and V13, conv §2 and §5, prompt §7-3; 00_summary VT-05
("回路雑音 = 損失（C = 28 dB で 0.158%）"). Tolerance: analytic identities 1e−12 relative (PLAN §6.3).
The elements are registered by ``detector_efficiency`` and applied through the same LossRegistry →
``offline.gaussian.pure_loss`` path as every other loss (prompt §7-3 implementation).
"""

import numpy as np
import pytest

import cranberry.blocks.det.testing as T
from cranberry.blocks.det import bhd
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline import gaussian as g
from cranberry.offline.loss_registry import LossElement, LossRegistry

pytestmark = pytest.mark.vt

REL = 1e-12
ONLY = {"qe": False, "circuit_noise": False, "imbalance_leak": False, "freq_response": False, "dip": False}


def registry(*blocks):
    r = LossRegistry(paths={"sqz": {"segments": ["fiber"], "observers": ["bhd", "gbs"]}})
    r.add(LossElement("fiber.il", 0.95, "fiber", {"gbs", "bhd"}, "[test]", True, segment="fiber"))
    for b in blocks:
        r.add_all(b.losses())
    return r


def mk(**over):
    return make_block(BlockSpec.from_obj(T.det_eff(**over)))


@pytest.mark.parametrize(("c_db", "loss", "tol"), [(28.0, 0.00158, 5e-6), (21.7, 0.0067, 5e-5)])
def test_circuit_noise_equivalent_loss(c_db, loss, tol):
    """1/(1 + C) = 0.158 % at 28 dB (Ref A '<0.2 %'), 0.67 % at 21.7 dB (Ref B) (T13 test 1, T12 V8)."""
    d = mk(clearance_lin=10 ** (c_db / 10), toggles={**ONLY, "circuit_noise": True})
    (el,) = [e for e in d.losses() if e.kind == "electronic"]
    assert 1 - el.eta == pytest.approx(loss, abs=tol)  # quoted values, last digit
    r = registry(d)
    eta_e = r.eta("sqz", "bhd") / 0.95
    assert eta_e == pytest.approx(1 / (1 + 10 ** (-c_db / 10)), rel=1e-12)
    # unsubtracted R_meas = (C R + 1)/(C + 1) equals the pure loss η_e for every input state
    c = 10 ** (c_db / 10)
    for cov in (g.squeezed_vacuum(1.2, 0.0), np.eye(2) * 3.1):
        theta = np.linspace(0, np.pi, 13)
        r_in = g.homodyne_variance(cov, theta) / 0.5
        r_out = g.homodyne_variance(g.pure_loss(cov, eta_e), theta) / 0.5
        np.testing.assert_allclose(r_out, bhd.r_meas_electronic(r_in, c), rtol=REL)
        fe = make_block(BlockSpec.from_obj(T.frontend(clearance_lin=c, toggles={
            "imbalance_leak": False, "freq_response": False, "dip": False})))
        np.testing.assert_allclose(fe.psd(0.0, r_in)["r_meas"], r_out, rtol=REL)


def test_qe_is_pure_loss_prompt_7_3():
    """η_q = 0.99: V_out = η_q V_in + (1 − η_q)/2 through the registry (T13 test 2)."""
    d = mk(qe=0.99, toggles={**ONLY, "qe": True})
    eta = registry(d).eta("sqz", "bhd") / 0.95
    assert eta == pytest.approx(0.99, rel=REL)
    for cov in (g.squeezed_vacuum(1.7, 0.3), np.eye(2) * 2.5):
        theta = np.linspace(0, np.pi, 11)
        v_in = g.homodyne_variance(cov, theta)
        np.testing.assert_allclose(g.homodyne_variance(g.pure_loss(cov, eta), theta), 0.99 * v_in + 0.01 * 0.5,
                                   rtol=REL)


def test_gbs_vs_bhd_separation_t12_v13():
    """BHD-only detector elements switched off ⇒ the BHD η equals the GBS-input η (T12 V13, D46)."""
    on = registry(mk(name="det", toggles={"freq_response": False, "dip": False}),
                  make_block(BlockSpec.from_obj(T.overlap(name="lo_mm"))))
    assert on.eta("sqz", "gbs") == pytest.approx(0.95, rel=REL)
    assert on.eta("sqz", "bhd") < on.eta("sqz", "gbs")
    off = registry(mk(name="det", on=False), make_block(BlockSpec.from_obj(T.overlap(name="lo_mm", on=False))))
    assert off.eta("sqz", "bhd") == off.eta("sqz", "gbs") == pytest.approx(0.95, rel=REL)
