"""opo_cw_spectrum block (WP-B6): contract, D55 frame, toggles, gain / coefficient tables (T02 §2, §8)."""

import numpy as np
import pytest

from cranberry.blocks.err import amplitude_gains_from_table, amplitude_gains_opo
from cranberry.blocks.source.opo_cw import make_opo_cw_spectrum, probe_phase_coef, squeezing_angle_coef
from cranberry.blocks.source.opo_testing import REF_B, REF_C, opo
from cranberry.core.block import MissingParameterError
from cranberry.offline import gaussian as g
from cranberry.path.coef import CoefKey
from cranberry.path.topology import as_topology
import cranberry.path.testing as pt


def test_contract_and_validation():
    b = opo(REF_B)
    assert b.losses() == [] and b.bands() == []  # intracavity loss is inside cov (via ρ)
    so = b.source_output
    assert so.basis == "cw_sideband" and so.cov.shape == (1, 2, 2)
    assert so.provenance["x"] == pytest.approx(np.sqrt(0.1 / 0.18))
    with pytest.raises(MissingParameterError):
        opo({k: v for k, v in REF_B.items() if k != "p_th_w"})
    with pytest.raises(MissingParameterError):  # every toggle must be given (D34)
        make_opo_cw_spectrum({**REF_B, "toggles": {"detuning": False}})
    with pytest.raises(ValueError, match="below threshold"):
        opo(REF_B, pump_w=0.18)
    with pytest.raises(ValueError, match="cavity_model"):
        opo(REF_B, cavity_model="bowtie")
    with pytest.raises(ValueError, match="loss_ic"):
        opo(REF_B, loss_ic=-0.1).source_output


@pytest.mark.parametrize("phi0", [0.0, 0.7, -2.0])
def test_d55_frame(phi0):
    """θ_min = θ_sq,0 = φ_OPA,0/2 + π/2 at Δ = 0 (θ_off = 0); rotation of the φ_OPA,0 = 0 result."""
    b = opo(REF_B, phi_opa0_rad=phi0, f_s_hz=[0.0, 1e6, 1e7])
    _, _, th = b.source_output.project(None)
    np.testing.assert_allclose(th, 0.0, atol=1e-13)
    ref = opo(REF_B, f_s_hz=[0.0, 1e6, 1e7]).source_output.cov
    np.testing.assert_allclose(b.source_output.cov, g.rotate(ref, phi0 / 2), atol=1e-14)


def test_toggles_off():
    fs = [0.0, 1e6, 2e7]
    # frequency_dependence off: every f_s equals the Ω = 0 value
    off = opo(REF_B, f_s_hz=fs, toggles={"frequency_dependence": False}).source_output.cov
    np.testing.assert_allclose(off, np.broadcast_to(off[0], off.shape), rtol=1e-15)
    on = opo(REF_B, f_s_hz=fs).source_output.cov
    np.testing.assert_allclose(on[0], off[0], rtol=1e-15)
    assert not np.allclose(on[2], off[2])
    # pump_dependent_loss off: L_ic(0)
    b = opo(REF_B, toggles={"pump_dependent_loss": False})
    assert b.loss_ic(b.x) == pytest.approx(0.00249)
    np.testing.assert_array_equal(b.source_output.cov, opo(REF_B, loss_ic=0.00249).source_output.cov)
    # detuning off: Δ0 ignored
    np.testing.assert_array_equal(opo(REF_B, detuning_rad_s=1e6).source_output.cov, opo(REF_B).source_output.cov)
    assert abs(float(opo(REF_B, detuning_rad_s=1e6, toggles={"detuning": True}).source_output.project(None)[2][0])) > 1e-3
    # x = 0: vacuum whatever else is on
    np.testing.assert_allclose(opo(REF_B, pump_w=0.0, detuning_rad_s=1e6, toggles={"detuning": True}).source_output.cov,
                               g.vacuum(1)[None], atol=1e-15)


def test_linewidth_model():
    b = opo(REF_C)
    kt, kl = b.rates(0.3)
    assert kt + kl == pytest.approx(np.pi * 1.75e9) and kt / (kt + kl) == pytest.approx(0.99)
    assert b.source_output.provenance["fwhm_hz"] == pytest.approx(1.75e9)
    p = opo(REF_B)
    assert p.kappa(p.x) == pytest.approx(299_792_458.0 * (0.123 + 0.00249 + 0.00222 * p.x) / (2 * 0.5))


def test_gain_table():
    b = opo(REF_B, gain_table_grid=[0.0, 0.05, 0.1, 0.15])
    t = b.source_output.gain_table
    assert t.axes == ("pump_w",) and t.value_shape == (2,)
    assert amplitude_gains_from_table(t, pump_w=0.1) == pytest.approx(amplitude_gains_opo(np.sqrt(0.1 / 0.18)))
    tx = opo(REF_C, gain_table_grid=[0.1, 0.5]).source_output.gain_table
    assert tx.axes == ("x",) and amplitude_gains_from_table(tx, x=0.5) == pytest.approx((2.0, 1 / 1.5))
    with pytest.raises(ValueError, match="below threshold"):
        opo(REF_B, gain_table_grid=[0.1, 0.2]).source_output


def test_coef_tables_and_topology_entries():
    b = opo(REF_B, coef_table_x=[0.5, 0.7, 0.8])
    tabs = b.source_output.coef_tables
    for x in (0.5, 0.7, 0.8):
        k = b.kappa(x)
        assert float(tabs["opo.detuning.squeezed"].interp(x=x)) == pytest.approx(squeezing_angle_coef(x, k), rel=1e-12)
        assert float(tabs["opo.detuning.probe"].interp(x=x)) == pytest.approx(probe_phase_coef(k), rel=1e-12)
    single = opo(REF_B)
    seg = pt.seg_spec("opo", 0.5, kind="opo", drivers={"detuning": "uniform"})
    beams = [pt.beam_spec(n, [pt.step_spec("opo")], "bhd") for n in ("sq", "pr")]
    topo = as_topology(pt.topo_spec([seg], beams, ["bhd"], single.coef_entries(["sq"], ["pr"])))
    k = single.kappa(single.x)
    assert topo.coef_value(CoefKey("opo", "detuning", "sq", "phase")) == pytest.approx(squeezing_angle_coef(single.x, k))
    assert topo.coef_value(CoefKey("opo", "detuning", "pr", "phase")) == pytest.approx(-0.5 / k)
    assert topo.coef_value(CoefKey("opo", "detuning", "pr", "group_delay")) == 0.0
