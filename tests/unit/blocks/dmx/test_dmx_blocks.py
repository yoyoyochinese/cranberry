"""DMX blocks (WP-B8): contracts (D34 required params, toggles off ⇒ identity, losses(), bands()), the
design formulas, compiled == reference (1e−12) and the Liu helpers."""

import math
import warnings

import numpy as np
import pytest

import cranberry.blocks.dmx.testing as dt
from cranberry.blocks.dmx.demux import stage_counts, transfer_matrix
from cranberry.blocks.dmx.liu import actuator_sub, liu_loop_delay_s, probe_in_switch_window
from cranberry.core.block import MissingParameterError
from cranberry.core.checks import ChkWarning
from cranberry.core.compile import run_reference
from cranberry.path.topology import C_M_PER_S


def blk(obj, fs=1e6):
    return dt.make(obj, fs)


def test_design_requires_every_choice():
    obj = dt.design()
    del obj["params"]["layout"]
    with pytest.raises(MissingParameterError):
        blk(obj)
    obj = dt.design()
    del obj["params"]["toggles"]["crosstalk"]
    with pytest.raises(MissingParameterError):
        blk(obj)
    passive = dt.design(switching="passive")
    for k in ("er_db", "t_rise_s", "sigma_j_s", "t_ring_s", "tau_mode_s", "f_tog_max_hz"):
        del passive["params"][k]
    blk(passive)  # passive needs no switching parameters


@pytest.mark.parametrize("layout", ["tree", "chain"])
@pytest.mark.parametrize("K", [1, 4, 5, 32])
def test_transfer_matrix_conserves_power(layout, K):
    t, eps = 0.97, 0.02
    T = transfer_matrix(K, layout, t, eps)
    S = stage_counts(K, layout)
    assert np.all(T >= 0)
    if layout == "chain":  # every stage either sends out or passes on: Σ_k′ T = t^{stages seen}
        assert np.all(T.sum(axis=0) <= 1 + 1e-15)
        np.testing.assert_allclose(np.diag(T), [(t * (1 - eps)) ** S[k] for k in range(K)], rtol=1e-14)
    elif K & (K - 1) == 0:
        np.testing.assert_allclose(T.sum(axis=0), t ** S[0], rtol=1e-14)
        np.testing.assert_allclose(np.diag(T), (t * (1 - eps)) ** S[0], rtol=1e-14)


def test_toggles_and_off_identity():
    b, _, _ = blk(dt.design(K=4))
    d = b.design()
    assert np.all(d["eta_sw"] < 1) and np.all(d["leak_in"] > 0)
    b, _, _ = blk(dt.design(K=4, toggles={"crosstalk": False}))
    np.testing.assert_allclose(b.design()["eta_sw"], 0.99**2, rtol=1e-15)
    assert np.all(b.design()["leak_in"] == 0)
    b, _, _ = blk(dt.design(K=4, toggles={"insertion_loss": False}))
    np.testing.assert_allclose(b.design()["eta_sw"], 0.995**2, rtol=1e-15)
    b, _, _ = blk(dt.design(K=4, toggles={"enabled": False}))
    np.testing.assert_array_equal(b.matrix(), np.eye(4))
    assert all(not e.enabled for e in b.losses())
    b, _, _ = blk(dt.design(K=4, toggles={"rise_time_clipping": True}))
    assert b.stage_eps() > 0.005  # edge 5 ns from a 10 ns edge adds leak


def test_design_table_and_losses():
    b, _, _ = blk(dt.design(K=4, f_rep_hz=20e6, N_b=4, n_guard=1, L0_m=1.0))
    tab = b.tables()["demux_design"]
    cols = tab.provenance["columns"]
    L = tab.values[:, cols.index("L_m")]
    np.testing.assert_allclose(np.diff(L), 5 * C_M_PER_S / (20e6 * dt.N_G), rtol=1e-14)
    assert L[0] == 1.0
    el = b.losses()
    assert [e.segment for e in el] == ["D0", "D1", "D2", "D3"]
    np.testing.assert_allclose([e.eta for e in el], tab.values[:, cols.index("eta_sw")], rtol=1e-15)


def test_feasibility_warnings():
    b, _, _ = blk(dt.design())  # PC (f_tog 1 MHz, t_r 10 ns) at 100 MHz, N_b = 1: both conditions fail
    with pytest.warns(ChkWarning):
        b.check()
    assert {e.code for e in b.chk_log.events} == {"demux_rise_time", "demux_toggle_rate"}
    ok, _, _ = blk(dt.design(f_rep_hz=1e6))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ok.check()
    guard, _, _ = blk(dt.design(N_b=50, n_guard=1))  # block routing + 1 guard slot: edge fits in 2 T_rep,
    guard.check()  # first stage toggles at 100 MHz / 102 = 0.98 MHz < 1 MHz
    assert not guard.chk_log.events


def test_switch_off_is_identity_and_bias_static():
    b, st, pr = blk(dt.switch(toggles={"enabled": False}, bias_offset_rad=0.3))
    run_reference(b, st, pr, 5)
    assert st["act.sw.dpsi"] == 0.0 and st["pow.sw.trans"] == 1.0 and st["leak.sw.frac"] == 0.0
    b, st, pr = blk(dt.switch(n_stages=2, toggles={"drive_phase_noise": False, "bias_drift": True},
                              bias_offset_rad=0.1, kappa_pm=0.5))
    run_reference(b, st, pr, 3)
    assert st["act.sw.dpsi"] == pytest.approx(2 * 0.5 * 0.1, rel=1e-15)
    assert st["pow.sw.trans"] == pytest.approx(math.cos(0.05) ** 4, rel=1e-15)
    np.testing.assert_allclose(st["drv.dmx.eom"], [0.1])


@pytest.mark.parametrize("mode", ["timeseries", "declared"])
def test_switch_compiled_matches_reference(mode):
    obj = dt.switch(mode=mode, toggles={"bias_drift": True}, bias_rw_rad_per_rts=0.3, bias_offset_rad=0.01)
    assert dt.compiled_matches_reference(obj, 1e6, 3000) <= 1e-12


def test_switch_bands():
    b, _, _ = blk(dt.switch(mode="declared"), fs=1e6)
    (d,) = b.bands()
    assert d.targets == ("dmx",) and d.variance_outside == pytest.approx(b.sigma_fast() ** 2, rel=1e-15)
    assert b.sigma_fast() == pytest.approx(0.0100, rel=1e-3)
    b, _, _ = blk(dt.switch(mode="timeseries"), fs=1e6)
    assert b.bands()[0].variance_outside == 0.0
    b, _, _ = blk(dt.switch(mode="timeseries", toggles={"bias_drift": True}, bias_rw_rad_per_rts=1e-3), fs=1e6)
    assert b.bands()[0].variance_outside == pytest.approx(5 * 1e-6 / (2 * math.pi**2 * 5e5), rel=1e-12)


def test_delay_line_timing_and_compiled():
    setup = {"block": "pathtest_table", "name": "T", "params": {"out": "drv.F31.dT", "shape": [4], "unit": "K",
                                                                "n_rows": 1}}
    import cranberry.path.testing  # noqa: F401  (registers pathtest_table)

    obj = [setup, dt.delay(in_dT="drv.F31.dT", dtau_s=1e-12, rep_rel_drift=1e-6, dL_m=0.01)]
    b, st, pr = blk(obj)
    st["drv.F31.dT_table"] = np.array([[1.0, 2.0, 3.0, 2.0]])
    run_reference(b, st, pr, 1)
    L = 31 * C_M_PER_S / (1e8 * dt.N_G)
    want = 1e-12 + 39e-15 * L * 2.0 + 31e-8 * 1e-6
    assert st["timing.dl.dt"] == pytest.approx(want, rel=1e-12)
    assert dt.compiled_matches_reference(obj, 1e6, 10) <= 1e-12
    off = [setup, dt.delay(in_dT="drv.F31.dT", dtau_s=1e-12, rep_rel_drift=1e-6,
                           toggles={"static_mismatch": False, "delay_thermal": False, "rep_rate_drift": False})]
    b, st, pr = blk(off)
    st["drv.F31.dT_table"] = np.array([[1.0, 2.0, 3.0, 2.0]])
    run_reference(b, st, pr, 1)
    assert st["timing.dl.dt"] == 0.0


def test_delay_line_offline():
    b, _, _ = blk(dt.delay(dL_m=0.01))
    assert b.L_m == pytest.approx(63.3077, rel=1e-5) and b.tau0_s == pytest.approx(31e-8, rel=1e-12)
    assert b.gdd_mismatch_s2 == pytest.approx(-21.7e-27 * 0.01)
    (e,) = b.losses()
    assert e.kind == "fiber" and e.eta == pytest.approx(10 ** (-0.2e-3 * b.L_m / 10)) and e.segment == "F31"
    b, _, _ = blk(dt.delay(toggles={"fiber_loss": False, "static_mismatch": False}, dL_m=0.01))
    assert not b.losses()[0].enabled and b.gdd_mismatch_s2 == 0.0


def test_liu_helpers():
    L, n_g = 63.0, dt.N_G
    assert liu_loop_delay_s(1e-6, "upstream", L, n_g) == 1e-6
    assert liu_loop_delay_s(1e-6, "downstream", L, n_g) == pytest.approx(1e-6 + n_g * L / C_M_PER_S, rel=1e-15)
    assert actuator_sub("upstream", L, 10) == 0 and actuator_sub("downstream", L, 10) == 9
    assert actuator_sub(31.5, L, 10) == 5
    with pytest.raises(ValueError):
        liu_loop_delay_s(0.0, 70.0, L, n_g)
    assert probe_in_switch_window(325e-9, 10e-9, 0.0, 40e-9, 320e-9)
    assert not probe_in_switch_window(355e-9, 10e-9, 0.0, 40e-9, 320e-9)
