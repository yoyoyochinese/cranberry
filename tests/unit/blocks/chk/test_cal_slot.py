"""``cal_slot_check`` (Phase 6 fix of reviews/phase6_wp6c_p.md C-3 / C-7): the F_liu20 numbers of the review (guard slot:
L_k × 1.25, data slot after the seed at 50 ns keeps 0.974 of the recovery, ≈ 2 % loss at click probability 0.8; F_wg100
10 ns: blind), the calibration-frame weighting, the small-count warning and the guard-slot requirement."""

import math

import pytest

from cranberry.core.registry import make_block
from cranberry.core.block import BlockSpec

LIU = dict(N_b=4, n_guard=0, K=32, t_rep_s=5e-8, cal_fraction=1.0, t_int_s=0.05, mu_per_gate=3.43, qe=0.93,
           bg_per_gate=1e-6, dead_time_s=1e-8, recovery_s=1.1e-8, sigma_target_rad=0.01, lock_ratio_R=1.0,
           min_gates=100, recovery_loss_max=0.01, toggles={"enabled": True})


def check(**kw):
    return make_block(BlockSpec.from_obj({"block": "cal_slot_check", "name": "c", "params": {**LIU, **kw}}))


def codes(b):
    return {e.code for e in b.chk_log.events}


def test_guard_slot_liu20():
    b = check(placement="guard", n_guard=1)
    r = b.result
    assert r["lk_factor"] == pytest.approx(1.25) and r["data_fraction"] == pytest.approx(0.8)
    assert r["dnu_max_hz"] == pytest.approx(r["dnu_max_without_slot_hz"] / 1.25)
    p = -math.expm1(-(0.5 * 0.93 * 3.43 + 1e-6))
    assert r["p_click"] == pytest.approx(p) and p == pytest.approx(0.8, abs=0.005)
    rec = -math.expm1(-(5e-8 - 1e-8) / 1.1e-8)
    assert rec == pytest.approx(0.974, abs=1e-3)
    assert r["eta_next_slot"] == pytest.approx(1 - p * (1 - rec))
    assert 0.02 < r["recovery_loss_mean"] < 0.022
    assert r["gates_per_window"] == pytest.approx(0.05 / (32 * 5 * 5e-8))
    assert codes(b) == {"dbs_cal_guard_slot", "dbs_recovery"}
    # F_wg100 (10 ns slots): the next slot is inside the dead time
    w = check(placement="guard", n_guard=1, N_b=50, t_rep_s=1e-8).result
    assert w["eta_next_slot"] == pytest.approx(1 - p)


def test_calibration_frames():
    b = check(placement="cal_frames", cal_fraction=0.2)
    r = b.result
    assert r["lk_factor"] == 1.0 and r["data_fraction"] == pytest.approx(0.8)
    assert r["gates_per_frame"] == 4 and r["f_gate_hz"] == pytest.approx(4 / (32 * 4 * 5e-8))
    assert r["gates_per_window"] == pytest.approx(4 * 0.2 * 0.05 / (32 * 4 * 5e-8))
    assert r["recovery_weight"] == pytest.approx(0.25)
    assert codes(b) == set()                      # 0.25 × 2.1 % < 1 %
    with pytest.raises(ValueError):
        check(placement="cal_frames", cal_fraction=1.0)


def test_ideal_and_few_counts():
    assert "dbs_cal_slot_ideal" in codes(check(placement="ideal_free"))
    b = check(placement="ideal_free", cal_fraction=0.01, t_int_s=0.01)   # the review's sweep corner: 15.6 gates
    assert b.result["gates_per_window"] == pytest.approx(15.6, abs=0.05)
    assert "dbs_few_counts" in codes(b)


def test_guard_needs_a_slot():
    with pytest.raises(ValueError, match="n_guard"):
        check(placement="guard", n_guard=0)
