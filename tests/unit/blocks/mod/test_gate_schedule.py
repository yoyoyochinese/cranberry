"""``gate_schedule`` (WP6-G): window timing, edges, leak, anti-sync detector gate, jitter, toggles, reference = kernel."""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.mod.gate_schedule import make_gate_schedule
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

pytestmark = pytest.mark.numba
FS = 1.0e6


def gate(name="gw", detector_gate="complement", write_phase=True, toggles=None, **kw):
    p = {"detector_gate": detector_gate, "write_phase": write_phase, "period_s": 100e-6, "lock_fraction": 0.655,
         "phase": 0.0, "edge_s": 0.0, "timing_jitter_s": 0.0, "leak": 0.0,
         "toggles": {"gating": True, "jitter": False, "edges": False, "leak": False, **(toggles or {})}}
    if detector_gate == "anti_sync":
        p.update(duty_det=0.25, settle_s=10e-6)
    p.update(kw)
    return {"block": "gate_schedule", "name": name, "params": p}


def series(g, n, compiled=False):
    slots = ["gate.gw.lock", "gate.gw.meas"] + (["gate.gw.phase"] if g["params"]["write_phase"] else [])
    rec, _, _ = pt.run_and_record([g, pt.record(slots, n)], FS, n, compiled=compiled, seed=3)
    return rec


def test_complement_windows_and_phase():
    r = series(gate(), 300)
    n = np.arange(300)
    x = (n % 100) / 100.0
    np.testing.assert_array_equal(r["gate.gw.lock"], np.where(x < 0.655, 1.0, 0.0))
    np.testing.assert_array_equal(r["gate.gw.meas"], 1.0 - r["gate.gw.lock"])
    np.testing.assert_allclose(r["gate.gw.phase"], 2 * math.pi * x, atol=1e-9)
    assert r["gate.gw.lock"][:100].sum() == 66  # ceil(65.5): Borealis 65 µs / 35 µs class


def test_edges_and_leak():
    r = series(gate(edge_s=4e-6, leak=1e-3, toggles={"edges": True, "leak": True}), 100)
    lock, meas = r["gate.gw.lock"], r["gate.gw.meas"]
    np.testing.assert_allclose(lock[:5], 1e-3 + (1 - 1e-3) * np.array([0, 0.25, 0.5, 0.75, 1.0]), rtol=1e-12)
    np.testing.assert_allclose(lock[66:71], 1e-3 + (1 - 1e-3) * np.array([0.875, 0.625, 0.375, 0.125, 0.0]), atol=1e-9)
    assert lock.min() == pytest.approx(1e-3) and lock.max() == 1.0
    assert meas[:70].sum() == 0 and meas[70:].all()  # measurement window after the falling edge (65.5 + 4 µs)
    off = series(gate(edge_s=4e-6, leak=1e-3, toggles={"edges": False, "leak": False}), 100)
    np.testing.assert_array_equal(off["gate.gw.lock"], np.where(np.arange(100) < 66, 1.0, 0.0))


def test_anti_sync_detector_window():
    r = series(gate(detector_gate="anti_sync", lock_fraction=0.5), 200)
    meas = r["gate.gw.meas"][:100]
    assert np.flatnonzero(meas).tolist() == list(range(60, 85))  # [D_pr T + settle, + duty_det T)
    assert np.all(r["gate.gw.lock"][60:85] == 0.0)


def test_gating_off_is_continuous_lock_and_measurement():
    r = series(gate(toggles={"gating": False, "jitter": True}, timing_jitter_s=1e-6), 150)
    assert np.all(r["gate.gw.lock"] == 1.0) and np.all(r["gate.gw.meas"] == 1.0)


def test_jitter_shifts_whole_periods():
    r = series(gate(timing_jitter_s=3e-6, toggles={"jitter": True}), 2000)
    lock = r["gate.gw.lock"]
    starts = np.flatnonzero(np.diff(lock) > 0.5) + 1
    falls = (np.flatnonzero(np.diff(lock) < -0.5) + 1)[1:]  # drop the end of the unjittered first window
    shift = starts - 100 * np.round(starts / 100)
    assert np.abs(shift).max() <= 15 and np.std(shift) > 1.0  # white per-period shifts of σ = 3 samples
    m = min(len(starts), len(falls))
    assert np.all(np.isin(falls[:m] - starts[:m], [65, 66, 67]))  # the whole window moves, its length stays


@pytest.mark.parametrize("mode", ["complement", "anti_sync"])
def test_reference_equals_kernel_and_chunk_invariance(mode):
    g = gate(detector_gate=mode, lock_fraction=0.4 if mode == "anti_sync" else 0.655, edge_s=2.5e-6, leak=1e-4,
             timing_jitter_s=2e-6, toggles={"jitter": True, "edges": True, "leak": True})
    a = series(g, 3000, compiled=False)
    b = series(g, 3000, compiled=True)
    for k in a:
        np.testing.assert_allclose(a[k], b[k], rtol=0, atol=1e-12)
    lay = SlotLayout(FS)
    blk = build(spec_from_obj([g]), lay, 3)
    st, pr = State(lay), Params(lay)
    for m in (7, 1000, 1, 1992):  # pieces of any length: same realisation (G11)
        run_compiled(blk, st, pr, m)
    assert float(st["gate.gw.lock"]) == a["gate.gw.lock"][-1]
    lay2 = SlotLayout(FS)
    blk2 = build(spec_from_obj([g]), lay2, 3)
    st2 = State(lay2)
    run_reference(blk2, st2, Params(lay2), 3000)
    np.testing.assert_array_equal(st["gate.gw.state"], st2["gate.gw.state"])


def test_validation():
    with pytest.raises(ValueError, match="lock_fraction"):
        make_gate_schedule(gate(lock_fraction=1.0)["params"])
    with pytest.raises(ValueError, match="phase"):
        make_gate_schedule(gate(phase=1.0)["params"])
    with pytest.raises(ValueError, match="timing_jitter"):
        make_gate_schedule(gate(timing_jitter_s=20e-6)["params"])
    with pytest.raises(ValueError, match="anti_sync"):
        make_gate_schedule(gate(detector_gate="anti_sync", lock_fraction=0.7, settle_s=10e-6, duty_det=0.25)["params"])
    with pytest.raises(ValueError, match="detector_gate"):
        make_gate_schedule(gate(detector_gate="sync")["params"])
    with pytest.raises(ValueError, match="2 block steps"):
        build(spec_from_obj([gate(period_s=1.5e-6)]), SlotLayout(FS), 0)
    blk = make_gate_schedule(gate(write_phase=False)["params"], name="gw")
    assert blk.writes == ("gate.gw.lock", "gate.gw.meas", "gate.gw.state")
    assert set(blk.runtime_params) == {"period_s", "lock_fraction", "phase", "edge_s", "timing_jitter_s", "leak"}
