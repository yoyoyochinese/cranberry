"""SlotLayout / State / Params contract (G02): flat arrays, named views, in-place updates."""

import numpy as np
import pytest

from cranberry.core.layout import LayoutFrozenError, SlotLayout
from cranberry.core.state import Params, State


def make_layout():
    lay = SlotLayout(f_sim_hz=1.0e6)
    lay.add("ctrl.lock1.u", unit="V")
    lay.add("J.arm1", (2, 2), np.complex128, init=0)
    lay.add("drv.arm1.dT", (4,), unit="K", ring_len=3)
    lay.add("lock.lockD.state", dtype=np.int64)
    lay.add_param("lock1.pid.kp", 0.25)
    lay.add_param("lock1.pid.sos", np.ones((2, 3)), (2, 3))
    return lay


def test_reserved_slots_and_sizes():
    lay = make_layout()
    assert lay.slot_names()[:2] == ("t", "n")
    assert lay.dt_s == pytest.approx(1e-6)
    sizes = lay.sizes
    assert sizes["S"] == 1 + 1 + 12  # t, ctrl.lock1.u, drv ring 3x4
    assert sizes["Z"] == 4
    assert sizes["I"] == 1 + 1 + 1  # n, ring head, lock state
    assert sizes["P"] == 7
    assert lay.slot("drv.arm1.dT").storage_shape == (3, 4)
    assert lay.slot("drv.arm1.dT").head >= 0 and lay.slot("ctrl.lock1.u").head == -1


def test_naming_rules():
    lay = SlotLayout(1.0)
    for bad in ["foo.x", "ctrl", "1ctrl.x", "ctrl..x", "ctrl.x-y"]:
        with pytest.raises(ValueError):
            lay.add(bad)
    lay.add("ctrl.a")
    with pytest.raises(ValueError):
        lay.add("ctrl.a")
    with pytest.raises(TypeError):
        lay.add("ctrl.b", dtype=np.float32)


def test_freeze_and_allocation_rules():
    lay = make_layout()
    with pytest.raises(LayoutFrozenError):
        State(lay)
    lay.freeze()
    with pytest.raises(LayoutFrozenError):
        lay.add("ctrl.late")


def test_named_views_are_in_place():
    lay = make_layout()
    lay.freeze()
    st = State(lay)
    st["ctrl.lock1.u"] = 1.5
    assert st.S[lay.slot("ctrl.lock1.u").offset] == 1.5
    assert float(st.ctrl.lock1.u) == 1.5
    st.ctrl.lock1.u = -2.0
    assert float(st["ctrl.lock1.u"]) == -2.0
    st["J.arm1"] = np.eye(2)
    np.testing.assert_array_equal(st.Z[lay.slot("J.arm1").offset : lay.slot("J.arm1").stop], [1, 0, 0, 1])
    v = st.view("drv.arm1.dT")
    v[1, 2] = 7.0
    assert st.drv.arm1.dT[1, 2] == 7.0
    assert st.ring_head("drv.arm1.dT") == 0
    assert st.n == 0 and st.t == 0.0
    snap = st.copy()
    st["ctrl.lock1.u"] = 9.0
    assert float(snap["ctrl.lock1.u"]) == -2.0


def test_params():
    lay = make_layout()
    lay.freeze()
    p = Params(lay)
    assert float(p["lock1.pid.kp"]) == 0.25
    p.set({"lock1.pid.kp": 0.5})
    assert p.P[lay.param("lock1.pid.kp").offset] == 0.5
    with pytest.raises(KeyError):
        p.set({"lock1.pid.kp": 1.0, "nope": 2.0})
    assert float(p["lock1.pid.kp"]) == 0.5  # nothing written when a name is invalid
    assert p.as_dict()["lock1.pid.sos"] == [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
