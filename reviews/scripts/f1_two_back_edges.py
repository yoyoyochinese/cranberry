"""Phase 5 M1 framework review, finding F1-1: a feedback loop with two back edges passes
check_pipeline silently and has a closure delay of 2/f_sim instead of the 1/f_sim of D55-1.

Same blocks, same tau_s = 0; only the position of the delay_buffer differs:
  A (correct) : plant -> sensor -> pid -> d -> act      (one back edge: act -> plant)
  B (misorder): plant -> sensor -> d -> pid -> act      (two back edges: pid -> d, act -> plant)
Both build without a PipelineCheckError. We record the sample at which the plant first moves after
a setpoint step at n = 0 (integrator plant: x[n+1] = x[n] + h*k*u).

Run: uv run python reviews/scripts/f1_two_back_edges.py
"""

import numpy as np

import cranberry.core.testing  # noqa: F401  (toy kinds)
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

FS = 1.0e6


def delay(name, inp, out, tau_s=0.0):
    return {"block": "delay_buffer", "name": name,
            "params": {"in": inp, "out": out, "tau_s": tau_s, "method": "linear", "unit": "V"}}


def blocks():
    plant = toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="rad", k=1.0)
    sensor = toy("toy_gain", "sensor", in_="psi.plant.x", out="v.y", unit="V", gain=1.0)
    pid = toy("toy_pid", "pid", in_="v.y", out="ctrl.pid.u", unit="V", kp=0.1, ki=0.0, kd=0.0, setpoint=1.0)
    d = delay("d", "ctrl.pid.u", "delay.d.out")
    act = toy("toy_gain", "act", in_="delay.d.out", out="act.u", unit="V", gain=1.0)
    return plant, sensor, pid, d, act


def first_move(order):
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(order), lay, 1)
    st, pa = State(lay), Params(lay)
    xs = []
    for _ in range(6):
        run_reference(blk, st, pa, 1)
        xs.append(float(st["psi.plant.x"]))
    return int(np.flatnonzero(np.abs(xs) > 0)[0]), xs


if __name__ == "__main__":
    plant, sensor, pid, d, act = blocks()
    nA, xA = first_move([plant, sensor, pid, d, act])
    nB, xB = first_move([plant, sensor, d, pid, act])
    print(f"A correct order : builds OK, plant first moves at sample {nA}  x={np.round(xA, 9)}")
    print(f"B misordered    : builds OK, plant first moves at sample {nB}  x={np.round(xB, 9)}")
    print(f"extra closure delay in B: {nB - nA} sample(s) = {(nB - nA) / FS * 1e6:.1f} us (D55-1 assumes 0)")
    assert nB - nA == 1
