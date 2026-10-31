"""G01 performance spike (WP-F1, PLAN §3 G01 condition): ~10 trivial-but-realistic blocks composed
through the real framework (build → compile_loop, njit closures per G02), 10⁶ steps.

Reports ns/step of the compiled loop (incl. per-piece ``prepare``: pre-drawing noise), µs/step of the
pure-Python reference, and compile time. G01 holds while the compiled loop is <= 1 µs/step.

Run: ``uv run pytest -m slow tests/bench -s`` or ``uv run python tests/bench/test_perf_spike.py``.
"""

import time

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.compile import clear_loop_cache, compile_loop, run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

pytestmark = [pytest.mark.slow, pytest.mark.numba]

FS = 1.0e6
N_COMPILED = 1_000_000
N_REFERENCE = 20_000
G01_LIMIT_NS = 1000.0


def bench_obj():
    """10 blocks: 2 noise IIRs, sine, sum, decimated PID, Thiran delay, saturating actuator,
    integrator plant, sensor gain, monitor sum."""
    return [
        toy("toy_noise_ar1", "n1", out="v.n1", unit="V", a=0.999, b=0.05),
        toy("toy_noise_ar1", "n2", out="v.n2", unit="V", a=0.9, b=0.2),
        toy("toy_sine", "dist", out="v.d", unit="V", amp=0.3, freq_hz=2.1e3),
        toy("toy_sum", "sens", ins=["v.x", "v.n1", "v.n2", "v.d"], gains=[1.0, 1.0, 0.5, 1.0], out="v.y", unit="V"),
        toy("toy_pid", "pid", 2, in_="v.y", out="ctrl.pid.u", unit="V", kp=0.4, ki=3e3, kd=2e-7, setpoint=0.0),
        {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out",
                                                           "tau_s": 3.37e-6, "method": "thiran1"}},
        toy("toy_sat", "act", in_="delay.dl.out", out="act.u", unit="V", gain=1.0, limit=0.8),
        toy("toy_integrator", "plant", in_="act.u", out="psi.plant.x", unit="V", k=2e4),
        toy("toy_gain", "pd", in_="psi.plant.x", out="v.x", unit="V", gain=1.0),
        toy("toy_sum", "mon", ins=["v.y", "act.u"], gains=[1.0, -1.0], out="metric.mon", unit="V"),
    ]


def make():
    lay = SlotLayout(FS)
    blk = build(spec_from_obj(bench_obj()), lay, 1)
    return blk, State(lay), Params(lay)


def measure():
    blk, st, p = make()
    assert len(list(iter_leaves(blk))) == 10
    clear_loop_cache()
    t0 = time.perf_counter()
    loop = compile_loop(blk)
    t_compile = time.perf_counter() - t0
    run_compiled(blk, st, p, 1000, loop)  # warm caches
    t0 = time.perf_counter()
    run_compiled(blk, st, p, N_COMPILED, loop)
    ns_step = (time.perf_counter() - t0) / N_COMPILED * 1e9
    blk_r, st_r, p_r = make()
    t0 = time.perf_counter()
    run_reference(blk_r, st_r, p_r, N_REFERENCE)
    us_step_ref = (time.perf_counter() - t0) / N_REFERENCE * 1e6
    assert np.isfinite(st.S).all()
    return {"compile_s": t_compile, "compiled_ns_per_step": ns_step, "reference_us_per_step": us_step_ref}


def test_g01_performance_spike():
    r = measure()
    print(f"\nG01 spike (10 blocks): compile {r['compile_s']:.2f} s, compiled {r['compiled_ns_per_step']:.1f} ns/step, "
          f"reference {r['reference_us_per_step']:.1f} us/step "
          f"(speed-up {r['reference_us_per_step'] * 1e3 / r['compiled_ns_per_step']:.0f}x)")
    assert r["compiled_ns_per_step"] <= G01_LIMIT_NS, "G01 condition violated: revisit G02 with the user"


if __name__ == "__main__":
    print(measure())
