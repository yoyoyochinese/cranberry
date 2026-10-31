"""core.rates step counting, core.rng.DrawBuffer (G11 sequential draws), runner Recorder."""

import numpy as np
import pytest

from cranberry.core.block import BaseBlock
from cranberry.core.compile import compile_loop, run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import Seq, build
from cranberry.core.rates import PIECE_SAMPLES, max_steps_in, n_steps_in
from cranberry.core.rng import DrawBuffer, stream
from cranberry.core.state import Params, State
from cranberry.runner.run import Recorder


def test_n_steps_in_brute_force():
    for r in (1, 2, 3, 7):
        for n0 in range(0, 15):
            for m in range(0, 20):
                expect = sum(1 for n in range(n0, n0 + m) if n % r == 0)
                assert n_steps_in(n0, m, r) == expect
                assert expect <= max_steps_in(m, r)
    with pytest.raises(ValueError):
        n_steps_in(0, 1, 0)


@pytest.mark.parametrize("rate_div", [1, 3])
def test_draw_buffer_sequential_for_any_pieces(rate_div):
    lay = SlotLayout(1.0)
    db = DrawBuffer("n1", 2, "standard_normal")
    db.init(lay, stream(5, "n1"), rate_div)
    lay.freeze()
    db.bind(lay)
    st = State(lay)
    assert db.capacity == max_steps_in(PIECE_SAMPLES, rate_div)
    got, n0 = [], 0
    for m in (5, 1, PIECE_SAMPLES, 13, 0, 400):
        db.refill(st, n0, m)
        for _ in range(n_steps_in(n0, m, rate_div) * 2):
            got.append(db.take(st))
        n0 += m
    expect = stream(5, "n1").standard_normal(len(got))
    np.testing.assert_array_equal(got, expect)
    with pytest.raises(ValueError):
        db.refill(st, 0, PIECE_SAMPLES + 1)


def test_draw_buffer_validation():
    with pytest.raises(ValueError):
        DrawBuffer("x", 1, "gumbel")
    with pytest.raises(ValueError):
        DrawBuffer("x", 0, "random")
    db = DrawBuffer("x", 1, "random", label="phase")
    assert db.slot_names == ("rng.x.phase.buf", "rng.x.phase.idx")


class Writer(BaseBlock):
    """Writes a float vector, an int and a complex 2×2 slot each step (Python and kernel)."""

    def __init__(self):
        super().__init__("writer", "w", {}, writes=("v.w", "lock.w.state", "J.w"))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add("v.w", (2,))
        layout.add("lock.w.state", (), np.int64)
        layout.add("J.w", (2, 2), np.complex128)

    def bind(self, layout):
        import numba

        o_v, o_i, o_z, o_n = (layout.slot(s).offset for s in ("v.w", "lock.w.state", "J.w", "n"))

        @numba.njit
        def k(S, Z, I, P):  # noqa: N803
            n = I[o_n]
            S[o_v] = n
            S[o_v + 1] = -n
            I[o_i] = 2 * n
            Z[o_z + 3] = n + 1j * n

        self.kernel = k
        self._o = (o_v, o_i, o_z, o_n)

    def __call__(self, state, params):
        o_v, o_i, o_z, o_n = self._o
        n = state.I[o_n]
        state.S[o_v], state.S[o_v + 1], state.I[o_i], state.Z[o_z + 3] = n, -n, 2 * n, n + 1j * n
        return state


@pytest.mark.parametrize("compiled", [False, True])
def test_recorder_float_int_complex(compiled):
    rec = Recorder(["J.w", "v.w", "lock.w.state"], decim=3, chunk_samples=10)
    lay = SlotLayout(1.0e3)
    blk = build(Seq([Writer(), rec]), lay, 0)
    st, p = State(lay), Params(lay)
    tables = []
    for _ in range(3):
        (run_compiled if compiled else run_reference)(blk, st, p, 10)
        tables.append(rec.harvest(st))
    out = rec.split(np.concatenate(tables), lay)
    n = np.arange(0, 30, 3)
    np.testing.assert_array_equal(out["v.w"], np.stack([n, -n], axis=1))
    np.testing.assert_array_equal(out["lock.w.state"], 2 * n)
    assert out["lock.w.state"].dtype == np.int64
    assert out["J.w"].shape == (10, 2, 2)
    np.testing.assert_array_equal(out["J.w"][:, 1, 1], n + 1j * n)
    if compiled:
        assert compile_loop(blk) is not None


def test_recorder_overflow_raises():
    rec = Recorder(["v.w"], decim=1, chunk_samples=4)
    lay = SlotLayout(1.0)
    blk = build(Seq([Writer(), rec]), lay, 0)
    st, p = State(lay), Params(lay)
    with pytest.raises(IndexError):
        run_reference(blk, st, p, 5)
