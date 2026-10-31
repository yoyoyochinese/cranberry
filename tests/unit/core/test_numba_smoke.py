"""numba installs and JIT-compiles on this platform, including the G02 kernel pattern: an @njit
closure with captured slot offsets, called from another @njit function (composition)."""

import numpy as np
import pytest

numba = pytest.importorskip("numba")

pytestmark = pytest.mark.numba


def test_njit_loop():
    @numba.njit
    def sumsq(x):
        s = 0.0
        for i in range(x.shape[0]):
            s += x[i] * x[i]
        return s

    x = np.arange(10.0)
    assert sumsq(x) == float((x * x).sum())


def test_closure_kernel_composition():
    def make_kernel(i_in, i_out, gain_idx):
        @numba.njit
        def kernel(S, Z, I, P):
            S[i_out] = P[gain_idx] * S[i_in]

        return kernel

    k1 = make_kernel(0, 1, 0)
    k2 = make_kernel(1, 2, 1)

    @numba.njit
    def loop(S, Z, I, P, n_steps):
        for _ in range(n_steps):
            k1(S, Z, I, P)
            k2(S, Z, I, P)
            I[0] += 1

    S = np.array([2.0, 0.0, 0.0])
    Z = np.zeros(1, dtype=np.complex128)
    I = np.zeros(1, dtype=np.int64)
    P = np.array([3.0, 0.5])
    loop(S, Z, I, P, 4)
    assert S[2] == 3.0 and I[0] == 4
