"""core.rng (G11): keyed streams are reproducible, order-independent and mutually independent."""

import numpy as np
import pytest

from cranberry.core import rng


def draws(root, name, n=8):
    return rng.stream(root, name).standard_normal(n)


def test_reproducible():
    np.testing.assert_array_equal(draws(20261031, "lock1.pid"), draws(20261031, "lock1.pid"))


def test_independent_of_creation_order_and_of_other_streams():
    a1 = rng.stream(7, "a")
    b1 = rng.stream(7, "b")
    xa1, xb1 = a1.standard_normal(5), b1.standard_normal(5)
    b2 = rng.stream(7, "b")  # created first this time, and "a" draws interleaved
    a2 = rng.stream(7, "a")
    xb2 = b2.standard_normal(5)
    _ = rng.stream(7, "c").standard_normal(100)  # an extra stream changes nothing
    xa2 = a2.standard_normal(5)
    np.testing.assert_array_equal(xa1, xa2)
    np.testing.assert_array_equal(xb1, xb2)


def test_different_names_and_roots_differ():
    assert not np.array_equal(draws(1, "dx1"), draws(1, "dx2"))
    assert not np.array_equal(draws(1, "dx1"), draws(2, "dx1"))


def test_streams_uncorrelated():
    x = rng.stream(3, "seg.thermal").standard_normal(200_000)
    y = rng.stream(3, "seg.mech").standard_normal(200_000)
    assert abs(np.corrcoef(x, y)[0, 1]) < 5 / np.sqrt(x.size)


def test_key_is_stable_across_processes():
    """SHA-256 based, not Python's salted hash(): pin the value so a change is noticed."""
    assert rng.stream_key("lock1.pid") == rng.stream_key("lock1.pid")
    assert rng.stream_key("lock1.pid") == (2809901333, 2587214443, 1906558330, 573695195)
    ss = rng.child_seed_sequence(5, "lock1.pid")
    assert ss.entropy == 5 and tuple(ss.spawn_key) == rng.stream_key("lock1.pid")


def test_root_seed_sequence_and_substream():
    root = np.random.SeedSequence(11, spawn_key=(3,))
    s = rng.stream(root, "x")
    assert tuple(s.bit_generator.seed_seq.spawn_key) == (3, *rng.stream_key("x"))
    parent = rng.stream(11, "noise")
    before = rng.stream(11, "noise").standard_normal(3)
    sub_a = rng.substream(parent, "ou0").standard_normal(4)
    np.testing.assert_array_equal(parent.standard_normal(3), before)  # parent not consumed
    np.testing.assert_array_equal(rng.substream(rng.stream(11, "noise"), "ou0").standard_normal(4), sub_a)
    assert not np.array_equal(rng.substream(parent, "ou1").standard_normal(4), sub_a)


@pytest.mark.parametrize("bad", [-1, 1.5, "3", True])
def test_bad_roots(bad):
    with pytest.raises((TypeError, ValueError)):
        rng.stream(bad, "x")


def test_bad_name():
    with pytest.raises(ValueError):
        rng.stream(1, "")
