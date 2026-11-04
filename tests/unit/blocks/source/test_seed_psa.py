"""``seed_psa`` (Phase 6 fix of reviews/phase6_wp6c_p.md C-1): a coherent seed through the phase-sensitive amplifier.
Against the independent α → α cosh r + α* sinh r of the review (p6_wp6c_counting.py G): the (π/2, 0) preparation keeps
the offset with intensity ratio e^{−4r}, the seed's own input phase is amplified e^{2r} (port k) / e^{−2r} (port 0), and
(±π/4) collapses the offset to 2 atan(e^{−2r}); amplification off is the identity; reference = njit."""

import math

import numpy as np
import pytest

from cranberry.blocks.source.seed_psa import psa_gain, psa_phase
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

R = 1.8


def psa_field(delta, r):
    a = np.exp(1j * delta)
    return a * math.cosh(r) + np.conj(a) * math.sinh(r)


def const(name, x):
    return {"block": "ramp", "name": name, "params": {"out": f"v.{name}", "waveform": "sawtooth", "amplitude_v": 0.0,
                                                     "period_s": 1.0, "offset_v": x, "toggles": {"enabled": True}}}


def build_psa(seed, sqz, delta, r=R, amp=True, w=1.0):
    objs = [const("s", seed), const("q", sqz),
            {"block": "seed_psa", "name": "psa", "params": {"in_seed": "v.s", "in_sqz": "v.q", "out": "v.o", "out_i": "v.i",
                                                            "r": r, "delta_rad": delta, "p_in_rel": w,
                                                            "toggles": {"amplification": amp}}}]
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(objs), lay, 1)
    return blk, State(lay), Params(lay)


def out(seed, sqz, delta, **kw):
    blk, st, pr = build_psa(seed, sqz, delta, **kw)
    run_reference(blk, st, pr, 1)
    return float(st["v.o"]), float(st["v.i"])


def test_against_the_field_model():
    for d in (0.0, 0.3, math.pi / 4, 1.2, math.pi / 2, 2.0, -0.7):
        a = psa_field(d, R)
        assert psa_phase(d, R) == pytest.approx(np.angle(a), abs=1e-12)
        assert psa_gain(d, R) == pytest.approx(abs(a) ** 2, rel=1e-12)


def test_review_table():
    # (π/2, 0): offset π/2 kept, intensity ratio e^{−4r} = 7.47e−4, input-phase gains 36.6 / 0.027
    assert psa_phase(math.pi / 2, R) - psa_phase(0.0, R) == pytest.approx(math.pi / 2, abs=1e-12)
    assert psa_gain(math.pi / 2, R) / psa_gain(0.0, R) == pytest.approx(7.47e-4, rel=1e-3)
    eps = 1e-6
    gk = (psa_phase(math.pi / 2 + eps, R) - psa_phase(math.pi / 2 - eps, R)) / (2 * eps)
    g0 = (psa_phase(eps, R) - psa_phase(-eps, R)) / (2 * eps)
    assert gk == pytest.approx(math.exp(2 * R), rel=1e-6) and gk == pytest.approx(36.6, abs=0.05)
    assert g0 == pytest.approx(math.exp(-2 * R), rel=1e-6)
    # (±π/4): equal intensities, offset 2 atan(e^{−2r}) = 0.055 rad
    assert psa_gain(math.pi / 4, R) == pytest.approx(psa_gain(-math.pi / 4, R), rel=1e-12)
    off = psa_phase(math.pi / 4, R) - psa_phase(-math.pi / 4, R)
    assert off == pytest.approx(2 * math.atan(math.exp(-2 * R)), rel=1e-12) and off == pytest.approx(0.055, abs=5e-4)


def test_block_output():
    # δ = delta + seed − sqz; out = sqz − π/2 + ψ(δ); out_i = p_in · g(δ)
    o, i = out(0.02, 0.3, math.pi / 2, w=0.5)
    d = math.pi / 2 + 0.02 - 0.3
    assert o == pytest.approx(0.3 - math.pi / 2 + psa_phase(d, R), abs=1e-12)
    assert i == pytest.approx(0.5 * psa_gain(d, R), rel=1e-12)
    # a seed on the squeezed axis follows its own phase × e^{2r} around the squeeze angle
    o1, _ = out(1e-4, 0.0, math.pi / 2)
    assert o1 == pytest.approx(math.exp(2 * R) * 1e-4, rel=1e-3)
    # amplification off (or r = 0): out = seed + delta − π/2, out_i = p_in
    for kw in ({"amp": False}, {"r": 0.0}):
        o, i = out(0.4, -1.1, 0.25, w=0.7, **kw)
        assert o == pytest.approx(0.4 + 0.25 - math.pi / 2, abs=1e-12) and i == pytest.approx(0.7, rel=1e-12)


@pytest.mark.numba
def test_reference_equals_compiled():
    a, sa, pa = build_psa(0.13, -0.4, 1.0)
    b, sb, pb = build_psa(0.13, -0.4, 1.0)
    run_reference(a, sa, pa, 3)
    run_compiled(b, sb, pb, 3)
    for k in ("v.o", "v.i"):
        assert float(sa[k]) == pytest.approx(float(sb[k]), rel=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        build_psa(0.0, 0.0, 0.0, r=-0.1)
