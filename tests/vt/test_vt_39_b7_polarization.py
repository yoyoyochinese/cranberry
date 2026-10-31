"""VT-39 (M4 part, WP-B7): polarization rows of the probe-cancellation tests (T10 §3.1, §3.4.2, §7
tests 8 and 9; T09 §1.6). The full VT-39 is WP-B8's ``test_vt_39_probe_cancellation.py``; the F6
groundwork is ``test_vt_39_f6_counterprop_lock.py``.

* T10 test 8 (SMF polarization of the backscatter): for a Haar-random M, ⟨|⟨p|MᵀM|s⟩|²⟩ = 1/3 ± 1 %
  (10⁵ samples of ``segment_jones``'s U0 generator); the same through ``segment_jones`` blocks (Haar U0
  per segment, statistical tolerance); M diagonal (PM) ⇒ 0.
* T10 test 9 (reciprocity, SMF orthogonal probe): a lossless M (Haar ``segment_jones`` chain + a static
  polarization controller) with M e_in ∝ p: the backward probe launched in s ⊥ p leaves in e_in^⊥
  (no power in e_in, 1e−24), the squeezed / probe phase sum is twice the common phase ψ and the
  difference is the retardance 2α, α = arg(M e_in)_p (1e−12).
"""

import math

import numpy as np
import pytest

import cranberry.path.testing as pt
from cranberry.blocks.pol.segment_jones import haar_su2
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.rng import stream
from cranberry.core.state import Params, State

ON = {"pol_birefringence_drift": True, "pol_crosstalk": True, "pol_smf_drift": True}
SMF = {"model": "smf", "axis_drivers": [], "axis_stokes": [1, 0, 0], "u0": "haar", "drift": "diffusion",
       "D_pol_rad2_per_s": 0.0, "k_T_rad_per_m_K": 0.0}
X, Y = np.array([1, 0]), np.array([0, 1])


def jseg(id, L=1.0, **kw):
    return pt.seg_spec(id, L, kind="smf", pol_model="jones", pol_params={"J_slot": f"J.{id}"}, **kw)


def sj(topo, ids, cfg=SMF):
    return {"block": "segment_jones", "name": "sj", "params": {"topology": topo, "segments": {s: dict(cfg) for s in ids},
                                                               "toggles": dict(ON)}}


def overlap_sq(M):
    return abs(X @ (M.T @ M) @ Y) ** 2  # |⟨p|MᵀM|s⟩|², p = x (real), s = y


def test_t10_8_haar_one_third():
    rng = stream(8, "vt39.haar")
    vals = np.array([overlap_sq(haar_su2(rng)) for _ in range(100_000)])
    assert np.mean(vals) == pytest.approx(1 / 3, rel=0.01)
    assert overlap_sq(np.diag([np.exp(0.7j), np.exp(-0.7j)])) == 0.0  # PM (diagonal)


def test_t10_8_through_segment_jones():
    ids = [f"s{i}" for i in range(400)]
    topo = pt.topo_spec([jseg(s) for s in ids], [pt.beam_spec("b", [pt.step_spec(ids[0])], "o", jones_out=(1, 0))],
                        ["o"], [])
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj([sj(topo, ids)]), lay, 39)
    s, p = State(lay), Params(lay)
    run_reference(blk, s, p, 1)
    vals = np.array([overlap_sq(np.array(s[f"J.{i}"])) for i in ids])  # one fiber per segment
    assert np.mean(vals) == pytest.approx(1 / 3, abs=4 * np.std(vals) / math.sqrt(len(vals)))
    pm = {"model": "pm", "axis_drivers": [], "dpsi_b0_rad": 0.9, "K_T_rad_per_m_K": -2.05, "theta_in_rad": 0.0,
          "theta_out_rad": 0.0}
    lay2 = SlotLayout(1e6)
    blk = build(spec_from_obj([sj(topo, ids[:1], pm)]), lay2, 0)
    s2, p2 = State(lay2), Params(lay2)
    run_reference(blk, s2, p2, 1)
    assert overlap_sq(np.array(s2["J.s0"])) == 0.0


def controller(c, beta):
    """SU(2) V with V·c = (e^{iβ}, 0) for a unit vector c (polarization controller + retardance β)."""
    return np.diag([np.exp(1j * beta), np.exp(-1j * beta)]) @ np.array([[np.conj(c[0]), np.conj(c[1])], [-c[1], c[0]]])


def test_t10_9_reciprocity_orthogonal_probe():
    ids = [f"s{i}" for i in range(5)]
    psi = 0.83  # common phase from segment s0: coefficient × constant dT

    def rig(V):
        segs = [jseg(s, drivers={"dT": "uniform"} if s == "s0" else None) for s in ids]
        segs.append(pt.seg_spec("pc", 0.0, pol_model="jones", pol_params={"J": V.tolist()}))
        fw = [pt.step_spec(s) for s in (*ids, "pc")]
        bw = [pt.step_spec(s, -1) for s in reversed((*ids, "pc"))]
        beams = [pt.beam_spec("sq", fw, "o", jones_in=X, jones_out=X),
                 pt.beam_spec("pr", bw, "o", jones_in=Y, jones_out=Y),
                 pt.beam_spec("pr_x", bw, "o", jones_in=Y, jones_out=X)]
        topo = pt.topo_spec(segs, beams, ["o"], pt.coef_spec("s0", "dT", "*", psi))
        pipe = [pt.pathtest("pathtest_table", "dT", "drv.s0.dT", (), "K", n_rows=3),
                pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=1), sj(topo, ids),
                *pt.path_pipeline(topo, {"las": "laser.las.phi"}, emission_time=False),
                pt.record(["beam.sq.phi", "beam.pr.phi", "beam.pr_x.power", "beam.pr.power"], 3)]
        return pt.run_and_record(pipe, 1e6, 3, values={"drv.s0.dT_table": [1.0, 1.0, 1.0]}, seed=91)

    _, st, _ = rig(np.eye(2))  # first pass: realise the fiber (fixed seed), then align it with the controller
    M = np.eye(2, dtype=complex)
    for s in ids:
        M = np.array(st[f"J.{s}"]) @ M
    c = M @ X
    rec, st, _ = rig(controller(c, 0.6))
    U = controller(c, 0.6) @ M
    alpha = np.angle(U[0, 0])
    assert abs(U[1, 0]) < 1e-15  # M e_in ∝ p
    assert rec["beam.pr_x.power"][-1] < 1e-24 * rec["beam.pr.power"][-1]  # probe leaves in e_in^⊥
    sq, pr = rec["beam.sq.phi"][-1], rec["beam.pr.phi"][-1]
    assert abs(np.angle(np.exp(1j * (sq + pr - 2 * psi)))) < 1e-12  # sum = 2ψ
    assert abs(np.angle(np.exp(1j * (sq - pr - 2 * alpha)))) < 1e-12  # difference = retardance 2α
    assert alpha == pytest.approx(0.6, abs=1e-12)
