"""Phase 5 M1 framework review: chunk-length bit-identity and compiled == reference for the F6 path
blocks with decimated writers / readers (rate_div 2 / 3 / 5 mixed), a per_sub white driver, a
counter-propagating beam, a jones segment with J_nr and a fractional delay_buffer.

Run: uv run python reviews/scripts/f1_chunk_multirate_path.py
"""

import numpy as np

import cranberry.path.testing as pt
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State

FS = 3.3e6
N = 20011


def pipeline():
    segs = [pt.seg_spec("c", 40.0, n_sub=7, spatial="distributed", drivers={"strain": "per_sub"}),
            pt.seg_spec("j", 3.0, drivers={"dT": "uniform"}, pol_model="jones",
                        pol_params={"J": [[0.6, 0.8], [-0.8, 0.6]], "J_nr_slot": "J.nr"})]
    beams = [pt.beam_spec("f", [pt.step_spec("c", 1), pt.step_spec("j", 1)], "o", jones_out=[1, 0]),
             pt.beam_spec("b", [pt.step_spec("j", -1), pt.step_spec("c", -1)], "o", jones_out=[1, 0],
                          freq_offset_hz=1.1e6)]
    topo = pt.topo_spec(segs, beams, ["o"], pt.coef_spec("c", "strain", "*", 1.0)
                        + pt.coef_spec("j", "dT", "*", 2.0, 1e-12))
    sp, ba = pt.path_pipeline(topo, {"las": "laser.las.phi"})
    sp["rate_div"] = 1
    ba["rate_div"] = 1
    white = pt.pathtest("pathtest_white", "w", "drv.c.strain", (7,), "1", sigma=1.0)
    white["rate_div"] = 2
    temp = pt.pathtest("pathtest_white", "t", "drv.j.dT", (), "K", sigma=0.1)
    temp["rate_div"] = 3
    las = pt.pathtest("pathtest_white", "las", "laser.las.phi", (), "rad", sigma=0.3)
    dly = {"block": "delay_buffer", "name": "d", "rate_div": 5,
           "params": {"in": "beam.f.phi", "out": "delay.d.out", "tau_s": 2.37e-6, "method": "thiran1", "unit": "rad"}}
    return [white, temp, las, sp, ba, dly,
            pt.record(["beam.f.phi", "beam.b.phi", "psi.c.f", "psi.c.b", "delay.d.out"], N, name="rec")]


def run(chunks, compiled):
    lay = SlotLayout(FS)

    lay.add("J.nr", (2, 2), np.complex128, "1")  # non-reciprocal factor N, written from outside
    blk = build(spec_from_obj(pipeline()), lay, 7, external=("J.nr",))
    st, pa = State(lay), Params(lay)
    c, s_ = np.cos(0.4), np.sin(0.4)
    st["J.nr"] = np.array([[c, -s_], [s_, c]]) * np.exp(0.2j)  # Faraday-like rotation
    done, i = 0, 0
    while done < N:
        m = min(chunks[i % len(chunks)], N - done)
        (run_compiled if compiled else run_reference)(blk, st, pa, m)
        done += m
        i += 1
    return np.array(st["rec.rec.buf"])


if __name__ == "__main__":
    base = run([N], True)
    for chunks in ([1], [7, 4096, 3], [4097], [5000, 11]):
        other = run(chunks, True)
        print(f"compiled chunks {chunks}: bit-identical = {np.array_equal(other, base)}")
    ref = run([4096, 333], False)
    rel = np.max(np.abs(ref - base)) / np.max(np.abs(base))
    print(f"reference vs compiled: max relative difference {rel:.2e}")
