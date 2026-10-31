"""VT-F2 (WP-B2 part): every WP-B2 block in one pipeline through ``runner.run`` — bit-identical
realisation for any ``chunk_samples`` (G11, PLAN §6.4), compiled = reference (1e−12, VT-44), and a
``RunSpec.toggles`` override that switches one noise source off leaves every other source's
realisation unchanged (G11, G13)."""

import numpy as np
import pytest

import cranberry.blocks.noise.testing as nt
import cranberry.path.testing as pt
from cranberry.config.spec import RunSettings, RunSpec
from cranberry.core.pipeline import spec_from_obj
from cranberry.runner.run import run

pytestmark = pytest.mark.numba
FS = 1.0e6
TOPO = nt.seg_topo({"a": (5.0, 1, "uniform"), "b": (9.0, 4, "per_sub"), "xtal": (0.045, 1, "uniform")},
                   beams=([pt.beam_spec("lo", [pt.step_spec("a", 1)], "det", source="dx1"),
                           pt.beam_spec("pr", [pt.step_spec("b", 1), pt.step_spec("fc", 1)], "det", source="dx1")],
                          ["det"]))
TOPO["segments"].append(pt.seg_spec("fc", 0.3))
REC = ["laser.dx1.phi", "laser.dx1.power", "drv.a.dT", "drv.a.strain", "opa.xt.dk", "v.voff", "psi.fc.pr",
       "beam.lo.phi", "beam.pr.phi"]


def pipeline():
    fc = {"block": "filter_cavity", "name": "fc1", "params": {
        "segment": "fc", "beam": "pr", "dnu_fc_hz": 1e5, "t0": 0.95, "source": "test", "in_phase": "laser.dx1.phi",
        "toggles": {"lock_residual": False, "laser_lowpass": True}}}
    seg_phase = pt.path_pipeline(TOPO, {"dx1": "laser.dx1.phi"})
    seg_phase[0]["params"]["segments"] = ["a", "b", "xtal"]
    return [nt.laser(), nt.thermal(TOPO, ["a", "b"], {"a": 3.0, "b": 30.0}, rate_div=7), nt.crystal(), nt.pd_offset(),
            nt.mech(TOPO, ["a", "b"], spatial_coherence={"b": 0.25}), nt.thermo(TOPO, ["a", "b"]), nt.gawbs(),
            fc, *seg_phase]


def spec(chunk, n=9000, toggles=None):
    return RunSpec(run=RunSettings(duration_s=n / FS, f_sim_hz=FS, chunk_samples=chunk, seed=20261031),
                   pipeline=spec_from_obj(pipeline()), toggles=toggles or {},
                   measure={"record": {"slots": REC, "decim": 3}})


def test_chunk_length_invariance_bitwise():
    base = run(spec(65536))
    assert base.mode == "compiled" and np.std(base.series["beam.pr.phi"]) > 0
    for chunk in (1, 997, 4096, 4097):
        res = run(spec(chunk))
        for k in REC:
            np.testing.assert_array_equal(res.series[k], base.series[k], err_msg=f"chunk={chunk} slot={k}")
        for name in base.final_state.layout.slot_names():
            if not name.endswith((".buf", ".idx")) and not name.startswith("rec."):
                np.testing.assert_array_equal(res.final_state[name], base.final_state[name], err_msg=name)


def test_reference_runner_matches_compiled():
    a = run(spec(1000, n=2000), mode="compiled")
    b = run(spec(999, n=2000), mode="reference")
    for k in REC:
        np.testing.assert_allclose(b.series[k], a.series[k], rtol=1e-12, atol=1e-300)


def test_toggle_override_leaves_other_sources_unchanged():
    on = run(spec(4096, n=4000))
    off = run(spec(4096, n=4000, toggles={"las.fm_flicker": False, "mn.peaks": False}))
    for k in ("drv.a.dT", "opa.xt.dk", "v.voff", "laser.dx1.power"):
        np.testing.assert_array_equal(off.series[k], on.series[k], err_msg=k)
    assert not np.array_equal(off.series["laser.dx1.phi"], on.series["laser.dx1.phi"])
    assert not np.array_equal(off.series["drv.a.strain"], on.series["drv.a.strain"])
