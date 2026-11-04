"""``viz.diagram``: stage / cluster classification, lock tracing, DOT output, rendering and the CLI (D85-13)."""

import pytest

from cranberry import __main__ as cli
from cranberry import config
from cranberry.config.overlay import REPO_ROOT
from cranberry.viz import diagram as dg

REFA = "presets/reference/refA.yaml"
FWG6 = "presets/future/F_wg6_nominal.yaml"


@pytest.fixture(scope="module")
def refa():
    return config.load(REPO_ROOT / REFA)


@pytest.fixture(scope="module")
def fwg6():
    return config.load(REPO_ROOT / FWG6, overlays=["demux/none"])


def test_leaves_stages_and_groups(fwg6):
    by = {lf.name: lf for lf in dg.pipeline_leaves(fwg6)}
    assert by["laser"].stage == "source"
    assert by["env.mech.ports"].group == "noise" and by["path.lock"].group == "bookkeeping"
    assert by["mi"].group == "bookkeeping" and by["liu0"].stage == "detection"
    assert by["path.phase"].stage == "optics" and by["pzt0"].stage == "control"
    assert by["sup0"].stack == "pipeline[23]" and by["tec1"].stack == "pipeline[27]"  # + the in_bw lpfs (G-3)


def test_lock_drives_trace_detector_and_lock(refa, fwg6):
    rows = {(r["actuator"], r["segment"], r["driver"]): r for r in dg.lock_drives(fwg6)}
    assert rows[("pzt0", "F0", "stretcher")]["lock"] == "L0" and rows[("tec0", "F0", "tec")]["detectors"] == ["liu0"]
    assert rows[("tec1", "F1", "tec")]["lock"] == "L1"
    assert rows[("dmxsw", "DMX0", "eom")]["lock"] is None  # switch drive: no detector upstream
    a = {r["actuator"]: r for r in dg.lock_drives(refa)}
    # through demodulator (err family, no beam read) to the beat detector
    assert a["eom1"]["lock"] == "lock1" and a["eom1"]["detectors"] == ["pd1"]
    assert a["eom2"]["lock"] == "lock2" and a["scan_act"]["lock"] is None


def test_signal_flow_collapse_expand_and_feedback(fwg6):
    dot = dg.signal_flow_dot(fwg6)
    assert '"grp:noise"' in dot and '"grp:bookkeeping"' in dot and '"b:env.mech.ports"' not in dot
    assert 'cluster_stack_pipeline[23]' in dot and "cluster_offline" in dot
    # the actuators write drv.F0.* read by segment_phase earlier in the pipeline: feedback edges
    assert '"b:pzt0" -> "b:path.phase"' in dot and "constraint=false" in dot
    full = dg.signal_flow_dot(fwg6, expand=["all"])
    assert '"b:env.mech.ports"' in full and '"grp:noise"' not in full and '"b:mi"' in full
    with pytest.raises(ValueError, match="unknown group"):
        dg.signal_flow_dot(fwg6, expand=["typo"])


def test_optical_layout_beams_couplers_locks(fwg6):
    dot = dg.optical_layout_dot(fwg6)
    assert '"src:laser" -> "seg:PUMP"' in dot and '"seg:G0" -> "obs:gbs"' in dot
    assert '"lock:L0" -> "act:pzt0"' in dot and '"act:pzt0" -> "seg:F0"' in dot
    assert '"obs:det0" -> "lock:L0"' in dot
    assert "peripheries=2" in dot  # combiner segments


def test_render_and_cli(tmp_path):
    try:
        dg.find_dot()
    except FileNotFoundError:
        pytest.skip("Graphviz dot not installed")
    rc = cli.main(["diagram", "refD_nominal", "--out", str(tmp_path), "-f", "svg", "-f", "dot"])
    assert rc == 0
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == ["refD_nominal.optical_layout.dot", "refD_nominal.optical_layout.svg",
                     "refD_nominal.signal_flow.dot", "refD_nominal.signal_flow.svg"]
    assert (tmp_path / "refD_nominal.signal_flow.svg").read_text().lstrip().startswith("<?xml")
    assert cli.main(["diagram", "no_such_preset", "--out", str(tmp_path)]) == 2
