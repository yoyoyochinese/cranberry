"""WP-I5 backward-compatible block additions (each checked against the behaviour it extends):

* ``demux_design``: ``er_lin`` (the name a YAML config loads ``er_db`` under, D56-1 / D61-5) and ``segments: null``
  + ``port_segments`` (register only the represented ports, G13);
* ``measure.metrics.*.lo_mode: optimal`` (the source's optimal LO mode; null for single-mode / CW bases);
* topology ``coef`` entries accept an informational ``id`` (keyed overlay merges).
"""

import math

import numpy as np
import pytest

from cranberry.blocks.dmx.testing import design
from cranberry.blocks.meas.metrics import _lo_mode
from cranberry.blocks.source import temporal_modes as tm
from cranberry.blocks.source.pulsed_testing import LIU, pulsed_opo
from cranberry.config.overlay import merge
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.offline.source_output import SourceOutput
from cranberry.path.topology import Topology
import cranberry.path.testing as pt


def _make(obj):
    return make_block(BlockSpec.from_obj(obj))


def test_demux_er_lin_equals_er_db():
    a = _make(design())
    p = dict(design()["params"])
    p["er_lin"] = 10 ** (p.pop("er_db") / 10)
    b = _make({"block": "demux_design", "name": "dmx", "params": p})
    np.testing.assert_allclose(a.matrix(), b.matrix(), rtol=1e-14, atol=0)
    with pytest.raises(ValueError, match="exactly one of er_db, er_lin"):
        _make({"block": "demux_design", "name": "dmx", "params": {**p, "er_db": 23.0}})


def test_demux_port_segments_register_selected_ports():
    full = _make(design(layout="chain"))
    p = dict(design(layout="chain")["params"], segments=None, port_segments={"DMX0": 0, "DMX1": 7})
    sel = _make({"block": "demux_design", "name": "dmx", "params": p})
    els = {e.segment: e.eta for e in sel.losses()}
    ref = {e.segment: e.eta for e in full.losses()}
    assert els == {"DMX0": ref["D0"], "DMX1": ref["D7"]}
    with pytest.raises(ValueError, match="port_segments"):
        _make({"block": "demux_design", "name": "dmx", "params": {**p, "port_segments": {"A": 0, "B": 0}}})


def test_lo_mode_optimal():
    so = pulsed_opo(LIU).source_output
    u, _ = tm.optimal_lo_mode(so)
    np.testing.assert_array_equal(_lo_mode("optimal", so), u)
    single = SourceOutput.single_mode_reference(loss=0.0, alpha_per_w=1.0, pump_w=1.0, phi_opa0_rad=0.0)
    assert _lo_mode("optimal", single) is None
    assert _lo_mode(None, single) is None
    with pytest.raises(ValueError, match="optimal"):
        _lo_mode("best", so)


def test_coef_id_is_ignored_and_merges_by_key():
    seg = pt.seg_spec("c", 1.0, drivers={"dT": "uniform"})
    beams = [pt.beam_spec("sq", [pt.step_spec("c")], "o")]
    coef = [{"id": "c.dT.all.phase", "segment": "c", "driver": "dT", "beam": "*", "effect": "phase", "value": 2.0},
            {"id": "c.dT.all.group_delay", "segment": "c", "driver": "dT", "beam": "*", "effect": "group_delay",
             "value": 0.0}]
    t1 = Topology.from_spec(pt.topo_spec([seg], beams, ["o"], coef))
    t2 = Topology.from_spec(pt.topo_spec([seg], beams, ["o"], [{k: v for k, v in e.items() if k != "id"} for e in coef]))
    assert t1.coef == t2.coef
    merged = merge({"coef": coef}, {"coef": [{"id": "c.dT.all.phase", "value": 3.0}]})["coef"]
    assert [e["value"] for e in merged] == [3.0, 0.0] and math.isclose(merged[0]["value"], 3.0)
