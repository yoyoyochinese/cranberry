"""Runner offline stage, topology injection and metric records (WP-F4)."""

import dataclasses
import json
import math

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
import cranberry.path.testing as pt
from cranberry.blocks.meas.testing import metric_cfg, path_blocks, refa_source, run_spec, single_port_topology, toy
from cranberry.config.spec import ConfigError
from cranberry.runner import record
from cranberry.runner.run import TOPOLOGY_KINDS, inject_topology, run

pytestmark = pytest.mark.numba


def tap_block(eta_tap=0.05):
    return {"block": "tap", "name": "tap1", "params": {
        "tap_ratio": eta_tap, "p_probe_w": 1e-3, "out": "v.tap", "segment": "sq_seg", "applies_to": ["gbs", "bhd"],
        "source": "[test]", "toggles": {"enabled": True}}}


def base_spec(explicit_topology=False, measure=None, toggles=None, sources=None):
    topo = single_port_topology()
    blocks = pt.path_pipeline(topo, {"las": "laser.las.phi"}) if explicit_topology else path_blocks({"las": "laser.las.phi"})
    pipe = [toy("toy_sine", "jit", out="drv.lo_seg.stretcher", unit="1", amp=0.01, freq_hz=1e3),
            toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
            tap_block(), *blocks,
            {"block": "metric_i", "name": "mi", "params": {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                                                          "setpoint_rad": 0.0}}]
    m = {"record": {"slots": ["beam.lo.phi"], "decim": 4}, "metrics": {"i": metric_cfg("refA.sqz")}}
    m.update(measure or {})
    return run_spec(pipe, f_sim_hz=1e6, n_samples=20_000, topology=topo,
                    offline={"blocks": sources if sources is not None else [refa_source()]}, measure=m,
                    toggles=toggles)


def test_offline_stage_registry_ledger_sources():
    res = run(base_spec())
    off = res.offline
    assert set(off.source_outputs) == {"refA.sqz"} and off.blocks[0].name == "refA.sqz"
    assert set(off.loss_registry.paths) == {"lo", "sqz"}
    assert off.loss_registry.eta("sqz", "bhd") == pytest.approx(0.95)  # the tap on sq_seg (G06)
    assert off.loss_registry.eta("lo", "bhd") == 1.0
    assert set(off.band_ledger.beams) == {"lo", "sqz"} and off.metric_names == ["i.port0"]
    assert off.band_ledger.metrics["i.port0"].f_max_hz == 1e12
    nl = res.metrics.noise_levels["i.port0"]
    assert nl["eta"] == pytest.approx(0.95) and nl["loss_path"] == "sqz" and nl["loss_target"] == "bhd"
    assert res.metrics.efficiency[0].eta_bhd == pytest.approx(0.95) and res.metrics.efficiency[0].eta_gbs is None
    assert res.metrics.metric_ii_status == "not_applicable"  # port-less config (G13)
    assert res.table_hashes.keys() == {"source.refA.sqz"}
    assert "metric.mi.port0" in res.series and res.series["metric.mi.port0"].shape == (5000,)
    assert res.metrics.metric_i[0].stats.n == 20_000  # full-rate moments despite decim 4
    assert res.metrics.sigma_rad[0] == pytest.approx(0.01 / math.sqrt(2), rel=1e-3)


def test_topology_injection_equals_explicit_route():
    a = run(base_spec(explicit_topology=False))
    b = run(base_spec(explicit_topology=True))
    assert a.series["metric.mi.port0"].tobytes() == b.series["metric.mi.port0"].tobytes()
    assert {"segment_phase", "beam_assembler", "lock_error", "metric_i", "metric_ii"} <= TOPOLOGY_KINDS
    spec = base_spec()
    inj = inject_topology(spec.pipeline, spec.topology)
    sp = [c for c in inj.children if c.name == "sp"][0]
    assert sp.params["topology"] is spec.topology
    assert "topology" not in [c for c in inj.children if c.name == "jit"][0].params
    assert "topology" not in [c for c in spec.pipeline.children if c.name == "sp"][0].params  # spec untouched
    assert inject_topology(spec.pipeline, {}) == spec.pipeline


def test_record_contains_metrics_loss_registry_band_ledger(tmp_path):
    res = run(base_spec(), record_root=tmp_path)
    names = {p.name for p in res.run_dir.iterdir()}
    assert {"metrics.json", "loss_registry.json", "band_ledger.json", "tables.json", "results.npz"} <= names
    m = json.loads((res.run_dir / "metrics.json").read_text())
    assert m["noise_levels"]["i.port0"]["eta"] == pytest.approx(0.95)
    assert m["metric_i"]["0"]["n"] == 20_000
    lr = json.loads((res.run_dir / "loss_registry.json").read_text())
    assert lr["eta"]["sqz"]["bhd"]["elements"] == ["tap1.tap"]
    bl = json.loads((res.run_dir / "band_ledger.json").read_text())
    assert bl["metrics"]["i.port0"]["f_max_hz"] == 1e12 and set(bl["beams"]) == {"lo", "sqz"}
    assert json.loads((res.run_dir / "tables.json").read_text()) == res.table_hashes
    stored = np.load(res.run_dir / "results.npz")
    assert stored["metric.mi.port0"].tobytes() == res.series["metric.mi.port0"].tobytes()
    # the record stores the spec as written (no injected topology) and reproduces the run
    again = run(record.load_record(res.run_dir))
    assert again.series["metric.mi.port0"].tobytes() == res.series["metric.mi.port0"].tobytes()


def test_offline_block_toggle_override_and_errors():
    spec = base_spec(sources=[{"block": "opa_waveguide", "name": "opa", "params": dict(
        model="lumped", alpha_per_w=9.06, pump_w=0.64, pump_ref="output", length_m=0.045, phi_opa0_rad=0.0,
        segment="sq_seg", psi_p_xtal_rad=0.0, gain_table_pump_w=None, waveguide_loss=0.02, loss_source="[test]",
        pump_rin={"f_lo_hz": 10.0, "f_hi_hz": 5e5, "variance_outside": 1e-8},
        toggles=dict(pump_loss=False, phase_mismatch=False, pump_depletion=False, pump_rin=True))}],
        measure={"metrics": {"i": metric_cfg("opa")}})
    on = run(spec).offline
    assert on.loss_registry.eta("sqz", "bhd") == pytest.approx(0.95 * 0.98)  # tap + waveguide (lumped)
    assert [d.quantity for d in on.band_ledger.decls] == ["rin"]
    off = run(dataclasses.replace(spec, toggles={"opa.pump_rin": False})).offline  # offline block via RunSpec.toggles
    assert off.band_ledger.decls == []
    with pytest.raises(ConfigError, match="unique"):
        run(base_spec(sources=[refa_source(), refa_source()]))
    with pytest.raises(ConfigError, match="offline.blocks"):
        run(base_spec(sources=[{"block": "single_mode_reference", "name": "bad", "params": {"loss": 0.1}}]))
    with pytest.raises(KeyError, match="no offline source"):
        run(base_spec(measure={"metrics": {"i": metric_cfg("nope")}}))
    with pytest.raises(KeyError, match="f_max_hz"):
        run(base_spec(measure={"metrics": {"i": {"source": "refA.sqz"}}}))


def test_run_without_metrics_or_topology_is_unchanged():
    spec = run_spec([toy("toy_sine", "s", out="v.s", unit="V", amp=1.0, freq_hz=1e3)], f_sim_hz=1e5, n_samples=100,
                    measure={"record": {"slots": ["v.s"], "decim": 1}})
    res = run(spec)
    assert res.metrics is None and res.table_hashes == {} and res.offline.record_parts() == {}
    assert res.summary()["v.s.rms"] == pytest.approx(np.sqrt(np.mean(res.series["v.s"] ** 2)))
