"""Optional ``mask_slot`` (WP6-G, D85-3) of ``metric_i`` / ``metric_ii``, ``measure.stray_light`` and
``runner.decompose``: absent = unchanged; present = evaluated where mask > 0.5 only."""

import warnings
from types import SimpleNamespace

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.blocks.meas.metrics import PhaseStats, evaluate_metrics
from cranberry.blocks.meas.testing import (
    metric_cfg,
    path_blocks,
    refa_source,
    run_spec,
    single_port_topology,
    toy,
)
from cranberry.runner.decompose import validity_violations
from cranberry.runner.run import run

pytestmark = pytest.mark.numba
N = 3000


def spec(mask: bool):
    gate = {"block": "gate_schedule", "name": "gw",
            "params": {"detector_gate": "complement", "write_phase": False, "period_s": 100e-6, "lock_fraction": 0.655,
                       "phase": 0.0, "edge_s": 0.0, "timing_jitter_s": 0.0, "leak": 0.0,
                       "toggles": {"gating": True, "jitter": False, "edges": False, "leak": False}}}
    mp = {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}}, "setpoint_rad": 0.0, "t_settle_s": 2e-4}
    if mask:
        mp["mask_slot"] = "gate.gw.meas"
    pipe = [toy("toy_sine", "jit", out="drv.lo_seg.stretcher", unit="1", amp=0.05, freq_hz=7.3e3),
            toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0), gate,
            *path_blocks({"las": "laser.las.phi"}), {"block": "metric_i", "name": "mi", "params": mp}]
    return run_spec(pipe, f_sim_hz=1e6, n_samples=N, topology=single_port_topology(), chunk_samples=1000,
                    offline={"blocks": [refa_source()]},
                    measure={"record": {"slots": ["gate.gw.meas"], "decim": 1, "antialias": "none"},
                             "metrics": {"i": metric_cfg("refA.sqz")}})


@pytest.mark.parametrize("mode", ["reference", "compiled"])
def test_metric_moments_and_series_inside_the_mask(mode):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full, masked = run(spec(False), mode=mode), run(spec(True), mode=mode)
    d = np.asarray(full.series["metric.mi.port0"])
    t = np.asarray(full.t_s)
    keep = (t >= 2e-4) & (np.asarray(full.series["gate.gw.meas"]) > 0.5)
    m = masked.metrics.metric_i[0]
    assert m.stats.n == keep.sum() == 34 * 28  # 34 samples per 100 µs window, 28 periods after the settle time
    np.testing.assert_allclose(m.delta_rad, d[keep], rtol=0, atol=1e-15)
    ref = PhaseStats.from_acc(np.array([keep.sum(), d[keep].sum(), (d[keep] ** 2).sum(), np.cos(2 * d[keep]).sum(),
                                        np.sin(2 * d[keep]).sum()]))
    assert m.stats.var == pytest.approx(ref.var, rel=1e-10)
    assert full.metrics.metric_i[0].stats.n == (t >= 2e-4).sum()  # no mask: every sample after t_settle_s
    assert masked.metrics.noise_levels["i.port0"]["eta"] == full.metrics.noise_levels["i.port0"]["eta"]


def test_stray_light_mask():
    series = {"leak.c": np.array([1.0, 2.0, 3.0, 4.0]), "gate.gw.meas": np.array([0.0, 1.0, 1.0, 0.0])}
    kw = dict(leaves=[], state=None, t_s=np.arange(4.0), series=series, topology=None, source_outputs={},
              loss_registry=None, band_ledger=None)
    plain = evaluate_metrics(**kw, measure={"stray_light": {"1": {"clicks_per_gate": "leak.c", "photons_per_mode": 0.5}}})
    masked = evaluate_metrics(**kw, measure={"stray_light": {"1": {"clicks_per_gate": "leak.c", "photons_per_mode": 0.5,
                                                                    "mask_slot": "gate.gw.meas"}}})
    assert plain.stray_light[1].clicks_per_gate == 2.5 and masked.stray_light[1].clicks_per_gate == 2.5 * 1.0
    assert masked.stray_light[1].photons_per_mode == 0.5
    series["leak.c"] = np.array([1.0, 2.0, 4.0, 8.0])
    assert evaluate_metrics(**kw, measure={"stray_light": {"1": {"clicks_per_gate": "leak.c", "photons_per_mode": 0.5,
                                                                  "mask_slot": "gate.gw.meas"}}}).stray_light[1].clicks_per_gate == 3.0


def test_decompose_validity_inside_the_mask():
    res = SimpleNamespace(series={"act.sat": np.array([0, 1, 0, 0]), "gate.gw.meas": np.array([0.0, 0.0, 1.0, 1.0])},
                          final_state=None)
    v = [{"slot": "act.sat", "abs_max": 0.5}]
    assert validity_violations(res, v) and not validity_violations(res, v, "gate.gw.meas")
    res.series["act.sat"] = np.array([0, 0, 0, 1])
    assert "inside gate.gw.meas" in validity_violations(res, v, "gate.gw.meas")[0]
