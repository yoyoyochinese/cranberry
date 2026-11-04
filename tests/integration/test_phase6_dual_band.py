"""Phase 6 WP6-C (D85-5): the two-layer lock overlay ``lock/dual_band_slow`` (P5) builds on the four future presets at
its first sweep point and its choice axes, and on F_liu20_nominal its slow layer (photon_counting_fringe on in-band
seeds → set point of the L-Liu loop, or a REF_k stretcher) closes with the right sign: it removes the static metric (ii)
offset the fast proxy leaves, windows close every T_int, the slow-loop stability check passes and the common-point
window row is suppressed by the slow loop's sensitivity."""

import importlib.util
import math
import warnings
from pathlib import Path

import pytest
import yaml

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.integration]

OV = "lock/dual_band_slow"
NOM = "lock/points/dual_band_slow_nominal"


def point(tmp, name, values):
    path = tmp / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"offline": {"params": values}}))
    return str(path)


@pytest.mark.parametrize("base", ["F_wg100", "F_wg6", "F_bor6", "F_liu20"])
def test_builds_on_every_future_preset(base):
    res = fs.build_check(fs.nominal(base), (OV,))
    assert set(res.metrics.metric_ii) == {1}
    assert "cnt.pcf1.window_index" in res.final_state.layout
    p = fs.params(res)
    assert p["dbs_f_gate_hz"] == pytest.approx(1.0 / p["t_frame_s"])


@pytest.mark.parametrize("values", [{"dbs_target": "ref_stretcher"}, {"dbs_seed_inject": "after_squeezer"}])
def test_choice_axes_build(tmp_path, values):
    res = fs.build_check(fs.nominal("F_liu20"), (OV, NOM, point(tmp_path, "pt", values)))
    p = fs.params(res)
    if "dbs_seed_inject" in values:
        assert p["dbs_in_a"] == "beam.dseedu1.phi"
    else:
        assert p["dbs_g_setpoint"] == 0.0 and p["dbs_ref_k_rad_per_v"] > 0.0


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("dbs")
    d = {"duration_s": 0.6}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {"base": run(fs.load(fs.nominal("F_liu20"), run=d)),
                "sp": run(fs.load(fs.nominal("F_liu20"), run=d, overlays=(OV, NOM))),
                "ref": run(fs.load(fs.nominal("F_liu20"), run=d,
                                   overlays=(OV, NOM, point(tmp, "ref", {"dbs_target": "ref_stretcher"}))))}


def test_slow_layer_removes_the_proxy_offset(runs):
    base = runs["base"].metrics.metric_ii[1].stats.mean
    assert abs(base) > 0.02                                    # the fast proxy alone leaves ≈ 35 mrad (nominal seed)
    for k in ("sp", "ref"):
        r = runs[k]
        assert r.valid and not [e for e in r.events if e["source"] == "chk.dbs_stability"]
        assert abs(r.metrics.metric_ii[1].stats.mean) < 0.2 * abs(base)
        st = r.final_state
        assert int(st["cnt.pcf1.window_index"]) == 12          # 0.6 s / 50 ms
        assert float(st["cnt.pcf1.x_hat"]) == pytest.approx(math.pi / 2, abs=0.05)
    p = fs.params(runs["sp"])
    assert p["dbs_cp_window_var_rad2"] < 0.05 * p["cp_window_var_rad2"]   # window row under the slow loop
    # the gates per window follow f_gate·T_int (one port-k slot per frame)
    assert float(runs["sp"].final_state["cnt.pcf1.n_gates"]) == pytest.approx(0.05 / p["t_frame_s"], abs=1.0)
