"""WP-I2: the refD preset (D01, D02, D35, D49, D51-8/9/11, D52, D53, D55-1; PLAN §2.9, §5.2).

* every quantity the documents leave open is a named sweep axis under ``offline.params`` (D34), and
  ``refD_nominal`` resolves all of them (single runs);
* device values reach the blocks as specified: DSC1 entries → canonical gains through the F5 mapping
  (D53, p.20 units), CV limits ±4 V, 915B 8669-1 6.49 rad/V, 914 two poles, PDB450C 10³ V/A, loop
  delay τ_L − 1/f_sim − explicit discretisation delays (D55-1), DAC load = oscilloscope ‖ 914 (D53);
* the sweep axes are live: sweeping refD itself (other axes pinned to the nominal point) changes the
  arm lengths, the laser decorrelation delay and γ_coh consistently;
* YAML round trip (VT-F1 form) and a compiled closed-loop run that locks at +π/2.
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import loads_yaml
from cranberry.hardware.devices import device_values, load_device
from cranberry.hardware.mapping import dsc1_to_canonical
from cranberry.runner.run import run
from cranberry.runner.sweep import expand

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parent / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = pytest.mark.numba

REFD = "presets/reference/refD.yaml"


def test_refD_unknowns_are_named_sweep_axes():
    spec = config.load(REFD)
    markers = config.sweep_markers(spec)
    axes = [k for k in markers if k.startswith("offline.params.")]
    deferred = [k for k in markers if not k.startswith("offline.params.")]
    assert len(axes) >= 30
    # everything else that is still open is an expression of those axes (deferred until a point is fixed)
    from cranberry.config.yaml_io import Expr

    assert all(isinstance(markers[k], Expr) for k in deferred)
    for name in ("dL_m", "h0_hz2_per_hz", "f_res_hz", "q_res", "latency_s", "r_914_ohm", "loop_sign",
                 "coupler_sign", "K_T_rad_per_m_K", "polarizer_er_db", "kappa2", "h_ref_rad2_hz"):
        assert f"offline.params.{name}" in axes, name
    nominal = config.load(rs.NOMINAL)
    assert config.sweep_markers(nominal) == {}
    assert set(axes) <= {f"offline.params.{k}" for k in nominal.offline["params"]} | {
        f"offline.params.{k}" for k in ("coupler_excess_db", "il_915b_db", "polarizer_er_db")}


def test_refD_device_values_and_mappings():
    spec = config.load(rs.NOMINAL)
    pid = rs.block_params(spec, "dsc1.pid")
    dev = device_values(load_device("dsc1"))
    canon, _ = dsc1_to_canonical({"p": 0.000241, "i": 0.001406, "d": 0.01108, "p_enabled": True,
                                  "i_enabled": True, "d_enabled": True}, dev)
    assert pid["kp"] == pytest.approx(canon.kp, rel=1e-15) == -2.41e-4
    assert pid["ki_per_s"] == pytest.approx(canon.ki_per_s, rel=1e-12) and pid["ki_per_s"] == pytest.approx(-703.0)
    assert pid["kd_s"] == pytest.approx(canon.kd_s, rel=1e-12) and pid["kd_s"] == pytest.approx(-2.216e-8)
    assert (pid["u_min_v"], pid["u_max_v"]) == (-4.0, 4.0)  # D35
    leaves, params = rs.built_leaves(spec)
    assert leaves["dsc1.pid"].values["discretisation"] == "dsc1" and leaves["dsc1.pid"].rate_div == 2
    assert leaves["dsc1.adc"].n_avg() == 2 and leaves["dsc1.adc"].lsb_v == pytest.approx(8 / 2**16)
    assert leaves["dsc1.dac"].lsb_v == pytest.approx(8 / 2**12) and leaves["dsc1.dac"].rate_div == 2
    assert leaves["dsc1.dac"].gain() == pytest.approx(5e5 / (5e5 + 220.0))  # 1 MΩ scope ‖ 1 MΩ 914 (D53)
    assert leaves["stretcher915b"].gain_rad_per_v(params) == pytest.approx(6.491, abs=5e-4)
    assert list(rs.block_params(spec, "amp914")["poles_hz"]) == [2.95e4, 1.12e5]  # D51-11
    tau = rs.block_params(spec, "dsc1.latency")["tau_s"]
    assert tau == pytest.approx(5.7e-6 - 1e-6 - 1.5e-6, abs=1e-15)  # D55-1 minus the explicit avg + ZOH
    dcf = rs.block_params(spec, "pdb450c")
    assert dcf["gain_v_per_a"] == 1e3 and dcf["v_sat_v"] == 4.6 and dcf["p_in_w"] == 3.33e-3
    assert (dcf["mu_abs"], dcf["mu_phase_rad"]) == (1.0, 0.0)  # μ from per_fringe (B7 contract)
    a1, a2 = rs.segment(spec, "a1"), rs.segment(spec, "a2")
    assert a1["L_m"] - a2["L_m"] == pytest.approx(9.0)  # ΔL centre (D53)
    assert a2["L_m"] == pytest.approx(2.1)


def _pinned_axes(**free):
    """Axes of refD pinned to the refD_nominal values, except ``free`` {name: values}."""
    nominal = config.load(rs.NOMINAL).offline["params"]
    refd = config.load(REFD)
    out = {}
    for path in config.sweep_markers(refd):
        if not path.startswith("offline.params."):
            continue
        name = path.split(".")[-1]
        lin = name[: -len("_db")] + "_lin" if name.endswith("_db") else name
        if name in free:
            out[path] = list(free[name])
        else:
            v = nominal[name] if name in nominal else 10 * math.log10(nominal[lin])
            out[path] = [v]
    return refd, out


def test_refD_sweep_axes_are_live():
    refd, axes = _pinned_axes(dL_m=[5.0, 15.0])
    points = expand(refd, axes)
    assert len(points) == 2
    for pt, dl in zip(points, (5.0, 15.0)):
        sp = pt.spec
        assert config.sweep_markers(sp) == {}
        a1, a2 = rs.segment(sp, "a1"), rs.segment(sp, "a2")
        assert a1["L_m"] - a2["L_m"] == pytest.approx(dl)
        tau = rs.block_params(sp, "dx1")["oob_pairs"][0]["tau_s"]
        assert tau == pytest.approx(1.468 * dl / rs.C_M_PER_S)
        h0 = sp.offline["params"]["h0_hz2_per_hz"]
        sig2 = 2 * math.pi**2 * h0 * tau - 4 * math.pi**2 * h0 * tau**2 * 5e5
        assert rs.block_params(sp, "pdb450c")["gamma_coh"] == pytest.approx(math.exp(-sig2 / 2), rel=1e-12)
        coef = next(c for c in sp.topology["coef"] if c["segment"] == "a1" and c["driver"] == "dT"
                    and c["effect"] == "phase")
        assert coef["value"] == pytest.approx(44.0 * a1["L_m"])
    # the pinned grid reproduces refD_nominal exactly at its own point
    _, axes9 = _pinned_axes(dL_m=[9.0])
    sp9 = expand(refd, axes9)[0].spec
    nominal = config.load(rs.NOMINAL)
    assert sp9.topology == nominal.topology
    assert [b.params for b in sp9.pipeline.children] == [b.params for b in nominal.pipeline.children]


def test_refD_yaml_round_trip():
    for path in (rs.NOMINAL, REFD):  # resolved point and the sweep spec with its deferred expressions
        spec = config.load(path)
        text = config.dump_yaml(spec)
        again = config.load(loads_yaml(text))
        assert again == spec
        assert config.dump_yaml(again) == text


def test_refD_nominal_locks_at_plus_half_pi():
    spec = rs.load(run={"duration_s": 0.1})
    res = run(spec)
    assert res.mode == "compiled"
    e = res.series["lock.lockD.err"]
    tail = e[e.size // 2:]
    st = res.final_state
    assert int(st["lock.dsc1.sup.n_loss"]) == 0 and int(st["lock.dsc1.sup.n_slip"]) == 0
    assert abs(float(np.mean(tail))) < 0.02  # PER lock-point shift ≤ arcsin(x/(1−x)) ≈ 8 mrad (x ≈ 7.9e−3)
    # in-band laser × ΔL dominates the residual: (2πτ)² h0 · f_N for f ≫ f_c (T15 §6.1: tens of mrad)
    tau = 1.468 * 9.0 / rs.C_M_PER_S
    assert 0.5 < float(np.std(tail)) / math.sqrt((2 * math.pi * tau) ** 2 * 2.546e4 * 5e5) < 1.2
    assert float(res.series["v.dsc1.cv"][-1]) == pytest.approx(-(math.pi / 2) / 6.488, abs=0.05)


def test_block_additions_lin_names_match_db_names():
    """WP-I2 block additions: ``er_lin`` / ``il_lin`` (polarizer_projection) and ``insertion_loss_lin``
    (B4 actuators) give exactly what ``er_db`` / ``il_db`` / ``insertion_loss_db`` give (D56-1: a YAML
    preset's ``*_db`` keys arrive renamed to ``*_lin``)."""
    from cranberry.core.compile import run_reference
    from cranberry.core.layout import SlotLayout
    from cranberry.core.pipeline import build, spec_from_obj
    from cranberry.core.state import Params, State

    def polarizer(**kw):
        p = {"in": "v.e", "out_jones": "v.oj", "out_dp": "v.dp", "segment": "s", "applies_to": ["probe"],
             "source": "test", "toggles": {"pol_drift_before_polarizer": True, "pol_extinction": True}, **kw}
        lay = SlotLayout(1e6)
        lay.add("v.e", (2,), np.complex128, "1")
        blk = build(spec_from_obj([{"block": "polarizer_projection", "name": "pol", "params": p}]), lay, 0,
                    external=("v.e",))
        st, pr = State(lay), Params(lay)
        st["v.e"] = np.array([np.cos(0.3), np.sin(0.3)], dtype=complex)
        run_reference(blk, st, pr, 1)
        return complex(st["v.oj"][1]), float(st["v.dp"]), blk.losses()[0].eta

    a = polarizer(er_db=28.0, il_db=0.3)
    b = polarizer(er_lin=10 ** 2.8, il_lin=10 ** 0.03)
    assert a[0] == pytest.approx(b[0], rel=1e-13) and a[1] == pytest.approx(b[1], rel=1e-13)
    assert a[2] == pytest.approx(b[2], rel=1e-13)

    from cranberry.core.registry import make_block
    from cranberry.core.block import BlockSpec

    def eta(**kw):
        p = {"in": "v.u", "gain_mode": "direct", "k_rad_per_v": 1.0, "polarity": 1.0, "modes": [],
             "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -5.0, "drive_max_v": 5.0,
             "hysteresis_frac": 0.0, "segment": "a1", "driver": "stretcher", "n_sub": 1, "sub": [0, 1],
             "ring_len": 0, "loss_applies_to": ["probe"], "loss_source": "test",
             "toggles": {"resonance": False, "hysteresis": False, "saturation": False, "insertion_loss": True}, **kw}
        els = make_block(BlockSpec("pzt_fiber_stretcher", "s", p)).losses()
        return [e.eta for e in els]

    assert eta(insertion_loss_db=1.0) == pytest.approx(eta(insertion_loss_lin=10 ** 0.1), rel=1e-13)
    assert eta(insertion_loss_db=None) == eta(insertion_loss_lin=None) == []
