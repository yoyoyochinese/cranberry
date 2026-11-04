"""Phase 6 P7 in-lock identification (D85-7) through the full runner.

* refD_nominal + ``ident/wire_refD`` + ``ident/in_lock_stepped``, every noise source off: L̂ from one locked run equals
  the exact multirate z-domain loop gain at the PID instants (``refD_support.loop_parts``; the loop is LTI at 500 kHz
  there) to 1e−4 / 0.01°, and the ``loop_stability_check`` model (D82 dsc1) within 3 % / 3° from 1 to 30 kHz.
* refD_nominal with its noise sources, through ``SimController.measure_loop_gain`` (the abstract-interface path):
  coherence ≥ 0.98 and the model within 6 % / 4° from 1 to 25 kHz at 0.1 V (the in-loop phase noise lowers the
  fringe slope by ⟨cos δ⟩ ≈ 2 %).
* F_wg6_nominal + ``ident/p7_literal`` (fast stage open, TEC rewired to the error): L̂ = −P̂·C₁ equals the model's
  fast path within 6 % / 5° up to 20 kHz.
* Every point overlay builds with both waveform overlays on its presets.
"""

import dataclasses
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import yaml

from cranberry import config
from cranberry.hardware.sim import SimController
from cranberry.core.registry import make_block
from cranberry.runner.run import inject_topology, offline_block_specs, run

_sp = importlib.util.spec_from_file_location("refD_support", Path(__file__).parent / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.integration, pytest.mark.numba]

NOISE_OFF = {**rs.NOISE_OFF, "pdb450c.per_crosscoupling": False, "pdb450c.visibility_hf_averaging": False}


def point(tmp, name, values):
    path = tmp / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"offline": {"params": values}}))
    return str(path)


def test_refD_overlay_noise_off_equals_z_domain_and_model(tmp_path):
    # 0.02 V: phase excursion ≤ 0.02 rad, fringe (sine) bias ≤ A²/8 ≈ 5e−5 (at 0.1 V it reaches 1e−3 at 18 kHz)
    pt = point(tmp_path, "pt", {"ident_n_freqs": 7, "ident_f_lo_hz": 1.0e3, "ident_f_hi_hz": 3.0e4,
                                "ident_amplitude": 0.02})
    spec = config.load(rs.NOMINAL, overlays=["ident/wire_refD", "ident/in_lock_stepped", pt])
    spec = type(spec)(**{**spec.__dict__, "toggles": {**spec.toggles, **NOISE_OFF}})
    res = run(spec)
    m = res.measurements["loop_identification"]["lockD"]
    f = m["f_hz"]
    assert f.size == 7 and f[0] == pytest.approx(1.0e3) and f[-1] == pytest.approx(3.0e4)
    fr = rs.fringe(spec)
    k_d = -fr["sign"] * fr["v0"] * fr["vis"]
    a, b, c, g = rs.loop_parts(spec, f, k_d=k_d)
    a2, b2, _, _ = rs.loop_parts(spec, f + spec.run.f_sim_hz / 2, k_d=k_d)
    lz = -g * c * (b * a + b2 * a2)
    np.testing.assert_allclose(np.abs(m["l_hat"] / lz), 1.0, atol=1e-4)
    assert np.max(np.abs(np.degrees(np.angle(m["l_hat"] / lz)))) < 0.01
    assert np.max(np.abs(m["mag_err"])) < 0.03 and np.max(np.abs(m["phase_err_deg"])) < 3.0
    assert m["f_u_hz"] == pytest.approx(m["model_f_u_hz"], rel=0.01)
    assert m["phase_margin_deg"] == pytest.approx(m["model_phase_margin_deg"], abs=1.0)


def test_refD_measure_loop_gain_with_noise():
    spec = config.load(rs.NOMINAL, overlays=["ident/wire_refD", "ident/in_lock_stepped"])
    ctl = SimController(spec)
    ctl.configure({"loops": {"lockD": {"pid_block": "dsc1.pid", "device": "dsc1"}},
                   "excitations": {"err": {"block": "ident.exc", "loop": "lockD", "relation": "minus_l",
                                           "reference": "v.ident.ref", "response": "v.dsc1.adc", "controller": None,
                                           "model": "chk.lockD_stability"}}})
    ctl.engage()
    out = ctl.measure_loop_gain("lockD", "stepped_sine", [1.0e3, 2.1e3, 4.3e3, 8.9e3, 18.35e3, 25e3], amplitude=0.1,
                                settle_s=5e-3, dwell_s=25e-3, t_start_s=0.02, amp_f_corner_hz=2.0e3, n_blocks=5)
    assert np.all(out["coherence"] >= 0.98)
    assert np.max(np.abs(out["mag_err"])) < 0.06 and np.max(np.abs(out["phase_err_deg"])) < 4.0
    assert out["f_u_hz"] == pytest.approx(out["model_f_u_hz"], rel=0.05)


def test_F_p7_literal_plant_matches_fast_path(tmp_path):
    pt = point(tmp_path, "pt", {"ident_n_freqs": 6, "ident_f_lo_hz": 500.0, "ident_f_hi_hz": 2.0e4})
    spec = config.load("presets/future/F_wg6_nominal.yaml", overlays=["ident/p7_literal", "ident/in_lock_stepped", pt])
    res = run(spec)
    m = res.measurements["loop_identification"]["lock0"]
    assert m["relation"] == "plant" and m["f_hz"].size == 6
    assert np.all(m["coherence"] > 0.9)
    assert np.max(np.abs(m["mag_err"])) < 0.06 and np.max(np.abs(m["phase_err_deg"])) < 5.0


def test_F_p7_literal_tec_hold_sign_and_settling(tmp_path):
    """Phase 5 review I-5: the TEC hold loop of p7_literal pulls a deliberate 0.3 rad error of the analyser's output
    offset back to the held point (negative feedback, ζ = 0.71 at f_hold = 0.5/(2π τ_th)); τ_th is shortened to 20 ms
    so that the run (≈ 10 hold time constants before the sweep starts) exercises the hold, then the sweep still
    matches the fast path."""
    pt = point(tmp_path, "pt", {"ident_n_freqs": 3, "ident_f_lo_hz": 1.0e3, "ident_f_hi_hz": 5.0e3,
                                "ident_t_start_s": 0.45, "ident_offset_err_rad": 0.3, "tec_tau_th_s": 0.02})
    spec = config.load("presets/future/F_wg6_nominal.yaml", overlays=["ident/p7_literal", "ident/in_lock_stepped", pt])
    p = spec.offline["params"]
    w_h = 2 * np.pi * p["ident_hold_f_hz"]
    assert 1 / (2 * np.sqrt(w_h * p["tec_tau_th_s"])) == pytest.approx(1 / np.sqrt(2))   # validity: ζ ≥ 0.5
    spec = dataclasses.replace(spec, measure={**spec.measure,
                                              "record": {**spec.measure["record"], "slots": ["lock.L0.err"]}})
    res = run(spec)
    t, e = np.asarray(res.t_s), np.asarray(res.series["lock.L0.err"])
    assert abs(np.mean(e[(t > 0.002) & (t < 0.005)])) > 0.2                # the offset error before the hold acts
    assert abs(np.mean(e[(t > 0.35) & (t < 0.45)])) < 0.02                 # held: back at the lock point
    m = res.measurements["loop_identification"]["lock0"]
    assert np.max(np.abs(m["mag_err"])) < 0.06 and np.max(np.abs(m["phase_err_deg"])) < 5.0


@pytest.mark.parametrize("preset, wire", [("presets/reference/refA.yaml", "ident/wire_refA"),
                                          ("presets/future/F_bor6_nominal.yaml", "ident/wire_F"),
                                          ("presets/future/F_liu20_nominal.yaml", "ident/p7_literal")])
@pytest.mark.parametrize("wave", ["ident/in_lock_stepped", "ident/in_lock_multisine"])
def test_point_overlays_build(preset, wire, wave):
    spec = config.load(preset, overlays=[wire, wave])
    config.validate(dataclasses.replace(spec, pipeline=inject_topology(spec.pipeline, spec.topology)))
    assert "loop_identification" in spec.measure
    ident = next(b for b in offline_block_specs(spec) if b.name == "ident")
    make_block(ident)


@pytest.mark.parametrize("ctrl, rdiv", [("controller/dsc1", 2), ("controller/fpga_generic", 1)])
def test_F_wire_F_runs_at_the_controller_rate(tmp_path, ctrl, rdiv):
    """Phase 5 review I-7 (D91-1 / D96-2): ident/wire_F's excitation runs at the fast stage's update rate (rate_div =
    pid_rate_div: 2 on controller/dsc1's 500 kHz) and the record sits on the update instants (decim a multiple of it,
    >= 2); with the DSC1 the in-lock L-hat still matches chk.lock_stability (dsc1 law, h = 2 us) through the runner."""
    from cranberry.config.spec import iter_block_specs

    pt = tmp_path / "pt.yaml"   # 0.1 rad of PZT phase: ~6 DSC1 DAC LSB (I-6); the DSC1 model with its 120 / 100 kHz poles
    pt.write_text("offline:\n  params:\n    ident_n_freqs: 4\n    ident_f_lo_hz: 500.0\n    ident_f_hi_hz: 1.0e+4\n"
                  "    ident_amplitude: !expr \"0.1 / offline.params.pzt_k_rad_per_v\"\n    dsc1_latency_s: 5.7e-6\n")
    spec = config.load("presets/future/F_wg6_nominal.yaml", overlays=[ctrl, "ident/wire_F", "ident/in_lock_stepped", str(pt)])
    blocks = {b.name: b for _, b in iter_block_specs(spec.pipeline)}
    assert blocks["ident.exc"].rate_div == blocks["pid0f"].rate_div == rdiv
    decim = spec.measure["record"]["decim"]
    assert decim % rdiv == 0 and decim >= 2
    m = run(spec).measurements["loop_identification"]["lock0"]
    assert np.all(m["coherence"] > 0.9)
    assert np.max(np.abs(m["mag_err"])) < 0.06 and np.max(np.abs(m["phase_err_deg"])) < 5.0
