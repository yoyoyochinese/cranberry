"""WP-I3 (M2): the refC preset (Verclas et al., Optica Quantum 4, 1 (2026); D03, D25, D34, D51-7, D55-1, D61;
PLAN §5.3, 00_summary §2.1 `refC` / `L-C`).

* every quantity Ref C leaves open is a named sweep axis under ``offline.params`` (D34; the fibre-noise
  correlation ``rho_conduit`` among them) and ``refC_nominal`` resolves all of them;
* the quantum_noise_power band integrals written in closed form in the preset equal ``qnl_band_integrals``
  of the opo_cw_spectrum source with the LossRegistry η (1e−9), and the stack gains follow the band
  assignment of [Ref C §2] (f_u1 > f_12 > f_23, positive offload);
* spectra vs [Ref C Fig. 4 / Fig. 5] through the runner's ESA hooks (cw_sideband metric read at each
  zero-span centre / swept bin): anti-squeezing 10.0 / 8.7 / 5.0 / 2.6 dB at 30 MHz / 0.25 / 0.69 / 1.21 GHz
  (T02 V-C1 ±0.1 dB; Fig. 4 read values 10.0 / 8.5 / 5.2 / 2.8 ±0.5 dB); squeezing monotonic in f and within
  ±0.7 dB of Fig. 4 (−4.7 / −4.7 / −3.0 / −2.2 dB; the model's −3.66 dB at 0.69 GHz is the largest gap,
  T02 §6.2 reports the same ≤ 0.5 dB class with θ̃ = 0.1 rad); the dark-noise-corrected swept spectrum
  equals R(f) of the metric, while the raw BHD_B output shows the interference dip at 0.83 GHz [Fig. 5];
* convergence (PLAN §6.4): f_sim × 2 and × 4 change the non-dither residual by ≤ 1/3 of 3 mrad, the chunk
  length is bit-identical, 10 seeds give a spread ≪ the 0.083 rad budget of VT-20;
* YAML round trip (VT-F1 form) of refC / refC_nominal and a recorded CLI run.
"""

import dataclasses
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from cranberry import config
from cranberry.blocks.err.qnl import qnl_band_integrals
from cranberry.config.yaml_io import Expr, loads_yaml
from cranberry.core.units import V_VAC
from cranberry.runner.convergence import convergence
from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parent / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]

AXES = {"x", "eta_fiber_link", "dL_m", "k_T_rad_per_m_K", "r_T_conduit_k_per_s", "r_T_spool_k_per_s",
        "dr_T_spool_k_per_s", "tau_th_conduit_s", "tau_th_spool_s", "h_ref_conduit_rad2_hz", "h_ref_spool_rad2_hz",
        "mech_gamma", "mech_f_c_hz", "rho_conduit", "rho_spool", "f_rho_hz", "h0_hz2_per_hz", "latency_s",
        "dac_bits", "f_u1_hz", "f_12_hz", "f_23_hz", "mirror_range_rad", "mirror_f_res_hz", "mirror_q",
        "stretcher_range_rad", "dip_r"}


@pytest.fixture(scope="module")
def full():
    """refC_nominal with the ESA hooks (0.62 s: 0.6 s zero spans after the 0.02 s settle)."""
    return run(rs.load(esa=True))


def test_refC_unknowns_are_named_sweep_axes():
    markers = config.sweep_markers(config.load(rs.PRESET))
    axes = {k.removeprefix("offline.params.") for k in markers if k.startswith("offline.params.")
            and not isinstance(markers[k], Expr)}
    assert axes == AXES
    assert all(isinstance(v, Expr) for k, v in markers.items() if k.removeprefix("offline.params.") not in AXES)
    assert config.sweep_markers(config.load(rs.NOMINAL)) == {}


def test_refC_qnl_integrals_and_stack_gains():
    res = run(rs.load(run={"duration_s": 0.025}))
    p = rs.params_of(res.spec)
    src = res.offline.source_outputs["opo"]
    eta = res.offline.loss_registry.eta("sqz", "bhd")
    f = np.linspace(p["f1_hz"], p["f2_hz"], 200_001)
    from cranberry.blocks.source.opo_cw import make_opo_cw_spectrum

    opo_p = dict(next(b for b in res.spec.offline["blocks"] if b["name"] == "opo")["params"])
    dense = make_opo_cw_spectrum({**opo_p, "f_s_hz": list(f)}).source_output
    vm, vp, _ = dense.project(None)
    r_m, r_p = (2 * (eta * v + (1 - eta) * V_VAC) for v in (vm, vp))
    num = qnl_band_integrals(f, np.ones_like(f), r_m, r_p, np.full_like(f, 10 ** -1.1), np.ones_like(f))
    for k, v in num.items():
        assert p[f"qnl_{k}"] == pytest.approx(v, rel=1e-9), k
    assert src.basis == "cw_sideband" and eta * p["escape_efficiency"] == pytest.approx(p["eta_total"], rel=1e-12)
    # stack (T16 §4.4): |K1 C1| K_d = 1 at f_u1, |K2 C2| = K1 at f_12, |K3 C3| = K2 at f_23; offload positive
    blocks = {b.name: b.params for b in _leaves(res.spec.pipeline)}
    k1 = math.pi / p["v_pi_v"]
    k2 = p["mirror_range_rad"] / p["out_range_v"]
    k3 = p["stretcher_range_rad"] / p["out_range_v"]
    assert blocks["pid1"]["ki_per_s"] == pytest.approx(-2 * math.pi * p["f_u1_hz"] / (p["qnl_slope_v_per_rad"] * k1))
    assert k2 * blocks["pid2"]["ki_per_s"] == pytest.approx(2 * math.pi * p["f_12_hz"] * k1)
    assert k3 * blocks["pid3"]["ki_per_s"] == pytest.approx(2 * math.pi * p["f_23_hz"] * k2)
    assert p["f_u1_hz"] > 10 * p["f_12_hz"] > 100 * p["f_23_hz"]
    assert blocks["dl1"]["tau_s"] == pytest.approx(p["latency_s"] - 1 / res.spec.run.f_sim_hz, abs=1e-18)   # D55-1
    # locked at Δ = 0 (mod π), the fast stage near the centre
    nl = res.metrics.noise_levels["i.port0"]
    assert abs(rs.at(nl, "mean_loop_rad", 3e7)) < 2e-3


def _leaves(spec):
    from cranberry.core.pipeline import BlockSpec, ParallelSpec, SeqSpec, StackSpec

    if isinstance(spec, BlockSpec):
        yield spec
    elif isinstance(spec, StackSpec):
        for c in (*spec.stages, *spec.actuators):
            yield from _leaves(c)
    else:
        for c in spec.children:
            yield from _leaves(c)


FIG4_F = (3e7, 2.5e8, 6.9e8, 1.21e9)
FIG4_ANTI = (10.0, 8.5, 5.2, 2.8)
FIG4_SQZ = (-4.7, -4.7, -3.0, -2.2)
VC1_ANTI = (10.0, 8.7, 5.0, 2.6)


def test_refC_zero_spans_vs_fig4(full):
    esa = full.measurements["esa"]
    tag = {3e7: "30m", 2.5e8: "250m", 6.9e8: "690m", 1.21e9: "1210m"}
    sq = [esa[f"sqz_{tag[f]}"]["mean_db"] for f in FIG4_F]
    an = [esa[f"anti_{tag[f]}"]["mean_db"] for f in FIG4_F]
    np.testing.assert_allclose(an, VC1_ANTI, atol=0.1)
    np.testing.assert_allclose(an, FIG4_ANTI, atol=0.5)
    np.testing.assert_allclose(sq, FIG4_SQZ, atol=0.7)
    assert np.all(np.diff(sq) > 0) and np.all(np.diff(an) < 0)
    # zero-span means agree with the metric's V′ at each centre (radiometer scatter of 0.6 s at VBW 300 Hz)
    nl = full.metrics.noise_levels["i.port0"]
    for f, s in zip(FIG4_F, sq):
        assert s == pytest.approx(rs.at(nl, "r_minus_db", f), abs=0.1)


def test_refC_swept_spectrum_and_dip(full):
    sp = full.measurements["esa"]["spectrum_sqz"]
    an = full.measurements["esa"]["spectrum_anti"]
    f = sp["f_hz"]
    nl = full.metrics.noise_levels["i.port0"]
    r_m = np.interp(f, nl["f_s_hz"], 10 ** (np.asarray(nl["r_minus_db"]) / 10))
    r_p = np.interp(f, nl["f_s_hz"], 10 ** (np.asarray(nl["r_plus_db"]) / 10))
    np.testing.assert_allclose(sp["r_lin"], r_m, rtol=1e-9)
    np.testing.assert_allclose(an["r_lin"], r_p, rtol=1e-9)
    # dark-noise corrected traces follow R(f). Radiometer scatter of (S − D)/(N − D) at RBW 300 kHz / VBW 300 Hz
    # with 11 dB clearance: ≈ 0.35 dB per squeezed bin at low f → 20-bin mean, 3σ = 0.25 dB; near the dip the
    # clearance drops to a few dB and the scatter grows (the "interference bump" of Fig. 5): median over all bins
    assert np.mean(sp["corrected_db"][:20]) == pytest.approx(np.mean(10 * np.log10(r_m[:20])), abs=0.25)
    assert np.nanmedian(sp["corrected_db"] - 10 * np.log10(r_m)) == pytest.approx(0.0, abs=0.3)
    for s in (slice(0, 20), slice(55, 75), slice(110, 130)):
        assert np.mean(an["corrected_db"][s]) == pytest.approx(np.mean(10 * np.log10(r_p[s])), abs=0.15)
    # > 3 dB squeezing below ≈ 0.7 GHz [Ref C §3]: the model crosses −3 dB between 0.69 and 1.21 GHz
    corr = 10 * np.log10(r_m)
    assert np.all(corr[f < 0.7e9] < -3.0) and corr[-1] > -3.0
    # raw (uncorrected) BHD_B output: the 0.83 GHz interference dip [Ref C Fig. 5 caption]
    fe = next(b for b in full.offline.blocks if b.name == "frontend")
    raw = fe.psd(f, 1.0)["shot"]
    assert f[np.argmin(raw)] == pytest.approx(0.83e9, abs=10e6)


def test_refC_convergence_fsim_chunk_seeds():
    spec = rs.load(run={"duration_s": 0.08})
    q = {"i.port0.std_rad": 3e-3}
    rep = convergence(spec, ("f_sim", "seed"), quantities=q, seeds=list(range(1, 11)), f_sim_factors=(2, 4))
    assert rep.checks["f_sim"]["passed"], rep.checks["f_sim"]
    seed = rep.checks["seed"]
    assert seed["passed"]
    assert seed["values"]["i.port0.std_rad"]["std"] < 1e-3                 # ≪ 0.083 rad (VT-20)
    # chunk length: bit-identical series and final physical state (refA pattern; rng.* bookkeeping excluded)
    a = run(spec)
    b = run(dataclasses.replace(spec, run=dataclasses.replace(spec.run, chunk_samples=100_003)))
    for k in a.series:
        assert a.series[k].tobytes() == b.series[k].tobytes(), k
    for n in a.final_state.layout.slot_names():
        if not n.startswith(("rec.", "rng.")):
            assert np.asarray(a.final_state[n]).tobytes() == np.asarray(b.final_state[n]).tobytes(), n


def test_refC_round_trip_and_cli_record(tmp_path):
    for p in (rs.PRESET, rs.NOMINAL):
        s = config.load(p)
        assert config.load(loads_yaml(config.dump_yaml(s))) == s
    from cranberry.__main__ import main

    nominal_short = tmp_path / "refC_short.yaml"
    nominal_short.write_text(f"extends: {rs.ROOT / rs.NOMINAL}\nrun: {{duration_s: 0.03}}\nmeasure: {{esa: null}}\n")
    out = tmp_path / "runs"
    assert main(["run", str(nominal_short), "--out", str(out)]) in (0, None)
    (run_dir,) = list(out.iterdir())
    m = json.loads((run_dir / "metrics.json").read_text())
    assert m["noise_levels"]["i.port0"]["f_s_hz"][1] == 3e7
    reg = json.loads((run_dir / "loss_registry.json").read_text())
    assert any(e["id"] == "optics.fiber_link" and "Ref C" in e["source"] for e in reg["elements"])


def test_refC_square_law_overlay_locks():
    """presets/axes/refC_detector/square.yaml (T14 §6.2): same loop gain, κ_det 0.97 → similar residual."""
    res_sq = run(rs.load(run={"duration_s": 0.08}, overlays=("presets/axes/refC_detector/square.yaml",)))
    res_log = run(rs.load(run={"duration_s": 0.08}))
    blocks = {b.name: b.params for b in _leaves(res_sq.spec.pipeline)}
    assert blocks["qnl"]["detector_law"] == "square" and blocks["qnl"]["kappa_det"] == 0.97
    nl_sq, nl_log = (r.metrics.noise_levels["i.port0"] for r in (res_sq, res_log))
    assert abs(rs.at(nl_sq, "mean_loop_rad", 3e7)) < 2e-3
    assert rs.at(nl_sq, "std_loop_rad", 3e7) == pytest.approx(rs.at(nl_log, "std_loop_rad", 3e7), rel=0.2)
