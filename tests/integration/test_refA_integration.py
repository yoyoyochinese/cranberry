"""Ref A integration (WP-I1, M1): measurement emulation through the runner hooks, convergence and run
records for presets/reference/refA.yaml.

* Fig. 3 phase-scanned trace (T04 §5.4, V9; C24): lock 2 off, the LO ramped at ≈ 70 rad/s (scan_ramp →
  scan actuator), zero span 3 MHz / RBW 1 MHz. The maxima reach the anti-squeezing level (±0.5 dB)
  whatever the VBW; the minima deepen monotonically with VBW (100 Hz → 1 kHz → 10 kHz). Their depth is
  NOT compared with the figure (C24: the figure's ≈ −10 dB minima at VBW 100 Hz are unexplained).
* Fig. 4 spectra (T04 §5.5, T13): swept 1–100 MHz, RBW 1 MHz, VBW 100 Hz through ``bhd_frontend``
  (35 MHz first-order TIA, 28 dB clearance at 3 MHz). (a) without circuit-noise subtraction the
  squeezing degrades towards 100 MHz as the clearance |H|²C₀ falls (analytic (R + 1/C)/(1 + 1/C)
  ±0.1 dB on 10-bin means); (b) with subtraction it stays flat at R (±0.1 dB); anti-squeezing is flat in
  both (≈ R′+). The 1 MHz probe peak of Fig. 4(a) is not emulated (line spectra are not in the hook).
* Convergence (PLAN §6.4, prompt §8; ``runner.convergence``): f_sim × 2, × 4 change R′∓ by ≤ 1/3 of
  ±0.10 dB, and R′− by ≤ 0.005 dB (statistical scale of a 0.05 s record) — the laser × FC ledger pair
  adds continuous − discrete expectation (M1-A-01), so the expected total does not depend on f_sim (the
  earlier "0.006 dB at × 2" of D63-1 was two O(h) errors cancelling); chunk length bit-identical;
  10 seeds: spread of R′− ≪ 0.10 dB.
* Single source (M1-A-02): the laser's FC pair (Δν_FC, τ_fc, τ_ref) equals ``fc.dnu_fc_hz``,
  ``fc.in_delay_s`` and ``Topology.timing`` of the beams lo / sqz, for refA, refA_tap and the Fig. 5
  overlays; the relaxation-oscillation peak option runs. ESA RBW (M1-A-05): −3 dB width = RBW.
* Run record with ``--out`` (CLI): config, seeds, git, env, metrics, loss registry, band ledger and the
  measurement arrays in results.npz; YAML round trip of the preset (VT-F1 form).
"""

import copy
import json
import math

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import load_yaml
from cranberry.runner.convergence import convergence
from cranberry.runner.run import run

pytestmark = [pytest.mark.numba, pytest.mark.integration]
PRESET = "presets/reference/refA.yaml"


def refa_tree():
    return load_yaml(PRESET)


@pytest.fixture(scope="module")
def base():
    return run(config.load(PRESET))


# --------------------------------------------------------------------------------------- Fig. 3


def test_fig3_phase_scan_minima_vs_vbw():
    t = refa_tree()
    t["toggles"] = {"pid2.feedback": False, "scan_ramp.enabled": True}  # lock 2 off, LO scanned
    zs = next(b for b in t["offline"]["blocks"] if b["name"] == "esa_zs")
    esa = {}
    for vbw in (100.0, 1e3, 1e4):
        blk = copy.deepcopy(zs)
        blk["name"] = f"esa_vbw{int(vbw)}"
        blk["params"]["vbw_hz"] = vbw
        blk["params"]["n_points"] = 3001  # resolve the ≈ 2.5 ms minima
        t["offline"]["blocks"].append(blk)
        esa[f"scan{int(vbw)}"] = {"block": blk["name"], "metric": "i.port0", "quadrature_rad": 0.0}
    t["measure"]["esa"] = esa
    res = run(config.load(t), keep_state=False)
    nl = res.metrics.noise_levels["i.port0"]
    traces = [res.measurements["esa"][k]["trace_db"] for k in ("scan100", "scan1000", "scan10000")]
    mins = [float(np.min(tr)) for tr in traces]
    maxs = [float(np.max(tr)) for tr in traces]
    assert mins[0] > mins[1] > mins[2]                       # monotonic in VBW (T04 V9)
    assert mins[0] > -6.0 and mins[2] < -9.0                 # 100 Hz fills the minima; 10 kHz nearly reaches R′−
    # maxima reach the anti-squeezing level at every VBW (the wider VBWs add radiometer scatter on top:
    # σ ≈ 0.36 dB per point at 10 kHz, so their maximum over 3001 points lies ≈ 1 dB higher)
    assert maxs[0] == pytest.approx(nl["v_plus_db"], abs=0.5)
    assert all(m > nl["v_plus_db"] - 0.5 for m in maxs[1:])
    # the scan really sweeps: ≈ 70 rad/s → ≥ 8 minima (period π) in 0.15 s
    tr = traces[2]
    below = tr < 0.0
    assert np.count_nonzero(below[1:] & ~below[:-1]) >= 8


# --------------------------------------------------------------------------------------- Fig. 4


def test_fig4_spectra_with_and_without_subtraction(base):
    sq = base.measurements["esa"]["spectrum_sqz"]
    anti = base.measurements["esa"]["spectrum_anti"]
    f = sq["f_hz"]
    c0 = 10**2.8 * (1 + (3e6 / 35e6) ** 2)
    c = c0 / (1 + (f / 35e6) ** 2)

    def meas(r):
        return 10 * np.log10((r + 1 / c) / (1 + 1 / c))

    r_sq = sq["r_lin"]
    blocks = [slice(0, 10), slice(45, 55), slice(90, 100)]
    for s in blocks:
        assert np.mean(sq["normalized_db"][s]) == pytest.approx(np.mean(meas(r_sq)[s]), abs=0.1)
        assert np.mean(sq["corrected_db"][s]) == pytest.approx(10 * math.log10(r_sq), abs=0.1)
        assert np.mean(anti["corrected_db"][s]) == pytest.approx(10 * math.log10(anti["r_lin"]), abs=0.1)
    # (a) degrades by ≈ 0.6 dB from 1 to 100 MHz, (b) does not
    assert np.mean(sq["normalized_db"][90:]) - np.mean(sq["normalized_db"][:10]) > 0.4
    assert abs(np.mean(sq["corrected_db"][90:]) - np.mean(sq["corrected_db"][:10])) < 0.1
    # the spectrum input is R before the detector's electronic loss (D58-1): η = 0.92 / η_e
    assert sq["eta"] == pytest.approx(0.92 * (1 + 10**-2.8), rel=1e-9)


# ---------------------------------------------------------------------------------- convergence


def test_convergence_fsim_chunk_seeds():
    t = refa_tree()
    t["run"]["duration_s"] = 0.05
    t["measure"].pop("esa")
    spec = config.load(t)
    q = {"i.port0.r_minus_db": 0.10, "i.port0.r_plus_db": 0.10}
    rep = convergence(spec, ("f_sim", "seed"), quantities=q, seeds=list(range(1, 11)), f_sim_factors=(2, 4))
    assert rep.checks["f_sim"]["passed"], rep.checks["f_sim"]
    for k in ("x2", "x4"):  # M1-A-01: no O(h) bias left (statistical scale only)
        assert abs(rep.checks["f_sim"]["values"][k]["i.port0.r_minus_db"]["change"]) <= 0.005, rep.checks["f_sim"]
    # chunk length: bit-identical series and final physical state. Checked here rather than with
    # convergence("chunk"), which also compares the DrawBuffer bookkeeping slots rng.<block>.buf/.idx —
    # they differ with the refill pattern although every drawn number is the same (false negative of
    # the framework check, reported by WP-I1).
    import dataclasses

    a = run(spec)
    b = run(dataclasses.replace(spec, run=dataclasses.replace(spec.run, chunk_samples=100_003)))
    assert set(a.series) == set(b.series)
    for k in a.series:
        assert a.series[k].tobytes() == b.series[k].tobytes(), k
    for n in a.final_state.layout.slot_names():
        if not n.startswith(("rec.", "rng.")):
            assert np.asarray(a.final_state[n]).tobytes() == np.asarray(b.final_state[n]).tobytes(), n
    seed = rep.checks["seed"]
    assert seed["passed"]
    assert seed["values"]["i.port0.r_minus_db"]["std"] < 0.01
    assert seed["values"]["i.port0.r_minus_db"]["mean"] == pytest.approx(-10.15, abs=0.10)


@pytest.mark.parametrize("preset, overlay", [(PRESET, None), (PRESET, "refA_curve/green"),
                                             (PRESET, "refA_curve/previous_work"),
                                             ("presets/reference/refA_tap.yaml", None)])
def test_fc_pair_single_source(preset, overlay):
    from cranberry.path.topology import Topology

    spec = config.load(preset, overlays=[overlay] if overlay else [])
    blocks = {c.name: c for c in spec.pipeline.children}
    (pair,) = blocks["koheras"].params["oob_pairs"]
    fc = blocks["fc"].params
    topo = Topology.from_spec(spec.topology)
    assert pair["target"] == fc["segment"] and fc["beam"] == "lo"
    assert pair["dnu_fc_hz"] == fc["dnu_fc_hz"] == (1.175e6 if overlay else 3.3e6)
    assert pair["tau_fc_s"] == fc["in_delay_s"] == pytest.approx(topo.timing("lo").tau_s, rel=1e-12)
    assert pair["tau_ref_s"] == pytest.approx(topo.timing("sqz").tau_s, rel=1e-12)


def test_relaxation_oscillation_peak_option():
    """`koheras.peaks` (UNVERIFIED stand-in, off in refA) adds S_ν in the 0.1–3 MHz band that sets θ."""
    t = refa_tree()
    t["run"]["duration_s"] = 0.02
    t["measure"].pop("esa")
    th = {}
    for on in (False, True):
        t["toggles"] = {"koheras.peaks": on}
        th[on] = run(config.load(copy.deepcopy(t)), keep_state=False).metrics.noise_levels["i.port0"]["sigma_total_rad"]
    assert th[True] > th[False] * 1.02  # ≈ +3.5 % (|1 − H_FC(0.7 MHz)|² · area / f0²)


def test_esa_gaussian_rbw_3db_width():
    """M1-A-05: refA's Gaussian RBW has −3 dB full width = RBW (enbw_factor = √(π/(4 ln 2)) = 1.0645)."""
    from cranberry.core.block import BlockSpec
    from cranberry.core.registry import make_block

    spec = config.load(PRESET)
    for b in spec.offline["blocks"]:
        if b["name"] in ("esa_zs", "esa_spec"):
            assert b["params"]["enbw_factor"] == pytest.approx(1.0645, abs=1e-4)
            blk = make_block(BlockSpec.from_obj(b))
            f = np.linspace(0, 1e6, 200_001)
            half = f[np.argmin(np.abs(blk._rbw_power(f) - 0.5))]
            assert 2 * half == pytest.approx(b["params"]["rbw_hz"], rel=1e-4)


# -------------------------------------------------------------------------------------- records


def test_cli_run_record_and_round_trip(tmp_path, base):
    from cranberry.__main__ import main

    out = tmp_path / "runs"
    assert main(["run", PRESET, "--out", str(out)]) in (0, None)
    (run_dir,) = list(out.iterdir())
    names = {p.name for p in run_dir.iterdir()}
    assert {"config.input.yaml", "config.resolved.yaml", "seeds.json", "git.json", "env.json", "results.npz",
            "metrics.json", "loss_registry.json", "band_ledger.json", "tables.json"} <= names
    m = json.loads((run_dir / "metrics.json").read_text())
    assert m["noise_levels"]["i.port0"]["r_minus_db"] == pytest.approx(
        base.metrics.noise_levels["i.port0"]["r_minus_db"], abs=1e-12)  # same seed → same result
    with np.load(run_dir / "results.npz") as z:
        assert {"measure.esa.squeezed.trace_db", "measure.esa.spectrum_sqz.corrected_db", "metric.mi.port0"} <= set(z.files)
    reg = json.loads((run_dir / "loss_registry.json").read_text())
    assert any(e["id"] == "sqz_opa.wg" and "Ref A" in e["source"] for e in reg["elements"])
    # VT-F1 form: load(dump(spec)) == spec for both presets
    for p in (PRESET, "presets/reference/refA_tap.yaml"):
        s = config.load(p)
        assert config.load(load_yaml_text(config.dump_yaml(s))) == s


def load_yaml_text(text):
    from cranberry.config.yaml_io import load_yaml as _ly
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.yaml"
        p.write_text(text)
        return _ly(p)
