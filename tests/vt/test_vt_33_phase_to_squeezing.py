"""VT-33 (WP-F4 part): phase → squeezing through the framework, metric invariance / decomposition, HF
add-on (S §4 VT-33 ← T04 V1, V5, V12, V13; G05, G13; PLAN §6.2).

* **T04 V1 via the full framework path**: topology (LO + squeezed beam at a BHD), path blocks with the
  topology injected by the runner, a deterministic LO jitter δ = m sin(2πft) with rms 9 mrad, the
  ``metric_i`` probe, the offline ``single_mode_reference`` source (Ref A: L = 8 %, α = 9.06 /W,
  P = 0.64 W) and ``phase_to_quadrature`` in the runner's measurement stage:
  R′− = −10.148 dB, R′+ = +20.556 dB (±1e−3 dB; the sinusoid's χ = J0(2m) differs from e^{−2σ²} by
  σ⁴ ≈ 7e−9), and the D48 window (±0.10 dB). ``exact`` and ``small_angle`` differ by < 1e−3 dB.
* **T04 V5 (Ref C dither, analytic mode)**: L = 21 %, V+ = 10 dB (e^{2r} = 12.4), m = 4.5°:
  −5.63 dB without jitter → −5.18 dB with the dither (±0.01 dB), and the time-series dither gives the
  same V′ as the analytic J0 (1e−6).
* **T04 V12**: metric (ii) is invariant under noise common to both ports (a segment shared just before
  the GBS; Θ_1 equal to 1e−12 rad) while metric (i) changes; Θ is additive in independent sources
  (sample by sample, 1e−12) in a linear loop; the four-component decomposition of a toy linear lock
  (Liu-type probe lock, PI + delay) sums to the all-on variance within 5 % (achieved ≈ 1 %), and is
  marked ``decomposition_invalid`` when the actuator saturates.
* **T04 V13**: out-of-band noise added through the ``BandLedger`` (coarse loop, f_N = 0.5 MHz, white
  phase noise above it declared up to f_max = 4 MHz) and the same noise in a fine-step time series
  (f_N = 4 MHz) give the same R′∓ within 0.01 dB (statistical 3σ ≈ 0.003 dB with 2^20 samples).
"""

import math

import numpy as np
import pytest
from scipy.special import j0

import cranberry.core.testing  # noqa: F401  (toy kinds)
import cranberry.path.testing as pt
from cranberry.blocks.meas.metrics import METRIC_II_COMPONENTS
from cranberry.blocks.meas.phase_to_quadrature import phase_to_quadrature
from cranberry.blocks.meas.testing import (
    metric_cfg,
    path_blocks,
    refa_source,
    run_spec,
    single_port_topology,
    toy,
    two_port_loop_spec,
)
from cranberry.core.registry import register_block
from cranberry.core.units import V_VAC, db_to_lin
from cranberry.offline import gaussian
from cranberry.offline.band_ledger import BandDecl
from cranberry.runner.decompose import decompose_metric_ii
from cranberry.runner.run import run

pytestmark = pytest.mark.numba

SIGMA = 0.009


def refa_jitter_spec(average="exact", f_sim=1.0e6, n=100_000, extra=(), metrics_extra=None, jitter_slot="drv.lo_seg.stretcher"):
    m = SIGMA * math.sqrt(2.0)
    pipe = [toy("toy_sine", "jit", out=jitter_slot, unit="1", amp=m, freq_hz=1.0e3),
            toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
            *extra,
            *path_blocks({"las": "laser.las.phi"}),
            {"block": "metric_i", "name": "mi", "params": {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                                                          "setpoint_rad": 0.0}}]
    return run_spec(pipe, f_sim_hz=f_sim, n_samples=n, topology=single_port_topology(),
                    offline={"blocks": [refa_source()]},
                    measure={"record": {"slots": ["beam.lo.phi"], "decim": 10},
                             "metrics": {"i": metric_cfg("refA.sqz", average=average, **(metrics_extra or {}))}})


@pytest.mark.parametrize("average", ["exact", "small_angle"])
def test_v1_refA_through_the_framework(average):
    res = run(refa_jitter_spec(average))
    assert res.mode == "compiled"
    nl = res.metrics.noise_levels["i.port0"]
    assert nl["sigma_loop_rad"] == pytest.approx(SIGMA, rel=1e-6)
    assert nl["r_minus_db"] == pytest.approx(-10.148, abs=1e-3)
    assert nl["r_plus_db"] == pytest.approx(20.556, abs=1e-3)
    assert abs(nl["r_minus_db"] + 10.15) <= 0.10 and abs(nl["r_plus_db"] - 20.56) <= 0.10  # D48 window
    assert nl["eta"] == 1.0 and nl["hf_variance_rad2"] == 0.0  # lumped §7-1 loss is inside R±
    if average == "exact":  # the measured χ is J0(2m) of the sinusoid (integer number of periods)
        r = math.sqrt(9.06 * 0.64)
        rm = gaussian.jitter_average(0.08 + 0.92 * math.exp(-2 * r), 0.08 + 0.92 * math.exp(2 * r),
                                     j0(2 * SIGMA * math.sqrt(2)), delta0_rad=0.0)
        assert nl["r_minus"] / V_VAC == pytest.approx(rm, rel=1e-6)  # path-delay interpolation ~1e-8


def test_v1_exact_vs_small_angle_framework():
    ex = run(refa_jitter_spec("exact")).metrics.noise_levels["i.port0"]
    sa = run(refa_jitter_spec("small_angle")).metrics.noise_levels["i.port0"]
    assert abs(ex["r_minus_db"] - sa["r_minus_db"]) < 1e-3 and abs(ex["r_plus_db"] - sa["r_plus_db"]) < 1e-3


def test_v5_refC_dither_analytic_and_timeseries():
    e2r = 12.4
    vm, vp = (0.21 + 0.79 / e2r) * V_VAC, (0.21 + 0.79 * e2r) * V_VAC
    m = math.radians(4.5)
    none = phase_to_quadrature(np.zeros(8), vm, vp, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0)
    ana = phase_to_quadrature(np.zeros(8), vm, vp, hf_variance_rad2=0.0, average="exact", setpoint_rad=0.0,
                              dither_m_rad=[m])
    assert none["v_plus_db"] == pytest.approx(10.0, abs=0.01)
    assert none["r_minus_db"] == pytest.approx(-5.63, abs=0.01)
    assert ana["r_minus_db"] == pytest.approx(-5.18, abs=0.01)
    assert ana["chi_dither"] == pytest.approx(0.99384, abs=1e-5)
    assert ana["sigma_dither_rad"] == pytest.approx(0.0555, abs=1e-4)
    t = np.arange(16 * 1000) / 16.0  # 1000 dither periods, 16 samples each
    ts = phase_to_quadrature(m * np.sin(2 * np.pi * t), vm, vp, hf_variance_rad2=0.0, average="exact",
                             setpoint_rad=0.0)
    assert ts["r_minus"] == pytest.approx(ana["r_minus"], rel=1e-6)


# ------------------------------------------------------------------------------------------------ V12


def test_v12_metric_ii_common_mode_invariance():
    on = run(two_port_loop_spec(com_noise=True, n_samples=1 << 14, record_decim=1))
    off = run(two_port_loop_spec(com_noise=False, n_samples=1 << 14, record_decim=1))
    th_on, th_off = on.series["metric.mii.port1"], off.series["metric.mii.port1"]
    assert np.max(np.abs(th_on - th_off)) < 1e-12  # Θ_1 does not see the shared segment
    d_on, d_off = on.series["metric.mi.port1"], off.series["metric.mi.port1"]
    assert np.std(d_on - d_off) > 1e-2  # metric (i) does (the LO does not traverse it)


def test_v12_additivity_of_independent_sources():
    def theta(**tg):
        toggles = {f"{n}.enabled": tg.get(n, False) for n in ("n_lock", "n_cp", "n_sq", "n_las", "n_com")}
        r = run(two_port_loop_spec(n_samples=1 << 14, record_decim=1, toggles=toggles))
        return r.series["metric.mii.port1"], r.metrics.metric_ii[1].stats.var

    both, v_both = theta(n_cp=True, n_sq=True)
    cp, v_cp = theta(n_cp=True)
    sq, v_sq = theta(n_sq=True)
    np.testing.assert_allclose(both, cp + sq, rtol=0, atol=1e-12)  # linear loop: superposition
    assert v_both == pytest.approx(v_cp + v_sq, rel=0.05)


def test_v12_decomposition_linearity_toy_loop():
    comps = decompose_metric_ii(two_port_loop_spec())
    c = comps[1]
    assert c.status == "ok" and c.reason is None
    assert set(c.variances_rad2) == set(METRIC_II_COMPONENTS)
    assert all(v > 0 for v in c.variances_rad2.values())
    assert abs(c.linearity_ratio - 1.0) <= 0.05
    assert c.linearity_ratio == pytest.approx(sum(c.variances_rad2.values()) / c.total_var_rad2)


def test_v12_decomposition_invalid_on_saturation():
    c = decompose_metric_ii(two_port_loop_spec(sat_limit=0.05))[1]
    assert c.status == "decomposition_invalid"
    assert "saturation" in c.reason and "drv.a1.stretcher" in c.reason


# ------------------------------------------------------------------------------------------------ V13


class _Flat:
    def __init__(self, s0):
        self.s0 = s0

    def __call__(self, f):
        return np.full_like(np.asarray(f, dtype=float), self.s0)


class BandWhite(pt.WhiteDriver):
    """White phase-noise driver (σ² = S0·f_hi per sample, flat one-sided PSD S0 up to f_hi = f_N) that
    declares the same PSD above its band to the ledger (G05)."""

    KIND = "vt33_white_band"

    def bands(self):
        return [BandDecl(source=self.name, quantity="phase", targets=(str(self.params["segment"]),),
                         f_lo_hz=float(self.params["f_lo_hz"]), f_hi_hz=float(self.params["f_hi_hz"]), unit="rad",
                         psd_outside=_Flat(float(self.params["psd_rad2_per_hz"])))]


register_block(BandWhite.KIND)(lambda params, *, name=None: BandWhite(params, name))

S0 = 1.0e-4 / 4.0e6  # rad²/Hz: σ_total = 10 mrad up to f_max = 4 MHz
F_MAX = 4.0e6


def v13_spec(f_sim):
    f_n = f_sim / 2
    n = int(round(0.131072 * f_sim))
    noise = {"block": "vt33_white_band", "name": "wn", "params": {
        "out": "drv.lo_seg.stretcher", "shape": [], "unit": "1", "sigma": math.sqrt(S0 * f_n),
        "segment": "lo_seg", "f_lo_hz": f_sim / n, "f_hi_hz": f_n, "psd_rad2_per_hz": S0}}
    pipe = [noise, toy("toy_sine", "las", out="laser.las.phi", unit="rad", amp=0.0, freq_hz=1.0),
            *path_blocks({"las": "laser.las.phi"}),
            {"block": "metric_i", "name": "mi", "params": {"ports": {"0": {"squeezed": "sqz", "lo": "lo"}},
                                                          "setpoint_rad": 0.0}}]
    # zero-length segments: no fractional-delay interpolation, which would low-pass white noise near
    # Nyquist (linear interpolation, 1 − 2μ(1 − μ) of the variance) and bias the comparison
    topo = single_port_topology(lo_len_m=0.0, sq_len_m=0.0)
    return run_spec(pipe, f_sim_hz=f_sim, n_samples=n, topology=topo, offline={"blocks": [refa_source()]},
                    measure={"metrics": {"i": metric_cfg("refA.sqz", f_max_hz=F_MAX)}}, seed=11, chunk_samples=1 << 16)


def test_v13_ledger_vs_fine_step_time_series():
    coarse = run(v13_spec(1.0e6))
    fine = run(v13_spec(2 * F_MAX))
    c, f = coarse.metrics.noise_levels["i.port0"], fine.metrics.noise_levels["i.port0"]
    assert c["hf_variance_rad2"] == pytest.approx(S0 * (F_MAX - 0.5e6), rel=1e-4)  # ledger: f_N → f_max
    assert f["hf_variance_rad2"] == pytest.approx(0.0, abs=1e-15)  # nothing above f_max
    assert c["sigma_total_rad"] == pytest.approx(0.01, rel=0.01)
    assert f["sigma_total_rad"] == pytest.approx(0.01, rel=0.01)
    assert abs(c["r_minus_db"] - f["r_minus_db"]) < 0.01
    assert abs(c["r_plus_db"] - f["r_plus_db"]) < 0.01
    # the HF part is what makes them agree: without it the coarse loop underestimates the jitter
    no_hf = phase_to_quadrature(coarse.metrics.metric_i[0].stats, c["v_minus"], c["v_plus"], hf_variance_rad2=0.0,
                                average="exact", setpoint_rad=0.0)
    assert f["r_minus_db"] - no_hf["r_minus_db"] > 0.3
    assert c["v_minus_db"] == pytest.approx(10 * math.log10(0.08 + 0.92 * math.exp(-2 * math.sqrt(9.06 * 0.64))),
                                            abs=1e-9)
    assert db_to_lin(c["r_minus_db"]) > db_to_lin(c["v_minus_db"])
