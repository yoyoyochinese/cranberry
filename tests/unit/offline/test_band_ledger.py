"""offline.band_ledger: out-of-band variance per metric and beam combination (G05, T04 §3.4, V13)."""

import json

import numpy as np
import pytest
from scipy.integrate import quad

from cranberry.core.fft import ifft_agrawal
from cranberry.core.rng import stream
from cranberry.core.units import noise_db
from cranberry.offline import gaussian as g
from cranberry.offline.band_ledger import BandDecl, BandLedger, integrate_psd
from cranberry.offline.source_output import SourceOutput


def lorentz(s0, fc):
    return lambda f: s0 / (1.0 + (np.asarray(f) / fc) ** 2)


def ledger_with_metric(f_max=1e6, weight=None, **kw):
    led = BandLedger(**kw)
    led.define_metric("i", f_max_hz=f_max, weight=weight)
    return led


def test_integrate_psd_closed_forms():
    s0, fc, a, b = 3e-10, 2e5, 5e4, 1e6
    exact = s0 * fc * (np.arctan(b / fc) - np.arctan(a / fc))
    assert integrate_psd(lorentz(s0, fc), a, b) == pytest.approx(exact, rel=1e-4)
    # power law over 8 decades
    assert integrate_psd(lambda f: 1.0 / f, 1.0, 1e8) == pytest.approx(np.log(1e8), rel=1e-4)
    # white FM S_ν = h0 → S_φ = h0/f² (conv §7)
    assert integrate_psd(lambda f: 2.0 + 0 * f, 1e3, 1e5, divide_f2=True) == pytest.approx(2.0 * (1e-3 - 1e-5), rel=1e-4)
    assert integrate_psd(lambda f: f, 5.0, 5.0) == 0.0
    with pytest.raises(ValueError):
        integrate_psd(lambda f: f, 0.0, 1.0)


def test_hf_variance_phase_psd_and_variance_rows():
    s0, fc = 3e-10, 2e5
    led = ledger_with_metric()
    led.declare(BandDecl("noise.a", "phase", ("sqz1",), 10.0, 5e4, "rad", psd_outside=lorentz(s0, fc)))
    led.declare(BandDecl("gawbs", "phase", ("sqz1",), 10.0, 5e4, "rad", variance_outside=2e-5))
    led.declare(BandDecl("rin", "rin", ("sqz1",), 10.0, 5e4, "1", variance_outside=1.0))  # not phase: ignored
    exact = s0 * fc * (np.arctan(1e6 / fc) - np.arctan(5e4 / fc))
    assert led.hf_variance("i", ["sqz1"]) == pytest.approx(exact + 2e-5, rel=1e-4)
    assert led.oob_variance(led.find(source="rin")[0], "i") == 1.0
    # f_hi above the metric's limit contributes nothing
    led2 = ledger_with_metric(f_max=1e4)
    led2.declare(BandDecl("noise.a", "phase", ("sqz1",), 10.0, 5e4, "rad", psd_outside=lorentz(s0, fc)))
    assert led2.hf_variance("i", ["sqz1"]) == 0.0
    with pytest.raises(KeyError):
        led.hf_variance("ii", ["sqz1"])


def test_freq_quantity_converted_to_phase():
    h0 = 50.0  # white FM [Hz²/Hz], Lorentzian FWHM π h0 (conv §7)
    led = ledger_with_metric(f_max=1e7)
    led.declare(BandDecl("laser", "freq", ("lo1",), 1.0, 1e5, "Hz", psd_outside=lambda f: h0 + 0 * f))
    assert led.hf_variance("i", {"lo1": 1.0}) == pytest.approx(h0 * (1e-5 - 1e-7), rel=1e-4)
    with pytest.raises(ValueError):
        BandDecl("laser", "freq", ("lo1",), 1.0, 1e5, "Hz", variance_outside=1e-3)


def test_metric_weighting():
    """Sideband weighting w(f) (T04 §3.2) multiplies the PSD inside the integral."""
    led = ledger_with_metric(f_max=1e6, weight=lambda f: 0.5 + 0 * f)
    led.declare(BandDecl("n", "phase", ("sqz1",), 1.0, 1e3, "rad", psd_outside=lambda f: 1e-12 + 0 * f))
    assert led.hf_variance("i", ["sqz1"]) == pytest.approx(0.5 * 1e-12 * (1e6 - 1e3), rel=1e-4)


def test_beam_combination_signs_and_segments():
    led = ledger_with_metric(beams={"sqz1": ["src", "arm1", "common"], "lo1": ["lo_arm", "common"], "sqz0": ["src", "arm0"]})
    led.declare(BandDecl("env.common", "phase", ("common",), 1.0, 1e3, "rad", variance_outside=4e-6))
    led.declare(BandDecl("env.arm1", "phase", ("arm1",), 1.0, 1e3, "rad", variance_outside=1e-6))
    led.declare(BandDecl("env.lo", "phase", ("lo_arm",), 1.0, 1e3, "rad", variance_outside=2e-6))
    led.declare(BandDecl("src.pump", "phase", ("src",), 1.0, 1e3, "rad", variance_outside=3e-6))
    # metric (i): δ = φ_LO − θ_sq — common segment cancels
    assert led.hf_variance("i", {"lo1": 1, "sqz1": -1}) == pytest.approx(1e-6 + 2e-6 + 3e-6)
    # metric (ii): Θ_1 = θ_sq,1 − θ_sq,0 — the shared source segment cancels
    assert led.hf_variance("i", {"sqz1": 1, "sqz0": -1}) == pytest.approx(4e-6 + 1e-6)
    # single beam, unsigned list
    assert led.hf_variance("i", ["sqz1"]) == pytest.approx(4e-6 + 1e-6 + 3e-6)
    # unsigned list with a declaration hitting two listed beams is ambiguous
    with pytest.raises(ValueError):
        led.hf_variance("i", ["sqz1", "lo1"])
    # direct beam targets work without define_beam
    led.declare(BandDecl("det", "phase", ("lo1",), 1.0, 1e3, "rad", variance_outside=5e-7))
    assert led.affects(led.find(source="det")[0], "lo1")


def test_find_and_json():
    led = ledger_with_metric()
    led.declare(BandDecl("dc_fringe.gamma_coh", "phase", ("pr1",), 1.0, 5e5, "rad", variance_outside=1e-4))
    led.declare(BandDecl("noise.a", "phase", ("sqz1",), 10.0, 5e4, "rad", psd_outside=lorentz(1e-10, 1e5)))
    rows = led.find(source="dc_fringe.gamma_coh")
    assert len(rows) == 1 and rows[0].variance_outside == 1e-4
    assert led.find(target="sqz1")[0].source == "noise.a"
    rec = led.to_json()
    json.dumps(rec)
    assert rec["decls"][1]["outside"] == "psd"
    assert rec["decls"][0]["oob_variance_by_metric"]["i"] == 1e-4
    assert rec["metrics"]["i"]["f_max_hz"] == 1e6


def test_band_decl_validation():
    with pytest.raises(ValueError):
        BandDecl("x", "phase", ("a",), 10.0, 5.0, "rad", variance_outside=0.0)
    with pytest.raises(ValueError):
        BandDecl("x", "phase", ("a",), 1.0, 5.0, "rad")
    with pytest.raises(ValueError):
        BandDecl("x", "phase", (), 1.0, 5.0, "rad", variance_outside=0.0)
    with pytest.raises(ValueError):
        BandDecl("x", "phase", ("a",), 1.0, 5.0, "rad", variance_outside=-1.0)
    led = ledger_with_metric()
    with pytest.raises(ValueError):
        led.define_metric("j", f_max_hz=np.inf, weight=None)


def test_t04_v13_style_ledger_vs_fine_time_series():
    """T04 V13 (as far as possible without F1): the measured V′ with (a) the in-loop band as a time
    series plus the ledger's χ_HF = e^{−2σ²_HF}, and (b) the full band as a fine-step time series,
    agree within 0.01 dB. Also: the ledger integral equals a direct numerical integration (quad).

    Fine series: Gaussian phase noise with one-sided PSD S(f) = S0/(1 + (f/fc)²) up to f_max, built
    by Fourier synthesis on the conv §8 wrapper (random phases, exact amplitudes √(2SΔf)). The loop
    band is f < f_hi (Nyquist of the loop rate); the ledger declares S outside it.
    """
    f_hi, f_max, fc = 5e4, 1e6, 2e5
    s0 = 7.3e-10  # σ²_in ≈ 3.6e−5, σ²_HF ≈ 1.6e−4 rad²
    psd = lorentz(s0, fc)
    n, dt = 2**18, 1.0 / 4e6
    df = 1.0 / (n * dt)
    k = np.arange(1, int(f_max / df) + 1)
    f = k * df
    rng = stream(2026, "test.band_ledger.v13")
    amp = np.sqrt(2.0 * psd(f) * df) * np.exp(1j * rng.uniform(0, 2 * np.pi, k.size))

    def synth(mask):
        x = np.zeros(n, dtype=complex)
        x[k[mask]] = amp[mask] * n * dt  # ifft_agrawal = fft(X)/(N dt) → Σ a e^{iφ} e^{−2πi k n/N}
        return ifft_agrawal(x, dt).real

    delta_fine = synth(np.ones(k.size, bool))
    delta_loop = synth(f < f_hi)

    led = BandLedger(beams={"sqz1": ["arm1"], "lo1": ["lo_arm"]})
    led.define_metric("i.port1", f_max_hz=f_max, weight=None)
    led.declare(BandDecl("noise.arm1", "phase", ("arm1",), df, f_hi, "rad", psd_outside=psd))
    var_hf = led.hf_variance("i.port1", {"lo1": 1.0, "sqz1": -1.0})
    var_quad = quad(psd, f_hi, f_max, limit=200)[0]
    assert var_hf == pytest.approx(var_quad, rel=1e-4)
    # time-series variance of the out-of-band part equals the ledger value (Parseval, Δf ≪ f_hi)
    assert np.var(delta_fine - delta_loop) == pytest.approx(var_hf, rel=2e-3)

    vm, vp, _ = SourceOutput.single_mode_reference(loss=0.08, alpha_per_w=9.06, pump_w=0.64, phi_opa0_rad=0.0).project(None)
    v_a = g.jitter_average(vm, vp, g.chi_timeseries(delta_loop) * g.chi_gaussian(np.sqrt(var_hf)), delta0_rad=0.0)
    v_b = g.jitter_average(vm, vp, g.chi_timeseries(delta_fine), delta0_rad=0.0)
    assert abs(float(noise_db(v_a) - noise_db(v_b))) < 0.01
    # without the ledger's out-of-band part the answer is visibly wrong (the test has teeth)
    v_no_hf = g.jitter_average(vm, vp, g.chi_timeseries(delta_loop), delta0_rad=0.0)
    assert abs(float(noise_db(v_no_hf) - noise_db(v_b))) > 0.3
