"""core.delay: ring + fractional delay (linear, thiran1), the delay_buffer block (PLAN §2.6, VT-44 part)."""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.delay import DelayDesign, DelayLine, make_delay_buffer
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.core.testing import toy

FS = 1.0e6
METHODS = ("linear", "thiran1")


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("k", [0, 1, 5, 17])
def test_integer_delay_is_exact(method, k):
    x = np.random.default_rng(0).standard_normal(64)
    y = DelayLine(k / FS, FS, method).process(x)
    np.testing.assert_array_equal(y[k:], x[: len(x) - k])
    assert np.all(y[:k] == 0.0)  # zero initial history


def test_linear_fractional_formula():
    x = np.random.default_rng(1).standard_normal(50)
    dl = DelayLine(3.25 / FS, FS, "linear")
    assert (dl.design.k, dl.design.mu) == (3, pytest.approx(0.25))
    y = dl.process(x)
    xp = np.concatenate([np.zeros(5), x])
    expect = 0.75 * xp[5 - 3 : 5 - 3 + 50] + 0.25 * xp[5 - 4 : 5 - 4 + 50]
    np.testing.assert_allclose(y, expect, rtol=0, atol=1e-15)


@pytest.mark.filterwarnings("ignore::cranberry.core.checks.ChkWarning")
@pytest.mark.parametrize("D", [0.3, 0.5, 1.0 + 1e-3, 2.7, 6.49, 6.5])
def test_thiran_design(D):
    d = DelayDesign(D, "thiran1")
    frac = D - d.k
    assert (0.5 <= frac < 1.5) or (D < 0.5 and d.k == 0)
    assert abs(d.a) < 1  # stable
    f = np.linspace(1e-6, 0.49, 200) * FS
    np.testing.assert_allclose(np.abs(d.freq_response(f, FS)), 1.0, atol=1e-13)  # all-pass


def test_thiran_short_delay_warns():
    """Phase 5 F1-6: thiran1 is documented for D >= 0.5; 0 < D < 0.5 emits a ChkWarning (linear does not,
    integer delays do not)."""
    import warnings

    from cranberry.core.checks import ChkWarning

    with pytest.warns(ChkWarning, match="thiran1_short_delay"):
        DelayDesign(0.3, "thiran1")
    with warnings.catch_warnings():
        warnings.simplefilter("error", ChkWarning)
        DelayDesign(0.3, "linear")
        DelayDesign(0.5, "thiran1")
        DelayDesign(0.0, "thiran1")
    lay = SlotLayout(FS)
    lay.add("v.x", (), np.float64, "V")
    blk = make_delay_buffer({"in": "v.x", "out": "delay.dd.out", "tau_s": 0.2 / FS, "method": "thiran1"}, name="dd")
    with pytest.warns(ChkWarning, match="dd"):
        blk.init(lay, None)


@pytest.mark.filterwarnings("ignore::cranberry.core.checks.ChkWarning")
@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("D", [2.3, 5.75, 0.4])
def test_low_frequency_phase_delay_equals_D(method, D):
    """Phase delay −arg H(e^{jω})/ω → D as ω → 0 (exact at DC for both; Thiran is maximally flat)."""
    d = DelayDesign(D, method)
    f = 1e-4 * FS
    w = 2 * np.pi * f / FS
    assert -np.angle(d.freq_response(f, FS)) / w == pytest.approx(D, abs=1e-6)
    # 5 % of Nyquist: phase delay within 5e-3 (Thiran) / 1e-2 (linear) sample of D
    f2 = 0.05 * FS / 2
    w2 = 2 * np.pi * f2 / FS
    tol = 5e-3 if method == "thiran1" else 1e-2
    assert -np.angle(d.freq_response(f2, FS)) / w2 == pytest.approx(D, abs=tol)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("D", [3.4, 7.85])
def test_simulated_phase_response_matches_analytic(method, D):
    """Steady-state response to a sinusoid equals the analytic H(e^{jω}) of the discrete delay."""
    N, m = 4096, 37  # 37 periods in the analysis window
    w = 2 * np.pi * m / N
    n = np.arange(3 * N)
    y = DelayLine(D / FS, FS, method).process(np.cos(w * n))
    seg = slice(2 * N, 3 * N)  # transients (|a|^n) long gone
    h_sim = 2.0 / N * np.sum(y[seg] * np.exp(-1j * w * n[seg]))
    h_ana = DelayDesign(D, method).freq_response(w * FS / (2 * np.pi), FS)
    assert abs(h_sim - h_ana) < 1e-10


def run_delay_block(method, tau_s, x, rate_div=1, compiled=False):
    obj = [{"block": "delay_buffer", "name": "d", "rate_div": rate_div,
            "params": {"in": "v.x", "out": "delay.d.out", "tau_s": tau_s, "method": method, "unit": "V"}}]
    lay = SlotLayout(FS)
    lay.add("v.x", unit="V")
    blk = build(spec_from_obj(obj), lay, 0, external=("v.x",))
    st, p = State(lay), Params(lay)
    out = []
    for xi in x:
        st["v.x"] = xi
        (run_compiled if compiled else run_reference)(blk, st, p, 1)
        out.append(float(st["delay.d.out"]))
    return np.array(out), blk


@pytest.mark.numba
@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("tau_s", [0.0, 4e-6, 4.6e-6])
def test_block_matches_delay_line_compiled_and_reference(method, tau_s):
    x = np.random.default_rng(2).standard_normal(40)
    ref = DelayLine(tau_s, FS, method).process(x)
    y_ref, blk = run_delay_block(method, tau_s, x)
    y_cmp, _ = run_delay_block(method, tau_s, x, compiled=True)
    np.testing.assert_array_equal(y_ref, ref)
    np.testing.assert_allclose(y_cmp, ref, rtol=1e-12, atol=1e-15)
    assert blk.children[0].slot_units == {"v.x": "V", "delay.d.out": "V"}


def test_decimated_block_delays_in_its_own_steps():
    """rate_div = 2: the ring advances every 2 base samples, delay in units of 2/f_sim; output held."""
    x = np.arange(1.0, 21.0)
    y, _ = run_delay_block("linear", 4e-6, x, rate_div=2)  # 2 own steps
    own_in = x[0::2]  # samples seen at n = 0, 2, 4, ...
    own_out = np.concatenate([[0.0, 0.0], own_in[:-2]])
    np.testing.assert_array_equal(y, np.repeat(own_out, 2))


def test_delay_buffer_params_validated():
    base = {"in": "v.x", "out": "delay.d.out", "tau_s": 1e-6, "method": "linear"}
    with pytest.raises(MissingParameterError):
        make_delay_buffer({k: v for k, v in base.items() if k != "method"}, name="d")
    with pytest.raises(ValueError):
        make_delay_buffer({**base, "method": "cubic"}, name="d")
    with pytest.raises(ValueError):
        make_delay_buffer({**base, "tau_s": -1.0}, name="d")
    lay = SlotLayout(FS)
    blk = make_delay_buffer(base, name="d")
    with pytest.raises(ValueError, match="not registered yet"):
        blk.init(lay, None)


def test_toy_import_side_effect_free_of_delay():
    assert toy("toy_gain", "g", in_="v.a", out="v.b", unit="V", gain=1.0)["params"]["in"] == "v.a"


@pytest.mark.parametrize("D", [0.3, 2.7, 5.0])
def test_interp_variance_deficit_closed_form(D):
    """D65-4 helper: one linear stage on a unit-variance white driver loses 2μ(1 − μ) of the variance
    (review F6-1 closed form); an integer delay loses nothing; an AR(1)-like 1/f² tail is integrable."""
    from cranberry.core.delay import interp_variance_deficit

    mu = D - np.floor(D)
    white = lambda f: np.full_like(f, 2.0 / FS)  # noqa: E731  one-sided, ∫ = 1
    assert interp_variance_deficit([[(D / FS, 1.0)]], FS, white) == pytest.approx(2 * mu * (1 - mu), abs=1e-6)
    assert interp_variance_deficit([[(D / FS, 1.0)], [(D / FS, -1.0)]], FS, white) == pytest.approx(
        4 * mu * (1 - mu), abs=1e-6)  # independent channels add
    assert np.isfinite(interp_variance_deficit([[(D / FS, 1.0)]], FS, lambda f: 1.0 / f**2))

