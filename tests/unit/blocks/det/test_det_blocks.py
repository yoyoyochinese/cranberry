"""DET blocks (WP-B5): contracts (D34, toggles, losses, device presets) and formulas of T13."""

import math

import numpy as np
import pytest

import cranberry.blocks.det.testing as T
from cranberry.blocks.det import bhd
from cranberry.blocks.det.pulsed import c1_first_order, crosstalk_coefficients, rho_adjacent
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.registry import make_block
from cranberry.core.state import Params, State
from cranberry.core.testing import toy
from cranberry.offline.loss_registry import LossRegistry


def mk(obj):
    return make_block(BlockSpec.from_obj(obj))


def reg(*blocks):
    r = LossRegistry(paths={"sqz": {"segments": [], "observers": ["bhd"]}})
    for b in blocks:
        r.add_all(b.losses())
    return r


def test_detector_efficiency_elements():
    d = mk(T.det_eff(on=True, **{"toggles": {"freq_response": False, "dip": False}}))
    els = {e.id.split(".")[-1]: e for e in d.losses()}
    assert set(els) == {"qe", "electronic", "imbalance"}
    assert els["qe"].kind == "detector" and els["electronic"].kind == "electronic"
    assert els["qe"].eta == 0.99
    assert els["electronic"].eta == pytest.approx(1 - 0.001582, abs=2e-6)  # C = 28 dB (T13 test 1)
    assert els["imbalance"].eta == pytest.approx(1 / (1 + 2.5e-4 * 61))
    assert reg(d).eta("sqz", "bhd") == pytest.approx(d.breakdown()["eta_det"], rel=1e-14)
    assert reg(d).eta("sqz", "gbs") == 1.0  # BHD-only elements never reach the GBS target (D46)


def test_toggles_off_identity_and_kept_in_record():
    d = mk(T.det_eff(on=False))
    r = reg(d)
    assert r.eta("sqz", "bhd") == 1.0 and len(r.breakdown("sqz", "bhd")) == 3
    assert d.breakdown()["eta_det"] == 1.0
    # subtract_dark removes the electronic loss only
    d2 = mk(T.det_eff(on=True, toggles={"subtract_dark": True, "freq_response": False, "dip": False}))
    assert not [e for e in d2.losses() if e.kind == "electronic"][0].enabled


def test_d34_no_defaults():
    obj = T.det_eff()
    del obj["params"]["clearance_lin"]
    with pytest.raises(MissingParameterError):
        mk(obj)
    obj = T.det_eff()
    del obj["params"]["toggles"]["dip"]
    with pytest.raises(MissingParameterError):
        mk(obj)
    with pytest.raises(MissingParameterError):  # both choices given
        mk(T.det_eff(i_n_a_rthz=1e-11, responsivity_a_per_w=1.0))


def test_clearance_scales_with_p_lo_and_from_noise():
    d = mk(T.det_eff(p_lo_w=0.0077, toggles={"freq_response": False, "dip": False}))
    assert d.clearance(3e6) == pytest.approx(10 ** 2.8 / 2)
    # C_0 = 2eℜP/i_n²: PDB450C at 3.33 mW, 69 pA/√Hz → ≈ 0.22 (T13 §4.1: −6.3 dB with the band-average NEP)
    obj = T.det_eff(device="pdb450c", p_lo_w=3.33e-3, wavelength_m=1.5488e-6,
                    toggles={"freq_response": False, "dip": False})
    for k in ("qe", "clearance_lin", "p_lo_ref_w"):
        del obj["params"][k]
    d = mk(obj)
    assert 10 * math.log10(d.clearance(0.0)) == pytest.approx(-6.3, abs=0.15)  # ℜ 1.0–1.03
    assert d.eta_q == pytest.approx(0.80, abs=0.005)  # ℜ = 1.0 A/W → η_q ≈ 0.80 (T13 §1.2)


def test_device_presets_load_and_spec_wins():
    d = mk(T.det_eff(device="bhd_generic_highqe", qe=0.97, toggles={"freq_response": False, "dip": False}))
    assert d.eta_q == 0.97
    obj = T.det_eff(device="bhd_generic_highqe", toggles={"freq_response": False, "dip": False})
    for k in ("qe", "clearance_lin", "p_lo_ref_w", "imbalance_abs2"):
        del obj["params"][k]
    d = mk(obj)
    assert d.eta_q == 0.99 and d.clearance(3e6) == pytest.approx(10 ** 2.8, rel=1e-12)
    assert d.delta_abs2 == pytest.approx(2.5e-4)  # CMRR 30 dB → |δ|² = 10^(−3)/4 (T13 §1.5)


def test_response_shapes():
    assert abs(bhd.h_tia(35e6, {"model": "first_order", "f3db_hz": 35e6})) ** 2 == pytest.approx(0.5)
    h2 = bhd.h_tia(np.array([0.0, 1e8]), {"model": "second_order", "f0_hz": 1e8, "q": 2.0})
    assert abs(h2[0]) == 1.0 and abs(h2[1]) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        bhd.h_tia(1.0, {"model": "third"})


def test_frontend_dc_sign_clip_load_and_monitor():
    fe = mk(T.frontend())
    assert fe.dc_voltage(3.33e-3, 0.0) == pytest.approx(3.33)
    assert mk(T.frontend(coupler_sign=-1)).dc_voltage(3.33e-3, 0.0) == pytest.approx(-3.33)
    assert fe.dc_voltage(5e-3, 0.0) == 4.6 and fe.dc_voltage(0.0, 5e-3) == -4.6
    assert mk(T.frontend(toggles={"saturation": False})).dc_voltage(5e-3, 0.0) == pytest.approx(5.0)
    half = mk(T.frontend(load="50ohm"))
    assert half.dc_voltage(3.33e-3, 0.0) == pytest.approx(3.33 / 2) and half.dc_voltage(9e-3, 0) == 2.3
    # MONITOR saturates at ~1 mW (PDB p.5, p.8): Ref D lock point 1.67 mW per PD → 10 V (T13 §4.1)
    assert fe.monitor(0.5e-3) == pytest.approx(5.0) and fe.monitor(1.67e-3) == 10.0


def test_pulsed_kernel():
    pk = mk(T.pulsed_kernel())
    c = pk.c
    assert c[0] == 1.0 and c[1] == pytest.approx(c1_first_order(50e6, 100e6), rel=1e-12)
    np.testing.assert_allclose(c, math.exp(-math.pi) ** np.arange(9), rtol=1e-10)
    assert pk.rho_adj == pytest.approx(c[1], rel=1e-9)
    assert pk.clearance_pulse == pytest.approx(1.3e7 / 730**2)  # Hansen 2001: ≈ 14 dB (T13 §4.3)
    assert 10 * math.log10(pk.clearance_pulse) == pytest.approx(13.9, abs=0.1)
    off = mk(T.pulsed_kernel(on=False))
    assert list(off.c) == [1.0] and off.clearance_pulse == math.inf and off.rho_adj == 0.0
    ring = crosstalk_coefficients({"model": "second_order", "f0_hz": 100e6, "q": 3.0}, 100e6, 1e-9, 6)
    assert np.any(ring < 0)  # ringing: sign-changing c_m (T13 §1.8)
    assert rho_adjacent([1.0, 0.0]) == 0.0
    with pytest.raises(ValueError):
        crosstalk_coefficients({"model": "first_order", "bandwidth_hz": 1e8}, 1e8, 2e-8, 3)  # W > T


@pytest.mark.numba
def test_offline_blocks_compiled_equals_reference():
    """All B5 blocks are offline / measurement stage: their kernels are no-ops; compiled == reference."""
    objs = [toy("toy_sine", "s", out="v.s", unit="V", amp=1.0, freq_hz=1e3), T.det_eff(), T.frontend(),
            T.pulsed_kernel(), T.zero_span(), T.swept(), T.pulse_sample(), T.overlap()]

    def run(fn):
        lay = SlotLayout(1e6)
        blk = build(spec_from_obj(objs), lay, 5)
        s, p = State(lay), Params(lay)
        fn(blk, s, p, 50)
        return s

    a, b = run(run_reference), run(run_compiled)
    np.testing.assert_allclose(a.S, b.S, rtol=1e-12, atol=0)
    assert a.n == b.n == 50


def test_pulsed_kernel_gate_shorter_than_mode_warns():
    """D46/D50 consistency check (ported from WP-B11's stand-in kernel when the two were merged)."""
    import warnings

    from cranberry.blocks.det.testing import pulsed_kernel
    from cranberry.core.checks import ChkWarning
    from cranberry.core.registry import make_block as _mk

    def build(**over):
        d = pulsed_kernel(**over)
        return _mk(BlockSpec(d["block"], d["name"], d["params"]))

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        short = build(window_s=5e-9, f_rep_hz=20e6, mode_length_s=25e-9)
        ok = build(window_s=30e-9, f_rep_hz=20e6, mode_length_s=25e-9)
        unchecked = build(window_s=5e-9, f_rep_hz=20e6)
    assert [e.code for e in short.chk_log.events] == ["gate_shorter_than_mode"]
    assert not ok.chk_log.events and not unchecked.chk_log.events
    assert any(issubclass(x.category, ChkWarning) for x in w)


def test_pulsed_kernel_charge_profile_ns_mode():
    """M5-11: δ pulses vs an extended charge profile. Closed form (first order, a = 2πB, profile inside the
    window, m ≥ 1): c_m = Σq e^{−a(mT − s)}(1 − e^{−aW}) / Σq (1 − e^{−a(W − s)}). The Borealis-type ns mode
    at 45 MHz (charge |u_opt|², 22 ns slot) gives a larger c_1 than the δ assumption."""
    from cranberry.blocks.source.pulsed_testing import BOREALIS, pulsed_opo

    resp, per, W = {"model": "first_order", "bandwidth_hz": 30e6}, 1 / 45e6, 1 / 45e6
    a = 2 * np.pi * 30e6
    s, q = np.array([2e-9, 9e-9]), np.array([0.7, 0.3])
    c = crosstalk_coefficients(resp, 45e6, W, 3, {"t_s": s.tolist(), "weight": q.tolist()})
    den = np.sum(q * (1 - np.exp(-a * (W - s))))
    for m in (1, 2, 3):
        assert c[m] == pytest.approx(np.sum(q * np.exp(-a * (m * per - s)) * (1 - np.exp(-a * W))) / den, rel=1e-12)
    delta = crosstalk_coefficients(resp, 45e6, W, 3)
    assert np.array_equal(delta, crosstalk_coefficients(resp, 45e6, W, 3, {"t_s": [0.0], "weight": [1.0]}))
    opo = pulsed_opo(BOREALIS, x_ss=2.4)
    u = opo.optimal_mode[0]
    keep = opo.grid < per
    prof = {"t_s": opo.grid[keep].tolist(), "weight": (np.abs(u[keep]) ** 2).tolist()}
    ns = crosstalk_coefficients(resp, 45e6, W, 3, prof)
    assert ns[1] > 1.5 * delta[1]  # the δ assumption underestimates c_1 for ns modes near 45 MHz
    d = T.pulsed_kernel(f_rep_hz=45e6, window_s=W, response=resp, pulse_charge_profile=prof)
    assert np.allclose(mk(d).c[:4], ns, rtol=1e-12)
    with pytest.raises(ValueError):
        crosstalk_coefficients(resp, 45e6, W, 3, {"t_s": [0.0, 1e-9], "weight": [1.0]})
