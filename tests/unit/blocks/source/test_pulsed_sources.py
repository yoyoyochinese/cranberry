"""pulsed_opo_source / opa_supermodes / opo_angle_map (WP-B11): contracts, toggles, frames.

Physics validation is in tests/vt/test_vt_45_pulsed_opo.py, test_vt_32_pulsed_opa.py, test_vt_03_pulsed_cw_limit.py.
"""

import numpy as np
import pytest

import cranberry.core.testing  # noqa: F401  (registers the toy_* driver blocks)
from cranberry.blocks.source.pulsed_testing import BOREALIS, pulsed_opo
from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.registry import make_block
from cranberry.core.state import Params, State

KC = np.pi * BOREALIS["fwhm_hz"]


def test_source_output_contract():
    b = pulsed_opo(BOREALIS)
    so = b.source_output
    n = b.grid.size
    assert so.basis == "temporal" and so.cov.shape == (2 * n, 2 * n) and so.modes.shape == (n, n)
    assert so.phi_opa0_rad == 0.0 and so.pump_path_phase_rad == {"opo": 0.0}
    assert set(so.coef_tables) == {"opo.detuning.squeezed", "opo.detuning.pump"}
    assert 0.0 <= so.tail_into_next_slot < 1e-12  # 6 MHz: the mode is long gone
    assert b(None, None) is None  # offline block: identity step
    assert b.losses() == [] and b.bands() == []


def test_detuning_toggle_off_ignores_detuning_and_pump_cavity_phase_reported():
    on = pulsed_opo(BOREALIS, delta_s_rad_s=0.02 * KC, delta_p_rad_s=0.04 * KC, toggles={"detuning": True})
    off = pulsed_opo(BOREALIS, delta_s_rad_s=0.02 * KC, delta_p_rad_s=0.04 * KC, toggles={"detuning": False})
    assert np.array_equal(off.cov_kernel, pulsed_opo(BOREALIS).cov_kernel)
    kp = np.pi * BOREALIS["pump_fwhm_hz"]
    assert on.source_output.pump_path_phase_rad["opo"] == pytest.approx(-np.arctan(0.04 * KC / kp), rel=1e-12)
    # cov has the static pump phase taken out: θ_off is the rest (signal detuning + chirp), not the full −1.06 d
    th = float(on.source_output.project(off.optimal_mode[0])[2])
    full = (on.angle_coefficients(1.297)[0] * 0.02 + on.angle_coefficients(1.297)[1] * 0.04) * KC
    static = 0.5 * on.source_output.pump_path_phase_rad["opo"]
    assert th == pytest.approx(full - static, abs=2e-4)


def test_pump_resonant_off_uses_input_envelope():
    b = pulsed_opo(BOREALIS, toggles={"pump_resonant": False})
    assert b.kappa_p is None and np.array_equal(b.pump_b(1e7), b.envelope.astype(complex))
    assert b.angle_coefficients(1.297)[1] == 0.0 and b.angle_coefficients(1.297)[2] == 0.0


def test_envelope_callable_equals_rect_sampled():
    tp = BOREALIS["pump_duration_s"]
    rect = pulsed_opo(BOREALIS, dt_kappa=1 / 40)
    dt = rect.dt_s
    env = pulsed_opo(BOREALIS, pump_shape="envelope", pump_envelope=lambda t_s: np.clip((tp - (t_s - 0.5 * dt)) / dt, 0, 1))
    assert np.allclose(env.envelope, rect.envelope, atol=1e-12)


def test_missing_and_bad_params():
    p = dict(BOREALIS)
    del p["x_ss"]
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("pulsed_opo_source", "o", p))
    with pytest.raises(MissingParameterError):
        make_block(BlockSpec("pulsed_opo_source", "o", dict(BOREALIS, toggles=dict(BOREALIS["toggles"], pulsed=False))))
    with pytest.raises(ValueError):
        pulsed_opo(BOREALIS, dt_kappa=0.1)
    with pytest.raises(ValueError):
        pulsed_opo(BOREALIS, toggles={"bogus": True})  # unknown toggle (D34 strict)


def test_opa_supermodes_toggles():
    base = dict(length_m=1e-3, alpha_per_w=1.0, peak_pump_w=0.25, pump_fwhm_s=1e-12, gvm_s_per_m=3e-10,
                beta2_signal_s2_per_m=1e-25, beta2_pump_s2_per_m=4e-25, n_points=64, time_window_s=8e-12, n_z=10,
                phi_opa0_rad=0.3, segment="xtal", psi_p_xtal_rad=0.2)
    b = make_block(BlockSpec("opa_supermodes", "a", dict(base, toggles={"multimode": True, "gvm": False, "gvd": False})))
    assert b._coeffs() == (0.0, 0.0, 0.0)
    so = b.source_output
    assert so.pump_path_phase_rad == {"xtal": 0.2}
    assert float(so.project(b.leading_mode)[2]) == pytest.approx(0.0, abs=1e-9)  # D55-4 frame with φ_OPA,0 = 0.3
    single = make_block(BlockSpec("opa_supermodes", "s", dict(base, toggles={"multimode": False, "gvm": True, "gvd": True})))
    assert single.source_output.basis == "single_mode" and single.schmidt_number == 1.0
    with pytest.raises(ValueError):
        make_block(BlockSpec("opa_supermodes", "x", dict(base, n_points=63, toggles={"multimode": True, "gvm": True, "gvd": True})))


# ---------------------------------------------------------------------------------------------- angle map


def _map_spec(on=True, ratio_mode=True, x_slot=False):
    coef = {"sqz": [[-0.6e-8, -0.2e-8], [-0.5e-8, -0.3e-8]], "pump": [[0.0, -0.1e-8], [0.0, -0.1e-8]]}
    toys = [{"block": "toy_sine", "name": "ds", "params": {"out": "opo.o.Delta_s", "unit": "rad/s", "amp": 1e6, "freq_hz": 1e3,
                                                         "toggles": {"enabled": True}}}]
    params = {"segment": "opo", "beams": ["sqz", "pump"], "in_delta_s": "opo.o.Delta_s", "in_delta_p": None,
              "delta_p_ratio": 2.0, "in_x_ss": None, "x_ss": 1.5, "x_grid": [1.0, 2.0], "coef": coef, "toggles": {"map": on}}
    if not ratio_mode:
        toys.append({"block": "toy_sine", "name": "dp", "params": {"out": "opo.o.Delta_p", "unit": "rad/s", "amp": 3e6,
                                                                 "freq_hz": 2e3, "toggles": {"enabled": True}}})
        params["in_delta_p"] = "opo.o.Delta_p"
        del params["delta_p_ratio"]
    if x_slot:
        toys.append({"block": "toy_sine", "name": "xs", "params": {"out": "pump.x_ss", "unit": "1", "amp": 3.0, "freq_hz": 5e2,
                                                                 "toggles": {"enabled": True}}})
        params["in_x_ss"] = "pump.x_ss"
        del params["x_ss"]
    return [*toys, {"block": "opo_angle_map", "name": "amap", "params": params}]


def _run(spec, n, compiled=False):
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(spec), lay, 1)
    s, p = State(lay), Params(lay)
    (run_compiled if compiled else run_reference)(blk, s, p, n)
    return s, p


def test_angle_map_values_ratio_and_interpolation():
    s, _ = _run(_map_spec(), 37)
    ds = float(s["opo.o.Delta_s"])
    cs, cp = -0.55e-8, -0.25e-8  # x = 1.5: midway
    assert float(s["psi.opo.sqz"]) == pytest.approx(cs * ds + cp * 2.0 * ds, rel=1e-12)
    assert float(s["psi.opo.pump"]) == pytest.approx(-0.1e-8 * 2.0 * ds, rel=1e-12)
    s, _ = _run(_map_spec(ratio_mode=False), 37)
    assert float(s["psi.opo.sqz"]) == pytest.approx(cs * float(s["opo.o.Delta_s"]) + cp * float(s["opo.o.Delta_p"]), rel=1e-12)


def test_angle_map_toggle_off_is_zero():
    s, _ = _run(_map_spec(on=False), 37)
    assert float(s["psi.opo.sqz"]) == 0.0 and float(s["psi.opo.pump"]) == 0.0


@pytest.mark.numba
@pytest.mark.parametrize("kw", [{}, {"ratio_mode": False}, {"x_slot": True}])
def test_angle_map_compiled_equals_reference(kw):
    a, _ = _run(_map_spec(**kw), 500)
    b, _ = _run(_map_spec(**kw), 500, compiled=True)
    for k in ("psi.opo.sqz", "psi.opo.pump"):
        assert abs(float(a[k]) - float(b[k])) <= 1e-12 * max(1.0, abs(float(a[k])))


def test_angle_map_from_pulsed_source_locked_residual():
    """V-P9 through the time loop: Δ_p = 2Δ_s, psi_sqz − ½ psi_pump = −0.543 Δ_s/κ_c (Liu-type pump lock)."""
    src = pulsed_opo(BOREALIS, x_ss=1.3, angle_table_x_ss=[1.2, 1.3, 1.4])
    mp = src.angle_map_params({"sqz": "squeezed", "pump": "pump"})
    spec = _map_spec()
    spec[-1]["params"].update(mp, x_ss=1.3)
    s, _ = _run(spec, 37)
    d = float(s["opo.o.Delta_s"]) / KC
    assert float(s["psi.opo.sqz"]) / d == pytest.approx(-1.063, abs=0.005)
    assert (float(s["psi.opo.sqz"]) - 0.5 * float(s["psi.opo.pump"])) / d == pytest.approx(-0.543, abs=0.005)


def test_angle_map_bad_coef_shape():
    spec = _map_spec()
    spec[-1]["params"]["coef"]["sqz"] = [[0.0, 0.0]]
    with pytest.raises(ValueError):
        build(spec_from_obj(spec), SlotLayout(1e6), 0)
