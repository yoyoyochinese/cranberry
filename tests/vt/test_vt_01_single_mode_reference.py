"""VT-01 (prompt §7-1, WP-B1 part): single-mode reference model.

Sources: S §4 VT-01; T01 V-T01-1; T04 V1; T12 V4; D48 (§7-1 reading); PLAN §6.2–§6.3.
Tolerances: analytic values R′− = −10.148 dB, R′+ = +20.556 dB to ±1e−3 dB (L = 8 %, α = 9.06 /W,
P = 0.64 W, θ = 9 mrad); small-angle vs exact Gaussian average < 1e−3 dB; D48 window
R′− = −10.15 ± 0.10 dB, R′+ = +20.56 ± 0.10 dB (lumped waveguide loss, as §7-1 defines it).
The ``phase_to_quadrature`` half (WP-F4) maps a δ(t) series; here the averages come from
``offline.gaussian`` (conv §6), which is what that block applies.
"""

import pytest

from cranberry.blocks.source.single_mode import r_prime
from cranberry.core.block import BlockSpec
from cranberry.core.registry import make_block
from cranberry.core.units import V_VAC, lin_to_db
from cranberry.offline.loss_registry import LossElement, LossRegistry

REF_A = dict(loss=0.08, alpha_per_w=9.06, pump_w=0.64, theta_rms_rad=0.009, phi_opa0_rad=0.0)
D48 = {"minus": (-10.15, 0.10), "plus": (20.56, 0.10)}


def in_window(rm_db, rp_db):
    return abs(rm_db - D48["minus"][0]) <= D48["minus"][1] and abs(rp_db - D48["plus"][0]) <= D48["plus"][1]


@pytest.mark.parametrize("average", ["exact", "small_angle"])
def test_analytic_values(average):
    b = make_block(BlockSpec("single_mode_reference", "refA.sqz", dict(REF_A, average=average)))
    rm, rp = b.r_prime_db()
    assert rm == pytest.approx(-10.148, abs=1e-3)
    assert rp == pytest.approx(20.556, abs=1e-3)
    assert in_window(rm, rp)


def test_small_angle_vs_exact():
    ex = make_block(BlockSpec("single_mode_reference", "a", dict(REF_A, average="exact"))).r_prime_db()
    sa = make_block(BlockSpec("single_mode_reference", "b", dict(REF_A, average="small_angle"))).r_prime_db()
    assert abs(ex[0] - sa[0]) < 1e-3 and abs(ex[1] - sa[1]) < 1e-3


def test_lumped_opa_with_registry_in_d48_window():
    """Same budget built the framework way: lumped opa_waveguide (L_wg 2 %) + registered downstream
    losses (optics 1 %, PD 1 %, circuit 0.2 %, mode mismatch 4 %; Ref A §4.2 → L = 7.98 %)."""
    toggles = dict(pump_loss=False, phase_mismatch=False, pump_depletion=False, pump_rin=False)
    opa = make_block(BlockSpec("opa_waveguide", "refA.opa", dict(
        model="lumped", alpha_per_w=9.06, pump_w=0.64, pump_ref="output", length_m=0.045, phi_opa0_rad=0.0,
        segment="xtal", psi_p_xtal_rad=0.0, gain_table_pump_w=None, waveguide_loss=0.02,
        loss_source="[Ref A §4.2]", toggles=toggles)))
    down = [("optics", 0.99, "component", "seg.optics", None), ("pd_qe", 0.99, "detector", None, "bhd1"),
            ("circuit", 0.998, "electronic", None, "bhd1"), ("mm", 0.96, "mode_overlap_spatial", None, "bhd1")]
    reg = LossRegistry(opa.losses(), paths={"sqz": {"segments": ["xtal", "seg.optics"], "observers": ["bhd1"]}})
    for i, eta, kind, seg, obs in down:
        reg.add(LossElement(id=i, eta=eta, kind=kind, applies_to={"bhd"}, source="[Ref A §4.2]", enabled=True,
                            segment=seg, observer=obs))
    eta = reg.eta("sqz", "bhd")
    assert 1 - eta == pytest.approx(0.0798, abs=1e-4)
    vm, vp, _ = opa.source_output.project(None)
    rm, rp = (eta * v / V_VAC + 1 - eta for v in (vm, vp))
    rpm, rpp = r_prime(rm, rp, 0.009, "exact")
    assert in_window(float(lin_to_db(rpm)), float(lin_to_db(rpp)))
