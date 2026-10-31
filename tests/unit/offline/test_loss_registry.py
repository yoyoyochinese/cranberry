"""offline.loss_registry: η per beam path and target (G06), Ref A/B/C budgets (T12 §3.1, §9 V4–V6, V13)."""

import json

import numpy as np
import pytest

from cranberry.core.units import noise_db
from cranberry.offline import MissingRegistrationError
from cranberry.offline import gaussian as g
from cranberry.offline.loss_registry import LossElement, LossRegistry
from cranberry.offline.source_output import SourceOutput

GB = frozenset({"gbs", "bhd"})


def el(id, eta, kind, applies_to, **loc):
    return LossElement(id=id, eta=eta, kind=kind, applies_to=applies_to, source=loc.pop("source", "[T12]"), enabled=loc.pop("enabled", True), **loc)


def ref_a_registry(bhd_only_enabled=True):
    """Ref A budget (T12 §3.1): waveguide 2 %, optics 1 %, QE 1 %, circuit noise 0.2 %, mode mismatch 4 %."""
    reg = LossRegistry(paths={"sqz1": {"segments": ["opa", "free_space"], "observers": ["bhd1"]}})
    reg.add(el("opa.waveguide", 0.98, "component", GB, segment="opa", source="[Ref A §4.2]"))
    reg.add(el("free_space.optics", 0.99, "component", GB, segment="free_space", source="[Ref A §4.2]"))
    reg.add(el("bhd1.qe", 0.99, "detector", {"bhd"}, observer="bhd1", source="[Ref A §3.1]", enabled=bhd_only_enabled))
    reg.add(el("bhd1.circuit", 0.998, "electronic", {"bhd"}, observer="bhd1", source="[Ref A §4.2]", enabled=bhd_only_enabled))
    reg.add(el("bhd1.overlap", 0.96, "mode_overlap_spatial", {"bhd"}, observer="bhd1", source="[Ref A §4.2]", enabled=bhd_only_enabled))
    return reg


def test_ref_a_budget():
    """T12 V4: 1 − 0.98·0.99·0.99·0.998·0.96 = 7.98 %; GBS-side (no BHD-only losses) 2.98 %."""
    reg = ref_a_registry()
    assert reg.loss("sqz1", "bhd") == pytest.approx(0.0798, abs=1e-4)
    assert reg.eta("sqz1", "bhd") == pytest.approx(0.98 * 0.99 * 0.99 * 0.998 * 0.96, rel=1e-14)
    assert reg.eta("sqz1", "gbs") == pytest.approx(0.98 * 0.99, rel=1e-14)
    # the lumped L = 8 % then gives §7-1 (T12 V4): −10.15 dB / +20.56 dB (±0.01 dB)
    vm, vp, _ = SourceOutput.single_mode_reference(loss=0.08, alpha_per_w=9.06, pump_w=0.64, phi_opa0_rad=0.0).project(None)
    chi = g.chi_gaussian(9e-3)
    assert float(noise_db(g.jitter_average(vm, vp, chi, delta0_rad=0.0))) == pytest.approx(-10.15, abs=0.01)
    assert float(noise_db(vp)) == pytest.approx(20.56, abs=0.01)
    # and the registry's own 7.98 % is within 0.01 dB of it
    vm2, vp2, _ = SourceOutput.single_mode_reference(
        loss=reg.loss("sqz1", "bhd"), alpha_per_w=9.06, pump_w=0.64, phi_opa0_rad=0.0
    ).project(None)
    assert float(noise_db(g.jitter_average(vm2, vp2, chi, delta0_rad=0.0))) == pytest.approx(-10.15, abs=0.015)


def test_ref_b_budget():
    """T12 V5: ρ = 0.123/0.1268, ζ = 0.99, ξ = 0.988 (η_mm = ξ²), η_PD = 0.998 → 1 − ρζξ²η = 0.0645."""
    reg = LossRegistry(paths={"sqz": {"segments": ["opo", "path"], "observers": ["bhd"]}})
    reg.add(el("opo.escape", 0.123 / (0.123 + 0.0038), "component", GB, segment="opo", source="[Ref B §4]"))
    reg.add(el("path.prop", 0.99, "component", GB, segment="path", source="[Ref B §4]"))
    reg.add(el("bhd.vis", 0.988**2, "mode_overlap_spatial", {"bhd"}, observer="bhd", source="[Ref B Eq. (2)]"))
    reg.add(el("bhd.qe", 0.998, "detector", {"bhd"}, observer="bhd", source="[Ref B §4]"))
    assert reg.loss("sqz", "bhd") == pytest.approx(0.0645, abs=1e-4)


def test_ref_c_budget():
    """T12 V6: 1 − 0.85 × 0.92 = 0.218 (fit 21 %); 1 km × 0.2 dB/km → 4.50 %."""
    fib = LossElement.from_il_db(
        "fiber.1km", 0.2, kind="fiber", applies_to=GB, source="[Ref C §3]", enabled=True, segment="fiber"
    )
    assert 1 - fib.eta == pytest.approx(0.0450, abs=1e-3)
    assert fib.il_db == pytest.approx(0.2, rel=1e-12)
    reg = LossRegistry(paths={"sqz": {"segments": ["table", "fiber"], "observers": ["bhd"]}})
    reg.add(el("table.pre_fiber", 0.85, "component", {"bhd"}, segment="table", source="[Ref C §3]"))
    reg.add(fib)
    reg.add(el("fiber.couple", 0.92 / fib.eta, "connector", GB, segment="fiber", source="[Ref C §3–4] [AE]"))
    assert reg.loss("sqz", "bhd") == pytest.approx(0.218, abs=1e-3)
    assert reg.eta("sqz", "gbs") == pytest.approx(0.92, rel=1e-12)


def test_gbs_vs_bhd_separation_v13():
    """T12 V13: turning off the BHD-only losses makes the BHD value equal the GBS-input value."""
    off = ref_a_registry(bhd_only_enabled=False)
    assert off.eta("sqz1", "bhd") == pytest.approx(off.eta("sqz1", "gbs"), rel=0, abs=0)
    # disabled elements stay in the record
    assert len(off.breakdown("sqz1", "bhd")) == 5
    assert [e["enabled"] for e in off.breakdown("sqz1", "bhd")].count(False) == 3


def test_paths_and_lo_target():
    reg = LossRegistry()
    reg.define_path("sqz1", ["arm1", "common"], ["bhd1"])
    reg.define_path("lo1", ["lo_arm", "common"], ["bhd1"])
    reg.add(el("arm1.conn", 0.95, "connector", GB, segment="arm1"))
    reg.add(el("lo_arm.conn", 0.9, "connector", {"lo"}, segment="lo_arm"))
    reg.add(el("common.splice", 0.99, "splice", {"gbs", "bhd", "lo"}, segment="common"))
    reg.add(el("bhd1.qe", 0.97, "detector", {"bhd", "lo"}, observer="bhd1"))
    assert reg.eta("sqz1", "bhd") == pytest.approx(0.95 * 0.99 * 0.97)
    assert reg.eta("sqz1", "gbs") == pytest.approx(0.95 * 0.99)
    assert reg.eta("lo1", "lo") == pytest.approx(0.9 * 0.99 * 0.97)
    np.testing.assert_allclose(reg.eta_array(["sqz1", "lo1"], "bhd"), [0.95 * 0.99 * 0.97, 0.99 * 0.97])
    assert reg.eta("sqz1", "probe") == 1.0  # nothing registered for the probe target
    with pytest.raises(KeyError):
        reg.eta("nope", "bhd")
    with pytest.raises(ValueError):
        reg.eta("sqz1", "esa")
    with pytest.raises(ValueError):
        reg.add(LossElement(id="floating", eta=0.9, kind="component", applies_to={"gbs"}, source="x", enabled=True))


def test_mandatory_registration_helper():
    """Modulators / actuators / probe couplers on the squeezed path must be registered (prompt §5.5, D10)."""
    reg = LossRegistry(paths={"sqz1": {"segments": ["arm1"], "observers": []}, "lo1": {"segments": ["lo"]}})
    required = ["arm1.915b", "arm1.probe_in"]
    assert reg.missing_registrations("sqz1", required) == required
    with pytest.raises(MissingRegistrationError):
        reg.require_registered("sqz1", required)
    reg.add(LossElement.from_il_db("arm1.915b.il", 0.05, kind="component", applies_to=GB, source="[T12 §3.3] [AE, UNVERIFIED]", enabled=True, segment="arm1"))
    # registered on the wrong path does not count
    reg.add(el("arm1.probe_in.il", 0.99, "component", GB, segment="lo"))
    assert reg.missing_registrations("sqz1", required) == ["arm1.probe_in"]
    reg2 = LossRegistry(paths={"sqz1": {"segments": ["arm1"]}})
    reg2.add(el("arm1.915b", 1.0, "component", GB, segment="arm1", enabled=False))  # η = 1 still counts
    reg2.add(el("arm1.probe_in.pbs", 0.89, "component", GB, segment="arm1"))
    reg2.require_registered("sqz1", required)
    # a prefix must end at a '.' boundary
    assert reg2.missing_registrations("sqz1", ["arm1.915"]) == ["arm1.915"]


def test_to_json_record():
    reg = ref_a_registry()
    rec = reg.to_json()
    text = json.dumps(rec, ensure_ascii=False)
    assert "[Ref A §4.2]" in text
    assert rec["eta"]["sqz1"]["bhd"]["eta"] == pytest.approx(reg.eta("sqz1", "bhd"))
    assert rec["eta"]["sqz1"]["gbs"]["elements"] == ["opa.waveguide", "free_space.optics"]
    assert rec["paths"]["sqz1"] == {"id": "sqz1", "segments": ["opa", "free_space"], "observers": ["bhd1"]}
    assert {e["id"] for e in rec["elements"]} == {e.id for e in reg.elements}
    zero = el("dead", 0.0, "component", {"gbs"}, segment="opa")
    assert zero.to_json()["il_db"] is None
