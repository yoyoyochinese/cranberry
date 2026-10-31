"""VT-F3 (WP-F3 / M0): double registration in BandLedger and LossRegistry overlap is a build-time error.

Sources: PLAN §6.2 VT-F3, G05, G06, T04 V13 / T12 V13 context. Canonical location is
tests/vt/test_vt_f3_double_registration.py (PLAN §6.1); kept here because WP-F3 owns
tests/unit/offline only (proposed move in the WP-F3 report).
"""

import pytest

from cranberry.offline import DoubleRegistrationError
from cranberry.offline.band_ledger import BandDecl, BandLedger
from cranberry.offline.loss_registry import LossElement, LossRegistry

pytestmark = pytest.mark.vt


def band(source, targets, f_lo=1.0, f_hi=5e5, quantity="phase", var=1e-6):
    return BandDecl(source=source, quantity=quantity, targets=targets, f_lo_hz=f_lo, f_hi_hz=f_hi, unit="rad", variance_outside=var)


def test_band_ledger_same_source_same_band_is_error():
    led = BandLedger([band("noise.arm1", ("arm1",))])
    with pytest.raises(DoubleRegistrationError):
        led.declare(band("noise.arm1", ("arm1",)))
    assert len(led.decls) == 1


def test_band_ledger_overlapping_and_disjoint_bands_of_one_source_are_errors():
    led = BandLedger([band("laser", ("lo1", "sqz1"), f_hi=1e3)])
    with pytest.raises(DoubleRegistrationError):
        led.declare(band("laser", ("sqz1",), f_lo=10.0, f_hi=2e3))  # overlapping band, shared target
    with pytest.raises(DoubleRegistrationError):
        led.declare(band("laser", ("sqz1",), f_lo=1e3, f_hi=5e5))  # disjoint band still double counts "outside"
    with pytest.raises(DoubleRegistrationError):
        BandLedger([band("x", ("a",)), band("x", ("a",))])  # constructor path too


def test_band_ledger_distinct_rows_are_fine():
    led = BandLedger()
    led.declare(band("noise.arm1", ("arm1",)))
    led.declare(band("noise.arm2", ("arm1",)))  # other source
    led.declare(band("noise.arm1", ("arm2",)))  # other target
    led.declare(band("noise.arm1", ("arm1",), quantity="dT"))  # other quantity
    assert len(led.decls) == 4


def _ov(id, kind, applies_to, **loc):
    return LossElement(id=id, eta=0.98, kind=kind, applies_to=applies_to, source="[T12]", enabled=True, **loc)


def test_loss_registry_duplicate_id():
    reg = LossRegistry(paths={"sqz": {"segments": ["s"]}})
    reg.add(_ov("s.conn", "connector", {"gbs"}, segment="s"))
    with pytest.raises(DoubleRegistrationError):
        reg.add(_ov("s.conn", "connector", {"gbs"}, segment="s"))


@pytest.mark.parametrize("kind", ["mode_overlap_spatial", "mode_overlap_temporal", "lo_shape"])
def test_loss_registry_overlap_counted_twice_on_one_path_pair(kind):
    # same observer → caught at add time without any path
    reg = LossRegistry()
    reg.add(_ov("bhd1.a", kind, {"bhd"}, observer="bhd1"))
    with pytest.raises(DoubleRegistrationError):
        reg.add(_ov("bhd1.b", kind, {"bhd", "gbs"}, observer="bhd1"))
    assert [e.id for e in reg.elements] == ["bhd1.a"]  # rolled back

    # different locations on the same path → caught against the defined path
    reg = LossRegistry(paths={"sqz1": {"segments": ["src", "arm1"], "observers": ["bhd1"]}})
    reg.add(_ov("src.overlap", kind, {"bhd"}, segment="src"))
    with pytest.raises(DoubleRegistrationError):
        reg.add(_ov("bhd1.overlap", kind, {"bhd"}, observer="bhd1"))
    assert len(reg.elements) == 1

    # path defined after the elements → caught at define_path
    reg = LossRegistry()
    reg.add(_ov("src.overlap", kind, {"bhd"}, segment="src"))
    reg.add(_ov("bhd1.overlap", kind, {"bhd"}, observer="bhd1"))
    with pytest.raises(DoubleRegistrationError):
        reg.define_path("sqz1", ["src"], ["bhd1"])
    assert "sqz1" not in reg.paths

    # disabled registrations count too ("a second registration is an error")
    reg = LossRegistry(paths={"sqz1": {"segments": ["src"], "observers": ["bhd1"]}})
    reg.add(_ov("src.overlap", kind, {"bhd"}, segment="src"))
    with pytest.raises(DoubleRegistrationError):
        reg.add(
            LossElement(id="bhd1.alt", eta=0.9, kind=kind, applies_to={"bhd"}, source="x", enabled=False, observer="bhd1")
        )


def test_loss_registry_overlaps_that_are_not_double_counts():
    reg = LossRegistry(
        paths={
            "sqz1": {"segments": ["src"], "observers": ["bhd1", "gbs1"]},
            "sqz2": {"segments": ["src2"], "observers": ["bhd2"]},
        }
    )
    reg.add(_ov("gbs1.overlap", "mode_overlap_spatial", {"gbs"}, observer="gbs1"))  # GBS pair
    reg.add(_ov("bhd1.overlap", "mode_overlap_spatial", {"bhd"}, observer="bhd1"))  # BHD pair
    reg.add(_ov("bhd1.temporal", "mode_overlap_temporal", {"bhd"}, observer="bhd1"))  # other kind
    reg.add(_ov("bhd1.lo_shape", "lo_shape", {"bhd"}, observer="bhd1"))  # other kind
    reg.add(_ov("bhd2.overlap", "mode_overlap_spatial", {"bhd"}, observer="bhd2"))  # other path
    reg.add(_ov("src.conn1", "connector", {"bhd"}, segment="src"))  # non-overlap kinds may repeat
    reg.add(_ov("src.conn2", "connector", {"bhd"}, segment="src"))
    assert reg.eta("sqz1", "bhd") == pytest.approx(0.98**5)
    assert reg.eta("sqz1", "gbs") == pytest.approx(0.98)
