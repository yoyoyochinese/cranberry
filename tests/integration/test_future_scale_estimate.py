"""WP-I5 (M6, PLAN §1.2, D13): 32-port scale estimate — memory and time per simulated port from the MARGINAL cost of a
real third port added to the 2-port F_liu20 configuration (reference port 0 + port k), extrapolated to 32 ports
(Phase 5 M6-09: the review measured 1.85 µs/sample per port → ≈ 59 µs/sample, D74-6 as corrected by D78). A report,
not a 32-port run (D13: build the 32-port configuration in Phase 6 when needed). The numbers are printed
(``pytest -s``) and sanity-checked."""

import importlib.util
from pathlib import Path

import pytest

_sp = importlib.util.spec_from_file_location("future_support", Path(__file__).parent / "future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)

pytestmark = [pytest.mark.numba, pytest.mark.integration, pytest.mark.slow]


def test_32_port_scale_estimate():
    est = fs.scale_estimate("F_liu20")
    print("\n[32-port scale estimate] " + fs.fmt_scale(est))
    assert est["metric_ii_ports"] == [1, 2]                         # the third port is real (own lock, metric (ii))
    assert 0 < est["us_per_sample_per_port"] < est["us_per_sample_2"]
    assert est["target_us_per_sample"] == pytest.approx(est["us_per_sample_2"] + 30 * est["us_per_sample_per_port"])
    # a port costs most of the 2-port loop (port 0 is short): the 32-port loop is ≳ 10× the 2-port one (review: ×17)
    assert est["target_us_per_sample"] > 10 * est["us_per_sample_2"]
    assert 0 < est["port_elems"] and est["target_state_bytes"] < 2e8   # memory is not the limit (draw buffers dominate)
