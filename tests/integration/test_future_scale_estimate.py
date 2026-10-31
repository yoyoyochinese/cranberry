"""WP-I5 (M6, PLAN §1.2, D13): 32-port scale estimate — memory and time per simulated port from the 2-port F_liu20
configuration (reference port 0 + port k), extrapolated to 32 ports. A report, not a 32-port run (D13: build the
32-port configuration in Phase 6 when needed). The numbers are printed (``pytest -s``) and sanity-checked."""

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
    assert 0 < est["port_k_elems"] < est["state_elems"]
    assert est["port_0_elems"] < est["port_k_elems"]                # port 0: short fibre, fewer subsegments
    assert 0 < est["us_per_sample_per_port"] < est["us_per_sample"]
    assert est["target_state_bytes"] < 1e8                           # memory is not the limit (≲ 100 MB)
    assert est["target_us_per_sample"] > est["us_per_sample"]
