"""mode_overlap_table (WP-B11): D51-6 axes, toggles, CW source rejected. Physics: tests/vt/test_vt_32_pulsed_opa.py."""

import numpy as np
import pytest

from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.registry import make_block

T0 = 0.6e-12
SRC = {"gaussian": dict(T0_s=T0, r=1.0, n_points=512, window_s=40 * T0, phi_opa0_rad=0.0)}


def table(**kw):
    p = dict(source=SRC, dt_grid_s=[-0.2e-12, 0.0, 0.2e-12], dL_grid_m=[0.0, 20.0], beta2_s2_per_m=-21.7e-27,
             beta3_s3_per_m=0.127e-39, pad=2, toggles={"dispersion": True, "jitter_fast": False})
    tg = dict(p["toggles"], **kw.pop("toggles", {}))
    return make_block(BlockSpec("mode_overlap_table", "ov", {**p, **kw, "toggles": tg}))


def test_axes_and_dispersion_toggle():
    on, off = table().table, table(toggles={"dispersion": False}).table
    assert on.axes == ("dt_s", "dL_m") and on.values.shape == (3, 2, 4)
    assert np.allclose(off.values[:, 0], off.values[:, 1], atol=1e-15)  # δL has no effect without dispersion
    assert on.values[1, 1, 3] < 1.0 and on.values[1, 0, 3] == pytest.approx(1.0, abs=1e-12)
    # timing and length are separate (D51-6): δτ alone does not rotate, δL alone gives arg κ = ½ arctan χ ≠ 0
    assert on.values[0, 0, 2] == pytest.approx(0.0, abs=1e-12) and abs(on.values[1, 1, 2]) > 1e-3


def test_cw_source_rejected_and_jitter_params_required():
    cw = {"block": "single_mode_reference", "name": "sm", "params": dict(loss=0.0, alpha_per_w=9.06, pump_w=0.64,
                                                                          theta_rms_rad=0.0, phi_opa0_rad=0.0, average="exact")}
    with pytest.raises(ValueError):
        table(source=cw).table
    with pytest.raises(MissingParameterError):
        table(toggles={"jitter_fast": True})
