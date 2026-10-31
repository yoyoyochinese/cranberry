"""MOD blocks: ``eom_phase`` (π u/V_π, bandwidth, insertion loss via losses()), ``dither_modulator``
(timeseries / analytic χ = J0(2m)), ``aom_shift`` (toggles)."""

import math

import numpy as np
import pytest
from scipy import special

import cranberry.core.testing  # noqa: F401
from cranberry.core.compile import run_reference
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.offline.loss_registry import BeamPath, LossRegistry

FS = 1.0e6
SRC = {"block": "ramp", "name": "u", "params": {"out": "v.u", "waveform": "triangle", "amplitude_v": 2.0,
                                               "period_s": 1e-4, "offset_v": 0.0, "toggles": {"enabled": True}}}


def make(obj, fs=FS):
    lay = SlotLayout(fs)
    blk = build(spec_from_obj(obj), lay, 6)
    return blk, State(lay), Params(lay)


def trace(obj, n, slots):
    blk, st, pr = make(obj)
    out = {s: np.empty(n) for s in slots}
    for i in range(n):
        run_reference(blk, st, pr, 1)
        for s in slots:
            out[s][i] = float(np.sum(st[s]))
    return blk, out


def eom(toggles, **kw):
    p = {"in": "v.u", "gain_mode": "vpi", "v_pi_v": 4.0, "polarity": 1.0, "poles_hz": [50e3],
         "discretisation": "zoh", "prewarp_hz": None, "drive_min_v": -1.0, "drive_max_v": 1.0,
         "segment": "sqz1", "driver": "eom", "n_sub": 1, "sub": [0, 1], "ring_len": 0,
         "insertion_loss_db": 3.0, "loss_applies_to": ["gbs", "bhd"], "loss_source": "UNVERIFIED (MPZ-LN-10, T16 §3.2)",
         "toggles": toggles}
    p.update(kw)
    return {"block": "eom_phase", "name": "eom", "params": p}


def test_eom_static_gain_and_toggles():
    off = {"bandwidth": False, "saturation": False, "insertion_loss": True}
    _, tr = trace([SRC, eom(off)], 200, ["v.u", "act.eom.dpsi", "drv.sqz1.eom"])
    np.testing.assert_allclose(tr["act.eom.dpsi"], math.pi * tr["v.u"] / 4.0, rtol=1e-15)
    np.testing.assert_array_equal(tr["drv.sqz1.eom"], tr["act.eom.dpsi"])
    # Ref C check: 200 mVpp at V_π = 4 V → 0.0785 rad amplitude (T16 §3.2)
    assert math.pi * 0.1 / 4.0 == pytest.approx(0.0785, abs=5e-5)
    _, tr = trace([SRC, eom({**off, "saturation": True})], 200, ["v.u", "act.eom.dpsi"])
    np.testing.assert_allclose(tr["act.eom.dpsi"], math.pi * np.clip(tr["v.u"], -1, 1) / 4.0, rtol=1e-15)
    _, tr = trace([SRC, eom({**off, "bandwidth": True})], 200, ["v.u", "act.eom.dpsi"])
    assert not np.allclose(tr["act.eom.dpsi"], math.pi * tr["v.u"] / 4.0)  # first-order lag acts


def test_eom_insertion_loss_on_squeezed_path_only():
    """The EOM registers its IL on its segment; the registry applies it to the squeezed path through that
    segment and not to a path elsewhere (G06, prompt §5.5): 3 dB caps η at 0.5."""
    blk, _, _ = make([SRC, eom({"bandwidth": False, "saturation": False, "insertion_loss": True})])
    reg = LossRegistry(paths={"sqz": BeamPath("sqz", ("src", "sqz1")), "lo": BeamPath("lo", ("lo_seg",))})
    reg.add_all(blk.losses())
    assert reg.eta("sqz", "gbs") == pytest.approx(10 ** -0.3)
    assert reg.eta("lo", "gbs") == 1.0
    blk, _, _ = make([SRC, eom({"bandwidth": False, "saturation": False, "insertion_loss": False})])
    assert not blk.losses()[0].enabled


def dither(mode, on=True):
    return {"block": "dither_modulator", "name": "dith", "params": {
        "dither_mode": mode, "m_rad": 0.0785, "f_m_hz": 1e5, "phase_rad": 0.0, "segment": "lo1",
        "driver": "dither", "n_sub": 2, "sub": [0, 2], "ring_len": 0, "toggles": {"enabled": on}}}


def test_dither_modes():
    blk, tr = trace([dither("timeseries")], 100, ["act.dith.dpsi", "drv.lo1.dither"])
    t = np.arange(100) / FS
    np.testing.assert_allclose(tr["act.dith.dpsi"], 0.0785 * np.sin(2 * np.pi * 1e5 * t), atol=1e-15)
    np.testing.assert_allclose(tr["drv.lo1.dither"] / 2, tr["act.dith.dpsi"], atol=1e-15)  # Σ density / N_sub = δψ (D61)
    assert np.sqrt(np.mean(tr["act.dith.dpsi"] ** 2)) == pytest.approx(0.0555, abs=5e-4)  # Ref C rms
    blk, tr = trace([dither("analytic")], 50, ["act.dith.dpsi"])
    np.testing.assert_array_equal(tr["act.dith.dpsi"], 0.0)
    assert blk.children[0].chi() == pytest.approx(special.j0(2 * 0.0785))
    blk, tr = trace([dither("timeseries", on=False)], 50, ["act.dith.dpsi"])
    np.testing.assert_array_equal(tr["act.dith.dpsi"], 0.0)
    assert blk.children[0].chi() == 1.0


def test_aom_toggle_and_loss():
    p = {"in": "v.u", "order": 1, "input": "frequency", "k_hz_per_v": 1e3, "polarity": 1.0, "segment": "pr1",
         "driver": "eom", "n_sub": 1, "sub": [0, 1], "ring_len": 0, "insertion_loss_db": 1.0,
         "loss_applies_to": ["probe"], "loss_source": "UNVERIFIED (AOM efficiency)",
         "toggles": {"enabled": False, "insertion_loss": True}}
    blk, tr = trace([SRC, {"block": "aom_shift", "name": "aom", "params": p}], 50, ["act.aom.dpsi"])
    np.testing.assert_array_equal(tr["act.aom.dpsi"], 0.0)
    assert blk.losses()[0].applies_to == {"probe"}
