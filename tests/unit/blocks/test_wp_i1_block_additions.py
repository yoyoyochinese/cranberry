"""Backward-compatible block additions made by WP-I1 (Ref A integration):

* ``loss_budget`` (new offline block): one LossElement per config entry; toggle; validation.
* ``loss_elements`` (``eom_phase`` / PZT actuators): ``insertion_loss_lin`` (the name a YAML
  ``insertion_loss_db`` gets on load, G16-i) accepted as the alternative of ``insertion_loss_db``.
* ``laser_source`` ``oob_pairs`` entry ``{target, dnu_fc_hz}``: |1 − H_FC|² S_ν above f_N (T05 §6).
* ``filter_cavity`` ``in_delay_s``: the low pass acts on φ_L(t − τ) so the LO is H_FC φ_L(t − τ_b).
"""

import math

import numpy as np
import pytest
from scipy.integrate import quad

import cranberry.blocks.noise.testing as nt
from cranberry import config
from cranberry.blocks.err.testing import rig
from cranberry.blocks.noise.laser import s_nu_model
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled, run_reference
from cranberry.offline.loss_registry import BeamPath, LossRegistry

pytestmark = pytest.mark.numba


# ------------------------------------------------------------------------------------- loss_budget


def budget(enabled=True, **over):
    e = {"id": "optics", "eta": 0.99, "kind": "component", "applies_to": ["gbs", "bhd"], "segment": "sq",
         "observer": None, "source": "[test]"}
    e.update(over)
    return {"block": "loss_budget", "name": "lb", "params": {"elements": [e], "toggles": {"enabled": enabled}}}


def test_loss_budget_registers_elements_and_toggle():
    root, _, _ = rig([budget()], 1e6, external={})
    (el,) = root.losses()
    assert (el.id, el.eta, el.kind, el.segment, el.observer) == ("lb.optics", 0.99, "component", "sq", None)
    reg = LossRegistry(paths={"sqz": BeamPath("sqz", ("sq",))})
    reg.add_all(root.losses())
    assert reg.eta("sqz", "bhd") == pytest.approx(0.99) and reg.eta("sqz", "gbs") == pytest.approx(0.99)
    root, _, _ = rig([budget(enabled=False)], 1e6, external={})
    reg = LossRegistry(paths={"sqz": BeamPath("sqz", ("sq",))})
    reg.add_all(root.losses())
    assert reg.eta("sqz", "bhd") == 1.0 and not root.losses()[0].enabled


@pytest.mark.parametrize("over, err", [({"observer": "bhd"}, ValueError), ({"segment": None}, ValueError),
                                       ({"eta": 1.2}, ValueError), ({"kind": "mirror"}, ValueError)])
def test_loss_budget_validation(over, err):
    with pytest.raises(err):
        rig([budget(**over)], 1e6, external={})
    bad = budget()
    del bad["params"]["elements"][0]["source"]
    with pytest.raises(MissingParameterError):
        rig([bad], 1e6, external={})


# ---------------------------------------------------------------------------- insertion_loss_lin


def eom(**loss):
    p = {"in": "v.u", "gain_mode": "vpi", "v_pi_v": 4.0, "polarity": 1.0, "poles_hz": [], "discretisation": "zoh",
         "prewarp_hz": None, "drive_min_v": -1.0, "drive_max_v": 1.0, "segment": "s", "driver": "eom", "n_sub": 1,
         "sub": [0, 1], "ring_len": 0, "loss_applies_to": ["gbs", "bhd"], "loss_source": "[test]",
         "toggles": {"bandwidth": False, "saturation": False, "insertion_loss": True}, **loss}
    return {"block": "eom_phase", "name": "eom", "params": p}


def test_insertion_loss_lin_equals_db():
    a, _, _ = rig([eom(insertion_loss_db=3.0)], 1e6, external={"v.u": "V"})
    b, _, _ = rig([eom(insertion_loss_lin=10 ** 0.3)], 1e6, external={"v.u": "V"})
    assert a.losses()[0].eta == pytest.approx(b.losses()[0].eta, rel=1e-15)
    c, _, _ = rig([eom(insertion_loss_lin=None)], 1e6, external={"v.u": "V"})
    assert c.losses() == []
    for bad, err in ((eom(), MissingParameterError), (eom(insertion_loss_db=1.0, insertion_loss_lin=1.2),
                                                       MissingParameterError), (eom(insertion_loss_lin=0.5), ValueError)):
        blk, _, _ = rig([bad], 1e6, external={"v.u": "V"})
        with pytest.raises(err):
            blk.losses()


def test_insertion_loss_db_from_yaml_is_accepted():
    """A loaded config renames insertion_loss_db → insertion_loss_lin (config.db); the block takes it."""
    tree = {"run": {"duration_s": 1e-5, "f_sim_hz": 1e6, "chunk_samples": 16, "seed": 1},
            "pipeline": [{"block": "ramp", "name": "u", "params": {"out": "v.u", "waveform": "triangle",
                                                                   "amplitude_v": 0.0, "period_s": 1.0,
                                                                   "offset_v": 0.0, "toggles": {"enabled": True}}},
                         eom(insertion_loss_db=3.0)]}
    spec = config.load(tree)
    leaf = [c for c in spec.pipeline.children if c.name == "eom"][0]
    assert "insertion_loss_lin" in leaf.params and "insertion_loss_db" not in leaf.params
    config.validate(spec)


# ----------------------------------------------------------------------------- laser oob transfer


def test_laser_oob_pair_filter_cavity_transfer():
    dnu = 3.3e6
    root, _, _ = rig([nt.laser(on=("fm_white",), oob_pairs=[{"target": "FC", "dnu_fc_hz": dnu},
                                                           {"target": "d", "tau_s": 1e-7}])], 2e7, external={})
    decls = {d.targets[0]: d for d in root.bands()}
    s_nu = s_nu_model(nt.laser(on=("fm_white",))["params"])
    f = np.array([1e7, 3e7, 1e9])
    x2 = (2 * f / dnu) ** 2
    np.testing.assert_allclose(decls["FC"].psd_outside(f), x2 / (1 + x2) * s_nu(f), rtol=1e-14)
    np.testing.assert_allclose(decls["d"].psd_outside(f), 4 * np.sin(np.pi * f * 1e-7) ** 2 * s_nu(f), rtol=1e-14)
    assert decls["FC"].f_hi_hz == 1e7
    # the full-band integral of |1 − H_FC|² h0/f² is π h0/Δν (T05 §6)
    h0 = 8e3
    tot = 2 * h0 / dnu * quad(lambda xx: 1 / (1 + xx * xx), 0, np.inf)[0]  # x = 2f/Δν
    assert tot == pytest.approx(math.pi * h0 / dnu, rel=1e-9)
    with pytest.raises(MissingParameterError):
        rig([nt.laser(oob_pairs=[{"target": "x"}])], 1e6, external={})[0].bands()


# --------------------------------------------------------------------------- filter cavity delay


def fc(**over):
    p = {"segment": "fc", "beam": "lo", "dnu_fc_hz": 1e6, "t0": 1.0, "source": "test", "in_phase": "laser.dx1.phi",
         "toggles": {"lock_residual": False, "laser_lowpass": True}}
    p.update(over)
    return {"block": "filter_cavity", "name": "fc1", "params": p}


def test_filter_cavity_input_delay():
    """ψ = H φ(t − τ) − φ(t − τ): equals the undelayed output shifted by τ (integer + fractional part,
    linear interpolation); absent / zero delay keeps the original behaviour; compiled = reference."""
    fs, n = 1e7, 300
    rng = np.random.default_rng(3)
    x = np.cumsum(rng.standard_normal(n)) * 1e-3

    def trace(compiled, **over):
        root, st, pa = rig([fc(**over)], fs, external={"laser.dx1.phi": "rad"})
        out = np.zeros(n)
        for k in range(n):
            st["laser.dx1.phi"] = x[k]
            (run_compiled if compiled else run_reference)(root, st, pa, 1)
            out[k] = float(st["psi.fc.lo"])
        return out

    base = trace(False)
    np.testing.assert_array_equal(trace(False, in_delay_s=0.0), base)
    np.testing.assert_array_equal(trace(False, in_delay_s=None), base)
    # integer delay of 3 samples: output = base shifted by 3 (zeros before)
    d3 = trace(False, in_delay_s=3 / fs)
    np.testing.assert_allclose(d3[3:], base[:-3], rtol=0, atol=1e-15)
    # fractional 2.25 samples: linear interpolation of the delayed input → filter is linear
    xi = np.concatenate([np.zeros(3), x])
    xd = 0.75 * xi[1:-2] + 0.25 * xi[:-3]  # x(n − 2.25) on the grid, zeros before t = 0
    a = math.exp(-math.pi * 1e6 / fs)
    y, prev, out = 0.0, 0.0, np.zeros(n)
    for k in range(n):
        y = a * y + (1 - a) * 0.5 * (xd[k] + prev)
        prev = xd[k]
        out[k] = y - xd[k]
    np.testing.assert_allclose(trace(False, in_delay_s=2.25 / fs), out, rtol=0, atol=1e-15)
    np.testing.assert_allclose(trace(True, in_delay_s=2.25 / fs), out, rtol=0, atol=1e-14)
