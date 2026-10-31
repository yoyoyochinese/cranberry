"""Backward-compatible block additions made by WP-I1 (Ref A integration):

* ``loss_budget`` (new offline block): one LossElement per config entry; toggle; validation.
* ``loss_elements`` (``eom_phase`` / PZT actuators): ``insertion_loss_lin`` (the name a YAML
  ``insertion_loss_db`` gets on load, G16-i) accepted as the alternative of ``insertion_loss_db``.
* ``laser_source`` ``oob_pairs`` entry ``{target, dnu_fc_hz, tau_fc_s, tau_ref_s}`` (T05 §6; corrected in
  Phase 5, M1-A-01 / D65): white FM → a ``phase`` row = continuous total − expected value of the discrete
  loop chain; flicker / peaks → a ``freq`` row |H_FC e^{−i2πfδ} − 1|² S_ν above f_N.
* ``filter_cavity`` ``in_delay_s``: the low pass acts on φ_L(t − τ) so the LO is H_FC φ_L(t − τ_b).
"""

import math

import numpy as np
import pytest

import cranberry.blocks.noise.testing as nt
from cranberry import config
from cranberry.blocks.err.testing import rig
from cranberry.blocks.noise.laser import fc_pair_white_variances, s_nu_model
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


C_LIGHT = 299792458.0
TAU_LO = (5.0 * 1.468 + 0.3 * 1.0003) / C_LIGHT                                   # refA: τ_b of lo
TAU_SQ = (1.0 * 1.0003 + 0.3 * 1.0003 + 0.2 * 1.0003 + 0.045 * 2.2 + 1.0 * 1.0003) / C_LIGHT  # τ_b of sqz


def test_laser_oob_pair_filter_cavity_rows():
    dnu = 3.3e6
    pair = {"target": "FC", "dnu_fc_hz": dnu, "tau_fc_s": TAU_LO, "tau_ref_s": TAU_SQ}
    root, _, _ = rig([nt.laser(on=("fm_white", "peaks"), oob_pairs=[pair, {"target": "d", "tau_s": 1e-7}])], 2e7,
                     external={})
    decls = {(d.targets[0], d.quantity): d for d in root.bands()}
    assert set(decls) == {("FC", "phase"), ("FC", "freq"), ("d", "freq")}
    # white FM: continuous − discrete expectation (no PSD above f_N: the loop already holds the aliased part)
    v_cont, v_disc = fc_pair_white_variances(8e3, dnu, TAU_LO, TAU_SQ, 5e-8)
    assert decls["FC", "phase"].variance_outside == pytest.approx(v_cont - v_disc, rel=1e-12) and v_cont > v_disc
    # band-limited generators (peaks): |H_FC e^{−i2πfδ} − 1|² S_ν above f_N
    s_pk = s_nu_model(nt.laser(on=("peaks",))["params"])
    f = np.array([1e7, 3e7, 1e9])
    h = np.exp(-2j * np.pi * f * (TAU_LO - TAU_SQ)) / (1 + 2j * f / dnu)
    np.testing.assert_allclose(decls["FC", "freq"].psd_outside(f), np.abs(h - 1) ** 2 * s_pk(f), rtol=1e-14)
    s_nu = s_nu_model(nt.laser(on=("fm_white", "peaks"))["params"])
    np.testing.assert_allclose(decls["d", "freq"].psd_outside(f), 4 * np.sin(np.pi * f * 1e-7) ** 2 * s_nu(f), rtol=1e-14)
    # white only → no freq row on FC; laser noise off → no row
    root, _, _ = rig([nt.laser(on=("fm_white",), oob_pairs=[pair])], 2e7, external={})
    assert [(d.quantity, d.targets) for d in root.bands()] == [("phase", ("FC",))]
    root, _, _ = rig([nt.laser(on=(), oob_pairs=[pair])], 2e7, external={})
    assert root.bands() == []
    with pytest.raises(MissingParameterError):
        rig([nt.laser(oob_pairs=[{"target": "x"}])], 1e6, external={})[0].bands()


def test_fc_pair_white_variances_closed_forms():
    h0, dnu = 200 / math.pi, 3.3e6
    # δ = 0: continuous π h0/Δν; discrete 2π²h0h(1 + a)/(4(1 − a)) (no interpolation, exact geometric sum)
    for fs in (1e7, 2e7, 1e9):
        a = math.exp(-math.pi * dnu / fs)
        vc, vd = fc_pair_white_variances(h0, dnu, 0.0, 0.0, 1 / fs)
        assert vc == pytest.approx(math.pi * h0 / dnu, rel=1e-14)
        assert vd == pytest.approx(2 * math.pi**2 * h0 / fs * (1 + a) / (4 * (1 - a)), rel=1e-12)
    # continuous with delay = full-band integral ∫|H e^{−iωδ} − 1|² h0/f² df, both signs of δ
    for d in (TAU_LO - TAU_SQ, TAU_SQ - TAU_LO):
        vc, _ = fc_pair_white_variances(h0, dnu, max(d, 0.0), max(-d, 0.0), 5e-8)
        x = np.geomspace(1e-2, 1e13, 2_000_001)
        g = np.abs(np.exp(-2j * np.pi * x * d) / (1 + 2j * x / dnu) - 1) ** 2 * h0 / x**2
        assert vc == pytest.approx(float(np.trapezoid(g * x, np.log(x))), rel=1e-4)
    # refA (M1-A-01 table): truth 8.173e-5 rad² (9.041 mrad); discrete chain at 10 / 20 / 40 MHz
    vals = [fc_pair_white_variances(h0, dnu, TAU_LO, TAU_SQ, 1 / fs) for fs in (1e7, 2e7, 4e7, 8e7, 1.6e8)]
    assert vals[0][0] == pytest.approx(8.1734e-5, rel=1e-4)
    np.testing.assert_allclose([v[1] for v in vals], [7.4655e-5, 7.3063e-5, 7.4887e-5, 7.8445e-5, 7.9873e-5], rtol=1e-4)


def test_fc_pair_white_variances_matches_direct_chain():
    """Brute-force sum (no geometric tail) and a Monte Carlo of the discrete chain itself: random walk with
    increments of variance 2π²h0h, linear interpolation at τ_fc / τ_ref, the cavity low pass."""
    import scipy.signal

    from cranberry.path.segments import frac_delay

    h0, dnu, h = 200 / math.pi, 1.0e6, 1 / 2e7
    _, vd = fc_pair_white_variances(h0, dnu, TAU_LO, TAU_SQ, h)
    a, n = math.exp(-math.pi * dnu * h), 4000
    kf, mf = frac_delay(TAU_LO, h)
    kr, mr = frac_delay(TAU_SQ, h)
    x, r = np.zeros(n), np.zeros(n)
    x[kf], x[kf + 1], r[kr], r[kr + 1] = 1 - mf, mf, 1 - mr, mr
    y = scipy.signal.lfilter([(1 - a) / 2] * 2, [1, -a], x)
    assert vd == pytest.approx(2 * math.pi**2 * h0 * h * float(np.sum(np.cumsum(y - r) ** 2)), rel=1e-10)
    rng = np.random.default_rng(5)
    phi = np.cumsum(rng.standard_normal(2_000_000)) * math.sqrt(2 * math.pi**2 * h0 * h)
    lo = scipy.signal.lfilter([(1 - a) / 2] * 2, [1, -a], np.convolve(phi, x[: kf + 2])[: phi.size])
    dl = lo - np.convolve(phi, r[: kr + 2])[: phi.size]
    assert float(np.var(dl[10_000:])) == pytest.approx(vd, rel=0.03)


def test_fc_pair_reference_cavity():
    """Phase 5 M3B-06 (Ref B: SHG on the pump vs MCC on the LO): ``dnu_ref_hz`` puts the reference beam through a
    cavity too. Continuous = full-band ∫|H_FC e^{−iωδ} − H_ref|² h0/f² df (both signs of δ, equal and unequal
    linewidths); Δν_ref → ∞ is the one-cavity pair; the discrete value is the plain sum of the chain; the freq row
    is |H_FC e^{−i2πfδ} − H_ref|² S_ν."""
    import scipy.signal

    from cranberry.path.segments import frac_delay

    h0, h = 200 / math.pi, 1 / 2e6
    x = np.geomspace(1e-2, 1e13, 2_000_001)
    for dnu, dnu_ref in ((1e6, 1e6), (1e6, 1e7), (3e6, 2e5)):
        for d in (TAU_LO - TAU_SQ, TAU_SQ - TAU_LO):
            tf, tr = max(d, 0.0), max(-d, 0.0)
            vc, _ = fc_pair_white_variances(h0, dnu, tf, tr, h, dnu_ref)
            g = np.abs(np.exp(-2j * np.pi * x * d) / (1 + 2j * x / dnu) - 1 / (1 + 2j * x / dnu_ref)) ** 2 * h0 / x**2
            assert vc == pytest.approx(float(np.trapezoid(g * x, np.log(x))), rel=1e-4), (dnu, dnu_ref, d)
    assert fc_pair_white_variances(h0, 1e6, TAU_LO, TAU_SQ, h, 1e15)[0] == pytest.approx(
        fc_pair_white_variances(h0, 1e6, TAU_LO, TAU_SQ, h)[0], rel=1e-6)
    dnu, dnu_ref, n = 1e6, 2e5, 4000
    _, vd = fc_pair_white_variances(h0, dnu, TAU_LO, TAU_SQ, h, dnu_ref)
    a, ar = math.exp(-math.pi * dnu * h), math.exp(-math.pi * dnu_ref * h)
    kf, mf = frac_delay(TAU_LO, h)
    kr, mr = frac_delay(TAU_SQ, h)
    xx, r = np.zeros(n), np.zeros(n)
    xx[kf], xx[kf + 1], r[kr], r[kr + 1] = 1 - mf, mf, 1 - mr, mr
    y = scipy.signal.lfilter([(1 - a) / 2] * 2, [1, -a], xx) - scipy.signal.lfilter([(1 - ar) / 2] * 2, [1, -ar], r)
    assert vd == pytest.approx(2 * math.pi**2 * h0 * h * float(np.sum(np.cumsum(y) ** 2)), rel=1e-10)
    pair = {"target": "FC", "dnu_fc_hz": dnu, "tau_fc_s": TAU_LO, "tau_ref_s": TAU_SQ, "dnu_ref_hz": dnu_ref}
    root, _, _ = rig([nt.laser(on=("fm_white", "peaks"), oob_pairs=[pair])], 2e6, external={})
    decls = {(d.targets[0], d.quantity): d for d in root.bands()}
    s_pk = s_nu_model(nt.laser(on=("peaks",))["params"])
    f = np.array([1e7, 3e7, 1e9])
    hh = np.exp(-2j * np.pi * f * (TAU_LO - TAU_SQ)) / (1 + 2j * f / dnu) - 1 / (1 + 2j * f / dnu_ref)
    np.testing.assert_allclose(decls["FC", "freq"].psd_outside(f), np.abs(hh) ** 2 * s_pk(f), rtol=1e-14)
    vc, vd = fc_pair_white_variances(8e3, dnu, TAU_LO, TAU_SQ, h, dnu_ref)
    assert decls["FC", "phase"].variance_outside == pytest.approx(max(vc - vd, 0.0), rel=1e-12)


def test_fc_pair_overcounting_chain_warns_and_declares_zero():
    """δ = 0 (no interpolation): the cavity discretisation over-represents π h0/Δν by ≈ ε²/12 — no negative
    ledger rows exist, so the row is 0 with a warning (ponytail ceiling in laser.bands)."""
    root, _, _ = rig([nt.laser(on=("fm_white",), oob_pairs=[{"target": "FC", "dnu_fc_hz": 3.3e6}])], 2e7,
                     external={})
    with pytest.warns(UserWarning, match="exceeds the continuous"):
        (d,) = root.bands()
    assert d.variance_outside == 0.0


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
