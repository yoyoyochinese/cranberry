"""WP-B7 LOSS blocks: per_budget (PER accumulation, loss registration, toggle) and coil_birefringence
(915B length, birefringent drift and its sweep, bend PER, bend loss). Values: T09 §1.3, §3.4, D51-8."""

import math

import numpy as np
import pytest

from cranberry.blocks.loss.coil_birefringence import make_coil_birefringence
from cranberry.blocks.loss.per_budget import (angle_from_per_db, make_per_budget, per_db_from_angle, per_db_from_x,
                                              x_from_per_db)
from cranberry.core.block import MissingParameterError
from cranberry.core.compile import run_compiled
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves, spec_from_obj
from cranberry.core.state import Params, State
from cranberry.offline.loss_registry import LossRegistry

REFD = [  # Ref D preset-like values (D51-8; connector PER [AE] T09 §3.5)
    {"id": "polz", "type": "polarizer", "per_db": 28.0, "source": "[D51-8] ER ≥ 28 dB"},
    {"id": "fc1", "type": "connector", "theta_deg": 3.0, "source": "[AE]"},
    {"id": "c1", "type": "coupler", "per_db": 25.0, "source": "[D51-8] PNH1550R5A2 PER ≥ 25 dB"},
    {"id": "coil", "type": "coil", "x": 7.8e-4, "source": "coil_birefringence"},
    {"id": "fib", "type": "fiber", "h_per_m": 1e-5, "L_m": 10.0, "source": "[Corning]"},
]


def pb(elements=REFD, loss=None, on=True):
    return make_per_budget({"elements": elements, "loss": loss, "toggles": {"pol_crosstalk": on}}, name="pb")


def test_conversions_round_trip():
    for per in (20.0, 25.0, 35.2):
        assert per_db_from_angle(angle_from_per_db(per)) == pytest.approx(per, rel=1e-12)
        assert per_db_from_x(x_from_per_db(per)) == pytest.approx(per, rel=1e-12)
    assert per_db_from_angle(math.radians(1.0)) == pytest.approx(35.2, abs=0.05)  # T09 §1.3 table
    assert per_db_from_angle(math.radians(5.7)) == pytest.approx(20.0, abs=0.05)


def test_power_sum_and_angles():
    b = pb()
    want = (x_from_per_db(28) + math.sin(math.radians(3)) ** 2 + x_from_per_db(25) + 7.8e-4 + math.tanh(1e-4))
    assert b.x_total == pytest.approx(want, rel=1e-14)
    assert b.per_total_db == pytest.approx(per_db_from_x(want), rel=1e-14)
    assert b.theta_rad("fc1") == pytest.approx(math.radians(3.0), rel=1e-12)
    assert b.losses() == []


def test_loss_registration_and_toggle():
    loss = {"segment": "pm_in", "applies_to": ["gbs"], "source": "[T09 §1.3]"}
    reg = LossRegistry(pb(loss=loss).losses(), paths={"sq": {"segments": ["pm_in"]}})
    assert reg.eta("sq", "gbs") == pytest.approx(1 - pb().x_total, rel=1e-14)
    off = pb(loss=loss, on=False)
    assert off.x_total == 0.0 and off.per_total_db == math.inf and off.theta_rad("fc1") == 0.0
    reg = LossRegistry(off.losses(), paths={"sq": {"segments": ["pm_in"]}})
    assert reg.eta("sq", "gbs") == 1.0 and len(reg.breakdown("sq", "gbs")) == 1


def test_per_budget_validation():
    with pytest.raises(ValueError, match="exactly one"):
        pb([{"id": "a", "type": "connector", "theta_deg": 1.0, "per_db": 30.0, "source": "x"}])
    with pytest.raises(ValueError, match="type"):
        pb([{"id": "a", "type": "magic", "x": 0.1, "source": "x"}])
    with pytest.raises(MissingParameterError):
        make_per_budget({"elements": [], "toggles": {"pol_crosstalk": True}})


COIL = {"n_turns": 38, "turn_length_m": 0.22, "bend_radius_m": 0.015, "cladding_radius_m": 62.5e-6,
        "bend_coeff": 0.136, "beat_length_m": 4.2e-3, "wavelength_m": 1548.8e-9, "K_T_rad_per_m_K": -2.05,
        "bend_loss_db": 0.1, "orientation": "random", "segment": "a1_915b", "applies_to": ["probe"], "source": "[T09 §3.4, UNVERIFIED]"}


def coil(on=(True, True), **kw):
    return make_coil_birefringence({**COIL, **kw, "toggles": {"pol_bend_coil": on[0], "bend_loss": on[1]}},
                                   name="c915b")


def test_915b_length_and_drift_sweep():
    c = coil()
    assert c.length_m == pytest.approx(8.36, rel=1e-12)
    assert c.dpsi_dT_rad_per_K == pytest.approx(-17.138, abs=1e-3)  # T09 §4.1: −17 rad/K
    assert coil(K_T_rad_per_m_K=-1.7).dpsi_dT_rad_per_K == pytest.approx(-14.2, abs=0.05)  # sweep ends
    assert coil(K_T_rad_per_m_K=-2.9).dpsi_dT_rad_per_K == pytest.approx(-24.2, abs=0.05)
    assert 2 * math.pi / abs(c.dpsi_dT_rad_per_K) == pytest.approx(0.37, abs=0.01)  # one fringe cycle per 0.37 K


@pytest.mark.parametrize("R,delta0,xt_db", [(0.015, 3.2e-3, -31.1), (0.020, 1.8e-3, -36.1)])
def test_bend_per(R, delta0, xt_db):
    c = coil(bend_radius_m=R)
    assert c.delta0 == pytest.approx(delta0, rel=0.02)
    assert 10 * math.log10(c.x_bend) == pytest.approx(xt_db, abs=0.1)  # T09 §3.4 [AE]
    assert c.per_bend_db == pytest.approx(per_db_from_x(c.x_bend))


def test_coil_toggles_and_loss():
    reg = LossRegistry(coil().losses(), paths={"pr": {"segments": ["a1_915b"]}})
    assert reg.eta("pr", "probe") == pytest.approx(10 ** (-0.01))
    off = coil(on=(False, False))
    assert off.x_bend == 0.0 and off.per_bend_db == math.inf
    reg = LossRegistry(off.losses(), paths={"pr": {"segments": ["a1_915b"]}})
    assert reg.eta("pr", "probe") == 1.0 and len(reg.breakdown("pr", "probe")) == 1


@pytest.mark.numba
def test_offline_blocks_run_in_a_compiled_pipeline():
    pipe = [{"block": "per_budget", "name": "pb", "params": {"elements": REFD, "loss": None,
                                                             "toggles": {"pol_crosstalk": True}}},
            {"block": "coil_birefringence", "name": "c915b",
             "params": {**COIL, "toggles": {"pol_bend_coil": True, "bend_loss": True}}}]
    lay = SlotLayout(1e6)
    blk = build(spec_from_obj(pipe), lay, 0)
    s, p = State(lay), Params(lay)
    run_compiled(blk, s, p, 10)
    assert s.n == 10 and [leaf.kind for leaf in iter_leaves(blk)] == ["per_budget", "coil_birefringence"]


def test_coil_aligned_orientation():
    """Coherent bends (T09 §3.4): mean over the birefringent phases 4·n·ε², peak 16·n²·ε² when all
    phases match (L_b = 2πR, straight = πR), rotating with temperature; random orientation unchanged."""
    alpha = 0.3
    c = coil(orientation="aligned", bend_angle_rad=alpha)
    eps2 = (c.delta0 * math.sin(2 * alpha)) ** 2
    assert c.x_bend == pytest.approx(c.x_bend_at(0.0), rel=1e-15)
    xs = [c.x_bend_at(dT) for dT in np.linspace(0.0, 400.0, 40001)]  # many cycles of every phase
    assert np.mean(xs) == pytest.approx(4 * 38 * eps2, rel=0.03)
    assert np.max(xs) <= 16 * 38 ** 2 * eps2 * (1 + 1e-12)
    R = 0.015
    peak = coil(orientation="aligned", bend_angle_rad=math.pi / 4, beat_length_m=2 * math.pi * R, turn_length_m=4 * math.pi * R)
    assert peak.x_bend == pytest.approx(16 * 38 ** 2 * peak.delta0 ** 2, rel=1e-9)
    assert coil(orientation="aligned", bend_angle_rad=alpha, on=(False, True)).x_bend_at(1.0) == 0.0
    assert coil(bend_angle_rad=alpha).x_bend == pytest.approx(2 * 38 * coil().delta0 ** 2, rel=1e-15)
    with pytest.raises(MissingParameterError):
        coil(orientation="aligned")
    with pytest.raises(ValueError, match="orientation"):
        coil(orientation="twisted")
    with pytest.raises(ValueError, match="2πR"):
        coil(orientation="aligned", bend_angle_rad=alpha, turn_length_m=0.05)


YAML_OFFLINE = """
run: {duration_s: 1.0e-5, f_sim_hz: 1.0e+6, chunk_samples: 16, seed: 1, seed_mode: common}
pipeline:
  - {block: ramp, name: u, params: {out: v.u, waveform: triangle, amplitude_v: 0.0, period_s: 1.0,
                                    offset_v: 0.0, toggles: {enabled: true}}}
offline:
  blocks:
    - block: per_budget
      name: pb
      params:
        elements:
          - {id: polz, type: polarizer, per_db: 28.0, source: "[D51-8]"}
          - {id: c1, type: coupler, per_db: 25.0, source: "[D51-8]"}
          - {id: fc1, type: connector, theta_deg: 3.0, source: "[AE]"}
        loss: {segment: pm_in, applies_to: [gbs], source: "[T09 §1.3]"}
        toggles: {pol_crosstalk: true}
    - block: coil_birefringence
      name: c915b
      params: {n_turns: 38, turn_length_m: 0.22, bend_radius_m: 0.015, cladding_radius_m: 6.25e-5,
               bend_coeff: 0.136, beat_length_m: 4.2e-3, wavelength_m: 1.5488e-6, K_T_rad_per_m_K: -2.05,
               bend_loss_db: 0.1, orientation: random, segment: a1_915b, applies_to: [probe],
               source: "[T09 §3.4, UNVERIFIED]", toggles: {pol_bend_coil: true, bend_loss: true}}
"""


def test_yaml_db_names_build_like_python(tmp_path):
    """M1D-03 / D61-5 (protected: YAML ↔ Python equivalence): config loading renames per_db → per_lin and
    bend_loss_db → bend_loss_lin (D56-1); both blocks build from the loaded spec and equal the Python ones,
    also after dump_yaml → reload."""
    import yaml

    from cranberry import config
    from cranberry.core.registry import ensure_loaded, make_block
    from cranberry.runner.run import offline_block_specs

    path = tmp_path / "loss.yaml"
    path.write_text(YAML_OFFLINE)
    spec = config.load(path)
    config.validate(spec)
    specs = {b.name: b for b in offline_block_specs(spec)}
    assert "bend_loss_lin" in specs["c915b"].params and "bend_loss_db" not in specs["c915b"].params
    assert "per_lin" in specs["pb"].params["elements"][0] and "per_db" not in specs["pb"].params["elements"][0]
    ensure_loaded()
    pb_y, c_y = make_block(specs["pb"]), make_block(specs["c915b"])
    els = [{"id": "polz", "type": "polarizer", "per_db": 28.0, "source": "s"},
           {"id": "c1", "type": "coupler", "per_db": 25.0, "source": "s"},
           {"id": "fc1", "type": "connector", "theta_deg": 3.0, "source": "s"}]
    pb_p = pb(els, loss={"segment": "pm_in", "applies_to": ["gbs"], "source": "s"})
    assert pb_y.x == pytest.approx(pb_p.x, rel=1e-13)
    assert pb_y.losses()[0].eta == pytest.approx(pb_p.losses()[0].eta, rel=1e-13)
    c_p = coil()
    assert c_y.x_bend == pytest.approx(c_p.x_bend, rel=1e-13)
    assert c_y.losses()[0].eta == pytest.approx(c_p.losses()[0].eta, rel=1e-13)
    assert c_y.losses()[0].eta == pytest.approx(10 ** (-0.01), rel=1e-13)
    # dump → reload: the linear names carry no dB suffix and reload unchanged
    again = {b.name: b for b in offline_block_specs(config.load(yaml.safe_load(config.dump_yaml(spec))))}
    assert make_block(again["pb"]).x == pb_y.x
    assert make_block(again["c915b"]).losses()[0].eta == c_y.losses()[0].eta


def test_lin_names_exactly_one_and_range():
    """D34: exactly one of the dB / linear names; the linear values must be physical."""
    tg = {"toggles": {"pol_bend_coil": True, "bend_loss": True}}
    no_db = {k: v for k, v in COIL.items() if k != "bend_loss_db"}
    with pytest.raises(MissingParameterError):
        coil(bend_loss_lin=1.02)  # both names
    with pytest.raises(MissingParameterError):
        make_coil_birefringence({**no_db, **tg})  # neither
    with pytest.raises(ValueError, match="loss factor"):
        make_coil_birefringence({**no_db, "bend_loss_lin": 0.9, **tg})
    with pytest.raises(ValueError, match="exactly one"):
        pb([{"id": "a", "type": "coupler", "per_db": 25.0, "per_lin": 316.2, "source": "x"}])
    with pytest.raises(ValueError, match="per_lin"):
        pb([{"id": "a", "type": "coupler", "per_lin": 0.0, "source": "x"}])
