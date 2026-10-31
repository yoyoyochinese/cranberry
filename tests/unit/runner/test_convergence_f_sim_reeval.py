"""``runner.convergence('f_sim')`` re-evaluates expressions of ``run.f_sim_hz`` and keeps the device rates
of ``keep_rate`` blocks (WP-F2c, D63-9). ``runner.run.run`` is replaced by a recorder: only the specs the
check builds are tested here (the numerics are covered by the integration convergence tests)."""

import pytest

from cranberry import config
from cranberry.config import FrozenExpressionsWarning
from cranberry.config.yaml_io import loads_yaml
from cranberry.runner.convergence import convergence

TEXT = """
run: {duration_s: 1.0e-4, f_sim_hz: 2.0e7, chunk_samples: 64, seed: 5, seed_mode: per_point}
pipeline:
  - {block: anything, name: dl, params: {tau_s: !expr "1.0e-6 - 1 / run.f_sim_hz"}}
  - {block: anything, name: pid, rate_div: 2, params: {kp: 1.0}}
  - {block: anything, name: adc, rate_div: 1, params: {}}
measure:
  convergence: {quantities: {q: 1.0}}
"""


class _Result:
    def __init__(self, spec):
        self.spec = spec

    def summary(self):
        return {"q": 0.0}


@pytest.fixture
def ran(monkeypatch):
    import cranberry.runner.run as runmod

    seen = []

    def fake_run(spec, **kw):
        seen.append(spec)
        return _Result(spec)

    monkeypatch.setattr(runmod, "run", fake_run)
    return seen


def _leaf(spec, name):
    return next(b for b in spec.pipeline.children if b.name == name)


def test_f_sim_check_re_evaluates_tau_s_at_x1_x2_x4(ran):
    spec = config.load(loads_yaml(TEXT))
    rep = convergence(spec, ["f_sim"], f_sim_factors=(2, 4))
    assert rep.checks["f_sim"]["passed"] is True
    base, x2, x4 = ran
    for s, k in ((base, 1), (x2, 2), (x4, 4)):
        assert s.run.f_sim_hz == 2e7 * k
        assert _leaf(s, "dl").params["tau_s"] == pytest.approx(1e-6 - 1 / (2e7 * k), abs=1e-18)
        assert s.run.seed == 5  # same root seed (common random numbers), not a per-point seed
        assert _leaf(s, "pid").rate_div == 2  # no keep_rate: rate_div unchanged


def test_keep_rate_scales_rate_div_argument_and_config(ran):
    spec = config.load(loads_yaml(TEXT))
    convergence(spec, ["f_sim"], f_sim_factors=(2, 4), keep_rate=["pid"])
    _, x2, x4 = ran
    assert (_leaf(x2, "pid").rate_div, _leaf(x4, "pid").rate_div) == (4, 8)
    assert (_leaf(x2, "adc").rate_div, _leaf(x4, "adc").rate_div) == (1, 1)
    assert _leaf(x4, "dl").params["tau_s"] == pytest.approx(1e-6 - 1 / 8e7, abs=1e-18)
    ran.clear()
    cfg = config.load(loads_yaml(TEXT.replace("{q: 1.0}}", "{q: 1.0}, keep_rate: [pid, adc]}")))
    convergence(cfg, ["f_sim"], f_sim_factors=(2,))
    assert (_leaf(ran[1], "pid").rate_div, _leaf(ran[1], "adc").rate_div) == (4, 2)
    with pytest.raises(KeyError, match="nope"):
        convergence(spec, ["f_sim"], f_sim_factors=(2,), keep_rate=["nope"])
    with pytest.raises(ValueError, match="integer factors"):
        convergence(spec, ["f_sim"], f_sim_factors=(1.5,), keep_rate=["pid"])


def test_spec_without_unresolved_tree_warns_and_freezes(ran):
    spec = config.load(loads_yaml(TEXT))
    bare = config.spec_from_tree(config.spec_to_tree(spec))
    with pytest.warns(FrozenExpressionsWarning):
        convergence(bare, ["f_sim"], f_sim_factors=(2,))
    assert _leaf(ran[1], "dl").params["tau_s"] == _leaf(spec, "dl").params["tau_s"]
