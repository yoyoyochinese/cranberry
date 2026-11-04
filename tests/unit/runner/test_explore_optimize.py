"""``runner.explore`` / ``runner.optimize`` and the CLI reports of ablate / decompose / convergence, on tiny
budgets (D85-13, prompt §8)."""

import dataclasses
import json
import math

import pytest

import cranberry.core.testing  # noqa: F401  (toy blocks)
import cranberry.path.testing  # noqa: F401
from cranberry import __main__ as cli
from cranberry import config
from cranberry.blocks.meas.testing import two_port_loop_spec
from cranberry.config.overlay import REPO_ROOT
from cranberry.runner import explore
from cranberry.runner.run import run

pytestmark = [pytest.mark.numba, pytest.mark.filterwarnings("ignore::UserWarning")]
SMALL = dict(n_samples=1 << 13, record_decim=8)


def test_combine_effects_keeps_what_mattered():
    spec = two_port_loop_spec(**SMALL)
    out = explore.combine_effects(spec, "ii.port1.std_rad",
                                  toggles=["n_sq.enabled", "n_cp.enabled", "n_com.enabled", "n_lock.enabled"],
                                  threshold=0.05)
    assert set(out["impact"]) == {"n_sq.enabled", "n_cp.enabled", "n_com.enabled", "n_lock.enabled"}
    # n_com drives a segment common to both ports: it cancels in metric (ii) and does not matter
    assert "n_com.enabled" not in out["mattered"] and "n_sq.enabled" in out["mattered"]
    assert out["mattered"] == sorted(out["mattered"], key=lambda t: -abs(out["impact"][t]))
    assert abs(out["residual"]) < 0.05 * out["baseline"]
    assert len(out["ladder"]) == len(out["mattered"]) and out["ladder"][-1][1] < out["baseline"]


def test_compare_seeds_and_differs():
    a = two_port_loop_spec(**SMALL)
    b = dataclasses.replace(a, toggles={"n_sq.enabled": False})
    res = explore.compare_seeds({"all": a, "no_sq": b}, [1, 2, 3], decompose=True)
    s = res["all"]["ii.port1.std_rad"]
    assert s["n"] == 3 and len(s["values"]) == 3 and s["min"] <= s["mean"] <= s["max"] and s["std"] > 0
    comp = res["all"]["ii.port1.squeezed_only_rad2"]
    # decompose sets the group toggles itself (runner.decompose), so the components ignore the toggle override
    assert comp["mean"] > 0 and comp["std"] > 0
    assert res["all"]["ii.port1.decomposition_ok"]["mean"] == 1.0 and "eff.port1.eta_gbs" in res["all"]
    assert explore.differs(res["all"]["ii.port1.std_rad"], res["no_sq"]["ii.port1.std_rad"])
    assert not explore.differs(s, s)
    only = explore.compare_seeds({"all": a}, [1], quantities=["i.port1.std_rad"])
    assert list(only["all"]) == ["i.port1.std_rad"] and only["all"]["i.port1.std_rad"]["std"] == 0.0


@pytest.mark.integration
def test_run_quantities_lock_and_actuator_range_on_refD():
    spec = config.load(REPO_ROOT / "presets/reference/refD_nominal.yaml")
    spec = dataclasses.replace(spec, run=dataclasses.replace(spec.run, duration_s=1e-3))
    ranges = explore.actuator_ranges(spec)
    assert ranges["stretcher915b"] == ("v.amp914", -5.0, 5.0)
    q = explore.run_quantities(run(explore.with_recorded(spec, [s for s, _, _ in ranges.values()])))
    assert q["lock.dsc1.sup.n_loss"] == 0 and q["lock.dsc1.sup.lock_duration_s"] == pytest.approx(1e-3)
    assert 0 < q["range.stretcher915b"] < 1 and q["valid"] == 1.0


# ------------------------------------------------------------------------------------------------
# optimize


def test_optimize_constrained_analytic():
    pytest.importorskip("botorch")
    from cranberry.runner.optimize import optimize

    def f(x):  # minimum of (a-0.8)² + (b-0.2)² subject to a <= 0.5 → (0.5, 0.2)
        return (x["a"] - 0.8) ** 2 + (x["b"] - 0.2) ** 2, {"a_max": x["a"] - 0.5}

    res = optimize(f, {"a": (0, 1), "b": (0, 1)}, n_init=6, n_iter=8, seed=1, raw_samples=32, num_restarts=2)
    assert len(res.history) == 14 and res.best is not None and res.best.feasible
    assert res.best.x["a"] <= 0.5 and res.best.objective < 0.15


def test_optimize_log_bounds_failures_and_maximize():
    pytest.importorskip("botorch")
    from cranberry.runner.optimize import optimize

    def f(x):
        if x["g"] > 50:
            raise RuntimeError("unstable")
        return -abs(math.log10(x["g"]) - 1.0), {}

    res = optimize(f, {"g": (0.1, 1000, "log")}, n_init=4, n_iter=4, minimize=False, seed=0, raw_samples=16,
                   num_restarts=2)
    assert any(e.error and "unstable" in e.error for e in res.history) or all(e.x["g"] <= 50 for e in res.history)
    assert res.best is not None and res.best.objective == max(e.objective for e in res.history if e.error is None)
    with pytest.raises(ValueError):
        optimize(f, {"g": (0, 1, "log")})


def test_runner_evaluator_and_lock_constraint():
    pytest.importorskip("botorch")
    from cranberry.runner.optimize import lock_maintained, optimize, runner_evaluator, servo_bump_limit

    spec = two_port_loop_spec(**SMALL)
    free = run(dataclasses.replace(spec, toggles={}))  # reference for the bump constraint (same loop here)
    ev = runner_evaluator(spec, "i.port1.std_rad",
                          {"lock": lock_maintained, "bump": servo_bump_limit(free, "lock.L1.err", 1e6, nperseg=128)},
                          seeds=[3])
    obj, cons = ev({"pipeline.pid.params.kp": 0.3})
    assert obj == pytest.approx(run(spec).summary()["i.port1.std_rad"])
    assert cons["lock"] < 0 and cons["bump"] < 0
    res = optimize(ev, {"pipeline.pid.params.kp": (0.05, 0.6)}, n_init=2, n_iter=1, raw_samples=8, num_restarts=1)
    assert len(res.history) == 3 and all(e.error is None for e in res.history)


def test_controller_evaluator_call_order():
    from cranberry.runner.optimize import controller_evaluator

    calls = []

    class Fake:
        def set_params(self, p):
            calls.append(("set", dict(p)))
            return dict(p)

        def reset_integrators(self):
            calls.append("reset")

        def engage(self):
            calls.append("engage")

        def disengage(self):
            calls.append("disengage")

        def acquire(self, duration_s, signals):
            calls.append(("acquire", duration_s, tuple(signals)))
            return {"e": [1.0, -1.0], "t_s": [0.0, 1.0]}

    ev = controller_evaluator(Fake(), 0.1, ["e"], lambda d: max(map(abs, d["e"])), {"c": lambda d: -1.0})
    assert ev({"lock1.pid.kp": 2.0}) == (1.0, {"c": -1.0})
    assert calls == [("set", {"lock1.pid.kp": 2.0}), "reset", "engage", ("acquire", 0.1, ("e",)), "disengage"]


# ------------------------------------------------------------------------------------------------
# CLI reports


def test_cli_ablate_decompose_convergence_print_and_save(tmp_path, capsys):
    spec = two_port_loop_spec(**SMALL)
    spec = dataclasses.replace(spec, measure={**spec.measure, "convergence": {"quantities": {"i.port1.std_rad": 1.0},
                                                                              "seeds": [1, 2]}})
    path = tmp_path / "toy.yaml"
    path.write_text(config.dump_yaml(spec))

    assert cli.main(["ablate", str(path), "--toggle", "n_sq.enabled", "-q", "ii.port1.std_rad",
                     "--save", str(tmp_path / "abl.json")]) == 0
    out = capsys.readouterr().out
    assert "n_sq.enabled" in out and "ii.port1.std_rad" in out
    saved = json.loads((tmp_path / "abl.json").read_text())
    assert saved["impact"]["n_sq.enabled"]["ii.port1.std_rad"]["delta"] < 0 and saved["valid"]["baseline"] is True

    assert cli.main(["decompose", str(path), "--save", str(tmp_path / "dec.json")]) == 0
    out = capsys.readouterr().out
    assert "port 1" in out and "squeezed_only" in out and "ok" in out
    assert json.loads((tmp_path / "dec.json").read_text())["1"]["status"] == "ok"

    assert cli.main(["convergence", str(path), "--what", "chunk", "--what", "seed",
                     "--save", str(tmp_path / "conv.json")]) == 0
    out = capsys.readouterr().out
    assert "chunk" in out and "passed=True" in out and "i.port1.std_rad" in out
    assert json.loads((tmp_path / "conv.json").read_text())["chunk"]["passed"] is True
