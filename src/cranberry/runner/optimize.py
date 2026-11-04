"""Constrained Bayesian optimisation of named parameters with BoTorch (prompt §8, Phase 6 step 5, D85-13).

The loop is backend-agnostic: it sees an evaluator ``evaluate(x) -> (objective, {constraint: c})`` with
``x`` = {parameter name: value} and every constraint feasible when ``c <= 0``. Two evaluators:

* :func:`runner_evaluator` — directly through the runner: each name is a sweep path (``runner.sweep``,
  e.g. ``pipeline.pid0f.params.kp``, ``offline.params.f_c_hz``); the point is expanded from the spec's
  unresolved tree, so dependent ``!expr`` values follow (WP-F2c), and run for each seed (objective: mean
  over seeds; constraints: worst seed).
* :func:`controller_evaluator` — through ``hardware.ControllerInterface`` only (``set_params`` →
  ``reset_integrators`` → ``engage`` → ``acquire`` → ``disengage``), so the same loop drives
  ``SimController`` or real hardware (Moku:Delta / DSC1). Names are the canonical runtime names
  (``lock1.pid.kp``); the applied (clipped) values are what is recorded.

Constraint helpers for the §8 constraints: :func:`lock_maintained` (no lock loss in any supervisor, run
valid) and :func:`servo_bump_limit` (peak in-loop / free-running |S|² ≤ limit; the free-running run
does not depend on the controller parameters and is computed once).

:func:`optimize` — Sobol initial design, then per iteration one point maximising the log constrained
expected improvement (``LogConstrainedExpectedImprovement``) of a multi-output ``SingleTaskGP``
(objective + constraints, inputs scaled to [0, 1]^d — ``"log"`` bounds in log10 —, outputs
standardised). A failed evaluation (exception or non-finite value) is kept in the history with its
error and enters the GP as the worst objective seen and a violated constraint. Needs the ``opt`` extra
(``uv sync --extra opt``; botorch, gpytorch, torch).
"""

from __future__ import annotations

import dataclasses
import math
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.hardware.interface import ControllerInterface
    from cranberry.runner.run import RunResult

__all__ = [
    "Evaluation",
    "OptimizeResult",
    "optimize",
    "runner_evaluator",
    "controller_evaluator",
    "lock_maintained",
    "servo_bump_limit",
]

Evaluator = Callable[[Mapping[str, float]], tuple[float, Mapping[str, float]]]


@dataclass
class Evaluation:
    """One evaluated point: ``x`` (requested values), ``objective``, ``constraints`` ({name: c}, feasible
    when c <= 0), ``feasible``, ``error`` (message of a failed evaluation, else None)."""

    x: dict[str, float]
    objective: float
    constraints: dict[str, float]
    feasible: bool
    error: str | None = None


@dataclass
class OptimizeResult:
    """``history`` (every evaluation in order), ``best`` (best feasible, None if none), ``minimize``."""

    history: list[Evaluation] = field(default_factory=list)
    best: Evaluation | None = None
    minimize: bool = True


def _parse_bounds(bounds: Mapping[str, Sequence[Any]]) -> list[tuple[str, float, float, bool]]:
    out = []
    for name, b in bounds.items():
        lo, hi, *scale = b
        log = bool(scale) and scale[0] == "log"
        if scale and scale[0] not in ("log", "linear"):
            raise ValueError(f"bounds[{name!r}]: scale must be 'linear' or 'log', got {scale[0]!r}")
        lo, hi = float(lo), float(hi)
        if not hi > lo or (log and lo <= 0):
            raise ValueError(f"bounds[{name!r}]: need hi > lo (> 0 for log), got {b!r}")
        out.append((name, math.log10(lo) if log else lo, math.log10(hi) if log else hi, log))
    return out


def _to_x(u: Sequence[float], spec: Sequence[tuple[str, float, float, bool]]) -> dict[str, float]:
    x = {}
    for ui, (name, lo, hi, log) in zip(u, spec):
        v = lo + float(ui) * (hi - lo)
        x[name] = 10.0**v if log else v
    return x


def optimize(
    evaluate: Evaluator,
    bounds: Mapping[str, Sequence[Any]],
    *,
    n_init: int = 5,
    n_iter: int = 10,
    minimize: bool = True,
    seed: int = 0,
    num_restarts: int = 4,
    raw_samples: int = 64,
    callback: Callable[[Evaluation], None] | None = None,
) -> OptimizeResult:
    """Constrained BO loop (module docstring). ``bounds``: {name: (lo, hi)} or (lo, hi, "log")."""
    import torch
    from botorch.acquisition.analytic import LogConstrainedExpectedImprovement, LogExpectedImprovement
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.optim import optimize_acqf
    from gpytorch.mlls import ExactMarginalLogLikelihood

    spec = _parse_bounds(bounds)
    if n_init < 2:
        raise ValueError("optimize needs n_init >= 2")
    d = len(spec)
    res = OptimizeResult(minimize=minimize)
    us: list[list[float]] = []
    cnames: list[str] | None = None

    def run_point(u: Sequence[float]) -> None:
        nonlocal cnames
        x = _to_x(u, spec)
        try:
            obj, cons = evaluate(x)
            obj, cons = float(obj), {k: float(v) for k, v in cons.items()}
            err = None if math.isfinite(obj) and all(math.isfinite(v) for v in cons.values()) else "non-finite value"
        except Exception as e:  # noqa: BLE001 - a failed point is data for the loop, kept with its error
            obj, cons, err = math.nan, {}, f"{type(e).__name__}: {e}"
        if cnames is None and err is None:
            cnames = list(cons)
        ev = Evaluation(x, obj, cons, err is None and all(v <= 0 for v in cons.values()), err)
        res.history.append(ev)
        us.append(list(map(float, u)))
        if ev.feasible and (res.best is None or (obj < res.best.objective if minimize else obj > res.best.objective)):
            res.best = ev
        if callback is not None:
            callback(ev)

    sobol = torch.quasirandom.SobolEngine(d, scramble=True, seed=seed)
    for u in sobol.draw(n_init, dtype=torch.double).tolist():
        run_point(u)
    torch.manual_seed(seed)
    unit = torch.stack([torch.zeros(d, dtype=torch.double), torch.ones(d, dtype=torch.double)])
    for _ in range(n_iter):
        names = cnames or []
        ok = [e for e in res.history if e.error is None]
        if len(ok) < 2:  # nothing to model yet: keep exploring
            run_point(sobol.draw(1, dtype=torch.double)[0].tolist())
            continue
        sign = -1.0 if minimize else 1.0
        worst = min(sign * e.objective for e in ok)
        cmax = {c: max([1.0] + [abs(e.constraints[c]) for e in ok if c in e.constraints]) for c in names}
        y = [[sign * e.objective if e.error is None else worst]
             + [e.constraints.get(c, cmax[c]) if e.error is None else cmax[c] for c in names] for e in res.history]
        X = torch.tensor(us, dtype=torch.double)
        Y = torch.tensor(y, dtype=torch.double)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = SingleTaskGP(X, Y, outcome_transform=Standardize(m=Y.shape[-1]))
            fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
            feas = [y_[0] for e, y_ in zip(res.history, y) if e.feasible]
            best_f = max(feas) if feas else float(Y[:, 0].min())
            if names:
                acq = LogConstrainedExpectedImprovement(model, best_f=best_f, objective_index=0,
                                                        constraints={i + 1: (None, 0.0) for i in range(len(names))})
            else:
                acq = LogExpectedImprovement(model, best_f=best_f)
            cand, _ = optimize_acqf(acq, bounds=unit, q=1, num_restarts=num_restarts, raw_samples=raw_samples)
        run_point(cand[0].tolist())
    return res


# ------------------------------------------------------------------------------------------------
# evaluators


def runner_evaluator(
    spec: RunSpec,
    objective: str | Callable[[RunResult], float],
    constraints: Mapping[str, Callable[[RunResult], float]] | None = None,
    *,
    seeds: Sequence[int] | None = None,
    mode: str = "auto",
) -> Evaluator:
    """Evaluator running ``spec`` with ``x`` set at its sweep paths (module docstring).

    ``objective``: a ``RunResult.summary()`` key or a callable on the result. ``constraints``: {name:
    callable(result) -> c}, feasible when c <= 0. ``seeds``: None → the spec's seed.
    """
    from cranberry.runner.run import _as_runspec, run
    from cranberry.runner.sweep import expand

    rs = _as_runspec(spec)
    rs = dataclasses.replace(rs, run=dataclasses.replace(rs.run, seed_mode="common"))
    seeds = [rs.run.seed] if seeds is None else [int(s) for s in seeds]
    obj_fn = (lambda r: r.summary()[objective]) if isinstance(objective, str) else objective
    constraints = dict(constraints or {})

    def evaluate(x: Mapping[str, float]) -> tuple[float, dict[str, float]]:
        point = expand(rs, {k: [float(v)] for k, v in x.items()})
        if len(point) != 1:
            raise ValueError(f"runner_evaluator: the spec has sweep axes of its own ({len(point)} points)")
        objs, cons = [], {c: -math.inf for c in constraints}
        for s in seeds:
            p = point[0].spec
            r = run(dataclasses.replace(p, run=dataclasses.replace(p.run, seed=s)), mode=mode)
            objs.append(float(obj_fn(r)))
            for c, fn in constraints.items():
                cons[c] = max(cons[c], float(fn(r)))
        return float(np.mean(objs)), cons

    return evaluate


def controller_evaluator(
    ctl: ControllerInterface,
    duration_s: float,
    signals: Sequence[str],
    objective: Callable[[Mapping[str, np.ndarray]], float],
    constraints: Mapping[str, Callable[[Mapping[str, np.ndarray]], float]] | None = None,
) -> Evaluator:
    """Evaluator through the abstract controller interface only (module docstring); ``objective`` and
    ``constraints`` act on the ``acquire`` data ({signal: samples, "t_s": time})."""
    constraints = dict(constraints or {})

    def evaluate(x: Mapping[str, float]) -> tuple[float, dict[str, float]]:
        ctl.set_params(dict(x))
        ctl.reset_integrators()
        ctl.engage()
        try:
            data = ctl.acquire(duration_s, list(signals))
        finally:
            ctl.disengage()
        return float(objective(data)), {c: float(fn(data)) for c, fn in constraints.items()}

    return evaluate


# ------------------------------------------------------------------------------------------------
# constraints


def lock_maintained(result: RunResult) -> float:
    """c <= 0 when every lock supervisor saw no lock loss and the run is valid (``RunResult.valid``):
    c = Σ n_loss + (1 if invalid) − 0.5."""
    from cranberry.runner.explore import run_quantities

    q = run_quantities(result)
    n = sum(v for k, v in q.items() if k.startswith("lock.") and k.endswith(".n_loss"))
    return float(n + (0.0 if result.valid else 1.0) - 0.5)


def servo_bump_limit(
    free: RunResult,
    slot: str,
    limit: float,
    *,
    nperseg: int | None = None,
    f_lo_hz: float = 0.0,
    f_hi_hz: float = math.inf,
) -> Callable[[RunResult], float]:
    """Constraint factory: c = max |S|² in [f_lo, f_hi] − ``limit`` with |S|² = PSD of the recorded
    ``slot`` in the candidate run / in ``free`` (the loop-open run, same seed and record)."""
    from cranberry.runner.explore import sensitivity_s2, servo_bump

    x_free = np.asarray(free.series[slot], dtype=float)
    fs = 1.0 / float(np.diff(free.t_s[:2])[0])

    def c(result: RunResult) -> float:
        f, s2 = sensitivity_s2(result.series[slot], x_free, fs, nperseg=nperseg)
        return servo_bump(f, s2, f_lo_hz, f_hi_hz)[0] - float(limit)

    return c
