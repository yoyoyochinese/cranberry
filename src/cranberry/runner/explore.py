"""Phase 6 exploration helpers on top of the runner (prompt §8, ``docs/phase6_prompt.md`` steps 2–4, D85-13).

* :func:`combine_effects` — §8 "then combine only effects that mattered": from an ablation
  (``runner.ablate``), the toggles whose removal changes a chosen quantity by more than a threshold;
  then one run with only those effects on (every other candidate toggle off) and a cumulative ladder
  that switches the mattering effects off one by one in order of impact (shows non-additivity).
* :func:`run_quantities` — the comparison quantities of one run as a flat {name: value}: metric (i) and
  (ii) moments and R′∓ (``RunResult.summary()``), η_k (``eff.port<k>.*``), stray light per port
  (``stray.port<k>.clicks_per_gate`` / ``.photons_per_mode``), per lock supervisor the slips, losses,
  relocks and lock duration (``lock.<sup>.n_slip`` ..., ``lock.<sup>.lock_duration_s`` = time of the
  first loss, or the run duration when the lock held), and per actuator the range use
  (``range.<actuator>`` = max |x − centre| / half range over the recorded series, 1 = at the limit).
* :func:`compare_seeds` — several configurations × the same seeds (common random numbers): mean, sample
  standard deviation, min / max of every quantity, optionally with the four metric-(ii) components per
  seed (``runner.decompose``, G13: ``ii.port<k>.<component>_rad2``). :func:`differs` applies the
  prompt's rule "claim only differences larger than the spread".
* :func:`sensitivity_s2` / :func:`servo_bump` — |S(f)|² as the ratio of the in-loop to the free-running
  error PSD (same seed, loop open vs closed; the Ref D comparison's estimator, D81-1) and its peak (the
  servo bump) in a band.

Actuator range use needs the actuator's input (PZT / EOM: ``params.in`` against ``drive_min_v`` /
``drive_max_v`` or ``input_min_v`` / ``input_max_v``) or temperature (thermal: ``act.<name>.dT`` against
``dt_min_k`` / ``dt_max_k``) to be recorded; :func:`with_recorded` adds those slots to
``measure.record`` (``compare_seeds`` does it). Without a record the final state is used.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = [
    "RANGE_KEYS",
    "combine_effects",
    "actuator_ranges",
    "with_recorded",
    "run_quantities",
    "compare_seeds",
    "differs",
    "sensitivity_s2",
    "servo_bump",
]

#: (min key, max key, slot template) of an actuator's range; ``{in}`` = ``params["in"]``, ``{name}`` = block name.
RANGE_KEYS: tuple[tuple[str, str, str], ...] = (
    ("drive_min_v", "drive_max_v", "{in}"),
    ("input_min_v", "input_max_v", "{in}"),
    ("dt_min_k", "dt_max_k", "act.{name}.dT"),
)


def combine_effects(
    spec: RunSpec,
    quantity: str,
    *,
    toggles: Sequence[str] | None = None,
    threshold: float = 0.05,
    relative: bool = True,
    ablation: Mapping[str, RunResult] | None = None,
    ladder: bool = True,
    mode: str = "auto",
) -> dict[str, Any]:
    """Ablate, keep the effects that mattered for ``quantity`` (a ``RunResult.summary()`` key), combine them.

    An effect (toggle) matters when switching it off changes ``quantity`` by more than ``threshold``
    (times |baseline| when ``relative``). ``ablation``: an earlier ``ablate`` result to reuse (its
    toggles are the candidates); else ``ablate(spec, toggles or list_toggles(spec))``.

    Returns {"baseline": value, "impact": {toggle: delta}, "mattered": [toggles by |delta|, largest
    first], "combined": value with only the mattering candidates on (every other candidate off),
    "residual": combined − baseline, "ladder": [(toggle switched off cumulatively, value), ...],
    "results": {"ablation": ..., "combined": RunResult, "ladder": [RunResult, ...]}}.
    """
    from cranberry.runner.ablate import ablate, list_toggles
    from cranberry.runner.run import _as_runspec, run

    rs = _as_runspec(spec)
    abl = dict(ablation) if ablation is not None else ablate(rs, list(toggles) if toggles is not None else list_toggles(rs), mode=mode)
    base = abl["baseline"].summary()[quantity]
    impact = {t: r.summary()[quantity] - base for t, r in abl.items() if t != "baseline"}
    scale = abs(base) if relative else 1.0
    mattered = sorted((t for t, d in impact.items() if abs(d) > threshold * scale), key=lambda t: -abs(impact[t]))
    rest = [t for t in impact if t not in mattered]

    def off(keys: Iterable[str]) -> RunSpec:
        return dataclasses.replace(rs, toggles={**dict(rs.toggles), **{k: False for k in keys}})

    combined = run(off(rest), mode=mode)
    cval = combined.summary()[quantity]
    steps, step_results = [], []
    if ladder:
        for i, t in enumerate(mattered):
            r = run(off(mattered[: i + 1]), mode=mode)
            step_results.append(r)
            steps.append((t, r.summary()[quantity]))
    return {"baseline": base, "impact": impact, "mattered": mattered, "combined": cval, "residual": cval - base,
            "ladder": steps, "results": {"ablation": abl, "combined": combined, "ladder": step_results}}


def actuator_ranges(spec: RunSpec) -> dict[str, tuple[str, float, float]]:
    """{actuator block name: (slot, lo, hi)} for every pipeline block with a range pair of :data:`RANGE_KEYS`."""
    from cranberry.config.spec import iter_block_specs

    out = {}
    for _, b in iter_block_specs(spec.pipeline):
        for lo_k, hi_k, tpl in RANGE_KEYS:
            lo, hi = b.params.get(lo_k), b.params.get(hi_k)
            if isinstance(lo, (int, float)) and isinstance(hi, (int, float)) and hi > lo:
                slot = tpl.format(name=b.name, **{"in": b.params.get("in", "")})
                if slot:
                    out[b.name] = (slot, float(lo), float(hi))
                    break
    return out


def with_recorded(spec: RunSpec, slots: Iterable[str]) -> RunSpec:
    """``spec`` with ``slots`` added to ``measure.record`` (a full-rate record of them when there is none)."""
    slots = list(slots)
    if not slots:
        return spec
    measure = dict(spec.measure or {})
    rec = dict(measure.get("record") or {"slots": [], "decim": 1})
    rec["slots"] = list(rec["slots"]) + [s for s in slots if s not in rec["slots"]]
    measure["record"] = rec
    return dataclasses.replace(spec, measure=measure)


def _supervisors(spec: RunSpec) -> list[str]:
    from cranberry.config.spec import iter_block_specs

    return [b.name for _, b in iter_block_specs(spec.pipeline) if b.kind == "lock_supervisor"]


def run_quantities(result: RunResult) -> dict[str, float]:
    """Comparison quantities of one run (module docstring)."""
    out = dict(result.summary())
    out["valid"] = float(result.valid)
    m = result.metrics
    if m is not None:
        for k, s in m.stray_light.items():
            out[f"stray.port{k}.clicks_per_gate"] = float(s.clicks_per_gate)
            out[f"stray.port{k}.photons_per_mode"] = float(s.photons_per_mode)
    st = result.final_state
    rs = result.spec

    def last(slot: str) -> float | None:
        if slot in result.series and np.size(result.series[slot]):
            return float(np.asarray(result.series[slot])[-1])
        if st is not None and slot in st.layout:
            return float(np.asarray(st[slot]).reshape(-1)[0])
        return None

    for sup in _supervisors(rs):
        for key in ("n_slip", "n_loss", "n_relock", "t_loss"):
            v = last(f"lock.{sup}.{key}")
            if v is not None:
                out[f"lock.{sup}.{key}"] = v
        t_loss = out.get(f"lock.{sup}.t_loss")
        if t_loss is not None:
            out[f"lock.{sup}.lock_duration_s"] = t_loss if t_loss >= 0 else float(rs.run.duration_s)
    for name, (slot, lo, hi) in actuator_ranges(rs).items():
        if slot in result.series and np.size(result.series[slot]):
            x = np.asarray(result.series[slot], dtype=float)
        elif st is not None and slot in st.layout:
            x = np.asarray(st[slot], dtype=float)
        else:
            continue
        out[f"range.{name}"] = float(np.max(np.abs(x - 0.5 * (lo + hi))) / (0.5 * (hi - lo)))
    return out


def _decomposed(spec: RunSpec, mode: str) -> dict[str, float]:
    from cranberry.runner.decompose import decompose_metric_ii

    out = {}
    for k, c in decompose_metric_ii(spec, mode=mode).items():
        for name, v in c.variances_rad2.items():
            out[f"ii.port{k}.{name}_rad2"] = float(v)
        out[f"ii.port{k}.total_rad2"] = float(c.total_var_rad2)
        out[f"ii.port{k}.linearity_ratio"] = float(c.linearity_ratio)
        out[f"ii.port{k}.decomposition_ok"] = float(c.status == "ok")
    return out


def compare_seeds(
    configs: Mapping[str, RunSpec],
    seeds: Sequence[int],
    *,
    decompose: bool = False,
    quantities: Sequence[str] | None = None,
    mode: str = "auto",
) -> dict[str, dict[str, dict[str, Any]]]:
    """{label: {quantity: {"mean", "std", "min", "max", "n", "values"}}} over ``seeds`` for each config.

    Every config runs with the same seeds (common random numbers). ``decompose``: also the metric-(ii)
    components per seed (1 + number of non-empty groups extra runs per seed; needs
    ``measure.decompose``; the decomposition sets the group toggles itself, so a config that differs
    only by a group toggle override has the same components). ``quantities``: keep only these keys (default: every key present in all
    seeds of that config). ``std`` is the sample standard deviation (ddof = 1; 0 for one seed).
    """
    from cranberry.runner.run import _as_runspec, run

    if not seeds:
        raise ValueError("compare_seeds needs at least one seed")
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for label, spec in configs.items():
        rs = _as_runspec(spec)
        rs = with_recorded(rs, [slot for slot, _, _ in actuator_ranges(rs).values()])
        rows = []
        for s in seeds:
            rss = dataclasses.replace(rs, run=dataclasses.replace(rs.run, seed=int(s)))
            row = run_quantities(run(rss, mode=mode))
            if decompose:
                row.update(_decomposed(rss, mode))
            rows.append(row)
        keys = [q for q in rows[0] if all(q in r for r in rows)] if quantities is None else list(quantities)
        stats = {}
        for q in keys:
            v = np.array([r.get(q, math.nan) for r in rows], dtype=float)
            stats[q] = {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)) if v.size > 1 else 0.0,
                        "min": float(np.min(v)), "max": float(np.max(v)), "n": int(v.size), "values": v.tolist()}
        out[label] = stats
    return out


def differs(a: Mapping[str, float], b: Mapping[str, float], k: float = 1.0) -> bool:
    """True when two :func:`compare_seeds` entries differ by more than their spread:
    |mean_a − mean_b| > k·(std_a + std_b) (prompt §8 / Phase 6 step 4)."""
    return abs(a["mean"] - b["mean"]) > k * (a["std"] + b["std"])


def sensitivity_s2(in_loop: np.ndarray, free: np.ndarray, fs_hz: float, *, nperseg: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(f [Hz], |S(f)|²) = PSD(in-loop error) / PSD(free-running error) (Welch, ``core.fft.psd_onesided``).
    Both series at ``fs_hz``, same disturbance realisation (same seed); the DC bin is dropped."""
    from cranberry.core.fft import psd_onesided

    f, s_in = psd_onesided(np.asarray(in_loop, dtype=float), fs_hz, nperseg=nperseg)
    _, s_free = psd_onesided(np.asarray(free, dtype=float), fs_hz, nperseg=nperseg)
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = s_in / s_free
    return f[1:], s2[1:]


def servo_bump(f_hz: np.ndarray, s2: np.ndarray, f_lo_hz: float = 0.0, f_hi_hz: float = math.inf) -> tuple[float, float]:
    """(peak |S|², its frequency [Hz]) of :func:`sensitivity_s2` in [f_lo, f_hi] (finite values only)."""
    f, s = np.asarray(f_hz), np.asarray(s2)
    sel = (f >= f_lo_hz) & (f <= f_hi_hz) & np.isfinite(s)
    if not np.any(sel):
        raise ValueError(f"servo_bump: no finite |S|² in [{f_lo_hz:g}, {f_hi_hz:g}] Hz")
    i = int(np.argmax(np.where(sel, s, -np.inf)))
    return float(s[i]), float(f[i])
