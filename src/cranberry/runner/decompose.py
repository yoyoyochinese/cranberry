"""``decompose_metric_ii(spec)``: four-component split of metric (ii) by reruns (G13). Owner: WP-F4.

Runs all-on plus one rerun per component noise-source set with the same RNG streams (G11), then checks
linearity: the component variances must sum to the all-on variance within 5 % (T04 V12). Saturation or
lock loss marks the decomposition invalid.

Procedure (WP-F4; T04 §4)
    * ``groups``: {component (``METRIC_II_COMPONENTS``, all four keys required): [toggle keys
      ``"<block>.<toggle>"`` of the noise sources of that component]}. A toggle belongs to one group;
      an empty list means the component has no source in this configuration (variance 0, no run).
    * all-on run: every group toggle True; component run c: the toggles of c True, of every other
      group False (all other toggles as in the spec — a noise source left out of the groups is then
      counted in every component run, which the linearity check flags when more than one group is
      non-empty; list every noise source). Same root seed → identical realisations of the
      sources that stay on (streams keyed by block name, G11), so in a linear loop
      Θ_all(t) = Σ_c Θ_c(t) sample by sample and Var(Θ_all) = Σ_c Var(Θ_c) up to the sample
      cross-covariances of independent sources.
    * variances are Var(Θ_k) (mean removed) from the metric probes' full-rate accumulators **plus the
      ledger part** (Phase 5 M6-06): the out-of-band variance the ``BandLedger`` adds to metric
      ``ii.port<k>`` for the run's beams (above the block Nyquist frequencies, and any explicit
      declaration such as a window row), evaluated in the same run. Every noise source declares its band
      only while its toggle is on, so the component run c carries exactly the ledger rows of its group,
      and the all-on total carries all of them (the f_N split of conv §6 is invisible to the sum).
      Toggle keys that are not noise sources but change the loop's own error (ADC / DAC quantisation:
      the DAC limit cycle, D65-3) belong to a group like any other key — the future presets put them in
      ``locked_path_residual`` (loop electronics); left out they would be counted in every run.
    * ``validity``: [{slot, abs_max}] — a run is invalid when |slot| ≥ abs_max at any recorded sample
      (saturation flags / counters, actuator outputs at their limit, lock-loss counters written by the
      CTRL / ACT / lock_supervisor blocks). The slots are added to ``measure.record`` (decim 1 when
      the spec records nothing). Invalid in any run (or a ``loop_unstable`` CHK event, which makes
      ``RunResult.valid`` False, M6-01), or linearity off by more than ``tolerance``
      → ``status = "decomposition_invalid"`` with the reason.

Config (``RunSpec.measure.decompose``, or the keyword arguments)::

    decompose:
      groups: {locked_path_residual: [...], common_points: [...], laser_decorrelation_oob: [...],
               squeezed_only: [...]}
      validity: [{slot: ctrl.lock1.sat, abs_max: 0.5}, ...]    # [] for none
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.blocks.meas.metrics import METRIC_II_COMPONENTS, MetricIIComponents
from cranberry.runner.run import INVALIDATING_CHK

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["decompose_metric_ii", "validity_violations"]


def validity_violations(result: RunResult, validity: Sequence[Mapping[str, Any]]) -> list[str]:
    """Messages for every ``{slot, abs_max}`` reached (recorded series, else the final state)."""
    out = []
    for v in validity:
        slot, lim = v["slot"], float(v["abs_max"])
        if slot in result.series:
            x = np.asarray(result.series[slot], dtype=float)
        elif result.final_state is not None and slot in result.final_state.layout:
            x = np.asarray(result.final_state[slot], dtype=float)
        else:
            raise KeyError(f"validity slot {slot!r} is neither recorded nor in the final state")
        if np.any(np.abs(x) >= lim):
            out.append(f"{slot} reached |x| >= {lim:g} (max {float(np.max(np.abs(x))):g})")
    return out


def _var_with_ledger(result: RunResult, k: int) -> float:
    """Var(Θ_k) of the time loop + the ledger's out-of-band variance of metric ii.port<k> (module docstring)."""
    m = result.metrics.metric_ii[k]
    name = f"ii.port{k}"
    led = result.offline.band_ledger if result.offline is not None else None
    hf = led.hf_variance(name, m.beams) if led is not None and name in result.offline.metric_names else 0.0
    return float(m.stats.var) + float(hf)


def decompose_metric_ii(
    spec: RunSpec,
    *,
    groups: Mapping[str, Sequence[str]] | None = None,
    validity: Sequence[Mapping[str, Any]] | None = None,
    tolerance: float = 0.05,
    mode: str = "auto",
) -> dict[int, MetricIIComponents]:
    """{port k: components of Θ_k}; ports without metric (ii) are omitted (module docstring).

    ``groups`` / ``validity`` default to ``spec.measure["decompose"]``; ``tolerance`` is the T04 V12
    linearity bound (5 %).
    """
    from cranberry.runner.run import _as_runspec, run

    rs = _as_runspec(spec)
    cfg = (rs.measure or {}).get("decompose") or {}
    groups = cfg.get("groups") if groups is None else groups
    validity = cfg.get("validity") if validity is None else validity
    if groups is None or validity is None:
        raise KeyError("decompose_metric_ii needs measure.decompose.groups and .validity (D34: no defaults)")
    if set(groups) != set(METRIC_II_COMPONENTS):
        raise ValueError(f"decompose groups must have exactly the keys {METRIC_II_COMPONENTS}, got {sorted(groups)}")
    owner: dict[str, str] = {}
    for c, ts in groups.items():
        for t in ts:
            if t in owner:
                raise ValueError(f"toggle {t!r} is in groups {owner[t]!r} and {c!r}")
            owner[t] = c
    for v in validity:
        if set(v) != {"slot", "abs_max"}:
            raise ValueError(f"validity entries are {{slot, abs_max}}, got {dict(v)!r}")

    measure = dict(rs.measure or {})
    vslots = [v["slot"] for v in validity]
    if vslots:
        rec = dict(measure.get("record") or {"slots": [], "decim": 1})
        rec["slots"] = list(rec["slots"]) + [s for s in vslots if s not in rec["slots"]]
        measure["record"] = rec

    def variant(on: str | None) -> RunSpec:
        tg = {**dict(rs.toggles), **{t: (on is None or owner[t] == on) for t in owner}}
        return dataclasses.replace(rs, toggles=tg, measure=measure)

    runs: dict[str | None, RunResult] = {None: run(variant(None), mode=mode)}
    all_on = runs[None]
    if all_on.metrics is None or not all_on.metrics.metric_ii:
        return {}
    for c in METRIC_II_COMPONENTS:
        if groups[c]:
            runs[c] = run(variant(c), mode=mode)
    reasons = []
    for key, r in runs.items():
        bad = validity_violations(r, validity)
        if bad:
            reasons.append(f"{'all-on' if key is None else key} run: saturation / lock loss: " + "; ".join(bad))
        unstable = sorted({e["source"] for e in r.events if e.get("code") in INVALIDATING_CHK})
        if unstable:
            reasons.append(f"{'all-on' if key is None else key} run: linearly unstable loop ({', '.join(unstable)})")
    out: dict[int, MetricIIComponents] = {}
    for k, m in all_on.metrics.metric_ii.items():
        total = _var_with_ledger(all_on, k)
        var = {c: (_var_with_ledger(runs[c], k) if c in runs else 0.0) for c in METRIC_II_COMPONENTS}
        s = sum(var.values())
        ratio = s / total if total > 0 else (1.0 if s == 0 else float("inf"))
        why = list(reasons)
        if not abs(ratio - 1.0) <= tolerance:
            why.append(f"linearity: Σ component variances / all-on = {ratio:.4f} (|·−1| > {tolerance:g})")
        out[k] = MetricIIComponents(var, total, ratio, "decomposition_invalid" if why else "ok",
                                    "; ".join(why) if why else None)
    return out
