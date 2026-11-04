"""Standard plots of the Phase 6 report (D84-3, D85-13). matplotlib only (``plot`` extra).

Each helper draws into ``ax`` (a new figure when None) and returns the Axes, so a report can lay out
several panels and save them itself (``ax.figure.savefig(...)``). Inputs are plain arrays / mappings,
e.g. from ``runner.explore`` (``compare_seeds``, ``sensitivity_s2``) and ``RunResult.series``:

* :func:`plot_psd` — one-sided PSD(s) on log–log axes with shaded frequency bands (e.g. the loop
  bandwidth, the metric band up to f_max, the ledger's out-of-band region).
* :func:`plot_metric_ii_components` — the four metric-(ii) components (G13) as grouped bars in mrad rms,
  one group per configuration, optional seed spread as error bars.
* :func:`plot_vs_length` — a metric against fiber length, mean ± spread over seeds per configuration.
* :func:`plot_sensitivity` — |S(f)|² (in-loop / free-running PSD ratio) with the 0 dB line and an
  optional servo-bump limit.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = ["plot_psd", "plot_metric_ii_components", "plot_vs_length", "plot_sensitivity", "COMPONENT_LABELS"]

#: Short labels of ``METRIC_II_COMPONENTS`` (G13).
COMPONENT_LABELS: dict[str, str] = {
    "locked_path_residual": "locked path",
    "common_points": "common points",
    "laser_decorrelation_oob": "laser decorrelation",
    "squeezed_only": "squeezed only",
}


def _ax(ax: Any) -> Any:
    if ax is not None:
        return ax
    import matplotlib.pyplot as plt

    return plt.subplots(figsize=(6.4, 4.0), layout="constrained")[1]


def plot_psd(
    curves: Mapping[str, tuple[np.ndarray, np.ndarray]],
    *,
    bands: Sequence[tuple[float, float, str]] = (),
    ylabel: str = "PSD [rad²/Hz]",
    title: str | None = None,
    ax: Any = None,
) -> Any:
    """``curves``: {label: (f [Hz], S)}; ``bands``: [(f_lo, f_hi, label)] shaded behind the curves."""
    ax = _ax(ax)
    for i, (lo, hi, label) in enumerate(bands):
        ax.axvspan(lo, hi, color=f"C{(i + 5) % 10}", alpha=0.12, lw=0, label=label)
    for label, (f, s) in curves.items():
        f, s = np.asarray(f), np.asarray(s)
        sel = f > 0
        ax.loglog(f[sel], s[sel], lw=1.0, label=label)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    if curves or bands:
        ax.legend(fontsize=8)
    if title:
        ax.set_title(title)
    return ax


def plot_metric_ii_components(
    variances_rad2: Mapping[str, Mapping[str, float]],
    *,
    spread_rad2: Mapping[str, Mapping[str, float]] | None = None,
    title: str | None = None,
    ax: Any = None,
) -> Any:
    """Grouped bars of √variance [mrad] per component; ``variances_rad2``: {config: {component: var}}
    (``MetricIIComponents.variances_rad2``, or seed means from ``compare_seeds``). ``spread_rad2``: the
    same shape with the seed standard deviation of the variance (error bar = its effect on √var)."""
    ax = _ax(ax)
    comps = list(COMPONENT_LABELS)
    comps += [c for v in variances_rad2.values() for c in v if c not in comps]
    n = max(len(variances_rad2), 1)
    w = 0.8 / n
    x = np.arange(len(comps))
    for i, (cfg, var) in enumerate(variances_rad2.items()):
        v = np.array([max(var.get(c, 0.0), 0.0) for c in comps])
        rms = np.sqrt(v) * 1e3
        err = None
        if spread_rad2 and cfg in spread_rad2:
            sd = np.array([spread_rad2[cfg].get(c, 0.0) for c in comps])
            lo = np.sqrt(np.clip(v - sd, 0.0, None)) * 1e3
            hi = np.sqrt(v + sd) * 1e3
            err = np.vstack([rms - lo, hi - rms])
        ax.bar(x + (i - (n - 1) / 2) * w, rms, w, yerr=err, capsize=2, label=cfg)
    ax.set_xticks(x, [COMPONENT_LABELS.get(c, c) for c in comps], rotation=15)
    ax.set_ylabel("metric (ii) component [mrad rms]")
    ax.grid(True, axis="y", alpha=0.25)
    if variances_rad2:
        ax.legend(fontsize=8)
    if title:
        ax.set_title(title)
    return ax


def plot_vs_length(
    lengths_m: Sequence[float],
    series: Mapping[str, Any],
    *,
    ylabel: str,
    logx: bool = True,
    logy: bool = False,
    title: str | None = None,
    ax: Any = None,
) -> Any:
    """Metric vs fiber length with seed spread. ``series``: {label: values} where values is an array of
    shape (n_lengths, n_seeds) (mean and sample std are taken here) or a pair (mean, std) of arrays of
    length n_lengths. Drawn as mean line with markers and a ±std band."""
    ax = _ax(ax)
    L = np.asarray(lengths_m, dtype=float)
    for label, v in series.items():
        if isinstance(v, tuple) and len(v) == 2:
            mean, std = (np.asarray(a, dtype=float) for a in v)
        else:
            a = np.asarray(v, dtype=float)
            if a.ndim != 2 or a.shape[0] != L.size:
                raise ValueError(f"plot_vs_length: {label!r} must be (n_lengths, n_seeds) or (mean, std)")
            mean = a.mean(axis=1)
            std = a.std(axis=1, ddof=1) if a.shape[1] > 1 else np.zeros(L.size)
        (line,) = ax.plot(L, mean, "o-", ms=3, lw=1.2, label=label)
        ax.fill_between(L, mean - std, mean + std, color=line.get_color(), alpha=0.2, lw=0)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("fiber length [m]")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    if series:
        ax.legend(fontsize=8)
    if title:
        ax.set_title(title)
    return ax


def plot_sensitivity(
    curves: Mapping[str, tuple[np.ndarray, np.ndarray]],
    *,
    limit: float | None = None,
    title: str | None = None,
    ax: Any = None,
) -> Any:
    """|S(f)|² curves {label: (f [Hz], |S|²)} on log–log axes, the |S|² = 1 line (in-loop = free) and an
    optional servo-bump ``limit`` (linear |S|²)."""
    ax = _ax(ax)
    for label, (f, s2) in curves.items():
        f, s2 = np.asarray(f), np.asarray(s2)
        sel = (f > 0) & np.isfinite(s2) & (s2 > 0)
        ax.loglog(f[sel], s2[sel], lw=1.0, label=label)
    ax.axhline(1.0, color="0.3", lw=0.8, ls="--", label="|S|² = 1 (in-loop = free)")
    if limit is not None:
        ax.axhline(limit, color="C3", lw=0.8, ls=":", label=f"limit {limit:g} ({10 * math.log10(limit):.1f} dB)")
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("|S|² = PSD in-loop / PSD free")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    if title:
        ax.set_title(title)
    return ax
