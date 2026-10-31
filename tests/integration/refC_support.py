"""Shared helpers of the Ref C tests (WP-I3): spec variants of ``refC_nominal`` and small analytic helpers.
Not a test module (loaded by path from the VT / integration files: ``tests`` is not a package under
``--import-mode=importlib``).

Spec variants go through ``config.load`` with an inline overlay (``extends`` the nominal preset), so every
expression (QNL band integrals, loop gains, τ_L − 1/f_sim, fibre coefficients) is re-evaluated.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import numpy as np

from cranberry import config
from cranberry.config.spec import RunSpec

ROOT = Path(__file__).resolve().parents[2]
PRESET = "presets/reference/refC.yaml"
NOMINAL = "presets/reference/refC_nominal.yaml"
TIMESERIES = "presets/axes/refC_dither/timeseries.yaml"

#: every stochastic / drifting source off (the dither and the loop remain)
NOISE_OFF: dict[str, bool] = {
    "laser.fm_white": False,
    "env.thermal.conduit.drift": False, "env.thermal.spool_sq.drift": False, "env.thermal.spool_lo.drift": False,
    "env.mech.conduit.mech": False, "env.mech.spool.mech": False,
    "qnl.radiometer_noise": False,
}

#: T04 §9.1 reading of Ref C: x such that the jitter-free anti-squeezing at 30 MHz is 10.0 dB with η_tot = 0.79
#: (e^{2r} = 12.4 in the single-mode form; −5.63 / −5.18 / −4.32 dB). brentq on the preset's source, 1e−12.
X_T04 = 0.5585017864288694


def load(params: dict[str, Any] | None = None, run: dict[str, Any] | None = None,
         toggles: dict[str, bool] | None = None, measure: dict[str, Any] | None = None,
         overlays: tuple[str, ...] = (), esa: bool = False) -> RunSpec:
    """``refC_nominal`` with ``offline.params`` / ``run`` / ``measure`` overrides; the ESA hooks are dropped
    unless ``esa`` (they need a 0.6 s record)."""
    tree: dict[str, Any] = {"extends": NOMINAL}
    if params:
        tree["offline"] = {"params": dict(params)}
    if run:
        tree["run"] = dict(run)
    m = dict(measure or {})
    if not esa:
        m["esa"] = None
    if m:
        tree["measure"] = m
    spec = config.load(tree, overlays=[ROOT / o for o in overlays])
    if toggles:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **toggles})
    return spec


def at(levels: dict[str, Any], key: str, f_hz: float) -> float:
    """Noise-level ``key`` of a cw_sideband metric at the grid frequency ``f_hz``."""
    f = np.asarray(levels["f_s_hz"], dtype=float)
    i = int(np.argmin(np.abs(f - f_hz)))
    assert abs(f[i] - f_hz) < 1.0, f"{f_hz} Hz is not on the source grid"
    return float(np.asarray(levels[key])[i])


def params_of(spec: RunSpec) -> dict[str, Any]:
    """Resolved ``offline.params`` of a spec."""
    return dict(spec.offline["params"])
