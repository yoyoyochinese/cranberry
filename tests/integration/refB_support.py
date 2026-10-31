"""Shared helpers of the Ref B tests (WP-I4): spec variants of ``refB_nominal`` and the Ref B Eq. (1) map.
Not a test module (loaded by path from the VT / integration files: ``tests`` is not a package under
``--import-mode=importlib``).

Spec variants go through ``config.load`` with an inline overlay (``extends`` the nominal preset), so every
expression (x, κ, L(x), A±, error slopes, loop gains, τ_L − 1/f_sim, OPO coefficients) is re-evaluated.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from cranberry import config
from cranberry.config.spec import RunSpec

ROOT = Path(__file__).resolve().parents[2]
PRESET = "presets/reference/refB.yaml"
NOMINAL = "presets/reference/refB_nominal.yaml"
SINGLE_LAYER = "presets/axes/refB_theta/single_layer_pzt.yaml"   # θ̃ = 3.9° (Fig. 3 dotted)
MEASURED_LOSS = "presets/axes/refB_loss/measured_100mw.yaml"     # L = 0.0038 [Ref B §3]
TIMESERIES = "presets/axes/refB_dither/timeseries.yaml"
MCC = "presets/axes/refB_mcc/lowpass.yaml"

#: index of f = 1 MHz (the ESA frequency) on the source grid [2.95e5, 1e6]
I_1MHZ = 1

#: every stochastic source and quantiser off (loops, dither and OPO length lock remain)
NOISE_OFF: dict[str, bool] = {
    "tisa.fm_white": False, "tisa.peaks": False, "env.mech.mech": False, "env.shared.mech": False,
    "opo_dn.length_noise": False, "pd1.shot_noise": False, "pd1.electronic_noise": False,
    "bhdbeat.shot_noise": False, "bhdbeat.electronic_noise": False,
    "adc1.quantization": False, "adc2.quantization": False, "dac1.quantization": False, "dac2.quantization": False,
}


def load(params: dict[str, Any] | None = None, run: dict[str, Any] | None = None,
         toggles: dict[str, bool] | None = None, measure: dict[str, Any] | None = None,
         pipeline: list[dict[str, Any]] | None = None, overlays: tuple[str, ...] = (), esa: bool = False) -> RunSpec:
    """``refB_nominal`` with ``offline.params`` / ``run`` / ``measure`` / pipeline-param overrides; the ESA hooks are
    dropped unless ``esa`` (they need ≥ 0.1 s after the settle)."""
    tree: dict[str, Any] = {"extends": NOMINAL}
    if params:
        tree["offline"] = {"params": dict(params)}
    if run:
        tree["run"] = dict(run)
    if pipeline:
        tree["pipeline"] = list(pipeline)
    if toggles:
        tree["toggles"] = dict(toggles)
    m = dict(measure or {})
    if not esa:
        m["esa"] = None
    if m:
        tree["measure"] = m
    return config.load(tree, overlays=[ROOT / o for o in overlays])


def params_of(spec: RunSpec) -> dict[str, Any]:
    """Resolved ``offline.params`` of a spec."""
    return dict(spec.offline["params"])


def offline_block_params(spec: RunSpec, name: str) -> dict[str, Any]:
    """Resolved params of the offline block ``name``."""
    return dict(next(b for b in spec.offline["blocks"] if b["name"] == name)["params"])


def db(r: Any) -> Any:
    return 10.0 * np.log10(r)


def eq1_db(v_minus_db: float, v_plus_db: float, theta_rad: float) -> tuple[float, float]:
    """Ref B Eq. (1): R′∓ = R∓ cos²θ̃ + R± sin²θ̃ [dB] from the jitter-free R∓ [dB]."""
    rm, rp = 10 ** (v_minus_db / 10), 10 ** (v_plus_db / 10)
    c2, s2 = np.cos(theta_rad) ** 2, np.sin(theta_rad) ** 2
    return float(db(rm * c2 + rp * s2)), float(db(rp * c2 + rm * s2))
