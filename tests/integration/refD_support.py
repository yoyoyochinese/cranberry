"""Shared helpers of the Ref D tests (WP-I2): spec variants of ``refD_nominal`` and the exact
multirate z-domain loop of the preset. Not a test module (loaded by path from the VT / integration
files: ``tests`` is not a package under ``--import-mode=importlib``).

Spec variants go through ``config.load`` with an inline overlay (``extends`` the nominal preset), so
every expression (τ_L − 1/f_sim − explicit delays, γ_coh, arm lengths, ...) is re-evaluated.

Loop model (``loop_response``): the DSC1 PID / DAC run every 2 base samples, everything else every
sample. With d the phase disturbance on φ_D, A(f) the actuator path at the base rate (hold 1 + z⁻¹,
output LPF, loop-closure z⁻¹, delay buffer, 914, 915B, segment-transit interpolation, sign of φ_D),
B(f) = ½·K_d·H_in(f)·(1 + z⁻¹)/2 the sensor path including the decimation factor ½, C the PID
(z₂ = z²) and g the DAC load divider, the component of φ_D at f is (f' = f + f_sim/2: the image that
the decimation folds back onto f)

    φ(f)/d(f) = [1 − g C B(f') A(f')] / [1 − g C (B(f) A(f) + B(f') A(f'))].
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import signal

from cranberry import config
from cranberry.blocks.ctrl.filters import butter_sos
from cranberry.config.spec import RunSpec, spec_from_tree, spec_to_tree
from cranberry.core.delay import DelayDesign
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import build, iter_leaves
from cranberry.core.state import Params

ROOT = Path(__file__).resolve().parents[2]
NOMINAL = "presets/reference/refD_nominal.yaml"
C_M_PER_S = 299_792_458.0

#: every stochastic / drifting source and both quantisers off (deterministic linear loop)
NOISE_OFF: dict[str, bool] = {
    "dx1.fm_white": False, "dx1.fm_flicker": False, "dx1.rin": False,
    "env.thermal.drift": False, "env.mech.mech": False,
    "pol.jones.pol_smf_drift": False,
    "pdb.offset.offset": False, "pdb.offset.drift": False,
    "pdb450c.shot_noise": False, "pdb450c.electronic_noise": False,
    "dsc1.adc.noise": False, "dsc1.adc.quantization": False, "dsc1.dac.quantization": False,
}


def load(params: dict[str, Any] | None = None, run: dict[str, Any] | None = None,
         blocks: list[dict[str, Any]] | None = None, measure: dict[str, Any] | None = None,
         toggles: dict[str, bool] | None = None) -> RunSpec:
    """``refD_nominal`` with ``offline.params`` / ``run`` / named-block overrides (keyed merge)."""
    tree: dict[str, Any] = {"extends": NOMINAL}
    if params:
        tree["offline"] = {"params": dict(params)}
    if run:
        tree["run"] = dict(run)
    if blocks:
        tree["pipeline"] = list(blocks)
    if measure is not None:
        tree["measure"] = measure
    spec = config.load(tree)
    if toggles:
        spec = dataclasses.replace(spec, toggles={**spec.toggles, **toggles})
    return spec


def replace_block(spec: RunSpec, name: str, new: list[dict[str, Any]]) -> RunSpec:
    """``spec`` with the top-level pipeline block ``name`` replaced by the block dicts ``new`` (same place)."""
    tree = spec_to_tree(spec)
    pipe = tree["pipeline"]
    i = next(k for k, b in enumerate(pipe) if b.get("name") == name)
    tree["pipeline"] = pipe[:i] + list(new) + pipe[i + 1:]
    tree["toggles"] = {k: v for k, v in tree["toggles"].items() if k.rpartition(".")[0] != name}
    return spec_from_tree(tree)


def block_params(spec: RunSpec, name: str) -> dict[str, Any]:
    return dict(next(b for b in spec.pipeline.children if b.name == name).params)


def segment(spec: RunSpec, sid: str) -> dict[str, Any]:
    return next(s for s in spec.topology["segments"] if s["id"] == sid)


def demod(x: np.ndarray, n: np.ndarray, f_hz: float, fs_hz: float) -> complex:
    """Σ x e^{−i2πfn/fs} (ratio of two of these = transfer function, Laplace sign s = +iω)."""
    return complex(np.sum(np.asarray(x, float) * np.exp(-2j * np.pi * f_hz * np.asarray(n) / fs_hz)))


def fringe(spec: RunSpec) -> dict[str, float]:
    """V_0, 𝒱 and P_det of the dc_fringe block (T15 eq. 3; δp = 0)."""
    p = block_params(spec, "pdb450c")
    p1 = p["kappa1"] * p["t_arm1"] * p["p_in_w"]
    p2 = (1 - p["kappa1"]) * p["t_arm2"] * p["p_in_w"]
    pdet = p1 + p2
    return {"v0": p["gain_v_per_a"] * p["responsivity_a_per_w"] * pdet, "pdet": pdet,
            "vis": 4 * math.sqrt(p["kappa2"] * (1 - p["kappa2"]) * p1 * p2) / pdet, "sign": p["coupler_sign"]}


def built_leaves(spec: RunSpec) -> tuple[dict[str, Any], Params]:
    """Leaf blocks of ``spec`` after build (dynamics discretised) and the runtime params."""
    from cranberry.config.spec import apply_toggles
    from cranberry.core.pipeline import from_spec
    from cranberry.runner.run import inject_topology

    pipe = apply_toggles(inject_topology(spec.pipeline, spec.topology), spec.toggles)
    lay = SlotLayout(spec.run.f_sim_hz)
    blk = build(from_spec(pipe), lay, spec.run.seed)
    return {b.name: b for b in iter_leaves(blk)}, Params(lay)


def _sos_resp(sos: np.ndarray, f: np.ndarray, fs: float) -> np.ndarray:
    return signal.sosfreqz(sos, worN=f, fs=fs)[1]


def loop_parts(spec: RunSpec, f_hz: np.ndarray, *, k_d: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """(A(f), B(f), C(f), g) of the module docstring; ``k_d`` = ∂v/∂φ_D at the lock point [V/rad]."""
    fs = spec.run.f_sim_hz
    h = 1.0 / fs
    f = np.atleast_1d(np.asarray(f_hz, float))
    zi = np.exp(-2j * np.pi * f * h)
    leaves, params = built_leaves(spec)
    lp = leaves["dsc1.latency"].params
    dd = DelayDesign(float(lp["tau_s"]) * fs, "linear")
    h_delay = zi ** dd.k * ((1 - dd.mu) + dd.mu * zi)
    amp = leaves["amp914"].freq_response(f)
    st = leaves["stretcher915b"]
    k_act = st.gain_rad_per_v(params) * float(params[st.param_name("polarity")])
    h_str = st.freq_response(f)
    a1 = segment(spec, "a1")
    from cranberry.path.segments import frac_delay

    k1, mu1 = frac_delay(0.5 * a1["n_g"] * a1["L_m"] / C_M_PER_S, h)
    h_tr = zi ** k1 * ((1 - mu1) + mu1 * zi)
    tg = {**{k: v for k, v in block_params(spec, "dsc1.out_bw")["toggles"].items()}}
    out_bw = block_params(spec, "dsc1.out_bw")
    in_bw = block_params(spec, "dsc1.in_bw")
    h_out = _sos_resp(butter_sos("lowpass", out_bw["order"], out_bw["f_c_hz"], fs, out_bw["method"]), f, fs) \
        if tg["enabled"] else np.ones_like(zi)
    h_in = _sos_resp(butter_sos("lowpass", in_bw["order"], in_bw["f_c_hz"], fs, in_bw["method"]), f, fs) \
        if in_bw["toggles"]["enabled"] else np.ones_like(zi)
    adc = leaves["dsc1.adc"]
    n_avg = adc.n_avg()
    h_avg = sum(zi ** j for j in range(n_avg)) / n_avg
    pid = leaves["dsc1.pid"]
    r = pid.rate_div
    h_hold = sum(zi ** j for j in range(r))
    a = h_hold * h_out * zi * h_delay * amp * k_act * h_str * (-1.0) * h_tr  # φ_D = φ_lo − φ_pr: −ψ_a1
    b = k_d * h_in * h_avg / r
    c = pid.freq_response(f, params)
    g = leaves["dsc1.dac"].gain()
    return a, b, c, g


def sensitivity_z(spec: RunSpec, f_hz: np.ndarray, *, k_d: float) -> np.ndarray:
    """Exact φ_D(f)/d(f) of the multirate loop (module docstring; PID rate_div r ∈ {1, 2})."""
    fs = spec.run.f_sim_hz
    f = np.atleast_1d(np.asarray(f_hz, float))
    a, b, c, g = loop_parts(spec, f, k_d=k_d)
    r = built_leaves(spec)[0]["dsc1.pid"].rate_div
    if r == 1:
        return 1.0 / (1.0 - g * c * b * a)
    if r != 2:
        raise NotImplementedError("rate_div 1 or 2 only")
    a2, b2, _, _ = loop_parts(spec, f + fs / 2, k_d=k_d)
    return (1.0 - g * c * b2 * a2) / (1.0 - g * c * (b * a + b2 * a2))
