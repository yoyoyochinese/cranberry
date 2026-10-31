"""Test support for the ERR blocks: a rig builder and parameter sets (values from T14/T15; test use only).

Values are the notes' reference / placeholder numbers, not presets: Ref A [Ref A]/[AE] (T14 §3.2–3.4),
Ref B with x = √(P/P_th) = √(100/180) (T02, D51-10), Ref C (T14 §6.3), Ref D (T15 §7), detector
placeholders (T14 §2.5, 仮置き). Presets proper live in ``presets/`` (WP-I*).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np

import cranberry.blocks.err  # noqa: F401  (registers the kinds)
from cranberry.blocks.err._common import amplitude_gains_opa, amplitude_gains_opo
from cranberry.core.layout import SlotLayout
from cranberry.core.compile import run_compiled
from cranberry.core.pipeline import Seq, build, from_spec, spec_from_obj
from cranberry.core.state import Params, State

__all__ = ["rig", "record", "blk", "pbd", "bhd", "liu", "qnl", "dcf", "dem", "rfr", "qnl_flat", "PHASES", "REF_A_SQ", "REF_A_DET", "REF_B", "REF_C", "X_REF_B", "DET"]

PHASES = ("beam.pr.phi", "beam.p.phi", "beam.lo.phi", "beam.ref.phi")

#: Ref A squeezer r = √(9.06 × 0.64) [Ref A], e^{2r} = 123.5; detection OPA e^{2r_d} = 10 (仮置き, T14 §3.2)
REF_A_SQ = amplitude_gains_opa(math.sqrt(9.06 * 0.64))
REF_A_DET = amplitude_gains_opa(0.5 * math.log(10.0))
#: Ref B: x = √(P/P_th), P = 100 mW, P_th = 180 mW (T02 §0, D51-10) ≈ 0.745
X_REF_B = math.sqrt(100.0 / 180.0)
REF_B = amplitude_gains_opo(X_REF_B)
#: detector placeholders (T14 §2.5): ℜ = 1 A/W, G = 1e3 V/A, i_n = 10 pA/√Hz
DET = {"responsivity_a_per_w": 1.0, "gain_v_per_a": 1.0e3, "i_n_a_rthz": 10e-12}


def _db(x: float) -> float:
    return 10.0 ** (x / 10.0)


#: Ref C (T14 §6.3): R+ = +10 dB, R− = −4.7 dB, dark noise −11 dB, B = 292 MHz, m = 0.0785 rad
REF_C = {"r_plus": _db(10.0), "r_minus": _db(-4.7), "el": _db(-11.0), "b_hz": 292e6, "m": 0.0785}


def qnl_flat(n_sn: float = 1.0) -> dict[str, float]:
    """quantum_noise_power runtime params for a flat band with the Ref C numbers (R± incl. dark noise)."""
    rp, rm, b = REF_C["r_plus"] + REF_C["el"], REF_C["r_minus"] + REF_C["el"], REF_C["b_hz"]
    pm, pd = n_sn * b * (rp + rm) / 2, n_sn * b * (rp - rm) / 2
    return {"p_mean": pm, "p_delta": pd, "p_dark": n_sn * b * REF_C["el"], "q0": pm * pm / b, "q1": -2 * pm * pd / b,
            "q2": pd * pd / b}


def blk(kind: str, name: str, *, toggles: dict[str, bool], rate_div: int = 1, **params: Any) -> dict[str, Any]:
    """Pipeline dict element (trailing ``_`` in a key is dropped: ``in_``)."""
    p = {k.rstrip("_") if k.endswith("_") else k: v for k, v in params.items()}
    p["toggles"] = dict(toggles)
    obj: dict[str, Any] = {"block": kind, "name": name, "params": p}
    if rate_div != 1:
        obj["rate_div"] = rate_div
    return obj


def rig(blocks: list[Any], f_sim_hz: float, *, seed: int = 1, external: Iterable[str] | Mapping[str, str] = PHASES):
    """Build ``blocks`` with externally written scalar slots (names → unit, or names with unit rad);
    returns (root, State, Params)."""
    lay = SlotLayout(f_sim_hz)
    units = dict(external) if isinstance(external, Mapping) else dict.fromkeys(external, "rad")
    ext = tuple(units)
    for s, u in units.items():
        lay.add(s, (), np.float64, u)
    root = build(spec_from_obj(blocks), lay, seed, external=ext)
    return root, State(lay), Params(lay)


def record(blocks: list[Any], f_sim_hz: float, n: int, slots: Iterable[str], *, seed: int = 1,
           external: Iterable[str] | Mapping[str, str] = PHASES, init: Mapping[str, float] | None = None,
           chunk: int = 1 << 16) -> dict[str, np.ndarray]:
    """Run ``blocks`` (compiled) for ``n`` samples and return every sample of ``slots`` (runner's Recorder)."""
    from cranberry.runner.run import Recorder

    lay = SlotLayout(f_sim_hz)
    units = dict(external) if isinstance(external, Mapping) else dict.fromkeys(external, "rad")
    for s, u in units.items():
        lay.add(s, (), np.float64, u)
    rec = Recorder(list(slots), 1, chunk)
    root = build(Seq([from_spec(spec_from_obj(blocks)), rec]), lay, seed, external=tuple(units))
    st, pa = State(lay), Params(lay)
    for k, v in (init or {}).items():
        st[k] = v
    rows = []
    done = 0
    while done < n:
        m = min(chunk, n - done)
        run_compiled(root, st, pa, m)
        rows.append(rec.harvest(st))
        done += m
    return rec.split(np.concatenate(rows), lay)


_DITHER = {"f_dither_hz": 1.0e6, "ram_a": 0.0, "ram_phase_rad": 0.0}


def _mk(kind, name, base, tg, over, tg_over, rate_div=1):
    base.update(over)
    tg.update(tg_over or {})
    return blk(kind, name, toggles=tg, rate_div=rate_div, **base)


def pbd(name="pbd", mode="heterodyne", *, rf=False, noise=False, toggles=None, rate_div=1, **over):
    """parametric_beat_detector: Ref A phase-detection OPA (heterodyne) or Ref B lock 1 (dither)."""
    p = {"mode": mode, "in_pr": "beam.pr.phi", "in_p": "beam.p.phi", "out": f"v.{name}", "phi_opa_rad": 0.0,
         "p_in_w": 1e-3, "eta": 0.8, **DET}
    tg = {"shot_noise": noise, "electronic_noise": noise, "rf_rate_mode": rf}
    if mode == "heterodyne":
        p.update(a_plus=REF_A_DET[0], a_minus=REF_A_DET[1], f_shift_hz=1.0e6, f_sh_leak_v=0.0)
        tg["f_sh_leak"] = False
    else:
        p.update(a_plus=REF_B[0], a_minus=REF_B[1], dither_m_rad=0.2, **_DITHER)
        p["f_dither_hz"] = 295e3
        tg["ram_offset"] = False
    return _mk("parametric_beat_detector", name, p, tg, over, toggles, rate_div)


def bhd(name="bhd", mode="heterodyne", *, rf=False, noise=False, toggles=None, rate_div=1, **over):
    """bhd_probe_beat: Ref A lock 2 (heterodyne, squeezer A±) or Ref B lock 2 (dither / dc)."""
    p = {"mode": mode, "in_pr": "beam.pr.phi", "in_p": "beam.p.phi", "in_lo": "beam.lo.phi", "out": f"v.{name}",
         "phi_opa_rad": 0.0, "p_lo_w": 1e-3, "p_pr_w": 1e-9, "noise_rel_shot": 1.0, **DET}
    tg = {"shot_noise": noise, "electronic_noise": noise}
    if mode != "dc":
        tg["rf_rate_mode"] = rf
    if mode == "heterodyne":
        p.update(a_plus=REF_A_SQ[0], a_minus=REF_A_SQ[1], f_shift_hz=1.0e6)
    else:
        p.update(a_plus=REF_B[0], a_minus=REF_B[1])
        if mode == "dither":
            p.update(dither_m_rad=0.1, **_DITHER)
            p["f_dither_hz"] = 295e3
            tg["ram_offset"] = False
    return _mk("bhd_probe_beat", name, p, tg, over, toggles, rate_div)


def liu(name="liu", mode="homodyne_dc", *, rf=False, noise=False, gating=True, toggles=None, rate_div=1, **over):
    """liu_probe_interferometer (T14 §8.2 F placeholders: P̄_pr 1 µW, P_ref 1 mW, G 1e4)."""
    p = {"mode": mode, "in_pr": "beam.pr.phi", "in_ref": "beam.ref.phi", "out": f"v.{name}", "p_pr_avg_w": 1e-6,
         "p_ref_w": 1e-3, "duty": 0.1, "visibility": 1.0, "coupler_sign": 1.0, **DET}
    p["gain_v_per_a"] = 1e4
    tg = {"shot_noise": noise, "electronic_noise": noise, "gating": gating}
    if mode != "homodyne_dc":
        tg["rf_rate_mode"] = rf
    if mode == "homodyne_dc":
        p["v_off_v"] = 0.0
    elif mode == "heterodyne":
        p["f_shift_hz"] = 1.0e6
    else:
        p.update(dither_m_rad=0.1, **_DITHER)
        tg["ram_offset"] = False
    return _mk("liu_probe_interferometer", name, p, tg, over, toggles, rate_div)


def qnl(name="qnl", law="square", *, rf=False, noise=False, toggles=None, rate_div=1, **over):
    """quantum_noise_power with Ref C flat-band numbers (k = 1 V per power unit)."""
    p = {"detector_law": law, "in_lo": "beam.lo.phi", "in_p": "beam.p.phi", "out": f"v.{name}", "phi_opa_rad": 0.0,
         **qnl_flat(), "k_det": 1.0, "kappa_det": 1.0, "dither_m_rad": REF_C["m"], **_DITHER}
    tg = {"radiometer_noise": noise, "ram_offset": False, "rf_rate_mode": rf}
    return _mk("quantum_noise_power", name, p, tg, over, toggles, rate_div)


def dcf(name="dcf", *, noise=False, toggles=None, rate_div=1, **over):
    """dc_fringe with Ref D numbers (T15 §7): P_in 3.33 mW, ℜ 1, G 1e3, ideal coupler, μ = 1."""
    p = {"in_lo": "beam.lo.phi", "in_pr": "beam.pr.phi", "in_dp": None, "in_voff": None, "out": f"v.{name}",
         "out_err": None, "p_in_w": 3.33e-3, "kappa1": 0.5, "kappa2": 0.5, "t_arm1": 1.0, "t_arm2": 1.0,
         "responsivity_a_per_w": 1.0, "eps_r": 0.0, "gain_v_per_a": 1e3, "v_off_v": 0.0, "mu_abs": 1.0,
         "mu_phase_rad": 0.0, "gamma_coh": 1.0, "per_x": 0.0, "per_dchi_rad": 0.0, "coupler_sign": 1.0,
         "i_n_a_rthz": 10e-12, "v_sat_v": 4.6, "setpoint_v": 0.0, "loop_sign": 1.0}
    tg = {"shot_noise": noise, "electronic_noise": noise, "rin": False, "offset_drift": False,
          "optical_imbalance": True, "per_crosscoupling": True, "visibility_hf_averaging": True, "saturation": True}
    return _mk("dc_fringe", name, p, tg, over, toggles, rate_div)


def dem(name="dem", inp="v.pbd", *, rf=False, order=0, toggles=None, rate_div=1, **over):
    """demodulator (default: 2 f_sh = 2 MHz reference, LPF off)."""
    p = {"in": inp, "in_phase": None, "out_i": f"dem.{name}.i", "out_q": f"dem.{name}.q", "lpf_order": order,
         "f_ref_hz": 2.0e6, "phi_d_rad": 0.0, "v_off_v": 0.0, "lpf_hz": 2.0e4}
    tg = {"rf_rate_mode": rf, "demod_offset": False, "rf_reference_phase_noise": False}
    return _mk("demodulator", name, p, tg, over, toggles, rate_div)


def rfr(name="rfr", *, on=True, rate_div=1, **over):
    """rf_reference (white phase 1e-12 rad²/Hz + white FM 1e-2 Hz²/Hz; arbitrary test values)."""
    p = {"out": f"rf.{name}.phi", "band_targets": ["beam.pr"], "s_phi_white_rad2_hz": 1e-12, "s_nu_white_hz2_hz": 1e-2}
    return _mk("rf_reference", name, p, {"phase_noise": on}, over, None, rate_div)
