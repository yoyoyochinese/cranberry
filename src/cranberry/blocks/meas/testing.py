"""Spec builders for metric / offline-stage tests (WP-F4). Not physics; registers no block kinds.

Presets must still give every value explicitly (D34); these helpers only shorten tests. Path blocks
are written without a ``topology`` param, so the runner injects ``RunSpec.topology`` (WP-F4).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

__all__ = ["REF_A", "path_blocks", "metric_cfg", "refa_source", "run_spec", "single_port_topology", "toy"]

#: Ref A single-mode reference (prompt §7-1, T04 §9.1): L = 8 %, α = 9.06 /W, P = 0.64 W, θ = 9 mrad.
REF_A: dict[str, Any] = dict(loss=0.08, alpha_per_w=9.06, pump_w=0.64, theta_rms_rad=0.009, phi_opa0_rad=0.0)


def toy(kind: str, name: str, **params: Any) -> dict[str, Any]:
    """``core.testing.toy`` (imported lazily so this module does not register the toy kinds)."""
    from cranberry.core.testing import toy as _toy

    return _toy(kind, name, **params)


def path_blocks(sources: Mapping[str, str], *, locks: bool = False, emission_time: bool = True) -> list[dict[str, Any]]:
    """``segment_phase`` → ``beam_assembler`` (→ ``lock_error``) without topology params (injected)."""
    out = [{"block": "segment_phase", "name": "sp", "params": {
                "segments": "all", "toggles": {"phase": True, "group_delay": True, "pm_axis_drift": False}}},
           {"block": "beam_assembler", "name": "ba", "params": {
                "beams": "all", "sources": dict(sources),
                "toggles": {"emission_time": emission_time, "polarization": False}}}]
    if locks:
        out.append({"block": "lock_error", "name": "le", "params": {"locks": "all"}})
    return out


def metric_cfg(source: Any, *, f_max_hz: Any = 1.0e12, weight: Any = None, lo_mode: Any = None,
               average: str = "exact", include_hf_phase_noise: bool = True,
               dither_m_rad: Sequence[float] = ()) -> dict[str, Any]:
    """One ``measure.metrics.<i|ii>`` entry (every key given, D34)."""
    return {"f_max_hz": f_max_hz, "weight": weight, "source": source, "lo_mode": lo_mode, "average": average,
            "include_hf_phase_noise": include_hf_phase_noise, "dither_m_rad": list(dither_m_rad)}


def refa_source(name: str = "refA.sqz", **over: Any) -> dict[str, Any]:
    """Offline ``single_mode_reference`` block spec with the Ref A values."""
    return {"block": "single_mode_reference", "name": name, "params": {**REF_A, "average": "exact", **over}}


def single_port_topology(*, lo_len_m: float = 1.0, sq_len_m: float = 1.0, sq_driver: bool = False) -> dict[str, Any]:
    """LO and squeezed beam meeting at a BHD observer ``bhd``; the LO segment ``lo_seg`` has a
    ``stretcher`` driver (1 rad per unit), the squeezed segment ``sq_seg`` a ``strain`` driver
    (1 rad per unit) when ``sq_driver``."""
    import cranberry.path.testing as pt  # lazy: importing it registers the pathtest kinds

    segs = [pt.seg_spec("lo_seg", lo_len_m, kind="free_space", drivers={"stretcher": "uniform"}),
            pt.seg_spec("sq_seg", sq_len_m, kind="free_space", drivers={"strain": "uniform"} if sq_driver else {})]
    beams = [pt.beam_spec("lo", [pt.step_spec("lo_seg")], "bhd"),
             pt.beam_spec("sqz", [pt.step_spec("sq_seg")], "bhd")]
    coef = pt.coef_spec("lo_seg", "stretcher", "lo", 1.0)
    if sq_driver:
        coef += pt.coef_spec("sq_seg", "strain", "sqz", 1.0)
    return pt.topo_spec(segs, beams, ["bhd"], coef)


def run_spec(pipeline: list, *, f_sim_hz: float, n_samples: int, topology: Mapping[str, Any] | None = None,
             offline: Mapping[str, Any] | None = None, measure: Mapping[str, Any] | None = None,
             toggles: Mapping[str, bool] | None = None, seed: int = 1, chunk_samples: int = 8192):
    """A ``RunSpec`` from a pipeline list and sections."""
    from cranberry.config.spec import RunSettings, RunSpec
    from cranberry.core.pipeline import spec_from_obj

    return RunSpec(run=RunSettings(duration_s=n_samples / f_sim_hz, f_sim_hz=f_sim_hz, chunk_samples=chunk_samples,
                                   seed=seed),
                   pipeline=spec_from_obj(pipeline), topology=dict(topology or {}), offline=dict(offline or {}),
                   measure=dict(measure or {}), toggles=dict(toggles or {}))


#: Noise sources of the two-port toy loop → metric (ii) component (``measure.decompose.groups``).
TWO_PORT_GROUPS: dict[str, list[str]] = {
    "locked_path_residual": ["n_lock.enabled"],
    "common_points": ["n_cp.enabled"],
    "laser_decorrelation_oob": ["n_las.enabled"],
    "squeezed_only": ["n_sq.enabled"],
}


def two_port_topology(*, l1_m: float = 200.0) -> dict[str, Any]:
    """Toy two-port layout (G13: reference port 0 + port 1) for metric (ii) tests.

    sqz0: a0 → com → gbs; sqz1: a1 → s1 → com → gbs; probes pr1: a1 → det1 and ref0: a0 → det1 (lock
    L1 holds φ_pr1 − φ_ref0 = 0 with the stretcher on a1, a Liu-type probe lock); lo1: l1 → gbs
    (metric (i) of port 1). Drivers (1 rad per unit):
    a0.strain (reference path), a1.dT (locked path) + a1.stretcher (actuator), s1.strain (squeezed
    only), com.eom (common to both squeezed beams just before the GBS: cancels in Θ_1 exactly).
    """
    import cranberry.path.testing as pt

    segs = [pt.seg_spec("a0", 1.0, kind="free_space", drivers={"strain": "uniform"}),
            pt.seg_spec("a1", l1_m, kind="smf", drivers={"dT": "uniform", "stretcher": "uniform"}),
            pt.seg_spec("s1", 1.0, kind="free_space", drivers={"strain": "uniform"}),
            pt.seg_spec("com", 1.0, kind="free_space", drivers={"eom": "uniform"}),
            pt.seg_spec("l1", 1.0, kind="free_space")]
    beams = [pt.beam_spec("sqz0", [pt.step_spec("a0"), pt.step_spec("com")], "gbs"),
             pt.beam_spec("sqz1", [pt.step_spec("a1"), pt.step_spec("s1"), pt.step_spec("com")], "gbs"),
             pt.beam_spec("pr1", [pt.step_spec("a1")], "det1"),
             pt.beam_spec("ref0", [pt.step_spec("a0")], "det1"),
             pt.beam_spec("lo1", [pt.step_spec("l1")], "gbs")]
    coef = (pt.coef_spec("a0", "strain", "*", 1.0) + pt.coef_spec("a1", "dT", "*", 1.0)
            + pt.coef_spec("a1", "stretcher", "*", 1.0) + pt.coef_spec("s1", "strain", "*", 1.0)
            + pt.coef_spec("com", "eom", "*", 1.0))
    obs = [{"id": "gbs", "kind": "gbs_input"}, {"id": "det1", "kind": "detector"}]
    locks = [{"id": "L1", "combo": {"pr1": 1.0, "ref0": -1.0}, "setpoint_rad": 0.0, "observer": "det1"}]
    ports = [{"index": 0, "squeezed": "sqz0", "lo": None, "probe": None, "gbs_observer": "gbs", "bhd_observer": None},
             {"index": 1, "squeezed": "sqz1", "lo": "lo1", "probe": "pr1", "gbs_observer": "gbs",
              "bhd_observer": "gbs"}]
    return pt.topo_spec(segs, beams, obs, coef, locks, ports)


def two_port_loop_spec(*, n_samples: int = 1 << 17, f_sim_hz: float = 1.0e6, seed: int = 3,
                       sat_limit: float | None = None, com_noise: bool = True, lock: bool = True,
                       record_decim: int = 16, toggles: Mapping[str, bool] | None = None):
    """Linear toy loop on :func:`two_port_topology`: AR(1) noises n_lock (a1.dT), n_cp (a0.strain),
    n_sq (s1.strain), n_com (com.eom), a near random walk n_las (laser phase); PI lock on pr1 through a
    1-sample delay buffer to the a1 stretcher (``toy_sat`` with ``sat_limit`` instead of a unit gain
    when given). Metric probes ``mi`` (port 1) and ``mii`` (port 1 vs 0); ``measure.decompose`` uses
    :data:`TWO_PORT_GROUPS` and, with ``sat_limit``, the actuator output as validity slot."""
    h = 1.0 / f_sim_hz
    act = (toy("toy_sat", "act", in_="delay.dl.out", out="drv.a1.stretcher", unit="1", gain=1.0, limit=sat_limit)
           if sat_limit is not None else
           toy("toy_gain", "act", in_="delay.dl.out", out="drv.a1.stretcher", unit="1", gain=1.0))
    pipe = [toy("toy_noise_ar1", "n_lock", out="drv.a1.dT", unit="K", a=0.999, b=3e-3),
            toy("toy_noise_ar1", "n_cp", out="drv.a0.strain", unit="1", a=0.9, b=4e-3),
            toy("toy_noise_ar1", "n_sq", out="drv.s1.strain", unit="1", a=0.8, b=5e-3),
            toy("toy_noise_ar1", "n_com", out="drv.com.eom", unit="1", a=0.95, b=2e-2, enabled=com_noise),
            toy("toy_noise_ar1", "n_las", out="laser.las.phi", unit="rad", a=1.0 - 1e-5, b=1e-2),
            *path_blocks({"las": "laser.las.phi"}, locks=True),
            toy("toy_pid", "pid", in_="lock.L1.err", out="ctrl.pid.u", unit="V", kp=0.3 if lock else 0.0,
                ki=2 * math.pi * 1e4 if lock else 0.0, kd=0.0, setpoint=0.0),
            {"block": "delay_buffer", "name": "dl", "params": {"in": "ctrl.pid.u", "out": "delay.dl.out",
                                                              "tau_s": h, "method": "linear"}},
            act,
            {"block": "metric_i", "name": "mi", "params": {"ports": "all", "setpoint_rad": 0.0}},
            {"block": "metric_ii", "name": "mii", "params": {"ports": "all", "setpoint_rad": 0.0}}]
    validity = [{"slot": "drv.a1.stretcher", "abs_max": sat_limit}] if sat_limit is not None else []
    measure = {"record": {"slots": ["lock.L1.err"], "decim": record_decim},
               "decompose": {"groups": TWO_PORT_GROUPS, "validity": validity}}
    return run_spec(pipe, f_sim_hz=f_sim_hz, n_samples=n_samples, topology=two_port_topology(), measure=measure,
                    seed=seed, toggles=toggles)
