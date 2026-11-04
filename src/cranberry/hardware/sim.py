"""``SimController``: the ControllerInterface implemented on the time-loop simulation (PLAN §2.10).
Owner: WP-F5. Runs through ``runner.run(spec) -> RunResult`` (WP-F1).

Semantics
    The controller holds a base ``RunSpec`` plus the state set through the interface (parameter
    overrides, engaged flag, ramp and injection settings). ``acquire(duration_s, signals)`` builds the
    effective spec (base spec, ``run.duration_s = duration_s``, overrides written into the block params,
    toggles written into ``RunSpec.toggles``) and runs it from its initial state. Nothing is carried
    over between calls: each ``acquire`` is a reproducible experiment with the current settings
    (integrators start at their initial state, so ``reset_integrators`` is recorded but changes
    nothing). The seed is the spec's ``run.seed`` unless changed with ``set_params({"run.seed": n})``.

Block contract used by SimController (names the blocks of WP-B4 must provide)
    * Loop controller (``loops.<id>.pid_block``): runtime params ``kp``, ``ki_per_s``, ``kd_s`` (and
      optionally ``omega_is_rad_per_s``, ``omega_ds_rad_per_s``, ``ki2_per_s2``; 0 = none), i.e.
      ``hardware.mapping.CANONICAL_KEYS``, and a toggle ``feedback`` (T15 §8.4) that opens the loop
      when False. ``engage``/``disengage`` set ``"<pid_block>.feedback"`` in ``RunSpec.toggles``.
    * Ramp generator (``ramps.<channel>.block``): runtime params ``amplitude_v``, ``period_s``,
      ``offset_v`` and a toggle ``enabled``.
    * Sine injection (``injections.<point>.block``): runtime params ``freq_hz``, ``amplitude`` (unit of
      the injected signal) and a toggle ``enabled``.
    * Excitation (``excitations.<point>.block``, D85-7): an ``excitation`` block (runtime ``amplitude``,
      ``t_start_s``; toggle ``enabled``). Its waveform keys are structural; ``measure_loop_gain`` writes
      them into the block spec of that one run.

``configure`` config (all names are block / slot names of the spec)::

    loops:      {<loop id>: {pid_block: <name>, device: <devices/*.yaml name> | null}}
    ramps:      {<channel>: {block: <name>, device: <name> | null}}                     # optional
    injections: {<point>: {block: <name>, reference: <slot>, response: <slot>,
                           amplitude: <float>, settle_s: <float>, periods: <int>}}     # optional
    excitations: {<point>: {block: <excitation name>, loop: <loop id>, relation: minus_l | plant,
                            reference: <slot> | null, response: <slot>,
                            controller: {kp, ki_per_s, kd_s[, discretisation, h_s]} | model | null,
                            model: <loop_stability_check name in offline.blocks> | null}}   # optional

With a device, PID gains and ramp settings are clipped exactly as that device would clip them
(``hardware.mapping.apply_device_limits``, device ``ramp`` ranges), so simulation and hardware see the
same applied values.

Loops, engage and injections (D85-7, backward compatible)
    ``engage(loop_id)`` / ``disengage(loop_id)`` switch one configured loop (None: all, the original
    call); configure each stage of a stack that must switch alone as its own loop (``pid_block`` of that
    stage). ``acquire(..., injections_on=True | [points])`` records with injections / excitations on at
    their current runtime parameters (default False: all off, as before).

``measure_loop_gain(loop_id, method, freqs_hz, *, band_hz, amplitude, **settings)`` (the excitation
point whose ``loop`` is ``loop_id``, or ``point=``): one run with the excitation on and the loops as
engaged; the waveform keys go into the excitation block spec, ``amplitude`` / ``t_start_s`` are runtime
overrides of that run only; the excitation's x and the reference / response slots are added to
``measure.record`` (which the spec must have: its ``decim`` must put the samples on the controller's
update instants, e.g. 2 for the DSC1 PID at rate_div 2). Settings: ``settle_s`` (required), ``t_start_s``
(default 0), ``dwell_s`` (stepped_sine), ``period_s`` + ``n_periods`` (multisine, chirp), ``prbs_order`` +
``bit_s`` + ``n_periods`` (prbs), ``duration_s`` + ``nperseg_s`` (white), ``n_freqs`` (with ``band_hz``
instead of ``freqs_hz``), ``amp_f_corner_hz`` (default inf: flat), ``n_blocks`` (default 4),
``coherence_min`` (default 0). The run duration follows from the schedule (t_start + steps × dwell, or
t_start + settle rounded up to whole periods + n_periods periods). Analysis: ``identify`` of
``meas.loop_identification`` (L̂, S, T, coherence, margins; with ``model`` also the model comparison).
"""

from __future__ import annotations

import dataclasses
import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from cranberry.core.block import BlockSpec, MissingParameterError
from cranberry.core.pipeline import ParallelSpec, PipelineSpec, SeqSpec, StackSpec, from_spec, iter_leaves
from cranberry.hardware.devices import device_values, load_device
from cranberry.hardware.interface import ControllerInterface
from cranberry.hardware.mapping import CANONICAL_KEYS, DeviceRangeWarning, PidCoefficients, apply_device_limits

if TYPE_CHECKING:
    from cranberry.config.spec import RunSpec
    from cranberry.runner.run import RunResult

__all__ = ["SimController", "demodulate", "LOOP_TOGGLE", "RAMP_PARAMS", "INJECTION_PARAMS", "ENABLE_TOGGLE",
           "EXCITATION_PARAMS"]

#: Toggle of a loop's controller block that closes the loop (T15 §8.4).
LOOP_TOGGLE = "feedback"
#: Toggle of ramp / injection blocks.
ENABLE_TOGGLE = "enabled"
#: Runtime params of a ramp block.
RAMP_PARAMS: tuple[str, ...] = ("amplitude_v", "period_s", "offset_v")
#: Runtime params of a sine-injection block.
INJECTION_PARAMS: tuple[str, ...] = ("freq_hz", "amplitude")
_PID_REQUIRED: tuple[str, ...] = ("kp", "ki_per_s", "kd_s")
_INJ_KEYS = frozenset({"block", "reference", "response", "amplitude", "settle_s", "periods"})
_EXC_KEYS = frozenset({"block", "loop", "relation", "reference", "response", "controller", "model"})
#: Runtime params of an excitation block (D85-7).
EXCITATION_PARAMS: tuple[str, ...] = ("amplitude", "t_start_s")


def demodulate(x: np.ndarray, t_s: np.ndarray, f_hz: float) -> complex:
    """Single-bin lock-in: (2/N) Σ x(t) e^{−i2πft}. For x = A cos(2πft + φ) over whole periods this is
    A e^{iφ}, so response/reference = H(i·2πf) in the Laplace convention (s = +iω)."""
    x = np.asarray(x, dtype=np.float64)
    return complex(2.0 / x.size * np.sum(x * np.exp(-2j * np.pi * f_hz * np.asarray(t_s, dtype=np.float64))))


def _replace_params(spec: PipelineSpec, per_block: Mapping[str, Mapping[str, Any]]) -> PipelineSpec:
    """Copy of ``spec`` with ``per_block[name]`` merged into the params of the leaf named ``name``."""
    if isinstance(spec, BlockSpec):
        upd = per_block.get(spec.name)
        return dataclasses.replace(spec, params={**spec.params, **upd}) if upd else spec
    if isinstance(spec, (SeqSpec, ParallelSpec)):
        return dataclasses.replace(spec, children=tuple(_replace_params(c, per_block) for c in spec.children))
    if isinstance(spec, StackSpec):
        return dataclasses.replace(
            spec,
            stages=tuple(_replace_params(c, per_block) for c in spec.stages),
            actuators=tuple(_replace_params(c, per_block) for c in spec.actuators),
        )
    raise TypeError(f"not a pipeline spec: {spec!r}")


def _jsonable(x: Any) -> Any:
    if isinstance(x, Mapping):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


class SimController(ControllerInterface):
    """Simulation-backed controller (module docstring).

    Args:
        spec: fully resolved ``RunSpec`` (no sweep markers) whose pipeline contains the controller,
            ramp and injection blocks named in ``configure``.
        run_fn: the runner, ``RunSpec -> RunResult``; None → ``cranberry.runner.run.run`` (WP-F1),
            looked up at call time. Tests pass a stand-in.
    """

    def __init__(self, spec: RunSpec, *, run_fn: Callable[[RunSpec], RunResult] | None = None):
        from cranberry.core import registry

        registry.ensure_loaded()
        self._spec = spec
        self._run_fn = run_fn
        from cranberry.runner.run import inject_topology

        # the runner injects RunSpec.topology into path blocks at build time; do the same so presets that rely on it
        # (refD, F_*) can be driven through the interface (D85-7 addition)
        pipe = inject_topology(spec.pipeline, spec.topology)
        self._leaves = {leaf.name: leaf for leaf in iter_leaves(from_spec(pipe))}
        #: full runtime-param name → (block name, key)
        self._runtime: dict[str, tuple[str, str]] = {
            f"{name}.{key}": (name, key) for name, leaf in self._leaves.items() for key in leaf.runtime_params
        }
        self._config: dict[str, Any] | None = None
        self._reset_state()

    # ------------------------------------------------------------ internal state

    def _reset_state(self) -> None:
        self._loops: dict[str, dict[str, Any]] = {}
        self._ramps: dict[str, dict[str, Any]] = {}
        self._injections: dict[str, dict[str, Any]] = {}
        self._excitations: dict[str, dict[str, Any]] = {}
        self._devices: dict[str, dict[str, Any]] = {}
        self._overrides: dict[str, float] = {}
        self._ramp_on: dict[str, bool] = {}
        self._engaged: dict[str, bool] = {}
        self._inj_on: set[str] = set()
        self._seed = self._spec.run.seed
        self._last: RunResult | None = None
        self._ops: list[dict[str, Any]] = []
        self._warnings: list[str] = []

    def _leaf(self, name: str, role: str, params: Sequence[str], toggle: str) -> None:
        leaf = self._leaves.get(name)
        if leaf is None:
            raise KeyError(f"{role}: no block named {name!r} in the pipeline; blocks {sorted(self._leaves)}")
        missing = [p for p in params if p not in leaf.runtime_params]
        if missing:
            raise ValueError(f"{role} block {name!r} lacks runtime param(s) {missing} (SimController contract)")
        if toggle not in leaf.toggles:
            raise ValueError(f"{role} block {name!r} lacks toggle {toggle!r} (SimController contract)")

    def _device(self, name: str | None) -> dict[str, Any] | None:
        if name is None:
            return None
        if name not in self._devices:
            self._devices[name] = load_device(name)
        return device_values(self._devices[name])

    def _require_configured(self) -> None:
        if self._config is None:
            raise RuntimeError("SimController.configure(config) must be called first")

    def _value(self, full: str) -> Any:
        if full in self._overrides:
            return self._overrides[full]
        block, key = self._runtime[full]
        v = np.asarray(self._leaves[block].params[key], dtype=np.float64)
        return float(v) if v.ndim == 0 else v.copy()

    def _toggle_state(self) -> dict[str, bool]:
        t = {f"{lp['pid_block']}.{LOOP_TOGGLE}": self._engaged[lid] for lid, lp in self._loops.items()}
        t.update({f"{r['block']}.{ENABLE_TOGGLE}": self._ramp_on[ch] for ch, r in self._ramps.items()})
        t.update({f"{i['block']}.{ENABLE_TOGGLE}": pt in self._inj_on for pt, i in self._injections.items()})
        t.update({f"{x['block']}.{ENABLE_TOGGLE}": pt in self._inj_on for pt, x in self._excitations.items()})
        return t

    def _effective_spec(self, duration_s: float, extra_params: Mapping[str, float] = {}, extra_toggles: Mapping[str, bool] = {},
                        structural: Mapping[str, Mapping[str, Any]] = {}, record_slots: Sequence[str] = ()) -> RunSpec:
        per_block: dict[str, dict[str, Any]] = {b: dict(v) for b, v in structural.items()}
        for full, v in {**self._overrides, **extra_params}.items():
            block, key = self._runtime[full]
            per_block.setdefault(block, {})[key] = v
        run = dataclasses.replace(self._spec.run, duration_s=float(duration_s), seed=self._seed)
        toggles = {**self._spec.toggles, **self._toggle_state(), **extra_toggles}
        measure = self._spec.measure
        if record_slots:
            rec = (measure or {}).get("record")
            if rec is None:
                raise KeyError("measure_loop_gain needs measure.record in the spec (its decim sets the sample instants; D34)")
            slots = list(rec["slots"]) + [s for s in record_slots if s not in rec["slots"]]
            measure = {**measure, "record": {**rec, "slots": slots}}
        return dataclasses.replace(self._spec, run=run, pipeline=_replace_params(self._spec.pipeline, per_block),
                                   toggles=toggles, measure=measure)

    def _loop_ids(self, loop_id: str | None) -> list[str]:
        if loop_id is None:
            return list(self._loops)
        if loop_id not in self._loops:
            raise KeyError(f"unknown loop {loop_id!r}; configured {sorted(self._loops)}")
        return [loop_id]

    def _run(self, spec: RunSpec) -> RunResult:
        fn = self._run_fn
        if fn is None:
            from cranberry.runner.run import run as fn
        self._last = fn(spec)
        return self._last

    def _note(self, msgs: Sequence[str]) -> None:
        self._warnings.extend(msgs)

    # ------------------------------------------------------------ ControllerInterface

    def configure(self, config: Mapping[str, Any]) -> None:
        unknown = set(config) - {"loops", "ramps", "injections", "excitations"}
        if unknown:
            raise ValueError(f"unknown SimController config keys {sorted(unknown)}")
        if "loops" not in config:
            raise MissingParameterError("SimController config needs 'loops' (may be {})")
        self._reset_state()
        for lid, lp in config["loops"].items():
            if set(lp) != {"pid_block", "device"}:
                raise ValueError(f"loops.{lid}: expected keys {{pid_block, device}}, got {sorted(lp)}")
            self._leaf(lp["pid_block"], f"loops.{lid}", _PID_REQUIRED, LOOP_TOGGLE)
            self._device(lp["device"])
            self._loops[lid] = dict(lp)
        for ch, r in config.get("ramps", {}).items():
            if set(r) != {"block", "device"}:
                raise ValueError(f"ramps.{ch}: expected keys {{block, device}}, got {sorted(r)}")
            self._leaf(r["block"], f"ramps.{ch}", RAMP_PARAMS, ENABLE_TOGGLE)
            self._device(r["device"])
            self._ramps[ch] = dict(r)
            self._ramp_on[ch] = False
        for pt, inj in config.get("injections", {}).items():
            if set(inj) != _INJ_KEYS:
                raise ValueError(f"injections.{pt}: expected keys {sorted(_INJ_KEYS)}, got {sorted(inj)}")
            self._leaf(inj["block"], f"injections.{pt}", INJECTION_PARAMS, ENABLE_TOGGLE)
            if not (inj["settle_s"] >= 0 and inj["periods"] >= 1):
                raise ValueError(f"injections.{pt}: need settle_s >= 0 and periods >= 1")
            self._injections[pt] = dict(inj)
        for pt, ex in config.get("excitations", {}).items():
            self._configure_excitation(pt, ex)
        self._engaged = {lid: False for lid in self._loops}
        self._config = _jsonable(config)
        self._ops.append({"op": "configure"})

    def _configure_excitation(self, pt: str, ex: Mapping[str, Any]) -> None:
        from cranberry.blocks.meas.loop_identification import RELATIONS

        if set(ex) != _EXC_KEYS:
            raise ValueError(f"excitations.{pt}: expected keys {sorted(_EXC_KEYS)}, got {sorted(ex)}")
        if pt in self._injections:
            raise ValueError(f"excitations.{pt}: point name already used by an injection")
        self._leaf(ex["block"], f"excitations.{pt}", EXCITATION_PARAMS, ENABLE_TOGGLE)
        if self._leaves[ex["block"]].kind != "excitation":
            raise ValueError(f"excitations.{pt}: block {ex['block']!r} is not an excitation block")
        if ex["loop"] not in self._loops:
            raise KeyError(f"excitations.{pt}: unknown loop {ex['loop']!r}; configured {sorted(self._loops)}")
        if ex["relation"] not in RELATIONS:
            raise ValueError(f"excitations.{pt}: relation must be one of {RELATIONS}")
        if ex["relation"] == "minus_l" and ex["reference"] is None:
            raise ValueError(f"excitations.{pt}: relation minus_l needs a reference slot")
        if ex["relation"] == "plant" and ex["controller"] is None:
            raise ValueError(f"excitations.{pt}: relation plant needs a controller (mapping or 'model')")
        if ex["model"] is not None:
            self._model_params(ex["model"])
        self._excitations[pt] = dict(ex)

    def set_params(self, params: Mapping[str, float]) -> dict[str, float]:
        self._require_configured()
        params = dict(params)
        applied: dict[str, float] = {}
        if "run.seed" in params:
            seed = params.pop("run.seed")
            if int(seed) != seed or seed < 0:
                raise ValueError(f"run.seed must be an int >= 0, got {seed!r}")
            self._seed = int(seed)
            applied["run.seed"] = float(self._seed)
        unknown = [k for k in params if k not in self._runtime]
        if unknown:
            raise KeyError(f"unknown runtime parameter(s) {unknown}; known {sorted(self._runtime)}")
        new = {k: float(v) for k, v in params.items()}
        for lp in self._loops.values():
            block, dev = lp["pid_block"], lp["device"]
            keys = [f"{block}.{k}" for k in CANONICAL_KEYS if f"{block}.{k}" in self._runtime]
            if dev is None or not any(k in new for k in keys):
                continue
            cur = {k.rsplit(".", 1)[1]: new.get(k, self._value(k)) for k in keys}
            cur = {k: v for k, v in cur.items() if not (k.startswith("omega") and v == 0.0)}  # 0 = no saturation
            clipped, msgs = apply_device_limits(PidCoefficients.from_params(cur), dev, self._device(dev))
            self._note(msgs)
            back = clipped.as_params()
            for k in keys:
                if k in new:
                    new[k] = back.get(k.rsplit(".", 1)[1], 0.0)
        self._overrides.update(new)
        applied.update(new)
        self._ops.append({"op": "set_params", "applied": dict(applied)})
        return applied

    def get_params(self) -> dict[str, float]:
        return {name: self._value(name) for name in self._runtime}

    def engage(self, loop_id: str | None = None) -> None:
        self._require_configured()
        for lid in self._loop_ids(loop_id):
            self._engaged[lid] = True
        self._ops.append({"op": "engage"} if loop_id is None else {"op": "engage", "loop": loop_id})

    def disengage(self, loop_id: str | None = None) -> None:
        self._require_configured()
        for lid in self._loop_ids(loop_id):
            self._engaged[lid] = False
        self._ops.append({"op": "disengage"} if loop_id is None else {"op": "disengage", "loop": loop_id})

    def reset_integrators(self) -> None:
        self._require_configured()
        self._ops.append({"op": "reset_integrators"})  # every acquire starts from the initial state

    def ramp(self, channel: str, *, amplitude: float, period_s: float, offset: float, enable: bool) -> None:
        self._require_configured()
        if channel not in self._ramps:
            raise KeyError(f"unknown ramp channel {channel!r}; configured {sorted(self._ramps)}")
        r = self._ramps[channel]
        vals = {"amplitude_v": float(amplitude), "period_s": float(period_s), "offset_v": float(offset)}
        dev = self._device(r["device"])
        if dev is not None and "ramp" in dev:
            limits = {"amplitude_v": "amplitude_range_v", "period_s": "period_range_s", "offset_v": "offset_range_v"}
            for k, lim in limits.items():
                if lim in dev["ramp"]:
                    lo, hi = dev["ramp"][lim]
                    y = min(max(vals[k], lo), hi)
                    if y != vals[k]:
                        msg = f"ramp {channel}.{k} = {vals[k]:.6g} outside [{lo:g}, {hi:g}] ({r['device']}); clipped to {y:.6g}"
                        warnings.warn(msg, DeviceRangeWarning, stacklevel=2)
                        self._note([msg])
                        vals[k] = y
        self._overrides.update({f"{r['block']}.{k}": v for k, v in vals.items()})
        self._ramp_on[channel] = bool(enable)
        self._ops.append({"op": "ramp", "channel": channel, **vals, "enable": bool(enable)})

    def acquire(self, duration_s: float, signals: Sequence[str],
                injections_on: bool | Sequence[str] = False) -> dict[str, np.ndarray]:
        self._require_configured()
        points = {**self._injections, **self._excitations}
        if isinstance(injections_on, bool):
            on = set(points) if injections_on else set()
        else:
            on = set(injections_on)
            unknown = sorted(on - set(points))
            if unknown:
                raise KeyError(f"unknown injection / excitation point(s) {unknown}; configured {sorted(points)}")
        self._inj_on = on
        try:
            res = self._run(self._effective_spec(duration_s))
        finally:
            self._inj_on = set()
        op = {"op": "acquire", "duration_s": float(duration_s), "signals": list(signals)}
        if on:
            op["injections_on"] = sorted(on)
        self._ops.append(op)
        return self._collect(res, signals)

    def _model_params(self, name: str) -> dict[str, Any]:
        from cranberry.runner.run import offline_block_specs

        for b in offline_block_specs(self._spec):
            if b.name == name and b.kind == "loop_stability_check":
                return dict(b.params)
        raise KeyError(f"no loop_stability_check block {name!r} in offline.blocks")

    def measure_loop_gain(self, loop_id: str, method: str, freqs_hz: Sequence[float] | None = None, *,
                          band_hz: Sequence[float] | None = None, amplitude: float, point: str | None = None,
                          **settings: Any) -> dict[str, Any]:
        """One locked run with the excitation of ``loop_id`` on (module docstring)."""
        import math

        from cranberry.blocks.ctrl.excitation import excitation_schedule
        from cranberry.blocks.meas.loop_identification import identify

        self._require_configured()
        if point is None:
            pts = [pt for pt, ex in self._excitations.items() if ex["loop"] == loop_id]
            if len(pts) != 1:
                raise KeyError(f"loop {loop_id!r} has {len(pts)} excitation points {pts}; pass point=")
            point = pts[0]
        if point not in self._excitations:
            raise KeyError(f"unknown excitation point {point!r}; configured {sorted(self._excitations)}")
        ex = self._excitations[point]
        known = {"settle_s", "t_start_s", "dwell_s", "period_s", "n_periods", "prbs_order", "bit_s", "duration_s",
                 "nperseg_s", "n_freqs", "amp_f_corner_hz", "n_blocks", "coherence_min"}
        unknown = sorted(set(settings) - known)
        if unknown:
            raise TypeError(f"measure_loop_gain: unknown setting(s) {unknown}; known {sorted(known)}")
        if "settle_s" not in settings:
            raise TypeError("measure_loop_gain needs settle_s (D34)")
        t_start = float(settings.get("t_start_s", 0.0))
        struct: dict[str, Any] = {"waveform": method, "settle_s": float(settings["settle_s"]),
                                  "amp_f_corner_hz": float(settings.get("amp_f_corner_hz", math.inf))}
        for k in ("dwell_s", "period_s", "prbs_order", "bit_s"):
            if k in settings:
                struct[k] = settings[k]
        if method in ("stepped_sine", "multisine"):
            if freqs_hz is not None:
                struct["freqs_hz"] = [float(f) for f in freqs_hz]
            elif band_hz is not None and "n_freqs" in settings:
                struct.update(freqs_hz=None, f_start_hz=float(band_hz[0]), f_stop_hz=float(band_hz[1]),
                              n_freqs=int(settings["n_freqs"]), spacing="log")
            else:
                raise TypeError(f"measure_loop_gain({method!r}) needs freqs_hz, or band_hz with n_freqs")
        if method == "chirp":
            if band_hz is None:
                raise TypeError("measure_loop_gain('chirp') needs band_hz")
            struct.update(f_start_hz=float(band_hz[0]), f_stop_hz=float(band_hz[1]), spacing="log")
        block = ex["block"]
        sc = excitation_schedule({**self._leaves[block].params, **struct, "t_start_s": t_start,
                                  "amplitude": float(amplitude)}, block)
        if method == "stepped_sine":
            duration = t_start + len(sc["freqs_hz"]) * sc["dwell_s"]
        elif method == "white":
            if "duration_s" not in settings:
                raise TypeError("measure_loop_gain('white') needs duration_s")
            duration = float(settings["duration_s"])
        else:
            if "n_periods" not in settings:
                raise TypeError(f"measure_loop_gain({method!r}) needs n_periods")
            per = sc["period_s"]
            duration = t_start + (math.ceil(sc["settle_s"] / per - 1e-9) + int(settings["n_periods"])) * per
        duration += 2.0 / self._spec.run.f_sim_hz
        x_slot = f"ctrl.{block}.x"
        slots = [s for s in (x_slot, ex["reference"], ex["response"]) if s is not None]
        spec = self._effective_spec(duration, {f"{block}.amplitude": float(amplitude), f"{block}.t_start_s": t_start},
                                    {f"{block}.{ENABLE_TOGGLE}": True}, structural={block: struct}, record_slots=slots)
        res = self._run(spec)
        data = self._collect(res, slots)
        model = None if ex["model"] is None else self._model_params(ex["model"])
        sig = {"x": data[x_slot], "response": data[ex["response"]]}
        if ex["reference"] is not None:
            sig["reference"] = data[ex["reference"]]
        out = identify(sig, data["t_s"], sc, self._spec.run.f_sim_hz, relation=ex["relation"],
                       controller=ex["controller"], model_params=model,
                       model_stages="all" if ex["relation"] == "minus_l" else "first",
                       coherence_min=float(settings.get("coherence_min", 0.0)),
                       n_blocks=int(settings.get("n_blocks", 4)), band_hz=band_hz, nperseg_s=settings.get("nperseg_s"))
        out.update(loop_id=loop_id, point=point, amplitude=float(amplitude), duration_s=duration)
        self._ops.append({"op": "measure_loop_gain", "loop": loop_id, "point": point, "method": method,
                          "amplitude": float(amplitude), "duration_s": duration,
                          "freqs_hz": [float(f) for f in np.atleast_1d(out["f_hz"])]})
        return out

    @staticmethod
    def _collect(res: RunResult, signals: Sequence[str]) -> dict[str, np.ndarray]:
        missing = [s for s in signals if s not in res.series]
        if missing:
            raise KeyError(f"signal(s) {missing} not recorded; recorded {sorted(res.series)}")
        out = {s: np.asarray(res.series[s], dtype=np.float64) for s in signals}
        out["t_s"] = np.asarray(res.t_s, dtype=np.float64)
        return out

    def measure_transfer_function(self, injection_point: str, freqs_hz: np.ndarray) -> np.ndarray:
        self._require_configured()
        if injection_point not in self._injections:
            raise KeyError(f"unknown injection point {injection_point!r}; configured {sorted(self._injections)}")
        inj = self._injections[injection_point]
        block, ref, resp = inj["block"], inj["reference"], inj["response"]
        freqs = np.atleast_1d(np.asarray(freqs_hz, dtype=np.float64))
        if np.any(freqs <= 0):
            raise ValueError("freqs_hz must be > 0")
        out = np.empty(freqs.size, dtype=np.complex128)
        for i, f in enumerate(freqs):
            window = inj["periods"] / f
            spec = self._effective_spec(
                inj["settle_s"] + window,
                {f"{block}.freq_hz": float(f), f"{block}.amplitude": float(inj["amplitude"])},
                {f"{block}.{ENABLE_TOGGLE}": True},
            )
            data = self._collect(self._run(spec), [ref, resp])
            t = data["t_s"]
            sel = (t >= t[0] + inj["settle_s"]) & (t < t[0] + inj["settle_s"] + window)
            x_ref = demodulate(data[ref][sel], t[sel], f)
            if x_ref == 0:
                raise ZeroDivisionError(f"reference {ref!r} has no component at {f:g} Hz")
            out[i] = demodulate(data[resp][sel], t[sel], f) / x_ref
        self._ops.append({"op": "measure_transfer_function", "point": injection_point, "freqs_hz": freqs.tolist()})
        return out

    def status(self) -> dict[str, Any]:
        last = None
        if self._last is not None:
            locks = {k: _jsonable(np.asarray(v)[-1]) for k, v in self._last.series.items() if k.startswith("lock.") and len(v)}
            last = {
                "locks": locks,
                "events": _jsonable(list(self._last.events)),
                "run_dir": None if self._last.run_dir is None else str(self._last.run_dir),
            }
        return {
            "backend": "sim",
            "configured": self._config is not None,
            "engaged": bool(self._engaged) and all(self._engaged.values()),
            "loops": {lid: {"engaged": self._engaged[lid], "device": lp["device"]} for lid, lp in self._loops.items()},
            "ramps": dict(self._ramp_on),
            "last_run": last,
            "warnings": list(self._warnings),
        }

    def log_record(self) -> dict[str, Any]:
        return _jsonable(
            {
                "controller": "SimController",
                "preset": self._spec.preset,
                "seed": self._seed,
                "config": self._config,
                "devices": self._devices,
                "params": dict(self._overrides),
                "toggles": self._toggle_state() if self._config is not None else {},
                "engaged": bool(self._engaged) and all(self._engaged.values()),
                "loops_engaged": dict(self._engaged),
                "ops": self._ops,
                "warnings": self._warnings,
            }
        )
