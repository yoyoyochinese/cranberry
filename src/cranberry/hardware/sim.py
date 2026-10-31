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

``configure`` config (all names are block / slot names of the spec)::

    loops:      {<loop id>: {pid_block: <name>, device: <devices/*.yaml name> | null}}
    ramps:      {<channel>: {block: <name>, device: <name> | null}}                     # optional
    injections: {<point>: {block: <name>, reference: <slot>, response: <slot>,
                           amplitude: <float>, settle_s: <float>, periods: <int>}}     # optional

With a device, PID gains and ramp settings are clipped exactly as that device would clip them
(``hardware.mapping.apply_device_limits``, device ``ramp`` ranges), so simulation and hardware see the
same applied values.
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

__all__ = ["SimController", "demodulate", "LOOP_TOGGLE", "RAMP_PARAMS", "INJECTION_PARAMS", "ENABLE_TOGGLE"]

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
        self._leaves = {leaf.name: leaf for leaf in iter_leaves(from_spec(spec.pipeline))}
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
        self._devices: dict[str, dict[str, Any]] = {}
        self._overrides: dict[str, float] = {}
        self._ramp_on: dict[str, bool] = {}
        self._engaged = False
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
        t = {f"{lp['pid_block']}.{LOOP_TOGGLE}": self._engaged for lp in self._loops.values()}
        t.update({f"{r['block']}.{ENABLE_TOGGLE}": self._ramp_on[ch] for ch, r in self._ramps.items()})
        t.update({f"{i['block']}.{ENABLE_TOGGLE}": False for i in self._injections.values()})
        return t

    def _effective_spec(self, duration_s: float, extra_params: Mapping[str, float] = {}, extra_toggles: Mapping[str, bool] = {}) -> RunSpec:
        per_block: dict[str, dict[str, Any]] = {}
        for full, v in {**self._overrides, **extra_params}.items():
            block, key = self._runtime[full]
            per_block.setdefault(block, {})[key] = v
        run = dataclasses.replace(self._spec.run, duration_s=float(duration_s), seed=self._seed)
        toggles = {**self._spec.toggles, **self._toggle_state(), **extra_toggles}
        return dataclasses.replace(self._spec, run=run, pipeline=_replace_params(self._spec.pipeline, per_block), toggles=toggles)

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
        unknown = set(config) - {"loops", "ramps", "injections"}
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
        self._config = _jsonable(config)
        self._ops.append({"op": "configure"})

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

    def engage(self) -> None:
        self._require_configured()
        self._engaged = True
        self._ops.append({"op": "engage"})

    def disengage(self) -> None:
        self._require_configured()
        self._engaged = False
        self._ops.append({"op": "disengage"})

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

    def acquire(self, duration_s: float, signals: Sequence[str]) -> dict[str, np.ndarray]:
        self._require_configured()
        res = self._run(self._effective_spec(duration_s))
        self._ops.append({"op": "acquire", "duration_s": float(duration_s), "signals": list(signals)})
        return self._collect(res, signals)

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
            "engaged": self._engaged,
            "loops": {lid: {"engaged": self._engaged, "device": lp["device"]} for lid, lp in self._loops.items()},
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
                "engaged": self._engaged,
                "ops": self._ops,
                "warnings": self._warnings,
            }
        )
