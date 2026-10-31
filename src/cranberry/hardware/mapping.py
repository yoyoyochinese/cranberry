"""Canonical PID coefficients ↔ device-native settings (protected list; PLAN §2.10, T16 §4.3, §4.6, G16-h).

Canonical form (device-agnostic, continuous time, independent terms; T16 §4.3)::

    C(s) = K_p + K_i / (s + ω_IS) + K_i2 / s²  + K_d s / (1 + s/ω_DS)
    u    = C(s) · e,      e = v − S   (input voltage minus setpoint voltage)

with K_p [1], K_i [1/s], K_i2 [1/s²] (only when ``double_integrator``), K_d [s], ω_IS, ω_DS [rad/s]
(None or 0 = no saturation). ``e = v − S`` is conv §9 ("quantity − setpoint") applied to the controller
input voltage. A device that defines its error the other way (DSC1: e = S − v, DSC1 manual p.18)
gets the sign flipped inside its mapper (conv §9); the device's ``pid.error_sign`` says which.

Mappers return ``(settings, warnings)``. Values a device cannot represent are clipped / rounded to the
nearest representable setting, a :class:`DeviceRangeWarning` is emitted and its message is also
returned; the returned settings are what the device would actually run. A combination that no clipping
can fix (mixed signs on a sign-locked device) raises ``ValueError``.

Device-native settings
----------------------

``dsc1`` (DSC1 manual Rev. B; printed pages)
    ``{"p", "i", "d": signed entries (0 when disabled), "p_enabled", "i_enabled", "d_enabled": bool}``.
    |entry| ∈ [1e-5, 1e5] (p.6), signs locked together (p.19). Entry units per p.20 (D53, user-verified):
    ``ki_entry: per_sample`` → I = K_i·h, ``kd_entry: per_sample`` → D = K_d/h, h = 2 µs. No integrator
    saturation, derivative filter or double integrator (p.19): those canonical terms are dropped with a
    warning. Sign: p = −K_p, i = −K_i·h, d = −K_d/h (e = S − v).

    D53 check (Ref D settings p = 0.000241, i = 0.001406, d = 0.01108) ⇒ canonical
    **K_p = −2.41e-4, K_i = −703 s⁻¹, K_d = −2.216e-8 s** (no saturation, no I+).

``moku_delta`` (Moku:Delta Specifications v26-0817, PID Controller p.65; T16 §4.3)
    ``{"sign": ±1, "p_enabled", "kp_db", "i_enabled", "f_i_hz", "i2_enabled", "f_i2_hz", "is_enabled",
    "is_level_db", "d_enabled", "f_d_hz", "ds_enabled", "ds_level_db"}`` (disabled values None).
    Crossover parameterisation with reference gain g = |K_p| (1 when P is off): K_i = g·2πf_I,
    K_d = g/(2πf_D), K_i2 = g·2πf_I·2πf_I2, IS level = K_i/ω_IS, DS level = K_d·ω_DS (UNVERIFIED: the
    vendor transfer function is not in the sheet; T16 §4.3). Ranges from p.65.

``fpga_generic`` (T16 §4.3: Tustin)
    ``{"fs_hz", "discretisation": "tustin", "kp", "i", "i2", "d"}`` where each of ``i``/``i2``/``d`` is
    ``{"b": [...], "a": [...]}`` (difference-equation coefficients in z⁻¹, a[0] = 1) or None. Tustin
    without prewarping: s → 2 f_s (1 − z⁻¹)/(1 + z⁻¹).
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, fields
from typing import Any

from cranberry.core.block import MissingParameterError

__all__ = [
    "PidCoefficients",
    "CANONICAL_KEYS",
    "DeviceRangeWarning",
    "DEVICE_MAPPERS",
    "DEVICE_INVERSE",
    "map_canonical",
    "map_to_canonical",
    "apply_device_limits",
    "canonical_to_dsc1",
    "dsc1_to_canonical",
    "canonical_to_moku",
    "moku_to_canonical",
    "canonical_to_fpga",
    "fpga_to_canonical",
]

TWO_PI = 2.0 * math.pi


class DeviceRangeWarning(UserWarning):
    """A canonical value was clipped / rounded / dropped to fit a device."""


@dataclass(frozen=True)
class PidCoefficients:
    """Canonical PID coefficients (module docstring). ``ki2_per_s2`` is required iff ``double_integrator``."""

    kp: float
    ki_per_s: float
    kd_s: float
    omega_is_rad_per_s: float | None
    omega_ds_rad_per_s: float | None
    double_integrator: bool
    ki2_per_s2: float | None = None

    def __post_init__(self) -> None:
        for name in ("omega_is_rad_per_s", "omega_ds_rad_per_s"):
            w = getattr(self, name)
            if w is not None and not w >= 0:
                raise ValueError(f"{name} must be >= 0 or None, got {w!r}")
        if self.double_integrator and self.ki2_per_s2 is None:
            raise ValueError("double_integrator=True needs ki2_per_s2")
        if not self.double_integrator and self.ki2_per_s2 not in (None, 0.0):
            raise ValueError("ki2_per_s2 given but double_integrator=False")

    def transfer(self, f_hz: Any) -> Any:
        """C(s) at s = i·2πf (Laplace convention, s = +iω; scalar or numpy array)."""
        s = 1j * TWO_PI * f_hz
        c = self.kp + self.ki_per_s / (s + (self.omega_is_rad_per_s or 0.0))
        if self.double_integrator:
            c = c + self.ki2_per_s2 / s**2
        d = s if not self.omega_ds_rad_per_s else s / (1.0 + s / self.omega_ds_rad_per_s)
        return c + self.kd_s * d

    def as_params(self) -> dict[str, float]:
        """Runtime-parameter form ({key: float}, None entries omitted); keys ⊂ :data:`CANONICAL_KEYS`."""
        return {k: float(getattr(self, k)) for k in CANONICAL_KEYS if getattr(self, k) is not None}

    @classmethod
    def from_params(cls, params: Mapping[str, float]) -> PidCoefficients:
        """Inverse of :meth:`as_params`: missing saturation keys → None, ``ki2_per_s2`` present → I+."""
        unknown = set(params) - set(CANONICAL_KEYS)
        if unknown:
            raise KeyError(f"unknown PID parameter(s) {sorted(unknown)}; known {CANONICAL_KEYS}")
        ki2 = params.get("ki2_per_s2")
        return cls(
            float(params["kp"]), float(params["ki_per_s"]), float(params["kd_s"]),
            params.get("omega_is_rad_per_s"), params.get("omega_ds_rad_per_s"),
            ki2 is not None, None if ki2 is None else float(ki2),
        )


#: Runtime-parameter keys of a canonical PID (the ``pid`` block's runtime params, PLAN §2.10).
CANONICAL_KEYS: tuple[str, ...] = tuple(f.name for f in fields(PidCoefficients) if f.name != "double_integrator")


# ---------------------------------------------------------------- helpers


def _warn(msgs: list[str], msg: str) -> None:
    msgs.append(msg)
    warnings.warn(msg, DeviceRangeWarning, stacklevel=3)


def _clip(msgs: list[str], name: str, x: float, lo: float, hi: float) -> float:
    y = min(max(x, lo), hi)
    if y != x:
        _warn(msgs, f"{name} = {x:.6g} outside [{lo:.6g}, {hi:.6g}]; clipped to {y:.6g}")
    return y


def _error_sign(device: Mapping[str, Any]) -> float:
    """+1 if the device error is v − S (canonical), −1 if S − v (conv §9)."""
    s = str(device["pid"]["error_sign"]).replace(" ", "")
    if s == "v-S":
        return 1.0
    if s == "S-v":
        return -1.0
    raise ValueError(f"unknown pid.error_sign {device['pid']['error_sign']!r}; expected 'v - S' or 'S - v'")


def _common_sign(values: Iterable[float]) -> float:
    signs = {math.copysign(1.0, v) for v in values if v != 0}
    if len(signs) > 1:
        raise ValueError("device has sign-locked gains (all nonzero terms must share one sign); got mixed signs")
    return signs.pop() if signs else 1.0


def _sat(w: float | None) -> float | None:
    return w if w else None


def _require(value: Any, what: str) -> Any:
    if value is None:
        raise MissingParameterError(f"{what} is null in the device preset; set it in the preset (D34)")
    return value


def _db(x: float) -> float:
    return 20.0 * math.log10(x)


def _undb(x: float) -> float:
    return 10.0 ** (x / 20.0)


# ---------------------------------------------------------------- DSC1


def _dsc1_scales(pid: Mapping[str, Any]) -> tuple[float, float]:
    """(entered I / K_i, entered D / K_d) per ``ki_entry`` / ``kd_entry`` (DSC1 p.20, D53)."""
    h = float(pid["h_s"])
    ki = {"per_sample": h, "per_second": 1.0}
    kd = {"per_sample": 1.0 / h, "seconds": 1.0}
    try:
        return ki[pid["ki_entry"]], kd[pid["kd_entry"]]
    except KeyError:
        raise ValueError(f"unknown ki_entry/kd_entry {pid['ki_entry']!r}/{pid['kd_entry']!r}") from None


def canonical_to_dsc1(c: PidCoefficients, device: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ DSC1 entries (module docstring; DSC1 p.6, p.18–20; D53)."""
    pid = device["pid"]
    msgs: list[str] = []
    if _sat(c.omega_is_rad_per_s):
        _warn(msgs, "DSC1 has no integrator saturation (p.18-19); omega_is dropped (pure integrator)")
    if _sat(c.omega_ds_rad_per_s):
        _warn(msgs, "DSC1 derivative has no filter (p.19); omega_ds dropped")
    if c.double_integrator and c.ki2_per_s2:
        _warn(msgs, "DSC1 has no double integrator (p.19); ki2 dropped")
    sigma = _error_sign(device)
    si, sd = _dsc1_scales(pid)
    raw = {"p": sigma * c.kp, "i": sigma * c.ki_per_s * si, "d": sigma * c.kd_s * sd}
    if pid["gain_sign_locked"]:
        _common_sign(raw.values())
    lo, hi = pid["gain_abs_range"]
    out: dict[str, Any] = {}
    for term, x in raw.items():
        a = abs(x)
        if a == 0:
            y = 0.0
        elif a < lo:
            y = 0.0 if a < lo / 2 else lo  # nearest representable of {0 (disabled)} ∪ [lo, hi]
            _warn(msgs, f"DSC1 {term} = {x:.6g} below |{lo:g}|; set to {math.copysign(y, x):.6g}")
        else:
            y = _clip(msgs, f"DSC1 {term}", a, lo, hi)
        out[term] = math.copysign(y, x) if y else 0.0
        out[f"{term}_enabled"] = bool(y)
    return out, msgs


def dsc1_to_canonical(settings: Mapping[str, Any], device: Mapping[str, Any]) -> tuple[PidCoefficients, list[str]]:
    """DSC1 entries → canonical. A disabled term counts as 0 whatever its entry (p.18)."""
    pid = device["pid"]
    sigma = _error_sign(device)
    si, sd = _dsc1_scales(pid)
    x = {t: float(settings[t]) if settings[f"{t}_enabled"] else 0.0 for t in ("p", "i", "d")}
    return PidCoefficients(sigma * x["p"], sigma * x["i"] / si, sigma * x["d"] / sd, None, None, False), []


# ---------------------------------------------------------------- Moku:Delta


def canonical_to_moku(c: PidCoefficients, device: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ Moku:Delta PID Controller settings (gain in dB, crossover frequencies; Moku p.65, T16 §4.3)."""
    pid = device["pid"]
    msgs: list[str] = []
    sigma = _error_sign(device)
    kp, ki, kd = sigma * c.kp, sigma * c.ki_per_s, sigma * c.kd_s
    ki2 = sigma * (c.ki2_per_s2 or 0.0)
    sign = _common_sign((kp, ki, kd, ki2))
    out: dict[str, Any] = {"sign": int(sign)}

    if kp:
        ref = _clip(msgs, "Moku |kp|", abs(kp), *pid["kp_abs_range"])
        out.update(p_enabled=True, kp_db=_db(ref))
    else:
        ref = 1.0
        out.update(p_enabled=False, kp_db=None)

    out.update(i_enabled=False, f_i_hz=None, i2_enabled=False, f_i2_hz=None, is_enabled=False, is_level_db=None)
    if ki:
        f_i = _clip(msgs, "Moku f_i_hz", abs(ki) / (TWO_PI * ref), *pid["f_i_hz_range"])
        ki_abs = ref * TWO_PI * f_i
        out.update(i_enabled=True, f_i_hz=f_i)
        if ki2:
            f_i2 = _clip(msgs, "Moku f_i2_hz", abs(ki2) / (ki_abs * TWO_PI), pid["f_i2_hz_min"], f_i)
            out.update(i2_enabled=True, f_i2_hz=f_i2)
        w = _sat(c.omega_is_rad_per_s)
        if w:
            hi = min(pid["is_level_max"], ki_abs / (TWO_PI * pid["f_is_hz_min"]))
            level = _clip(msgs, "Moku IS level", ki_abs / w, ref, max(hi, ref))
            out.update(is_enabled=True, is_level_db=_db(level))
    else:
        if ki2:
            _warn(msgs, "Moku double integrator needs the integrator (f_i2 <= f_i, p.65); ki2 dropped")
        if _sat(c.omega_is_rad_per_s):
            _warn(msgs, "integrator saturation without integrator; omega_is dropped")

    out.update(d_enabled=False, f_d_hz=None, ds_enabled=False, ds_level_db=None)
    if kd:
        f_d = _clip(msgs, "Moku f_d_hz", ref / (TWO_PI * abs(kd)), *pid["f_d_hz_range"])
        kd_abs = ref / (TWO_PI * f_d)
        out.update(d_enabled=True, f_d_hz=f_d)
        w = _sat(c.omega_ds_rad_per_s)
        if w:
            hi = min(pid["ds_level_max"], kd_abs * TWO_PI * pid["f_ds_hz_max"])
            level = _clip(msgs, "Moku DS level", kd_abs * w, ref, max(hi, ref))
            out.update(ds_enabled=True, ds_level_db=_db(level))
    elif _sat(c.omega_ds_rad_per_s):
        _warn(msgs, "differentiator saturation without differentiator; omega_ds dropped")
    return out, msgs


def moku_to_canonical(settings: Mapping[str, Any], device: Mapping[str, Any]) -> tuple[PidCoefficients, list[str]]:
    """Moku:Delta settings → canonical (inverse of :func:`canonical_to_moku`)."""
    sigma = _error_sign(device) * float(settings["sign"])
    ref = _undb(settings["kp_db"]) if settings["p_enabled"] else 1.0
    kp = ref if settings["p_enabled"] else 0.0
    ki = ki2 = w_is = kd = w_ds = 0.0
    if settings["i_enabled"]:
        wi = TWO_PI * settings["f_i_hz"]
        ki = ref * wi
        if settings["i2_enabled"]:
            ki2 = ki * TWO_PI * settings["f_i2_hz"]
        if settings["is_enabled"]:
            w_is = ki / _undb(settings["is_level_db"])
    if settings["d_enabled"]:
        kd = ref / (TWO_PI * settings["f_d_hz"])
        if settings["ds_enabled"]:
            w_ds = _undb(settings["ds_level_db"]) / kd
    i2 = bool(settings["i_enabled"] and settings["i2_enabled"])
    return PidCoefficients(sigma * kp, sigma * ki, sigma * kd, _sat(w_is), _sat(w_ds), i2, sigma * ki2 if i2 else None), []


# ---------------------------------------------------------------- generic FPGA (Tustin)


def _tustin_first_order(b1: float, b0: float, a0: float, k: float) -> dict[str, list[float]]:
    """(b1 s + b0)/(s + a0) with s → k(1 − z⁻¹)/(1 + z⁻¹)."""
    den = k + a0
    return {"b": [(b1 * k + b0) / den, (b0 - b1 * k) / den], "a": [1.0, (a0 - k) / den]}


def canonical_to_fpga(c: PidCoefficients, device: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ generic FPGA Tustin difference equations (module docstring; T16 §4.3)."""
    fs = float(_require(device["pid"]["update_hz"], "fpga pid.update_hz"))
    k = 2.0 * fs
    msgs: list[str] = []
    sigma = _error_sign(device)
    w_is, w_ds = c.omega_is_rad_per_s or 0.0, c.omega_ds_rad_per_s or 0.0
    for name, w in (("omega_is", w_is), ("omega_ds", w_ds)):
        if w > math.pi * fs:
            _warn(msgs, f"{name} = {w:.6g} rad/s above Nyquist (pi*fs); Tustin warps it strongly")
    out: dict[str, Any] = {"fs_hz": fs, "discretisation": "tustin", "kp": sigma * c.kp, "i": None, "i2": None, "d": None}
    if c.ki_per_s:
        out["i"] = _tustin_first_order(0.0, sigma * c.ki_per_s, w_is, k)
    if c.double_integrator and c.ki2_per_s2:
        g = sigma * c.ki2_per_s2 / k**2
        out["i2"] = {"b": [g, 2.0 * g, g], "a": [1.0, -2.0, 1.0]}
    if c.kd_s:
        if w_ds:
            out["d"] = _tustin_first_order(sigma * c.kd_s * w_ds, 0.0, w_ds, k)
        else:
            _warn(msgs, "pure derivative under Tustin has a pole at z = -1 (Nyquist ringing); set omega_ds")
            g = sigma * c.kd_s * k
            out["d"] = {"b": [g, -g], "a": [1.0, 1.0]}
    return out, msgs


def fpga_to_canonical(settings: Mapping[str, Any], device: Mapping[str, Any]) -> tuple[PidCoefficients, list[str]]:
    """Tustin difference equations → canonical (exact inverse of :func:`canonical_to_fpga`)."""
    k = 2.0 * float(settings["fs_hz"])
    sigma = _error_sign(device)
    ki = w_is = kd = w_ds = 0.0
    ki2 = None
    if settings.get("i"):
        (b0, _), (_, a1) = settings["i"]["b"], settings["i"]["a"]
        w_is = k * (1.0 + a1) / (1.0 - a1)
        ki = b0 * (k + w_is)
    if settings.get("i2"):
        ki2 = sigma * settings["i2"]["b"][0] * k**2
    if settings.get("d"):
        (b0, _), (_, a1) = settings["d"]["b"], settings["d"]["a"]
        if a1 == 1.0:
            kd = b0 / k
        else:
            w_ds = k * (1.0 + a1) / (1.0 - a1)
            kd = b0 * (k + w_ds) / (w_ds * k)
    return PidCoefficients(sigma * settings["kp"], sigma * ki, sigma * kd, _sat(w_is), _sat(w_ds), ki2 is not None, ki2), []


# ---------------------------------------------------------------- dispatch

Mapper = Callable[[PidCoefficients, Mapping[str, Any]], tuple[dict[str, Any], list[str]]]
Inverse = Callable[[Mapping[str, Any], Mapping[str, Any]], tuple[PidCoefficients, list[str]]]

#: device preset name → canonical → device mapper
DEVICE_MAPPERS: dict[str, Mapper] = {"dsc1": canonical_to_dsc1, "moku_delta": canonical_to_moku, "fpga_generic": canonical_to_fpga}
#: device preset name → device → canonical mapper
DEVICE_INVERSE: dict[str, Inverse] = {"dsc1": dsc1_to_canonical, "moku_delta": moku_to_canonical, "fpga_generic": fpga_to_canonical}


def _device(device_name: str, device: Mapping[str, Any] | None) -> Mapping[str, Any]:
    if device is not None:
        return device
    from cranberry.hardware.devices import device_values, load_device

    return device_values(load_device(device_name))


def map_canonical(c: PidCoefficients, device_name: str, device: Mapping[str, Any] | None = None) -> tuple[dict[str, Any], list[str]]:
    """Dispatch to the mapper of ``device_name``; ``device`` = plain values (``devices.device_values``),
    None → load ``devices/<device_name>.yaml``. Returns (device settings, warnings)."""
    try:
        mapper = DEVICE_MAPPERS[device_name]
    except KeyError:
        raise KeyError(f"no PID mapping for device {device_name!r}; known {sorted(DEVICE_MAPPERS)}") from None
    return mapper(c, _device(device_name, device))


def map_to_canonical(settings: Mapping[str, Any], device_name: str, device: Mapping[str, Any] | None = None) -> tuple[PidCoefficients, list[str]]:
    """Device settings → canonical (inverse dispatch)."""
    try:
        inverse = DEVICE_INVERSE[device_name]
    except KeyError:
        raise KeyError(f"no PID mapping for device {device_name!r}; known {sorted(DEVICE_INVERSE)}") from None
    return inverse(settings, _device(device_name, device))


def apply_device_limits(c: PidCoefficients, device_name: str, device: Mapping[str, Any] | None = None) -> tuple[PidCoefficients, list[str]]:
    """Canonical → device → canonical: the coefficients the device would actually run (used by
    ``SimController.set_params`` so simulation and hardware clip identically)."""
    dev = _device(device_name, device)
    settings, msgs = map_canonical(c, device_name, dev)
    back, more = map_to_canonical(settings, device_name, dev)
    return back, msgs + more
