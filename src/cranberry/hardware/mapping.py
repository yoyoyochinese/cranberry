"""Canonical PID coefficients → device representations (PLAN §2.10, G16-h). Owner: WP-F5.

Canonical form (device-agnostic): K_p [1], K_i [1/s], K_d [s], integrator / differentiator saturation
corners ω_IS, ω_DS [rad/s] (None = none), optional double integrator. Device presets live in
``devices/*.yaml`` with ``source:`` page references. Out-of-range values warn and are clipped; the
clipped values are returned.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

__all__ = ["PidCoefficients", "DEVICE_MAPPERS", "map_canonical", "canonical_to_dsc1", "canonical_to_moku", "canonical_to_fpga"]


@dataclass(frozen=True)
class PidCoefficients:
    kp: float
    ki_per_s: float
    kd_s: float
    omega_is_rad_per_s: float | None
    omega_ds_rad_per_s: float | None
    double_integrator: bool


def canonical_to_dsc1(c: PidCoefficients, device: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ DSC1 entry values (K_i unit switch per ``ki_entry``, sign e = S − v, h). WP-F5."""
    raise NotImplementedError("hardware.mapping.canonical_to_dsc1 is implemented by WP-F5")


def canonical_to_moku(c: PidCoefficients, device: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ Moku representation (gain in dB and crossover frequencies). WP-F5."""
    raise NotImplementedError("hardware.mapping.canonical_to_moku is implemented by WP-F5")


def canonical_to_fpga(c: PidCoefficients, device: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """→ generic FPGA (Tustin) representation. WP-F5."""
    raise NotImplementedError("hardware.mapping.canonical_to_fpga is implemented by WP-F5")


#: device preset name → mapper
DEVICE_MAPPERS = {"dsc1": canonical_to_dsc1, "moku_delta": canonical_to_moku, "fpga_generic": canonical_to_fpga}


def map_canonical(c: PidCoefficients, device_name: str, device: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Dispatch to the mapper of ``device_name``; returns (device values, warnings)."""
    try:
        mapper = DEVICE_MAPPERS[device_name]
    except KeyError:
        raise KeyError(f"no PID mapping for device {device_name!r}; known {sorted(DEVICE_MAPPERS)}") from None
    return mapper(c, device)
