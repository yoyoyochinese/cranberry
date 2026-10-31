"""FIB offline tables from the classical fiber model (WP-B9, T06 §10.1, PLAN §4.4). Cached via ``TableCache``.

* ``fiber_linear_table`` — per length L: [η = e^{−αL}, τ_g = n_g L/c, GDD = β2 L, TOD = β3 L]
  (toggles ``loss``, ``beta2``, ``beta3``; off → η = 1 / 0).
* ``lo_nonlinear_table`` — RK4IP of a LO / probe pulse train per average power, saved along z:
  table ``lo_nonlinear`` (axes ``p_avg_w`` × ``z_m``) with columns :data:`LO_COLUMNS` and table
  ``gnlse_output_spectrum`` (axes ``p_avg_w`` × ``z_m`` × ``omega_rad_s`` ascending) = |Ã(z, Ω)|².
  ov = ⟨A_lin|A_NL⟩/(‖A_lin‖‖A_NL‖) with A_lin the linear-only propagation of the same input (the
  "linear mode", T06 §8.4). ∂(arg ov)/∂ln P (LO RIN → BHD phase, ``spm_rin_to_phase``) is exact
  from the tangent-linear RK4IP with δA(0) = A(0)/2: ∂ arg ov/∂ln P = Im(⟨A_lin|δA_NL⟩/⟨A_lin|A_NL⟩).

Pulse energy E = P_avg/f_rep; sech: T0 = τ_FWHM/(2 acosh √2), P0 = E/(2T0); gaussian
(|A|² ∝ e^{−T²/T0²}): T0 = τ_FWHM/(2√ln 2), P0 = E/(√π T0).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from cranberry.blocks.fib.gnlse import (
    C_LIGHT,
    TOGGLES,
    _check_fiber,
    grid_edge_monitor,
    model_from_params,
    propagate,
    pulse_field,
    time_axis,
)
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.fft import fftshift
from cranberry.core.registry import register_block
from cranberry.offline.tables import Table, TableCache, code_version, spec_hash

__all__ = ["LO_COLUMNS", "LINEAR_COLUMNS", "FiberLinearTable", "LoNonlinearTable", "make_fiber_linear_table",
           "make_lo_nonlinear_table", "pulse_t0_p0"]

LINEAR_COLUMNS = ("eta", "tau_g_s", "gdd_s2", "tod_s3")
LO_COLUMNS = ("arg_ov_rad", "ov_abs2", "darg_ov_dlnp_rad", "phi_nl_peak_rad", "centroid_shift_rad_s",
              "rms_width_ratio")
_REV = 1  # bump when the numerics change (cache invalidation, G16-g)


def pulse_t0_p0(shape: str, t_fwhm_s: float, energy_j: float) -> tuple[float, float]:
    """(T0, P0) of a ``sech`` / ``gaussian`` pulse from intensity FWHM and energy (module docstring)."""
    if shape == "sech":
        t0 = t_fwhm_s / (2 * np.arccosh(np.sqrt(2.0)))
        return t0, energy_j / (2 * t0)
    if shape == "gaussian":
        t0 = t_fwhm_s / (2 * np.sqrt(np.log(2.0)))
        return t0, energy_j / (np.sqrt(np.pi) * t0)
    raise ValueError(f"pulse shape must be 'sech' or 'gaussian', got {shape!r}")


def _ascending(name: str, key: str, xs) -> np.ndarray:
    a = np.asarray(xs, dtype=float)
    if a.ndim != 1 or a.size < 1 or (a.size > 1 and not np.all(np.diff(a) > 0)):
        raise ValueError(f"block {name!r}: {key} must be a non-empty strictly ascending list")
    return a


def _cached(kind: str, name: str, params: Mapping[str, Any], cache: TableCache | None, names, build) -> dict[str, Table]:
    """Each table of ``names`` from the cache, else from ``build()`` (one call builds all of them)."""
    cache = cache or TableCache()
    ver = code_version(f"fib.{kind}", _REV)
    memo: dict[str, Table] = {}

    def build_all():
        if not memo:
            memo.update(build())
        return memo

    out = {}
    for tname in names:
        key = spec_hash({"kind": kind, "table": tname, "params": dict(params)}, ver)
        out[tname] = cache.get_or_build(key, lambda t=tname, k=key: _with_prov(build_all()[t], k, ver, name))
    return out


def _with_prov(t: Table, key: str, ver: str, name: str) -> Table:
    return Table(t.name, t.grid, t.values, t.unit, {**t.provenance, "spec_hash": key, "code_version": ver, "block": name})


class FiberLinearTable(BaseBlock):
    """``fiber_linear_table`` (module docstring). Params (required, D34): ``toggles`` {loss, beta2, beta3},
    ``lengths_m`` (ascending), ``alpha_per_m``, ``n_g``, ``beta2_s2_per_m``, ``beta3_s3_per_m``."""

    toggle_names = ("loss", "beta2", "beta3")
    kind_name = "fiber_linear_table"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.kind_name
        require(params, "lengths_m", "alpha_per_m", "n_g", "beta2_s2_per_m", "beta3_s3_per_m", block=name)
        _ascending(name, "lengths_m", params["lengths_m"])
        super().__init__(self.kind_name, name, params)

    def __call__(self, state, params):
        return state

    def tables(self, cache: TableCache | None = None) -> dict[str, Table]:
        """{``fiber_linear``: Table(axis ``L_m``, columns :data:`LINEAR_COLUMNS`)}."""
        p, tg = self.params, self._toggles

        def build():
            L = np.asarray(p["lengths_m"], dtype=float)
            vals = np.stack([np.exp(-float(p["alpha_per_m"]) * L) if tg["loss"] else np.ones_like(L),
                             float(p["n_g"]) * L / C_LIGHT,
                             float(p["beta2_s2_per_m"]) * L * tg["beta2"],
                             float(p["beta3_s3_per_m"]) * L * tg["beta3"]], axis=-1)
            return {"fiber_linear": Table("fiber_linear", {"L_m": L}, vals, "mixed",
                                          {"columns": list(LINEAR_COLUMNS), "sources": ["T06 §10.1"]})}

        return _cached(self.kind_name, self.name, p, cache, ("fiber_linear",), build)


class LoNonlinearTable(BaseBlock):
    """``lo_nonlinear_table`` (module docstring). Params (required, D34): ``toggles`` (``fiber_gnlse``'s
    without ``cw_shortcut``), the ``fiber_gnlse`` fiber keys, ``pulse: {shape, t_fwhm_s}``, ``f_rep_hz``,
    ``p_avg_w`` (ascending), ``grid: {n_t, t_window_s}``, ``steps: {n_steps, n_save, control[, rtol]}``."""

    toggle_names = tuple(t for t in TOGGLES if t != "cw_shortcut")
    kind_name = "lo_nonlinear_table"

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.kind_name
        _check_fiber(params, name)
        require(params, "pulse", "f_rep_hz", "p_avg_w", "grid", "steps", block=name)
        require(params["pulse"], "shape", "t_fwhm_s", block=f"{name}.pulse")
        require(params["grid"], "n_t", "t_window_s", block=f"{name}.grid")
        require(params["steps"], "n_steps", "n_save", "control", block=f"{name}.steps")
        _ascending(name, "p_avg_w", params["p_avg_w"])
        super().__init__(self.kind_name, name, params)
        self.chk_log = ChkLog()

    def __call__(self, state, params):
        return state

    def _build(self) -> dict[str, Table]:
        p, g, s = self.params, self.params["grid"], self.params["steps"]
        t, dt = time_axis(g["n_t"], g["t_window_s"])
        model = model_from_params(p, {**self._toggles, "cw_shortcut": False}, int(g["n_t"]), dt)
        shift = model.n_t // 2
        p_avg = np.asarray(p["p_avg_w"], dtype=float)
        rows, spectra = [], []
        for pw in p_avg:
            t0, p0 = pulse_t0_p0(p["pulse"]["shape"], float(p["pulse"]["t_fwhm_s"]), pw / float(p["f_rep_hz"]))
            at0 = model.to_w(np.roll(pulse_field(p["pulse"]["shape"], p0, t0, t), -shift).astype(complex))
            pr = propagate(model, at0, float(p["L_m"]), n_steps=int(s["n_steps"]), n_save=int(s["n_save"]),
                           control=s["control"], rtol=s.get("rtol"), dat0=at0 / 2)
            grid_edge_monitor(model, pr.at_z[-1], self.name, self.chk_log)
            w0, c0 = model.rms_width(at0), model.centroid(at0)
            row = []
            for z, at, dat, phi in zip(pr.z, pr.at_z, pr.dat_z, pr.phi_nl_peak_rad, strict=True):
                lin = model.linear(at0, z)
                ov = np.vdot(lin, at)
                ovn = ov / np.sqrt(np.vdot(lin, lin).real * np.vdot(at, at).real)
                row.append([np.angle(ov), abs(ovn) ** 2, (np.vdot(lin, dat) / ov).imag, phi,
                            model.centroid(at) - c0, model.rms_width(at) / w0])
            rows.append(row)
            spectra.append(fftshift(np.abs(pr.at_z) ** 2, axes=-1))
        z = np.asarray(pr.z)
        prov = {"columns": list(LO_COLUMNS), "sources": ["T06 §8.4", "T06 §10.1"]}
        return {
            "lo_nonlinear": Table("lo_nonlinear", {"p_avg_w": p_avg, "z_m": z}, np.array(rows), "mixed", prov),
            "gnlse_output_spectrum": Table("gnlse_output_spectrum",
                                           {"p_avg_w": p_avg, "z_m": z, "omega_rad_s": fftshift(model.omega)},
                                           np.array(spectra), "W*s^2", {"sources": ["T06 §10.4"]}),
        }

    def tables(self, cache: TableCache | None = None) -> dict[str, Table]:
        """{``lo_nonlinear``, ``gnlse_output_spectrum``} (module docstring), cached by spec hash."""
        return _cached(self.kind_name, self.name, self.params, cache, ("lo_nonlinear", "gnlse_output_spectrum"),
                       self._build)


@register_block("fiber_linear_table")
def make_fiber_linear_table(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``fiber_linear_table`` (:class:`FiberLinearTable`)."""
    return FiberLinearTable(params, name)


@register_block("lo_nonlinear_table")
def make_lo_nonlinear_table(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``lo_nonlinear_table`` (:class:`LoNonlinearTable`)."""
    return LoNonlinearTable(params, name)
