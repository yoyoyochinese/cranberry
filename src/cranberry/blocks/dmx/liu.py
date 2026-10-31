"""Liu-type lock specifics not covered by B3 (detector) / F6 (path model): WP-B8 (T10 §3.2, T11 §10.1, D30, D40).

* Actuator position (T10 §3.2 "ループ遅延"): the probe enters at the GBS end (x = L) and is detected
  at the upstream end (x = 0). An actuator at x_a acts on the probe n_g x_a/c before detection, so the
  loop delay is τ_loop + n_g x_a/c: ``upstream`` (x_a = 0) adds nothing, ``downstream`` (x_a = L) adds
  the probe transit n_g L/c. The residual of distributed noise, eq. (6), is the same at both ends.
  The F6 path model produces this delay by itself (actuator subsegment ``actuator_sub``); use
  :func:`liu_loop_delay_s` for the ``delay_buffer`` ``tau_s`` (D55-1: tau_s = loop delay − 1/f_sim).
* Detection before / after the demux (D40, T11 §10.1 ``shared_with_probe``) is a topology choice: a
  probe beam path that includes the demux segment makes the slow switch phase part of the locked path
  (cancelled, computed by F6); otherwise ψ_sw stays in metric (ii) component 4 (squeezed-only).
  After the demux the probe passes only while port k's switch is open: :func:`probe_in_switch_window`.
"""

from __future__ import annotations

import math

from cranberry.path.topology import C_M_PER_S

__all__ = ["ACTUATOR_POSITIONS", "actuator_x_m", "actuator_sub", "liu_loop_delay_s", "probe_in_switch_window"]

ACTUATOR_POSITIONS = ("upstream", "downstream")


def actuator_x_m(position: str | float, L_m: float) -> float:  # noqa: N803
    """x_a [m] from the upstream (probe-detection) end: ``upstream`` 0, ``downstream`` L, or a number."""
    if position == "upstream":
        return 0.0
    if position == "downstream":
        return float(L_m)
    x = float(position)
    if not 0.0 <= x <= L_m:
        raise ValueError(f"actuator position must be in [0, L] or one of {ACTUATOR_POSITIONS}, got {position!r}")
    return x


def actuator_sub(position: str | float, L_m: float, n_sub: int) -> int:  # noqa: N803
    """Subsegment index of the actuator on an F6 segment of ``n_sub`` subsegments (0 = upstream end)."""
    return min(n_sub - 1, int(math.floor(actuator_x_m(position, L_m) / L_m * n_sub)))


def liu_loop_delay_s(tau_loop_s: float, position: str | float, L_m: float, n_g: float) -> float:  # noqa: N803
    """Total loop delay τ_loop + n_g x_a/c (module docstring)."""
    return float(tau_loop_s) + n_g * actuator_x_m(position, L_m) / C_M_PER_S


def probe_in_switch_window(t_arrive_s: float, tau_pr_s: float, window_start_s: float, window_s: float,
                           period_s: float) -> bool:
    """True if the probe pulse [t_arrive, t_arrive + τ_pr) at the switch lies inside port k's open window
    [window_start, window_start + window) of every frame (all times modulo ``period_s``)."""
    a = (t_arrive_s - window_start_s) % period_s
    return a + tau_pr_s <= window_s
