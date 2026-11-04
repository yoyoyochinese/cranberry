"""``cal_slot_check``: where the calibration gates of a slow in-band layer sit in the GBS frame, and what they cost
(offline CHK, G16-j, D34: warns, never stops; Phase 6 fix of reviews/phase6_wp6c_p.md C-3 / C-7; lit-pittaluga2021
§3.3 conditions 1 and 3). Owner: WP6-C fix.

Frame: K ports × B = N_b + n_guard slots, T_f = K·B·T_rep. Placements (``placement``):

* ``cal_frames`` (no extra slot): a fraction f = ``cal_fraction`` (< 1) of the frames carries no squeezed data; in such a
  frame the N_b slots of port k's block and of port 0's block hold the seeds (they meet at the GBS input), so N_b gates
  per calibration frame. L_k unchanged; the GBS keeps a fraction 1 − f of the frames. Recovery as
  ``forward_ref_stray_table``'s calibration frames: the data frame after a calibration frame (weight min(1, f/(1 − f)) of
  the data frames) starts (1 + n_guard)·T_rep after the last seed of the worst (last) port.
* ``guard`` (one slot per block, every frame with f = 1): the seed sits in a guard slot, which must exist
  (n_guard ≥ 1, overlay demux/guard_slot or lock/points/dbs_cal_guard: raises otherwise, as ``forward_ref_stray_table``
  guard_slot). L_k and the D78 linewidth requirement Δν ≤ σ²/(2πRΔτ), Δτ = (K − 1)·B·T_rep, change by B/(B − 1) against
  the frame without that slot (CHK ``dbs_cal_guard_slot``); the GBS keeps N_b/B of the slots. The demux edge must then
  fit before the seed (the guard no longer buys switching time). The seed is one slot from the next data slot (weight f).
* ``ideal_free`` (BAD, the pre-review assumption kept for comparison, D34): a free slot next to every block that costs
  neither a slot nor L_k (CHK ``dbs_cal_slot_ideal``).

SNSPD recovery into the next data slot (Δt after the seed gate, efficiency r(Δt) = 0 for Δt < ``dead_time_s``, else
1 − exp(−(Δt − dead)/``recovery_s``), as ``photon_counting_fringe``): the seed gate clicks with p = 1 − e^{−λ},
λ = qe·μ/2 + ``bg_per_gate`` (the fringe side, X = π/2), so that data gate loses p·(1 − r) of its efficiency;
frame-weighted mean > ``recovery_loss_max`` → CHK ``dbs_recovery``. Gates per counting window
G = (gates per frame)·f·T_int/T_f < ``min_gates`` → CHK ``dbs_few_counts`` (small-count regime: clamps, all-click windows,
the estimator leaves its asymptotic MLE / Cramér-Rao behaviour).

Params (all required, D34): ``placement``, ``N_b``, ``n_guard``, ``K``, ``t_rep_s``, ``cal_fraction``, ``t_int_s``,
``mu_per_gate``, ``qe``, ``bg_per_gate`` (known background per output: dark + squeezed vacuum b'), ``dead_time_s``,
``recovery_s``, ``sigma_target_rad``, ``lock_ratio_R``, ``min_gates``, ``recovery_loss_max``. Toggle ``enabled``.
:attr:`result`: gates_per_frame, f_gate_hz (gates per second before ``cal_fraction``, the ``photon_counting_fringe``
``f_gate_hz``), gates_per_window, data_fraction, lk_factor, dtau_s, dnu_max_hz, dnu_max_without_slot_hz, p_click,
recovery_gap_s, eta_next_slot, recovery_weight, recovery_loss_mean.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.checks import ChkLog
from cranberry.core.registry import register_block

__all__ = ["CalSlotCheck", "make_cal_slot_check", "cal_slot", "PLACEMENTS"]

KIND = "cal_slot_check"
PLACEMENTS = ("cal_frames", "guard", "ideal_free")
_KEYS = ("placement", "N_b", "n_guard", "K", "t_rep_s", "cal_fraction", "t_int_s", "mu_per_gate", "qe", "bg_per_gate",
         "dead_time_s", "recovery_s", "sigma_target_rad", "lock_ratio_R", "min_gates", "recovery_loss_max")


def cal_slot(p: Mapping[str, Any]) -> dict[str, float]:
    """The numbers of :attr:`CalSlotCheck.result` (module docstring)."""
    nb, ng, K, tr, f = int(p["N_b"]), int(p["n_guard"]), int(p["K"]), float(p["t_rep_s"]), float(p["cal_fraction"])
    b = nb + ng
    tf = K * b * tr
    place = p["placement"]
    per_frame = nb if place == "cal_frames" else 1
    base_b = b - 1 if place == "guard" else b
    dtau = (K - 1) * b * tr
    sig2r = float(p["sigma_target_rad"]) ** 2 / (2 * math.pi * float(p["lock_ratio_R"]))
    data = (1.0 - f) if place == "cal_frames" else (nb / b if place == "guard" else 1.0)
    gap = (1 + ng) * tr if place == "cal_frames" else tr
    weight = min(1.0, f / (1.0 - f)) if place == "cal_frames" else f
    dead, rec = float(p["dead_time_s"]), float(p["recovery_s"])
    r_next = 0.0 if gap < dead else -math.expm1(-(gap - dead) / rec)
    p_click = -math.expm1(-(0.5 * float(p["qe"]) * float(p["mu_per_gate"]) + float(p["bg_per_gate"])))
    loss = p_click * (1.0 - r_next)
    return {"gates_per_frame": float(per_frame), "f_gate_hz": per_frame / tf,
            "gates_per_window": per_frame * f * float(p["t_int_s"]) / tf, "data_fraction": data,
            "lk_factor": b / base_b, "dtau_s": dtau, "dnu_max_hz": sig2r / dtau,
            "dnu_max_without_slot_hz": sig2r / ((K - 1) * base_b * tr), "p_click": p_click, "recovery_gap_s": gap,
            "eta_next_slot": 1.0 - loss, "recovery_weight": weight, "recovery_loss_mean": weight * loss}


class CalSlotCheck(BaseBlock):
    """The ``cal_slot_check`` block (module docstring)."""

    toggle_names = ("enabled",)

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, *_KEYS, block=name)
        if params["placement"] not in PLACEMENTS:
            raise ValueError(f"block {name!r}: placement must be one of {PLACEMENTS}, got {params['placement']!r}")
        if params["placement"] == "guard" and int(params["n_guard"]) < 1:
            raise ValueError(f"block {name!r}: placement guard needs a guard slot, n_guard >= 1 (apply demux/guard_slot "
                             "or lock/points/dbs_cal_guard)")
        if not 0.0 < float(params["cal_fraction"]) <= 1.0:
            raise ValueError(f"block {name!r}: need 0 < cal_fraction <= 1")
        if params["placement"] == "cal_frames" and float(params["cal_fraction"]) >= 1.0:
            raise ValueError(f"block {name!r}: cal_frames with cal_fraction 1 leaves no data frame")
        if int(params["N_b"]) < 1 or int(params["K"]) < 2:
            raise ValueError(f"block {name!r}: need N_b >= 1 and K >= 2")
        super().__init__(KIND, name, params)
        self.chk_log = ChkLog()
        self.result: dict[str, float] | None = None
        if self._toggles["enabled"]:
            self.result = self._check()

    def __call__(self, state, params):
        return state  # offline block: not a time-loop stage

    def _check(self) -> dict[str, float]:
        p, log = self.params, self.chk_log
        r = cal_slot(p)
        if p["placement"] == "ideal_free":
            log.warn("dbs_cal_slot_ideal", self.name, "calibration gate in a free slot next to the block that costs no "
                     "slot and no L_k: the frame has none (n_guard = 0 means every slot is data); use cal_frames or guard",
                     **r)
        if p["placement"] == "guard":
            log.warn("dbs_cal_guard_slot", self.name,
                     f"seed in the guard slot: L_k × {r['lk_factor']:.3g}, GBS keeps {r['data_fraction']:.3g} of the "
                     f"slots; D78 linewidth requirement at σ = {float(p['sigma_target_rad']):.3g} rad: "
                     f"{r['dnu_max_hz']:.3g} Hz instead of {r['dnu_max_without_slot_hz']:.3g} Hz", **r)
        if r["recovery_loss_mean"] > float(p["recovery_loss_max"]):
            log.warn("dbs_recovery", self.name,
                     f"SNSPD recovery: the data slot {r['recovery_gap_s']:.3g} s after a seed gate (click p = "
                     f"{r['p_click']:.3g}) keeps η × {r['eta_next_slot']:.4f}; frame-weighted loss "
                     f"{r['recovery_loss_mean']:.3g} > {float(p['recovery_loss_max']):.3g}", **r)
        if r["gates_per_window"] < float(p["min_gates"]):
            log.warn("dbs_few_counts", self.name,
                     f"{r['gates_per_window']:.3g} calibration gates per counting window < {float(p['min_gates']):g}: "
                     "small-count regime (clamped / all-click windows; the estimator is not at its Cramér-Rao bound)",
                     **r)
        return r


@register_block(KIND)
def make_cal_slot_check(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``cal_slot_check`` (module docstring)."""
    return CalSlotCheck(params, name)
