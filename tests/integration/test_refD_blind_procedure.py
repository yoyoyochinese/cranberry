"""Ref D blind comparison (D04, PLAN §6.5): PLACEHOLDER — skipped until the end of Phase 5.

Nothing here reads, references or encodes `blind/` or any measured Ref D value. The test documents the
procedure and the inputs it will need; the coordinator fills it in at the end of Phase 5.

Procedure (PLAN §6.5):

1. Phase 4–5: Ref D is checked only analytically (VT-13 test_vt_13_closed_loop_tf.py, VT-14
   test_vt_14_refD_free_running.py, VT-21 test_vt_21_refD_thermal_ramp.py, VT-24
   test_vt_24_refD_analytic.py). Workers never open `blind/`.
2. Static quantities the user may measure before the end of Phase 5 (D36) — ΔL from fringe counting,
   fibre lengths, insertion losses — may replace the corresponding axes in ``presets/reference/refD.yaml``
   (``offline.params.dL_m``, ``L_arm2_m``, ``il_915b_db``, ``coupler_excess_db``) before the freeze.
   Everything else stays a sweep axis: 915B resonance (``f_res_hz``, ``q_res``), 914 input impedance
   (``r_914_ohm``), DSC1 latency (``latency_s``, ``in_out_poles``), anti-windup, ``h0_hz2_per_hz``,
   ``h_m1_hz2``, ``D_pol_rad2_per_s``, ``K_T_rad_per_m_K``, environmental models (D27), loop / coupler sign.
3. Freeze: when every Phase-5 VT passes, the coordinator fixes ``refD.yaml`` + the prediction config
   and tags ``refD-blind-freeze``. The prediction is the band over the axes (``runner.sweep``; seeds per
   point, ``run.seed_mode: per_point``) of: true φ_D(t) and its PSD, e_DSC1 = S − v PSD (v = ``v.pdb``),
   residual rms, lock duration (``lock.dsc1.sup.t_loss``), lock-point offset (mean of
   ``lock.lockD.err``). Result-file hashes go to PROCESS_LOG.md.
4. Only the coordinator opens the blind targets and compares them with the frozen prediction. No
   parameter is tuned.
5. Deviations are attributed by ablation (``runner.ablate``, one toggle at a time) and by the position
   inside the sweep ranges; a model fix is a new freeze and a separate comparison (the first is kept).

Inputs the filled-in test will need (all provided by the coordinator, none exist yet):
    * the frozen preset path / tag and the prediction run directory (``runs/<id>/``, hashes logged);
    * the measured quantities of step 3, in the same units / slot names, with their uncertainty;
    * the comparison metric per quantity (inside / outside the predicted band).
"""

import pytest

pytestmark = pytest.mark.skip(reason="Ref D blind comparison: end of Phase 5, coordinator only (D04, PLAN §6.5)")


def test_refD_blind_comparison_placeholder():
    raise AssertionError("not to be run before the refD-blind-freeze tag (PLAN §6.5)")
