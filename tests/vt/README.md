# Validation tests (VT) — checklist

Validation tests are on the protected list (AGENTS.md): never delete or weaken them. Source of truth for
the list, tolerances and targets: `PLAN.md` §6.2–§6.3 and `notes/00_summary.md` §4 (with the D48 / D51
readings of prompt §7). This file maps each VT ID to its future test file and owner so progress can be
ticked off.

Rules (PLAN §6.1, §6.3):

- One file per VT: `tests/vt/test_vt_XX_<slug>.py`. The module docstring cites the source test numbers
  (S §4, T-notes) and the tolerance.
- Files under `tests/vt/` are auto-marked `vt` (see `tests/conftest.py`); add `@pytest.mark.slow` for
  anything over ~60 s and `@pytest.mark.numba` when compiled kernels are exercised.
- Stochastic checks use fixed seeds (`core.rng.stream`) and a sample size whose 3σ statistical error is
  below the tolerance.
- Never read or reference `blind/` (D04).

| Done | VT | prompt §7 | Main modules | Kind | WP | M | File |
|---|---|---|---|---|---|---|---|
| [x] | VT-01 | 1 | `single_mode_reference`, `phase_to_quadrature` | unit | B1, F4 | M1 | `test_vt_01_single_mode_reference.py` |
| [x] | VT-02 | 1 | `refA` full (loss budget, lock to 9 mrad; ±0.10 dB, D48) | integ | I1 | M1 | `test_vt_02_refA_full.py` |
| [~] | VT-03 | 1 | CW limit (pulsed blocks identity), GNLSE CW shortcut | unit + integ | B9, B11 | M5 | `test_vt_03_cw_limit.py` |
| [x] | VT-04 | 2 | Fig. 5 (4 curves × 8 points, P_opt) | unit + integ | B1, I1 | M1 | `test_vt_04_refA_fig5.py` |
| [~] | VT-05 | 3 | `offline.gaussian`, `LossRegistry`, `pure_loss`, circuit noise = loss, OPO / pulsed OPO parts | unit | F3 (M0), B5, B6, B11 | M0–M5 | `test_vt_05_pure_loss.py` |
| [~] | VT-06 | 4 | linearized fiber, OPA g → 0, OPO x = 0, pulsed OPO (initial vacuum term) | unit | B10, B1, B6, B11 | M1–M5 | `test_vt_06_vacuum_limits.py` |
| [x] | VT-07 | 5 | `fiber_gnlse` soliton | unit | B9 | M5 | `test_vt_07_gnlse_soliton.py` |
| [x] | VT-08 | 5 | `fiber_gnlse` recurrence | unit | B9 | M5 | `test_vt_08_gnlse_recurrence.py` |
| [x] | VT-09 | 5 | `fiber_gnlse` step halving (error 1/16) | unit | B9 | M5 | `test_vt_09_gnlse_step_halving.py` |
| [x] | VT-10 | 5 | `fiber_gnlse` photon number with Raman | unit | B9 | M5 | `test_vt_10_gnlse_photon_number.py` |
| [x] | VT-11 | 6 | `raman_noise_table` (thermal part → 0 only, D48) | unit | B10 | M5 | `test_vt_11_raman_noise.py` |
| [x] | VT-12 | 6 | `gawbs_table` | unit | B10 | M5 | `test_vt_12_gawbs.py` |
| [~] | VT-13 | 7 | generic linear loop (M0, `SimController.measure_transfer_function`), Ref D z-domain (M1) | integ | F1, F5, B4, I2 | M0, M1 | `test_vt_13_closed_loop_tf.py` (Ref D part done, WP-I2) |
| [x] | VT-14 | 7 | feedback off → free-running (Ref D) | integ | I2, F6 | M1 | `test_vt_14_refD_free_running.py` |
| [x] | VT-15 | 7 | integrator loop residual, lock-2 reduction (A−/A+) | integ | B4, I1 | M1 | `test_vt_15_integrator_residual.py` |
| [x] | VT-16 | 8 | tap trade-off, both definitions (a)(b) (D51-4) | integ | I1 | M1 | `test_vt_16_tap_tradeoff.py` |
| [~] | VT-17 | 9 | Ref B Fig. 3 (x = √(P/P_th), Ω_B = 2πf/γ_B, T02 V-B3 table) | unit + integ | B6, I4 | M3 | `test_vt_17_refB_fig3.py` |
| [~] | VT-18 | 9 | Ref B definition discrimination | unit + integ | B6, I4 | M3 | `test_vt_18_refB_definitions.py` |
| [x] | VT-19 | 10 | dither contribution | unit + integ | F4, B3, I3 | M2 | `test_vt_19_dither_contribution.py`; Ref C preset `test_vt_19_refC_dither.py` (WP-I3) |
| [x] | VT-20 | 10 | total ≈ 0.1 rad, insensitive to controller model (10 %) | integ | I3 | M2 | `test_vt_20_refC_total_residual.py` |
| [~] | VT-21 | 10 | lock duration set by slow-actuator range | integ | B4, I2, I3, I5 | M1, M2, M6 | `test_vt_21_lock_duration.py`; Ref D thermal ramp `test_vt_21_refD_thermal_ramp.py` (WP-I2); Ref C stack `test_vt_21_refC_lock_duration.py` (WP-I3) |
| [~] | VT-22 | 11 | shared-segment cancellation (equal coefficients), passes with polarization off | integ | F6 (M0), I4 (M3), B7 | M0, M3, M4 | `test_vt_22_shared_segment_cancellation.py` |
| [~] | VT-23 | 11 (added) | no cancellation on OPO segment with unequal coefficients | integ | B6, I4, B11 | M3, M5 | `test_vt_23_unequal_coef_no_cancellation.py` |
| [x] | VT-24 | D04 | Ref D pre-blind analytic checks | unit + integ | I2, B3, B4, B2, B7 | M1 | `test_vt_24_refD_analytic.py` |
| [x] | VT-30 | — | OPA | unit | B1 | M1 | `test_vt_30_opa.py` |
| [~] | VT-31 | — | CW OPO | unit | B6 | M2 | `test_vt_31_cw_opo.py` |
| [x] | VT-32 | — | pulsed OPA (table interpolation 1e−3) | unit | B11 | M5 | `test_vt_32_pulsed_opa.py` |
| [~] | VT-33 | — | phase → squeezing, ESA, metric invariance / decomposition (V11, V12), HF add-on (V13) | unit + integ | F4, B5 | M1 (V1–V10, V13), M4 (V11, V12) | `test_vt_33_phase_to_squeezing.py` |
| [~] | VT-34 | — | laser (delay, coherence, fractional delay, Liu static compensation, filter cavity, ASE, CMRR) | unit + integ | B2, F6 | M1, M4 | `test_vt_34_laser.py` |
| [x] | VT-35 | — | GNLSE extras | unit | B9 | M5 | `test_vt_35_gnlse_extra.py` |
| [x] | VT-36 | — | fiber quantum extras | unit | B10 | M5 | `test_vt_36_fiber_quantum_extra.py` |
| [~] | VT-37 | — | environmental noise (coefficients and signs, PSD, correlation, counter-propagation, axis coefficients) | unit | B2, F6 | M1, M2, M4 | `test_vt_37_environmental_noise.py` |
| [x] | VT-38 | — | polarization | unit | B7 | M4 | `test_vt_38_polarization.py` |
| [~] | VT-39 | — | probe cancellation breakdown | integ | B8, F6, B7 | M4 | `test_vt_39_probe_cancellation.py` |
| [x] | VT-40 | — | demux | unit + integ | B8 | M4 | `test_vt_40_demux.py` |
| [~] | VT-41 | — | losses (Ref B / C budgets, 32-port budget, GBS vs BHD separation) | unit | F3, I3, I4, I5 | M2–M6 | `test_vt_41_loss_budgets.py` (Ref C part, WP-I3) |
| [x] | VT-42 | — | detection | unit | B5 | M1 | `test_vt_42_detection.py` |
| [x] | VT-43 | — | error signals (incl. baseband ↔ RF-rate equivalence) | unit | B3 | M1–M4 | `test_vt_43_error_signals.py` |
| [~] | VT-44 | — | actuators, delay buffer, stack = parallel, compiled = pure Python (1e−12), f_sim convergence | unit | F1, B4 | M0, M1 | `test_vt_44_delay_stack_compile.py` |
| [x] | VT-45 | — | pulsed OPO | unit | B11 | M5 | `test_vt_45_pulsed_opo.py` |
| [ ] | VT-F1 | — | YAML ↔ Python round trip (`spec == load(dump(spec))`, bit-identical results at equal seed) | integ | F2 | M0 | `test_vt_f1_yaml_round_trip.py` |
| [ ] | VT-F2 | — | chunk-length invariance (bit-identical, G11), f_sim doubling (§6.4) | integ | F1, B2 | M0, each M | `test_vt_f2_chunk_invariance.py` |
| [x] | VT-F3 | — | double registration in BandLedger / LossRegistry overlap is an error | unit | F3 | M0 | `test_vt_f3_double_registration.py` |
| [ ] | VT-F4 | — | run record (config, seed, commit, env) complete and reproduces the result | integ | F2 | M0 | `test_vt_f4_run_record.py` |

Already covered by WP-F0 unit tests (not VTs, framework properties of PLAN §6.1): conv §8 FFT
(`tests/unit/core/test_fft.py`), RNG stream reproducibility / independence (`test_rng.py`), numba smoke
(`test_numba_smoke.py`), slot layout and pipeline composition basics (`test_layout_state.py`,
`test_block_pipeline.py`).
