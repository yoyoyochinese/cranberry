# Phase 6 prompt: Exploration, including fusion candidates from the literature

Start only after D81-7 (the user trusts the simulator). Read `AGENTS.md`, `docs/squeezing_sim_prompt.md` §8,
`DECISIONS.md` and `PLAN.md` first. All rules there still apply: protected list, planning interrogation,
`PROCESS_LOG.md` logging, notes in Japanese (D19), never fabricate citations (mark UNVERIFIED).

## Goal

1. Run the §8 exploration on the future presets (`F_wg100`, `F_wg6`, `F_bor6`, `F_liu20`): ablation, combination
   of the effects that matter, BoTorch tuning per fiber length.
2. Settle the decisions left for Phase 6: LOCK_TOPOLOGY (D14), the LO delivery lock (D12), the phase budget
   allocation σ_t (D49), and whether 32 simultaneous ports are needed (D13).
3. Evaluate whether the seven schemes below can be built into or fused with our future system
   (single squeezer → demux → per-port PM fiber delays → GBS; Liu-type counter-propagating probe as the baseline).

## Candidate schemes (PDFs in `refs/`, file names start with first author and year)

| ID | Paper | Idea to import | Fusion with our system / what to test |
|---|---|---|---|
| P1 | Aghaee Rad et al. 2025, Nature 638, 912 (Aurora) | Classical reference pulses interleaved between quantum pulses in phase- and polarization-stabilized fiber delay modules; TEC (slow) + piezo fiber phase shifter (fast) | Same-path reference pulses in the demux guard slots of each port, as a true-path sensor or as slow calibration of the Liu proxy. Compare metric (ii) and stray clicks per gate with the Liu baseline |
| P2 | Madsen et al. 2022, Nature 606, 75 (Borealis; `refs/xanadu_2022.pdf`) | Counter-propagating lock beam with alternating measurement and stabilization windows | Liu probe switched off during measurement windows, controller held (sample-and-hold). Sweep window period and duty ratio against drift between windows, aliasing (Dick effect) and stray light |
| P3 | Nardelli et al. 2026, Optica Quantum 4, 138 | CW reference 8 nm away from the quantum band, choppers anti-synchronized with the detector gate, digital low-pass below the chopper rate, path-length actuators | WDM co-propagating probe (`probe/sep_wdm` with an offset overlay) plus time gating. Find the gate-rate window between residual drift (too slow) and gating noise (too fast) |
| P4 | Švarc et al. 2023, Opt. Express 31, 12562 | Different-wavelength co-propagating reference, set point adapted to the summed outputs, 0.05° over 15 h at 8.5 m | (a) Bias of a different-wavelength proxy from the dispersion of the thermo-optic coefficient (see below): check that the model contains it, then size it per fiber length and temperature drift. (b) Intensity-normalized set point for the DC-fringe lock (Ref D, `notes/dc-fringe-lock-ref-d.md`) under RIN and power drift |
| P5 | Pittaluga et al. 2021, Nat. Photon. 15, 530 | Dual-band stabilization: bright off-band reference for the fast loop, dim in-band reference counted by SNSPDs for the residual | Two-layer lock: fast proxy (Liu or WDM probe) + slow in-band correction from SNSPD clicks of a weak seed in calibration slots. Trade photon budget and calibration rate against proxy bias removal |
| P6 | Zhong et al. 2021, PRL 127, 180502 (Jiuzhang 2.0) | 1450 nm reference co-propagating with the squeezed light in GBS; phases programmed by electrical delays on the reference signal | Program GBS port phases through the lock set point (reference RF phase) instead of extra modulators on the squeezed path. Test set-point step response, settling time and cross-port consistency |
| P7 | Kashiwazaki et al. 2023, APL, DOI 10.1063/5.0144385 | Open-loop response measured while a slow loop holds the lock; demodulation frequency chosen from the phase margin; EOM (fast) + fiber stretcher (slow) cascade | Emulate the in-lock identification through the abstract controller interface so the same procedure runs on Moku:Delta. Use the identified model to set BoTorch bounds and the actuator split |

For P4 (a), with a path-length actuator holding the reference phase, the signal phase drifts by

$$
\delta\phi_s = \frac{2\pi L\,n_s}{\lambda_s}\left(\frac{\beta_s}{n_s}-\frac{\beta_r}{n_r}\right)\delta T,
\qquad \beta \equiv \frac{\partial n}{\partial T},
$$

where $L$ is the fiber length, $n_{s,r}$ the indices and $\lambda_{s,r}$ the wavelengths of signal and reference, and
$\delta T$ the temperature change. A phase-only actuator (e.g. an AOM frequency) instead leaves
$\left(1-\lambda_s/\lambda_r\right)\delta\phi_s^{\mathrm{free}}$, where $\delta\phi_s^{\mathrm{free}}$ is the free-running
signal phase change. Verify both against the code.

## Procedure

1. **Check `refs/`.** If a PDF is missing, stop and tell me which one. For each paper write a short
   `notes/lit-<citekey>.md`: what the paper does (with figure and page numbers), parameters we can transfer,
   the mapping to blocks and overlays, and items marked UNVERIFIED. Confirm from the PDFs the Methods numbers of
   P1 and P2 (window lengths, pulse lengths, actuator models, phase noise per loop); the survey took them from
   secondary extraction.
2. **Baseline (§8).** Ablation on the four nominal future presets with the existing lock presets
   (`presets/locks/`), sweeping fiber length; 1-port and 2-port pair configurations.
3. **Express P1–P7 as configurations.** Use existing blocks, overlays and presets first; keep existing presets
   unchanged (D34). If a scheme needs a new block, list all such blocks for P1–P7 first, run one grill-me
   round (coordinating agent only), record the decisions, and wait for my approval before implementing.
   A new block gets a fresh-context Phase 5 review before its results are used.
4. **Compare and combine.** For each scheme alone and for promising fusions (e.g. P2 gating on the Liu probe;
   P3 + P5 on a WDM probe; P1 as the slow layer under the Liu probe), report: metric (i), metric (ii) with its
   four components (G13), stray clicks per gate and photons per mode, η_k, lock duration and slips, actuator range
   use, per fiber length. Use 10–30 seeds; claim only differences larger than the spread.
5. **Optimize.** BoTorch per fiber length (constraints: lock maintained, servo bump within limits). Variables
   include gain, filter cutoffs, actuator split, gate period and duty ratio, demodulation frequency.
6. **Report** in `notes/findings.md`. Every insight needs a physical explanation and a check against step size,
   chunk size and seed. End with:
   - a recommendation for D14, D12 and D49, with the evidence;
   - a ranking of P1–P7 by gain in metric (ii), stray light, hardware cost and complexity, and Moku:Delta feasibility;
   - bench experiments that would test the top two schemes on our hardware.

Stop and ask me when a needed parameter is missing (placeholder rule) or when a design choice would be
expensive to reverse. Do not modify the frozen Ref D prediction (`predictions/refD/`).
