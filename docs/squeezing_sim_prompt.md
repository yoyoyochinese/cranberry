# Project Prompt: Modular Simulation of Phase-Stabilized Squeezed Light (OPA → Fiber → BHD)

> Placeholders are written as `{{NAME}}`. Do not invent values for them.
> If a placeholder is still unfilled when you need it, use the paper's value listed next to it (marked "ref:")
> and record that substitution in `PROCESS_LOG.md`. If no reference value exists, stop and ask.

---

## 0. Role and goal

You are part of a small team of coding/research agents building a **device-agnostic, modular simulation**
of an experiment that generates squeezed light with an OPA, propagates it through optical fiber, and
measures it by balanced homodyne detection (BHD) under active phase stabilization.

The simulation's purposes are:

1. Predict the measured squeezing / anti-squeezing level (dB relative to shot noise) as time series and
   as frequency spectra, for a given set of enabled physical effects and control configurations.
2. Derive the required phase-locking performance (residual phase jitter) for each fiber length to deliver
   a target squeezing level to a downstream GBS circuit.
3. Serve as a test bed where agents later explore combinations of noise sources, filters, actuators and
   controller parameters to find good configurations and physically interesting insights.

The final deliverable also includes a **reusable English skill** describing the workflow used to build it
(see §10).

---

## 1. Reference experiment

Primary reference (read it fully; it is in `refs/`):
K. Hirota, T. Kashiwazaki et al., "Generation of 10-dB squeezed light from a broadband waveguide optical
parametric amplifier with improved phase locking method" (PPLN waveguide OPA, 10.1 dB squeezing).

Key elements of that setup the simulation must be able to represent:

- CW seed laser (1545.32 nm) + EDFA + SHG (772.66 nm pump).
- Fundamental split into LO path and probe path.
- Probe path: two AOMs (net frequency shift 1 MHz) + EOM (actuator for lock 1).
- LO path: EOM (actuator for lock 2) + filter cavity to suppress ASE.
- Probe and pump combined, then split into two branches:
  - Squeezer OPA (PPLN waveguide, 45 mm) → squeezed light → BHD.
  - Phase-detection OPA (separate PPLN waveguide) → PD → error signal for **lock 1** (pump–probe phase).
- **Lock 2** (LO–probe phase): error signal from the BHD output (1 MHz beat between residual probe and LO).
- Measurement: ESA in zero-span mode, center 3 MHz, RBW 1 MHz, VBW 100 Hz; also spectra up to 100 MHz.
- Conventional alternative topology: tap the squeezed light after the squeezer OPA for lock 1.
  Both topologies must be expressible as different pipeline configurations.

### Additional reference topologies (must also be expressible as configurations)

**Reference B — co-propagating probe, OPO source** (in `refs/`):
Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with
improvement of phase stability in homodyne measurement", Opt. Express 15, 4321 (2007).
- Source: sub-threshold bow-tie OPO with PPKTP (output coupler T = 0.123, round trip ~500 mm,
  P_th = 180 mW), cavity spectrum R± = 1 ± ηξ²ζρ·4x/((1∓x)² + 4Ω²).
- Pump-dependent intracavity loss L(x) = 0.00249 + 0.00222x → block parameters may be functions of
  other parameters.
- Probe **co-propagates** with the squeezed light through the OPO (same path), phase-modulated at
  295 kHz by a PZT (dither).
- Lock 1 (pump–probe): DC intensity of the transmitted probe on a PD, locked at max/min parametric gain
  via PZT1 (mechanical resonance ~14 kHz).
- Lock 2 (LO–probe): BHD output demodulated → PZT3 (resonance ~22 kHz).
- Laser excess-noise peaks at 89 kHz and 177 kHz (reason for the 295 kHz modulation frequency).
- Note for the research agent: check Eq. (2) against the paper's own numbers (G = 18.7, ρ, ζ, ξ, η)
  for both squeezing and anti-squeezing, and document any definitional differences (e.g. of G or x)
  you find before using it as a validation target.

**Reference C — quantum noise locking (no probe), two fibers** (in `refs/`):
S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized
GHz-bandwidth squeezed vacuum states of light between two buildings", Optica Quantum 4, 1 (2026).
- Source: monolithic ppKTP OPO (linewidth fitted ~1.75 GHz), squeezed vacuum without coherent
  displacement.
- Squeezed light and LO each travel through **separate 1 km fibers** (same conduit → partially
  correlated environmental phase noise).
- Error signal: 1 MHz phase modulation on the **LO** by a fiber EOM (Vπ = 4 V, 200 mVpp → ±4.5°);
  BHD high-pass output → gain → envelope detector (log amplifier, 8–300 MHz band) → demodulation at
  1 MHz. No probe beam.
- Three **stacked** PI(D) controllers, each driving one actuator with a different band: fiber EOM
  (kHz), piezo mirror (below), fiber stretcher (Hz to sub-Hz); plus a ramp generator for acquisition.
- Lock duration limited by the slow actuator's range (minutes).
- The dither itself contributes to residual phase noise (total fitted ~0.1 rad).

**Differences in our experiment** (the simulation must support both the reference and ours):

- Our laser is **pulsed**: {{LASER_MODE}} (CW / pulsed), {{PULSE_WIDTH}}, {{REP_RATE}}, {{AVG_POWER}}.
- Our squeezed light (and possibly LO/probe) propagates through **fibers of several lengths**:
  {{FIBER_LENGTHS}}. This breaks the reference paper's assumption that pump and probe see a nearly common
  path between the split point and the two OPAs, and adds fiber-induced phase noise to the squeezed path.
- Our hardware: Moku:Delta (as one device preset only), EOM(s), high-QE balanced detector, probe PD.
- Our lock topology: {{LOCK_TOPOLOGY}} (ref A phase-detection OPA / conventional tap / ref B
  co-propagating probe / ref C quantum noise lock / other).

---

## 2. Parameters (fill in; reference values given where available)

### Source and OPA
- Wavelengths: {{LAMBDA_FUND}} (ref: 1545.32 nm), {{LAMBDA_PUMP}} (ref: 772.66 nm)
- Laser mode: {{LASER_MODE}} (ref: CW); pulse width {{PULSE_WIDTH}}; repetition rate {{REP_RATE}}
- Laser phase/intensity noise spectra: {{LASER_NOISE_DATA}} (file path or model)
- Source type: {{SOURCE_TYPE}} (single-pass waveguide OPA / OPO cavity; ref A: waveguide OPA)
- OPO parameters (if OPO): output coupler T {{OPO_T}}, intracavity loss {{OPO_LOSS}} (constant or
  function of pump), round-trip length {{OPO_ROUNDTRIP}}, threshold {{OPO_PTH}}, linewidth {{OPO_LINEWIDTH}}
- OPA: type {{OPA_TYPE}} (ref: PPLN waveguide), length {{OPA_LENGTH}} (ref: 45 mm),
  propagation loss {{OPA_PROP_LOSS}} (ref: 0.1 dB/cm), SHG efficiency α {{OPA_ALPHA}} (ref: 906 %/W),
  pump power {{PUMP_POWER}} (ref: 640 mW)
- Phase-detection OPA (if used): {{PD_OPA_PARAMS}}; probe/pump powers into it: {{PD_OPA_POWERS}} (ref: probe ~mW)

### Fiber (per length)
- Lengths: {{FIBER_LENGTHS}}; fiber type: {{FIBER_TYPE}}
- β2 {{BETA2}}, β3 {{BETA3}}, γ {{GAMMA_NL}}, loss {{FIBER_LOSS}}, Raman fraction f_R {{F_RAMAN}},
  cladding diameter {{CLAD_DIAM}} (for GAWBS), temperature {{TEMPERATURE}}
- Which beams travel in fiber: {{BEAMS_IN_FIBER}} (squeezed / LO / probe / pump)
- Environmental (thermal/acoustic) phase-noise model or measured PSD: {{FIBER_PHASE_NOISE_DATA}}

### Detection
- Detector QE {{DET_QE}} (ref: ~99%), bandwidth {{DET_BW}} (ref: 35 MHz),
  shot-noise-to-circuit-noise clearance {{CIRCUIT_CLEARANCE}} (ref: 28 dB at 3 MHz, 15.4 mW LO)
- LO power {{LO_POWER}} (ref: 15.4 mW); LO filter cavity linewidth {{FILTER_CAVITY_LW}}; ASE level {{ASE_LEVEL}}
- Mode-matching efficiency (squeezed–LO) {{MODE_MATCH}} (ref: ~96%, i.e. ~4% loss)
- Other optical loss {{OTHER_LOSS}} (ref: 1% optics + 1% PD conversion)
- Measurement emulation: ESA center {{ESA_CENTER}} (ref: 3 MHz), RBW {{ESA_RBW}} (ref: 1 MHz),
  VBW {{ESA_VBW}} (ref: 100 Hz), spectrum span {{ESA_SPAN}} (ref: 0–100 MHz)

### Modulation and actuators
- Error-signal scheme per lock: {{MOD_SCHEME}} (ref: AOM frequency shift 1 MHz, heterodyne demodulation)
- EOM: Vπ {{EOM_VPI}}, bandwidth {{EOM_BW}}, range {{EOM_RANGE}}
- Fiber stretcher (PZT): gain {{PZT_GAIN}} (rad/V), resonances {{PZT_RESONANCES}}, range {{PZT_RANGE}}

### Controller (device-agnostic; Moku:Delta is one preset)
- ADC/DAC: sample rate {{ADC_FS}}, bits {{ADC_BITS}}, input range {{ADC_RANGE}}, DAC bits/range {{DAC_SPEC}}
- Total loop delay {{LOOP_DELAY}}, controller update rate {{CTRL_RATE}}
- Filter/PID structure and initial gains: {{PID_INIT}}
- Moku:Delta preset values: {{MOKU_PRESET}} (look up from official documentation; cite the page)

### Targets
- Target squeezing at GBS input: {{TARGET_SQZ_DB}}; GBS mode count / requirements: {{GBS_REQUIREMENTS}}

---

## 3. Conventions

- ħ = 1, quadratures (q̂, p̂), vacuum variance = 1/2. All reported noise levels are in dB relative to
  shot noise: S_dB = 10 log10(V / V_SN).
- Phase definitions: write down, in `notes/conventions.md`, the sign convention for every phase
  (pump, probe, LO, squeezing angle) before implementing anything. Every module must refer to it.
- SI units internally. Every equation in notes defines every symbol.

---

## 4. Workflow (phases)

Keep `PROCESS_LOG.md` updated throughout: for each phase, what was done, what worked, what was
corrected and why. This log is the source for the skill in §10.

### Phase 0 — Clarify (coordinating agent + me)

Before research, the coordinating agent clarifies unfilled placeholders and priorities with me,
following the [planning interrogation rules](../AGENTS.md#planning-interrogation). Decisions go to
`DECISIONS.md`.

### Phase 1 — Research (parallel agents, one effect/topic per worktree)

Survey first (one agent): list every effect that could affect squeezing or the locks in this system.
At minimum consider: OPA squeezing with distributed loss in the waveguide; pump depletion / pump noise;
laser phase noise and ASE (and the LO filter cavity); fiber dispersion, SPM/XPM, Raman (stimulated and
spontaneous, thermal phonons), GAWBS, Brillouin, thermal/acoustic phase noise, polarization drift;
detector QE, bandwidth, circuit noise; mode mismatch; probe leakage into the BHD; heterodyne error-signal
generation and demodulation; actuator dynamics; ADC/DAC quantization; loop delay.

Deep dive (one agent per topic) writes `notes/<topic>.md` containing:
- Governing equations with every symbol defined
- Assumptions and validity range
- Frequency band where the effect acts
- Scaling with fiber length, power, temperature, pulse width
- **Order-of-magnitude estimate with our parameters** (use reference values if placeholders are unfilled,
  and say so)
- Validation test(s) that the implementation must pass
- Primary references (author, title, journal, year, DOI). Never fabricate a citation; mark anything you
  could not verify as UNVERIFIED.

Notes are written in {{NOTES_LANGUAGE}}.

### Phase 2 — Synthesis and triage

One agent merges the notes into `notes/00_summary.md`: which effects are kept, in which simulation
stage, and why. Effects estimated to be negligible for all planned fiber lengths are kept as optional
blocks but disabled by default.

### Phase 3 — Plan (approval gate)

Write `PLAN.md`: architecture (see §5), module list with interfaces, data flow, test plan, and the
split of implementation work across worktrees. Then review the plan with me following the
[planning interrogation rules](../AGENTS.md#planning-interrogation), update `PLAN.md` and
`DECISIONS.md`, and **stop and wait for my approval before implementing.**

### Phase 4 — Implement (parallel worktrees per module)

### Phase 5 — Review loop

A reviewer agent with fresh context checks each module against its note and cited sources:
conventions, signs, factors of 2, units, numerical convergence. Findings go back to the implementing
agent. Repeat until all validation tests in §7 pass.

### Phase 6 — Exploration (after I confirm the simulator is trustworthy)

See §8.

---

## 5. Architecture requirements

### 5.1 Two stages

- **Offline stage** (per fiber length / per configuration): classical pulse propagation by GNLSE solved
  with RK4IP, plus a quantum-noise layer (linearized fluctuations or truncated-Wigner stochastic GNLSE;
  justify the choice in the plan). Output: lookup tables (e.g. squeezing/anti-squeezing transfer, excess
  noise spectra, LO–signal temporal mode overlap, nonlinear phase vs power).
  - Effects inside the GNLSE (dispersion, SPM, Raman, loss, noise terms) are toggled as **terms of the
    equation**, not as sequential pipeline steps.
  - `LASER_MODE = CW` must also work (GNLSE reduces or is skipped).
- **Time-loop stage**: time-domain simulation of slow phase noise, error-signal generation, and control,
  using the offline tables.

### 5.2 Time-loop pipeline

- Each block is created by a factory function that takes parameters and returns a step function
  `step(state, params) -> state`.
- The experiment runner receives a list of blocks and applies them in order at every time step
  (or every chunk; see 5.4).
- Blocks can contain their own list of sub-blocks with the same interface (nested pipelines), e.g.
  a controller made of `hpf >> pid >> eom_actuator` and `lpf >> pid >> fiber_stretcher` in parallel.
- Feedback: controllers write outputs to `state`; actuators read them through explicit **delay buffers**
  (this is how loop delay is modeled).
- **The state must track phases per optical path** (pump, probe to squeezer OPA, probe to phase-detection
  OPA, squeezed path including each fiber, LO path), not a single phase error. Locks are defined as
  relations between path phases. This is required to model the loss of common-path cancellation when
  fibers of different lengths are inserted.
- Pipelines can be built both from Python lists in `if __name__ == "__main__":` and from YAML config
  files (the latter for automated exploration). Both produce the same objects.

### 5.3 Device-agnostic parameters

Only generic parameters (sample rates, bit depths, ranges, bandwidths, delays, actuator transfer
functions) appear in the code. Device presets (Moku:Delta, others) are config files. Provide an
abstract controller interface that a later real-hardware wrapper (Moku Python API) can implement, so the
same BoTorch optimization loop can drive either the simulation or the hardware.

### 5.4 Performance

Heterodyne error signals (ref: 1 MHz) and ESA emulation (ref: 3 MHz center) require MHz-rate signals,
while lock dynamics need ≥ 0.1 s of simulated time. Choose and justify one of: (a) baseband/complex-
envelope modeling of the RF signals, (b) chunked vectorized processing with loop delay ≥ chunk length,
(c) compiled per-sample loop (e.g. numba). The choice must not prevent modeling loop delay accurately.

### 5.5 Block families needed to cover references A, B, C and our setup

Each family is a set of interchangeable blocks with a common interface. Switching between reference
experiments must require only configuration changes, never code changes.

- **Source**: single-pass waveguide OPA (ref A), OPO with cavity linewidth (refs B, C), pulsed source
  (ours). Parameters may be callables of other parameters (e.g. pump-dependent loss).
- **Topology / paths**: a probe may be absent (ref C), co-propagating with the squeezed light (ref B),
  or in a separate branch with its own detection OPA (ref A). Implemented through the per-path phase
  model (§5.2): co-propagation = shared path segments, so common-mode noise cancels automatically.
- **Noise sources**: may be independent or **partially correlated between paths** (correlation
  coefficient as a parameter), e.g. two fibers in one conduit.
- **Modulators**: placeable on any path (probe, LO, pump, squeezed). Every modulator on the squeezed
  path must apply its insertion loss to the squeezed state.
- **Error-signal generators**:
  1. Coherent probe, heterodyne (AOM frequency shift) + demodulation (ref A).
  2. Dither + DC intensity extremum lock on parametric gain (ref B, lock 1).
  3. Probe–LO interference in the BHD + dither demodulation (ref B, lock 2).
  4. Quantum noise lock: band-limited noise power of the BHD output → envelope/log detector →
     dither demodulation (ref C). This requires the **statistical fluctuation of a finite-bandwidth,
     finite-time noise-power estimate** (radiometer-type noise) as the sensor noise, and the noise
     variance as a function of the instantaneous phase, integrated over the detection band using the
     source spectrum.
- **Controllers**: parallel band-split (HPF/LPF per actuator) **and** stacked cascades where each stage
  both drives its actuator and feeds the next stage (ref C). Signal adders, ramp generators.
- **Actuators**: transfer function with mechanical resonances (ref B: 14 kHz, 22 kHz), range limits with
  saturation, and lock-loss detection with optional re-acquisition (ramp) logic.
- **Detector**: frequency response including non-flat features (ref C: interference dip near 0.8–1 GHz),
  circuit noise, QE.

---

## 6. Outputs

- Time series (emulating ESA zero span): shot noise, squeezed noise, anti-squeezed noise in dB, plus the
  phase-scanned trace (as in the reference Fig. 3).
- Spectra: noise power vs frequency (absolute and normalized to shot noise), with and without circuit
  noise subtraction (as in reference Fig. 4).
- Residual phase error: time series and PSD, per lock.
- Squeezing / anti-squeezing vs pump power (as in reference Fig. 5).
- Every run saves: config, random seed, git commit hash, and results.

---

## 7. Validation tests (must pass before Phase 6)

1. **Reference reproduction (single-mode model).** With
   R± = L + (1 − L)·exp(±2√(αP)) and R′± ≈ R± cos²θ + R∓ sin²θ,
   where L is total effective loss, α is the SHG efficiency, P is pump power and θ is RMS phase jitter,
   the parameters L = 8%, α = 906 %/W, P = 640 mW, θ = 9 mrad must give about −10.1 dB squeezing and
   about +20.5 dB anti-squeezing (the paper reports 20.2 ± 0.2 dB measured). The full simulation with the
   reference loss budget and a lock tuned to 9 mrad must agree with this within a tolerance you justify.
2. Reproduce the shape of reference Fig. 5 (L = 8%/θ = 9 mrad and L = 12%/θ = 14 mrad curves).
3. Pure loss: V_out = η V_in + (1 − η)·V_vac exactly.
4. No nonlinearity, no loss: vacuum stays vacuum.
5. GNLSE/RK4IP: fundamental soliton invariance; 2nd-order soliton recurrence at the soliton period;
   photon number conserved with Raman on; step-halving reduces global error ~16×.
6. Raman noise → 0 as temperature → 0; GAWBS peak frequencies consistent with cladding diameter.
7. Control: closed-loop transfer function from simulation matches the analytic one for a linear
   test case; with all feedback disabled, residual phase error equals the free-running noise.
8. Topology test: conventional tap topology shows the tap-ratio trade-off (SNR vs added loss) described
   in the reference; phase-detection-OPA topology removes the tap loss.
9. Reference B: reproduce the trends of its Fig. 3 (squeezing / anti-squeezing vs normalized pump power
   for θ = 1.5° and 3.9°, with pump-dependent loss L(x)), including the optimum pump power x_opt where
   squeezing starts to degrade. Report any mismatch with Eq. (2) found in the research phase.
10. Reference C: the dither alone (±4.5° sinusoid, i.e. ~0.056 rad RMS) must appear as a residual phase
    noise contribution; total residual phase noise with the reference's noise sources should be of order
    0.1 rad. Lock duration must be limited by the slow actuator's range when drift exceeds it.
11. Co-propagation test: when probe and squeezed light share all path segments, path noise in shared
    segments must not appear in the corresponding lock's residual error.

---

## 8. Exploration (Phase 6)

- First ablation: disable one effect/block at a time and record its impact on squeezing at the GBS input.
- Then combine only effects that mattered.
- Optimize controller parameters (filter cutoffs, PID gains, actuator split frequency) with BoTorch
  (constrained: lock maintained, servo bump within limits), per fiber length.
- Report findings in `notes/findings.md`. Every claimed insight must include a physical explanation and a
  check that it is not a numerical artifact (step size, chunk size, seed).

---

## 9. Skill write-up (after development is finished)

From `PROCESS_LOG.md`, write an **English** skill (`SKILL.md`) describing the generic workflow:
research → synthesis/triage → plan (approval gate) → implement in parallel → review loop → exploration.
Include the note template, validation-first rules, placeholder handling, and lessons learned.
Keep project-specific physics out of the skill body (put it in an example section).
