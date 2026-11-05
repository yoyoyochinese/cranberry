---
name: literature-to-simulator
description: Multi-agent workflow for building a validated research simulator from papers: research → synthesis/triage → plan (approval gate) → parallel implementation → fresh-context review → exploration. Use when starting or running a project that turns literature into a simulation of a physical experiment, coordinating parallel worker agents in worktrees, or planning a blind validation or a configuration sweep.
---

# Literature to simulator

One **coordinator** (the agent in the user's session) drives the phases, owns the shared records, and talks to
the user. **Workers** (subagents, one task each, usually in their own git worktree) research, implement, review
and explore. The user decides; the coordinator never answers a decision on the user's behalf.

## Shared records (single sources of truth)

| File | Owner | Holds |
|---|---|---|
| `DECISIONS.md` | coordinator | Numbered decisions (decision, reason, date, phase) and an **Open** list. Every agent reads it before a task. |
| `PROCESS_LOG.md` | coordinator | Per phase: what was done, what worked, what was corrected and why, every skill invocation (phase, skill, why, outcome, helped?). Source for the skill write-up. |
| `notes/conventions.md` | coordinator, user-approved | Signs, units, FFT and normalisation conventions, symbol table. Every reused library is checked against it. |
| `notes/<topic>.md` | one worker each | Research notes; template in [NOTE-TEMPLATE.md](NOTE-TEMPLATE.md). |
| `PLAN.md` | planner worker, user-approved | Architecture, interfaces, data flow, test plan, work split. |
| `reviews/<phase>_<module>.md` + `reviews/scripts/` | reviewer workers | Findings, each backed by a reproduction script. |

Workers write only their own files. A change a worker wants in a coordinator-owned file goes into its report as
a **proposal**; the coordinator applies proposals in one batch. With this rule, 16 parallel workers committed with
zero file conflicts.

## Global rules

- **Placeholders.** Inputs the user has not filled are written `{{NAME}}`. Use the reference paper's value when
  one exists, mark the substitution in the output and `PROCESS_LOG.md`; when none exists, stop and ask.
- **Presets, not defaults.** A recommended configuration is one named preset among switchable alternatives.
  Keep bad configurations as presets too; the simulator exists to compare them.
- **Validation first.** Each note names the validation test the implementation must pass; tests are written
  before or with the code, and exploration starts only after they pass and the user trusts the simulator.
- **Citations.** Every value carries its provenance tag (paper, manual, own estimate). A citation that could not
  be verified is marked UNVERIFIED, never filled in from memory.
- **Protected list.** Framework pieces the user names as requirements (runner, state, config loading,
  device-agnostic interfaces, reference-reproduction blocks, toggles, validation tests, run logging) are kept
  even when currently unused. Minimal-code tools apply only inside single block bodies.
- **Commits.** Work on a development branch; commit at phase steps; push and PRs only when asked. Check the
  current branch before every commit.

## Phases

Each phase ends on its completion criterion. Record outcomes in `DECISIONS.md` and `PROCESS_LOG.md` before
moving on.

### 0. Clarify

Interview the user in **frontier** rounds: ask every question whose prerequisites are settled, number them, give
a recommended answer worded so "yes" accepts it; hold dependent questions for a later round.
Before asking, look facts up yourself (dispatch a worker); only decisions go to the user. Open each reference the
user cites and confirm its content matches its label.

Done when: the frontier is empty, every answer is a numbered decision, and the user has confirmed the summary.

### 1. Research

1. One worker writes a **survey**: every effect that could matter, with an order-of-magnitude estimate and a
   proposed topic split.
2. One worker per topic writes a deep-dive note in parallel (template). Workers that reach the same number
   independently have cross-validated it; record that.
3. The coordinator collects every worker's open questions, deduplicates them, settles the ones a convention or
   earlier decision already answers, and puts the rest to the user in frontier rounds.
4. The coordinator writes `notes/conventions.md` and gets user approval.

Done when: every survey topic has a merged note, every question from the notes is a decision or an Open entry,
and conventions are approved.

### 2. Synthesis and triage

One worker merges the notes into `notes/00_summary.md`: each effect kept, made optional, or dropped, with the
stage it belongs to, plus a table of every discrepancy between notes. Give it the coordinator's running memo of
worker reports too, so fixes mentioned only in reports reach the summary. Negligible effects stay as optional,
switched-off blocks.

Done when: every survey row has a triage verdict and every discrepancy has a resolution or an Open entry.

### 3. Plan (approval gate)

A planner worker drafts `PLAN.md` with explicit gate decisions. The coordinator checks it against the notes and
conventions in frontier rounds with the user, then **stops until the user approves the plan.**

Done when: the user has approved `PLAN.md` in writing.

### 4. Implement in parallel

Split by module into waves (framework first, then blocks, then integration of each reference setup). Brief every
worker with: the base branch to merge first, the names of blocks other workers are building concurrently, and the
rule that a substitute implementation gets a different name. After each merge, run the full fast suite and tell
running workers about any newly required configuration key before it breaks their presets.

Done when: every module is merged, the full suite passes, and each reference setup reproduces its validation
targets.

### 5. Review loop

A **fresh-context** reviewer checks each module against its note and sources: conventions, signs, factors of 2,
units, double counting, numerical convergence. Every finding ships with a reproduction script. Findings go back to
an implementing worker; repeat until the validation tests pass. When two reviewers disagree, adopt the one that
reproduces numerically and check whether the conclusion depends on the choice. Every new block gets this review
before its results are used.

For a setup with measured data, make a **blind** prediction: a fresh worker produces it without opening the
blinded targets, the coordinator freezes it (tag + result hashes in `PROCESS_LOG.md`), and only then are the
measurements opened and compared without tuning. Start the comparison with quantities that cancel the unknowns
(e.g. closed-loop ÷ open-loop ratios separate the loop model from the disturbance model).

Done when: every module's review is resolved and the user states they trust the simulator.

### 6. Exploration

1. **Ablation:** remove one effect at a time; record its impact on the target metric.
2. **Combine** only the effects that mattered.
3. **Compare** configurations and fusions with 10–30 seeds on common random numbers; claim a difference only when
   it exceeds the seed spread.
4. **Optimise** controller parameters (e.g. Bayesian optimisation) under stability constraints, per operating point.
5. **Report** in `notes/findings.md`: each insight with its physical explanation and an artefact check (step
   size halved, chunk size changed, seed). End with recommendations for the open decisions, a ranking of the
   candidate schemes, and bench experiments for the top two.

Evaluate each claim at the timescale it is about: a slow drift is invisible in a short run, and an optimiser
whose objective uses short runs cannot see it.

Done when: every open decision assigned to exploration has a recommendation with evidence, and every number in
the report traces to a script and seed set.

## Running parallel workers

- Worktrees may start from the wrong branch; every brief begins with "merge the development branch".
- A worker stopped by a usage limit resumes intact when messaged; respawning loses its context.
- A worker that cannot write a report file returns the text; the coordinator saves it verbatim and says so.
- Tests asserting "X is not implemented yet" break when another worker implements X; the integrating coordinator
  fixes them on merge.
- Ordering bugs between stages appear only in the full suite, never in a block's unit tests; run it after every merge.
- Keep at most a handful of heavy simulation processes running at once.

## Lessons learned

- Decide what quantity a control loop must hold fixed (at the measurement device, or at the downstream user)
  before choosing where the actuator goes.
- Text extracted from PDFs drops math symbols; confirm equation definitions on the rendered page.
- Mark provisional values when a worker cannot wait for another's result; marked values make discrepancies
  mechanical to trace in synthesis.
- Integration is where convention mismatches surface (a drive value as density vs sum, a delay counted twice).
- Verify the contents of data files before assuming two of them hold the same thing.
- Composable configuration overlays need conflict detection: an overlay that silently discards another's setting
  produces plausible wrong runs. Make the intended discard an explicit declaration.
- A numerical floor (discretisation, quantisation limit cycle) can masquerade as physics; when the step-halving
  check fails, report the value as an upper bound.
- Log each skill invocation as it happens; reconstructing the log later loses the "why".

## Example: phase-stabilised squeezed light (OPA → fiber → homodyne)

The project this skill was distilled from simulated squeezed light delivered through km-scale fiber delays to a
32-port boson-sampling circuit, with active phase locks.

- References A/B/C (published waveguide-OPA, OPO, and noise-lock experiments) became named presets and validation
  targets; the lab's own coherent-light interferometer (Ref D) was the blind-validation setup.
- Phase 0 decision that shaped the rest: locking the LO only fixes the phase at the detector, not at the circuit
  input, so the actuator moved to the squeezed path.
- Phase 1: an 84-effect survey, 16 deep-dive notes, ~60 worker questions merged into 14 user questions.
- Phase 5: the blind comparison matched the loop sensitivity function within 30 % and exposed two disturbance-model
  gaps (acoustic level ~10×, an unmodelled actuator side mode), later calibrated as overlays.
- Phase 6: two axes (bench squeezing experiment; circuit integration), eight literature schemes expressed as
  configurations. Finding: the remaining phase error is set by laser linewidth over the inter-port delay, which no
  lock removes; an in-band slow layer measuring the true path removed the drift terms (121 → 94 mrad at 30 s).
