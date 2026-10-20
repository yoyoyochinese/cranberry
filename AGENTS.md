# Agent instructions

The skill rules below apply when executing the squeezed-light simulation project described in
[docs/squeezing_sim_prompt.md](docs/squeezing_sim_prompt.md). Read that prompt in full before starting
its workflow. The prompt defines the physics, architecture, phases, validation, and deliverables;
this file defines skill usage for the coordinating agent and all workers in every worktree.

## Skill usage

**Agents invoke the skills below on their own, at the points specified. I will not invoke them
manually.** If an agent cannot invoke a skill (it is not installed for that agent, or it is configured
so that only the user can run it), the agent stops at that point, tells me exactly which skill to run
and with what input, and waits.

### Planning interrogation

- When:
  - Phase 0: grill-me, to resolve unfilled placeholders, lock topology priority, and which reference
    experiment (A/B/C) to reproduce first.
  - Phase 3: grill-with-docs on `PLAN.md`, so the plan is checked against `notes/` and
    `notes/conventions.md`. Fall back to grill-me if grill-with-docs is not installed.
  - Any later point where a design decision would be expensive to reverse (e.g. changing the `State`
    structure or the block interface).
- Who runs it: only the coordinating agent in my main session. Parallel worker agents never run it.
- Who answers: me. If I do not respond, do not guess the answers. Write the open questions to
  `DECISIONS.md` under "Open" and continue only with work that does not depend on them.
- Output: every resolved question goes to `DECISIONS.md` (decision, reason, date, phase). Every agent
  reads `DECISIONS.md` before starting a task.

### ponytail / ponytail-review (minimal code)

- When: Phase 4 only, and only while implementing the **internals of a single block** (the body of a
  step function, a numerical routine, a filter). Optionally ponytail-review on block internals after
  that block's tests pass.
- Not during: Phases 0–3, Phase 5, and any work on the framework layer.
- **Protected list** — these are requirements, not speculative flexibility. They must never be removed,
  merged, or simplified on the grounds of being unused, over-general, or speculative:
  - The experiment runner, `State`, the per-path phase model, delay buffers
  - The block factory interface and nested pipeline composition (parallel and stacked)
  - YAML config loading and its equivalence with Python-list pipelines
  - Device-agnostic parameters, device presets, the abstract controller interface (for real hardware)
  - Blocks that exist to reproduce references A, B, C, even if my current configuration does not use them
  - Toggle flags for GNLSE terms and noise sources
  - Logging of config, seed, and commit hash
  - Validation tests
- Reuse of libraries (numpy, scipy.signal, scipy.fft, existing GNLSE implementations) is encouraged,
  but every reused function must be checked against `notes/conventions.md` (e.g. FFT sign and
  normalization conventions, which differ between GNLSE codes) and wrapped behind the block interface.
- If ponytail is installed in a form that applies to every turn (e.g. as a plugin with hooks), the agent
  must still apply the protected list strictly, and should turn it off (return to normal mode) when
  moving from block internals to framework work.
- ponytail-review does not check correctness. The Phase 5 physics reviewer remains mandatory.

### Other skills

- An agent may use other installed skills (e.g. a research-papers skill for literature search) when
  they clearly help the current phase. Prefer them over general web search for finding and citing papers.

### Logging

Every skill invocation is logged in `PROCESS_LOG.md`: phase, skill, why it was invoked, outcome, and
whether it helped. This feeds the skill write-up in §10 of the project prompt.

