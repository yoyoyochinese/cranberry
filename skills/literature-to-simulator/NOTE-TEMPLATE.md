# Note templates

## Deep-dive note (`notes/<topic>.md`, phase 1)

```markdown
# <Topic ID> <Title> (<file-stem>)

- Phase / topic: <phase>, <topic ID> (survey rows <E..>, decisions <D..>)
- Date:
- Inputs: <prompt, DECISIONS range, conventions, survey, papers and manuals read>
- Value tags: [Ref X] paper, [manual] datasheet, [AE] own estimate (unverified), [UNVERIFIED] not checked at source

## 0. Summary
<conclusions first, each with its number and tag>

## 1. Symbols
<every symbol defined, with units; matches notes/conventions.md>

## 2. Governing equations
## 3. Assumptions and validity range
## 4. Frequency band where the effect acts
## 5. Scaling (fiber length, power, temperature, pulse width, ...)
## 6. Order-of-magnitude estimate with our parameters
<state which values are placeholders replaced by reference values>
## 7. Mapping to the simulator (blocks, parameters, toggles)
## 8. Validation tests the implementation must pass
<numbered, each with target and tolerance>
## 9. Questions for the user / proposals for coordinator-owned files
## 10. References
<author, title, journal, year, DOI; UNVERIFIED where not checked>
```

## Literature note (`notes/lit-<citekey>.md`, exploration)

```markdown
# <Citekey> — <one-line what the paper does>

## 0. Summary
## 1. What the paper does (figure and page numbers)
## 2. Numbers (parameters, budgets, measured residuals)
## 3. What transfers to our system, and what does not
## 4. Mapping to existing blocks/overlays; new blocks needed
## 5. Proposal
## 6. UNVERIFIED items and open questions
```

## Review report (`reviews/<phase>_<module>.md`)

```markdown
# Review: <module> (<work packages>)

## 1. Findings
| ID | Severity (blocker/major/minor) | Summary | Script |
## 2. <ID> (<severity>): <title>
### What happens
### Numbers (script, configuration, seed, duration)
### Proposed fix
```

## Decision entry (`DECISIONS.md`)

```markdown
### D<n> <title>
- Decision: <numbered items>
- Reason / evidence: <notes, scripts, figures>
- Date: <YYYY-MM-DD> / Phase <p>
```
