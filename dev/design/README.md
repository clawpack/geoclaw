# `dev/design/` — durable design docs

This directory holds **tracked, durable design documents** for GeoClaw
subsystems: roadmaps, refactor design records, and implementation plans that are
meant to stay current and be reviewed alongside the code they describe. It
follows the Clawpack `dev/` convention (cf. `amrclaw/dev/`, `clawutil/dev/`,
`pyclaw/development/`) — developer-facing material that lives with the code but is
distinct from the user-facing Sphinx docs in the separate `clawpack/doc` repo.

## Two tiers — why this exists

A "source of truth" living in a gitignored path is not one: it can't be reviewed,
can't be updated in the same PR as the code, can't be `git blame`d, and can be
lost. Several of these docs *were* lost that way (they lived in a gitignored
`.plans/`). The split below is the fix.

| Tier | Location | Tracked? | Purpose |
|---|---|---|---|
| **Durable design** | `dev/design/` | ✅ yes | Few, curated, owned design docs meant to stay current. |
| **Ephemeral scratch** | `.plans/` | ❌ gitignored | Per-task, multi-author working notes; allowed to go stale; never authority. |

`.plans/` is ignored via the global git excludesfile, not this repo's
`.gitignore` — leave it that way. Do **not** move ephemeral scratch here, and do
**not** un-ignore `.plans/`.

## Keeping these current

1. **Same-PR rule.** Any PR that changes a subsystem with a design doc here
   updates that doc's status/progress in the same PR.
2. **Verifiable status.** Prefer a "verified against commit X / PR #N" note over a
   bare "done", so drift is visible.
3. **Archive, don't delete.** When a doc is superseded, keep the reasoning and
   mark it superseded rather than removing history.

## Contents

- `met\_forcing\_refactor.md` — design record for the surge→met object-model + Fortran
  rename refactor (**\*\*complete\*\***; recovered verbatim from git).
- `met\_forcing\_docs\_plan.md` — documentation/transition-guide/tutorial plan
  (approved, execution deferred; recovered verbatim from git).
- `met\_forcing\_roadmap.md` — living post-refactor roadmap (S/M/L work items;
  reconstructed from the vault mirror, git-verified).
- `met\_wind\_field\_generator.md` — parametric TC wind-field generator design &
  verification plan (**\*\*partial reconstruction\*\***; original lost, no full mirror).
- `friction-file-plan.md` — file-based friction field input implementation plan
  (relocated from `.plans/`).
- `topo\_input\_roadmap.md` — living roadmap for the topo/dtopo input rework:
  merge order for the linear PR stack, superseded branches and what replaced
  them, and the deferred items (antimeridian seam gap, dtopo Fortran parity,
  1D preprocessing scope, build-flag hygiene).
