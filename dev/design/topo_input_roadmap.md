# Topography / Input Handling — Roadmap

*Status: living document. Records the merge order for the topo & dtopo input
rework, what has landed, what was superseded and why, and what is deliberately
deferred.*

**Provenance.** Written 2026-09-07 from the branch and PR state at that date,
and amended later the same day after test-merging the open PRs in order: §1
(the stack forks at #741; #745 rebased onto it), §2.3 (layout guards, and the
`make .output` trap), §4 (#745's verification re-anchored to its post-rebase
commits) and §5 (new — why the sequence broke). Every claim below was checked
against the repositories rather than from notes; where a claim is *not* verified
it says so. Companion to `met_forcing_roadmap.md`, which covers the
meteorological side of the same input-handling effort.

---

## 1. Where things stand

### Landed

| Change | PR | Merge |
|---|---|---|
| `fetch_remote_topo`; deprecate `read_netcdf`; refactor `etopotools` | clawpack/geoclaw#726 | `a3ee4b27` |
| User docs for the above | rjleveque/doc#5 | open, commit `c17e0eff` |

### Open, in intended merge order

The topo stack is linear through #741 and then **forks**: both #745 and #749
sit directly on #741 and neither contains the other. So #739 and #741 must
merge first, in order; after that #745 and #749 may merge in either order.

| Order | PR | Branch | Contents |
|---|---|---|---|
| 1 | clawpack/geoclaw#739 | `topo-crop-silent-failures` | 13 silent-failure fixes; remote + cross-seam paths reachable from `setrun` |
| 2 | clawpack/geoclaw#741 | `units-policy-conformance` | `UNITS_POLICY` registry, `dev/design/units_policy.md`, 4 conformance fixes |
| 3a | clawpack/geoclaw#745 | `fix-1d-topo-data-format` | 1D `topo.data` / `dtopo.data` layouts (2 commits on #741) |
| 3b | clawpack/geoclaw#749 | `topo-input-parity` | `coordinate_tools`, `gridded_input`, dtopo crop/coarsen parity, `coordinate_system` wrap gate (7 commits on #741) |

`topo-input-parity` is 12 commits from `master` and carries #739 and #741 with
it. Run `pytest tests/ -m "not remote"` on it before merging.

**#745 was rebased onto #741** (2026-09-07). It was originally cut from `master`
and conflicted with #739's `TopographyData.write()` refactor; resolving that at
merge time silently defeated the PR (see §5). Rebasing puts the resolution under
#745's own CI instead. Two consequences:

- Until #741 lands, #745's GitHub diff shows **20 files, +2864/−128** — it
  includes all of #739 and #741, because a cross-fork PR cannot be retargeted at
  a branch that exists only in the fork. Review the last two commits only
  (`d18e1733..db0d1ad2`). The diff collapses to its own 2 commits automatically
  once `master` contains #741.
- #745 and #749 merge cleanly against each other, **and correctly**: verified
  that the 1D `override_order` write still lands between `ntopofiles` and the
  record loop in the combined tree, which is the ordering the 1D Fortran reader
  requires. A clean textual merge is not sufficient evidence here — that is
  exactly the failure mode §5 records.

Independent of the stack, all off `master`:

| PR | Branch | Fixes |
|---|---|---|
| clawpack/geoclaw#742 | `fix-topo0save-dtopo-index` | `topo0save` slot indexing (a released v5.14.0 regression) |
| clawpack/geoclaw#743 | `fix-multilayer-openmp-link` | `$(FFLAGS)` dropped from the multilayer link line |
| clawpack/geoclaw#744 | `fix-isaac-setplot-storm-format` | isaac `setplot` hard-coded the storm file format |
| clawpack/clawutil#208 | `fix-1d-topo-data-format` | passes `num_dim` to `TopographyData` / `DTopoData` |

**Merge geoclaw#745 before clawutil#208.** #745 defaults to `num_dim=2` and is
inert alone; clawutil landing first against an older geoclaw would break every
geoclaw `setrun`. #208 guards that with a `TypeError` fallback, which should be
removed once the geoclaw requirement catches up.

### Superseded — close, do not merge

Both verified by comparing branch contents, not from memory.

**clawpack/geoclaw#727** (`topo-buffer-wrap`, commit `80d50c79`) — "Combine
`buffer` and anti-meridional wrapping". Its mechanism was **reimplemented**, not
dropped:

| #727 | Replaced by | Where |
|---|---|---|
| `data.py::_expand_wrapped_topo` | `_resolve_topo_records` | #739 |
| `TopoInspector(buffer=…)` baking the buffer into descriptor `crop_bounds` | Fortran `nbuf4 = topo_buffer(topo_idx)` on the `nc_has_crop` branch | #739 (`c893f238`) |
| `test_topo_entries_buffer_no_wrap` / `_with_wrap` | `test_topo_entries_no_wrap`, `test_topo_entries_wrap_required`, `test_cross_seam_entries_keep_buffer_and_coarsen` | #739 |
| `test_topo_data_netcdf_wrap_and_buffer` | `test_descriptor_crop_honors_buffer`, `test_descriptor_crop_without_buffer_does_not_cover_ring` — end-to-end against the real `xgeoclaw` | #739 |

Coverage downstream is strictly stronger. Note that `c893f238` also *fixed* a
bug in this area: leaving `nbuf4 = 0` on the `nc_has_crop` branch silently
dropped `topo_buffer` for every file `topo_entries()` writes. `master` still has
the single assignment; #739 has both.

**`topo-input-unify`** (branch only, no PR) — delete. It is #727's commit plus a
`topo-data-access-cleanup` merge plus six phase commits, and all six have 1:1
rebased equivalents on `topo-input-parity`:

`19573e4e→994e6efb` (coordinate_tools), `aeff31ad→16538959` (Phase 1),
`8e20b144→8a631fa1` (Phase 2), `7494f69f→de4cdb8c` (3a), `da0f3331→6f9969bf`
(3b), `6622c535→64988a92` (3c). Parity adds `286ad84d` on top and is rebased
onto #741, which `topo-input-unify` predates. Nothing unique remains.

Per this directory's "archive, don't delete" rule, close #727 and delete the
branch **with a note naming the replacement** — that note is the archive.

---

## 2. Deferred

### 2.1 Antimeridian seam gap — design decision needed

**Owner: unassigned. Blocks the coverage half of clawpack/geoclaw#740.**

A global cell-centered DEM cropped across ±180 is written to `topo.data` as two
entries. GEBCO spans `-179.99791666 … 179.99791666`, so for a domain of
`[-190, -60]` the split yields `[-190, -180.0021]` and `[-179.9979, -60]` with a
**one-cell hole** straddling the antimeridian, and GeoClaw reports that topo
does not cover the domain.

`buffer` cannot close it, and this is not a bug in the buffer handling: after
#739 the buffer *is* applied, but it is clamped to the file's own extent and
neither entry has data past its own edge. The two cells are physically adjacent
in *cell coverage*; the gap is an artifact of point registration, where cell
centres never reach ±180.

Closing it means letting the seam entry draw its buffer from the **other end of
the file** — the cell at `-179.9979` relabelled `+180.0021`. That is more than
#727 implemented, and it needs answers first:

- When is a grid "global enough" to wrap? Exact `360°` span is too strict for
  cell-centered data; what tolerance, expressed in cells or degrees?
- Does the rule apply to `coordinate_system == 1` (Cartesian)? Presumably not —
  the wrap gate added in `topo-input-parity` already answers this and should be
  reused rather than re-derived.
- Should the synthesized column be written into the descriptor, or handled in
  the Fortran reader?

### 2.2 dtopo crop/wrap parity in the Fortran reader

**Deferred to mandli/geoclaw#15** ("Complete dtopo crop/wrap parity with topo
(ASCII + netCDF) in the Fortran reader"), open. `topo-input-parity` brings the
Python side to parity; the Fortran dtopo reader still lags.

### 2.3 What preprocessing should 1D support?

**Owner: unassigned. Opened by clawpack/geoclaw#745.**

#745 unblocks `examples/1d_classic` by writing the pre-#726 layouts for
`num_dim == 1`, and warns when a preprocessing attribute is requested that 1D
will ignore. That is deliberate non-parity, not a stopgap — but it has not been
*decided*, only defaulted.

The question is which attributes are meaningful in 1D. `crop_extent` and
`coarsen` plausibly are, on `x` only; `y_shift` is meaningless; `z_shift` and
`negate_z` would be cheap. Answering it means either teaching
`src/1d_classic/shallow/topo_module.f90` the block format or agreeing that
preprocessing belongs in Python for 1D.

Whatever is decided, the layout itself is now pinned **by field name**:
`test_topo_data_1d_uses_legacy_layout` and
`test_topo_data_1d_override_order_written_once_for_multiple_files` assert the
ordered list of `topo.data` field labels, so a future `write()` refactor that
drops or moves `override_order` fails naming the field rather than reporting a
line-count mismatch. The two-file case exists because with a single topo file a
per-file `override_order` write is indistinguishable from the correct one. This
seam has broken once already (§5) and sits between two independently evolving
concerns — #739's record resolution and the 1D layout branch — so it is worth
the redundancy.

Note the coverage situation, which is why the #726 breakage went unnoticed for
a full release cycle. #745 adds unit tests pinning the 1D *file layouts*
(`tests/test_data.py`, marked `python`), but **nothing builds or runs a
`1d_classic` example anywhere in CI**.

CI selects entirely by marker, across five jobs:

| Workflow | Job | Selection |
|---|---|---|
| `testing.yml` | `python-tests` | `python and not remote` |
| `testing.yml` | `regression-tests` | `regression and not slow and not adjoint and not remote` |
| `slow-tests.yml` | `slow-tests` | `slow` |
| `slow-tests.yml` | `regression-tests` (matrix) | `regression and not slow and not adjoint`, one leg adding `remote` |
| `slow-tests.yml` | `python-tests` (matrix) | `python`, one leg `python and remote` |

No test in any of them touches 1D — `git grep -l 1d_classic -- tests/ .github/`
is empty. So the layout tests would have caught #726, but a change to
`src/1d_classic/shallow/topo_module.f90` still would not be: nothing compiles or
runs the 1D code. An end-to-end 1D case should be added alongside whatever is
decided here; they are cheap, the whole suite of six runs in seconds.

The `regression` marker is the one to use for that — it is a first-class CI
selector (it is what runs #742's `topo0save` case and the `topo_crop`
end-to-end suite), not a local-only convention.

**Prerequisite for that test: `make .output` is not enough.** All six
`examples/1d_classic` cases generate their own input files from Python, via
`.PHONY` `topo` (and, for `okada_dtopo`, `dtopo`) targets that **nothing
depends on** — only `make all` invokes them, and `all` also runs `.plots` and
`.htmls`. So a harness that calls `make .output` gets a missing or empty
`celledges.data` / `dtopo_okada.dtt1` and the Fortran aborts in the *file*
reader, not the settings reader — a confusingly similar-looking failure to the
one #745 fixes. Observed while verifying #745: a `0`-byte `dtopo_okada.dtt1`
produced `Fortran runtime error: End of file` at
`src/1d_classic/shallow/topo_module.f90:253`, well past the `dtopo.data` read
that was actually under test.

Either wire the generated inputs into `.output`'s prerequisites (the better
fix — it makes the examples work the way every other GeoClaw example does), or
have the test call `make topo dtopo` first. Prefer the former; the current
arrangement is a trap for any automated runner, not just CI.

### 2.4 Build-flag hygiene

Two separate items, both surfaced by clawpack/geoclaw#743.

**`examples/bouss/radial_flat/Makefile:76`** hardcodes `-fopenmp` rather than
`$(FFLAGS)`. It works today only because `-fopenmp` happens to be the flag that
matters; any other user `FFLAGS` entry needed at link time is dropped.

**`ALL_LFLAGS` should arguably always include `FFLAGS`** in
`clawutil/src/Makefile.common`. Today `LFLAGS ?= $(FFLAGS)` is cancelled by any
downstream Makefile that touches `LFLAGS` before `Makefile.common` is included —
the trap #743 fixes locally. Fixing it centrally would remove the trap for every
package at once, at the cost of double-adding `FFLAGS` for the Makefiles that
already include it explicitly (harmless when linking). It is a semantics change
across all of Clawpack and wants its own PR and discussion.

### 2.5 `-DNETCDF` requires `make new` — documented, not fixed

Long-standing clawutil behaviour, newly visible now that examples differ in
`FFLAGS`. `Makefile.common:66` derives `OBJECTS` from `SOURCES`, so objects are
written **next to the sources** in `$(CLAW)/geoclaw/src/…` and
`$(CLAW)/amrclaw/src/…`, shared by every example, with nothing tying an object
to the flags it was built with. `$(EXE)` depends on `$(MAKEFILE_LIST)`, which
forces a **relink** but never a **recompile** — hence a NetCDF-built object
silently reused by a non-NetCDF example.

**Decision: document it** (`make new` when toggling NetCDF) and open a clawutil
issue. Two real fixes were considered and both rejected for now:

- A *flags-stamp dependency* is correct but makes two examples with different
  flags thrash the whole shared library on every `make`.
- *Per-flags object directories* are the better fix and too large for this pass.

---

## 3. Correction on record: clawpack/geoclaw#740

The `topo0save` bug fixed in #742 was initially believed to explain #740's
missing tsunami. **It does not**, and the reasoning is recorded here so it is
not re-tried.

`topo0save` only suppresses deformation when a real topo file outranks
GeoClaw's internal `topo_for_dtopo` grid, which is unconditionally marked for
update. #740 coarsens GEBCO 20× to `0.167°` while `dtopo_usgs100227.tt3` is
`0.101°` (verified from the file header). The dtopo is finer, so
`topo_for_dtopo` wins priority over the topo entries across the source region
and the seafloor does move there.

#742 is still a genuine released-v5.14.0 regression with a demonstrated
behavioural change; it simply is not this one. The cause of #740's missing
tsunami is **still unknown** and needs reproduction against the original GEBCO
file. Its second symptom, the coverage warning, is §2.1 above.

---

## 4. Verification notes

Claims in §1 were checked as follows, at 2026-09-07:

- Merge of #726: `git log main | grep a3ee4b27`.
- Shape of the stack: `git merge-base --is-ancestor` between each pair. Linear
  through #741; #745 and #749 are siblings on it (neither is an ancestor of the
  other), and their pairwise merge was simulated with `git merge-tree` and the
  resulting `data.py` inspected, not just checked for conflicts.
- #727 / `topo-input-unify` supersession: commit-by-commit subject comparison
  against `topo-input-parity`, plus a name-by-name check that every test #727
  added has a named equivalent downstream.
- #742: regression test fails before / passes after, with the rebuild confirmed
  fresh by comparing object and source mtimes; `chile2010` (6 tests) unchanged.
- #743: `LFLAGS` inspected via `make -p` before and after; both multi-layer
  examples build under `FFLAGS=-fopenmp` and `plane_wave` runs under
  `OMP_NUM_THREADS=2`.
- #744: both forcing families driven through the example's own `setrun`.
- #745 / clawutil#208, **at the original commit `31c31937`**: all six
  `examples/1d_classic` cases build and run; the pre-fix failure reproduces the
  reported error exactly.
- #745 **after the rebase onto #741** (`d18e1733`, `db0d1ad2`) — the code under
  the note above changed, so it was re-checked rather than carried over:
  `tests/` at `-m python` is 426 passed / 5 skipped / 3 xfailed; the two layout
  guards fail with `'topo_path' != 'override_order'` when the fix is reverted;
  `examples/1d_classic/okada_dtopo` runs to completion (60 frames, `t = 3600`)
  through both 1D Fortran readers, with `num_dim = 1` arriving via the real
  `ClawRunData` path rather than a hand-built object. Only that one example was
  re-run, not all six.

**Not verified:** that `topo-input-parity` passes CI on a machine other than the
author's, and the remote-fetch tests throughout (marked `remote`, skipped).
Also not re-run after #745's rebase: the other five `examples/1d_classic` cases,
and any Fortran regression suite (`-m regression`) on the combined tree.

---

## 5. Correction on record: the #745 merge resolution

Recorded because the failure mode is not specific to these PRs and cost a
debugging cycle.

All six open PRs reported `MERGEABLE` on GitHub while the sequence did not
merge. GitHub tests each PR against `master` only, so **mergeability is not
transitive**: #745 was cut from `master` and conflicts only appear once #739 has
landed. Nothing in the PR UI shows this.

The conflict is in `TopographyData.write()`. #739 moved the write loop behind
`_resolve_topo_records()` (a cross-seam type-4 crop expands into two entries, so
`ntopofiles` is not known until every file is resolved); #745 adds `num_dim`
branching to that same loop. A hand-resolution taken at merge time kept #745's
per-file early `continue` but **dropped its `override_order` write**. Since
`src/1d_classic/shallow/topo_module.f90:read_topo_settings` reads
`topo_missing / test_topography / ntopofiles / override_topo_order / topofname`
positionally, Fortran then reads the quoted path into a `logical` — reproducing
the exact bug #745 exists to fix, and failing #745's own
`test_topo_data_1d_uses_legacy_layout` with `assert 5 == 6`.

Two lessons, both cheap:

- **Resolve stacked conflicts in the PR, not in the merge.** A resolution made
  in an integration merge is exercised by no CI run anywhere. Rebasing #745 onto
  #741 put it under #745's own CI, which is what would have caught this.
- **A clean `git merge-tree` is not evidence of a correct merge.** Verify the
  merged artifact, not just the absence of conflict markers — for positional
  file formats especially, where a dropped line is silent until Fortran runtime.

Also confirmed while investigating, so it is not re-derived: #739/#741 and #742
touch `src/2d/shallow/topo_module.f90` in disjoint routines (the NetCDF
descriptor-crop branch around line 1481 versus the `topo0save`
dtopo-intersection loop around line 439), and #741 adds nothing to that file
beyond #739. The
`topo_type` loop variable versus `topo.topo_type` in the record loop is *not* a
behavioural difference — `_resolve_topo_records` back-assigns `topo.topo_type`
after inferring it — though the loop variable is the correct thing to write.
