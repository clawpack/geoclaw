# Topography / Input Handling — Roadmap

*Status: living document. Records the merge order for the topo & dtopo input
rework, what has landed, what was superseded and why, and what is deliberately
deferred.*

**Status (2026-10-13): the stack has landed.** Every geoclaw and clawutil PR
in §1 merged to `master` on 2026-10-08, in the order below (rjleveque/doc#5,
the user docs, is still open); the supersession clean-up is done;
§2.6 is closed by clawpack/geoclaw#755. The rest of §2 and §3 are still open.
§1 below the "Landed" table is kept as the record of the plan the merges
followed.

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
| This roadmap | clawpack/geoclaw#748 | `15bcbc30` |
| Silent-failure fixes (order 1) | clawpack/geoclaw#739 | `658738fe` |
| Units policy (order 2) | clawpack/geoclaw#741 | `96fcfeb9` |
| `topo0save` slot indexing | clawpack/geoclaw#742 | `947f5bde` |
| Multilayer OpenMP link | clawpack/geoclaw#743 | `3b94335b` |
| isaac `setplot` storm format | clawpack/geoclaw#744 | `4bfae834` |
| 1D `topo.data` / `dtopo.data` layouts (order 3a) | clawpack/geoclaw#745 | `ab5916e7` |
| `num_dim` passed to the data classes | clawpack/clawutil#208 | `3236167e` |
| Input parity (order 3b) | clawpack/geoclaw#749 | `1687a021` |
| Units policy, ASCII coverage (§2.6) | clawpack/geoclaw#755 | `58d0290c` |

Merged 2026-10-08 in Randy LeVeque's order (#748, #739, #741, #742, #743,
#745, then clawutil#208 after #745 as required below), followed by #749 and
#755. #755 re-lands #753, which had merged into its stacked base branch
`topo-input-parity` rather than `master`. `Test GeoClaw` on `master` passed at
each merge commit from `947f5bde` through `58d0290c`. Before merging, the
combined tree was tested locally (Python unit tests, the regression suite under
gcc debug and optimized flags, and `examples/1d_classic/okada_dtopo`); see the
update to §4's "not verified" note.

### In review / open

| Change | Tracking | Status |
|---|---|---|
| coarsen + `align` with no aligned point in the crop: Python `crop_indices` raises `EmptyAlignedCropWindow`; Fortran `apply_align_buffer_coarsen` stops naming the file and axis | clawpack/geoclaw#761 | clawpack/geoclaw#762 (open, branch `crop-align-overshoot`, head `c1ac4691`). Verified 2026-10-19: the new Python tests and the `topo_crop` Fortran regression (types 3 and 4) fail on `0ba4db7b` and pass with the fix; the two sweeps comparing against the old code are unchanged |
| One crop-coverage rule for every topo type (a crop the file doesn't cover is an error, with half a cell allowed at the edges), plus a tile-selection helper | clawpack/geoclaw#760 | coverage rule + Fortran backstop on branch `topo-crop-coverage` (stacked on #762), ahead of the `crop_bounds` PR; tile helper next. See §2.7 |

### Merge order as planned (historical)

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
| clawpack/geoclaw#756 | `fix-geotiff-reader` | type 5 (GeoTIFF) reader: dead `import gdal` replaced by rasterio → `osgeo` fallback; swapped row/column counts; pixel-corner origin read as a grid point; nodata not mapped to NaN; `read_header` type 5 branch raised `NameError` |

**Merge geoclaw#745 before clawutil#208.** #745 defaults to `num_dim=2` and is
inert alone; clawutil landing first against an older geoclaw would break every
geoclaw `setrun`. #208 guards that with a `TypeError` fallback, which should be
removed once the geoclaw requirement catches up.

### Superseded — close, do not merge

*Done:* #727 is closed and the `topo-input-unify` and `topo-input-parity`
branches are deleted (checked 2026-10-13).

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
  cell-centered data; what tolerance, expressed in cells or degrees? §2.8's
  PR A2 (Python wrap-on-read, required to match Fortran) needs the same
  answer.
- Does the rule apply to `coordinate_system == 1` (Cartesian)? Presumably not —
  the wrap gate added in `topo-input-parity` already answers this and should be
  reused rather than re-derived.
- Should the synthesized column be written into the descriptor, or handled in
  the Fortran reader?

### 2.2 dtopo crop/wrap parity in the Fortran reader

**Deferred to mandli/geoclaw#15** ("Complete dtopo crop/wrap parity with topo
(ASCII + netCDF) in the Fortran reader"), open. `topo-input-parity` brings the
Python side to parity; the Fortran dtopo reader still lags. Since §2.8's PR A
dropped `crop_bounds` from `dtopo.data`, #15's plan to parse it no longer
applies; see §2.8.

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

### 2.6 Units policy: ASCII coverage and override unification

**Done: clawpack/geoclaw#755 (`58d0290c`, re-landing #753).** `topo_ascii` is
now `assume+check` and `dtopo_ascii` is `warn+assume` (no deformation bound
exists for rule 5); `assume_units` is a `read()` argument on both readers; `ft`
converts; `MetInspector.assume_units` takes a dict, with the bool deprecated.
`units_policy.md` and the `UNITS_POLICY` registry are the record. The text
below is the plan as written before #755.

*Originally:* **Owner: unassigned. Opened by clawpack/geoclaw#741.** Policy itself lives in
`units_policy.md`; this entry records only what that PR deliberately left open.

#741 states GeoClaw's units policy and makes it enforceable — a `UNITS_POLICY`
registry in `units.py` that both `tests/test_units_policy.py` and the table in
`units_policy.md` are generated from. It closed four violations and left **two
rows open on purpose**, `topo_ascii` and `dtopo_ascii`:

> ASCII carries no units and has no override; elevation in cm or feet is read
> as meters with no message and no sanity check.

Unlike the other deferrals in this section, half of this one is already
machine-checked: both rows carry a `gap=` string and their conformance tests are
marked `xfail(strict=True)`, so they are reported on every run and closing
either without updating the registry fails the build. What is *not* recorded
anywhere tracked is the rest of the work, which is why it is written down here.

**The design question.** Rule 1 of the policy is "units must be declared in the
file; never silently assumed". An ASCII topo header is `ncols / nrows / xlower /
ylower / cellsize / nodata_value` — there is no field to declare them in, and
none is proposed. So `on_missing='raise'` is not available: it would refuse
every ASCII file GeoClaw has ever read. The honest answer is a distinct
conformance state — *assume the contract unit, and magnitude-check it* — which
has to be visibly different in the table from today's `silent-assume`, or the
policy quietly reads as though ASCII were exempt. Rule 5 is doing the whole job
alone on this path.

**What the work involves**, beyond the two rows:

- `assume_units` promoted from the netCDF-only `nc_params` channel to a
  first-class `read()` argument on both `Topography` and `DTopography` — the
  name is meaningless inside `nc_params` for a `.tt3`.
- `_check_magnitude` extended to ASCII reads, raising as it does for NetCDF,
  with `skip_sanity_check` as the documented opt-out. **This can break runs that
  work today** and wants its own release note.
- Unit-string matching widened: it is currently exact-string and case-sensitive,
  so `Meters` and `KM` raise, and `feet`/`ft` have no conversion at all — a file
  that correctly declares `units="ft"` cannot be read. Note the CF alias tables
  deliberately contain British spellings (`metre`, `metres per second`); those
  are *data* that real files declare and must not be "corrected".
- `MetInspector.assume_units` accepting `{role: unit}`. Its `bool` is
  defensible — Met covers several variables with different contract units — and
  `units_policy.md` says so; the dict form is for parity, not correction.

**Base it on #749, not #741.** `topo-input-parity` moves the format dispatch
into `gridded_input.py`'s reader registry, which is where the ASCII read path
now lives; work written against #741's `if/elif` ladder would be rewritten. The
same refactor also makes the job smaller: `gridded_input`'s docstring already
states that shared unit handling belongs in the product rather than the reader,
so the ASCII case falls out of hoisting it instead of being added twice.

### 2.7 Crops that run off a file (clawpack/geoclaw#760)

**Decision (2026-10-19): a crop the file doesn't cover is an error, for every
topo type, checked in Python at `make data`.** "Covered" allows half a cell
beyond the outermost grid points.

**Status (2026-10-19): implemented on branch `topo-crop-coverage`**, a PR of its
own ahead of the `crop_bounds` PR rather than folded into it. Python: one check,
`data._check_crop_covered`, for types 2/3 (in domain coordinates, after
`x/y_shift`) and type 4 (in file coordinates until the `crop_bounds` PR, with
the ±360 shifts allowed when `resolve_wrap` permits). One error message for all
cases, naming the file, its cell-edge extent and the usable overlap. Fortran:
type 4 now stops like types 2/3 when a requested crop misses its file. Verified:
the 22 rejection and tile-edge cases in `tests/test_topo_crop_coverage.py`
and `test_crop_missing_the_file_stops[4]` fail on `dc6b65d8`, and all 31 cases
plus both backstop runs pass with the change; full `tests/` passes (669). One example regression,
`bowl-slosh-netcdf[coarsened]`, caught a round-off in the first version (clipping
to a point rebuilt as `(x0 - dx/2) + dx/2` dropped the first column). It is now
pinned by `test_netcdf_crop_at_the_points_keeps_the_outermost_points`.

**Behavior changes:** an ASCII (type 2/3) crop past the file, including across
the antimeridian, now errors instead of being clipped by Fortran. Type 4's three
error messages become one.

Behavior on `0ba4db7b`, from a 0.25° cell-centered test tile written as type 3
and type 4:

| crop vs. file | `Topography.read` | `make data` t3 | `make data` t4 | Fortran t2/3 | Fortran t4 |
|---|---|---|---|---|---|
| partly outside | clips and warns | writes it | raises | clips silently | clips silently |
| entirely outside | returns the whole file | writes it | raises | `stop 1` | 0×0 grid kept in the list |
| exactly the tile's cell edges | ok | ok | raises | ok | ok |

Points to carry into the implementation:

- The type-4 check (`_validate_crop_bounds`, and the latitude check in
  `_resolve_topo_records`) compares the crop against the outermost cell
  *centers*. That is why a crop at a tile's own edges is rejected. The new rule's
  half-cell allowance fixes it.
- Type 4 currently raises one of three different errors, depending on which
  axis is off. When latitude is off, the message blames the antimeridian, even
  for a regional tile. Use one message that names the file extent and the
  intersection.
- Fortran type 4 turns an entry that misses the crop into `mx = my = 0` with
  extents `(0, 0)` and keeps it in the topo list. Nothing has checked how
  `topoarea`, `rectintegral` or the `topo0save` overlap test treat such an
  entry, so stop with a message naming the file instead (done; see Status
  above). In the one run before the fix
  (`test_crop_missing_the_file_stops[4]`, domain on the real file) it went
  through `topoarea` without crashing and only tripped the coverage check.
  Type 2/3 already stops with "topo_crop_extent does not overlap topo file",
  printing the file and the crop.
- The coverage check itself (`topo arrays do not cover domain`,
  `topo_module.f90` after `topoarea`) ended with a bare `stop`, so the run
  exited with status 0 and scripts and CI couldn't tell it from success. It is
  now `stop 1` (2026-10-19, same branch). Only the fatal branch changed: an
  uncovered area over 1% of the domain stops, while one between 1e-12 and 1%
  still only warns ("do not quite cover") and continues. GeoClaw never
  documented an exit status, so nothing is broken by the change.
  `test_descriptor_crop_without_buffer_does_not_cover_ring` now also asserts
  the nonzero status.
- Tile selection: Randy's use case is a long list of tiles, of which only those
  that intersect the crop should be used. Give it a helper (working name
  `topotools.tiles_covering(paths, crop_extent, coarsen=, align=, buffer=)`)
  that sets each `crop_extent` to the intersection and drops tiles whose kept
  window is empty. It should decide that with `crop_indices`, so that it agrees
  with Fortran. That relies on #761 (§1, "In review / open"). `kmltools` can
  then draw the extents actually used.
- Python and Fortran search for the aligned start from opposite ends of an axis
  stored descending (e.g. latitude north-to-south). Python works on an
  ascending view and searches from the south. Fortran works in file-index order
  and searches from the north. By reasoning, not by test, both pick the same
  coarse lattice: the `coarsen` candidates cover every residue class whichever
  end they start from. So they should agree on whether a crop is empty and on
  which points are kept, except perhaps when `align` sits exactly half-way
  between lattice points (a tie). The `topo_crop` regressions don't cover align
  on a descending axis.

### 2.8 Crop windows, `crop_bounds`, and the `topo.data` notation

**Owner: Kyle Mandli; Randy LeVeque reviews every PR, and A is the important
one. Planned 2026-10-15. PR A implemented 2026-10-21 on branch
`topo-crop-frame`, stacked on #763; B, A2 and C not started.**

Randy asked two questions about `topo.data`. First, what is the difference
between the per-file `crop_extent` line and the netCDF descriptor's
`crop_bounds`? Second, why does one file use three notations
(`value =: key`, `value # key`, `key = value`)? He also wants to plot
topography across the antimeridian from Python alone.

**`crop_bounds` is not a duplicate of `crop_extent`.** It is the crop for one
`topo.data` entry, in that entry's *file* coordinates. For a crop that crosses
the seam, `topo_entries` writes two entries for one file. Each has its own
`lon_wrap_offset` and `crop_bounds`, while both repeat the user's
`crop_extent`. Fortran uses `nc_crop_bounds` first for type 4
(`topo_module.f90`, "Priority: nc_crop_bounds …"). The dtopo descriptor's
`crop_bounds` has no parser case, so it only produces an "unknown key"
warning. The met descriptor's copy is ignored, because met carries its crop on
its own `crop_extent` line.

**Confirmed bug: netCDF and ASCII crop in different frames when `x_shift` is
set.** ASCII (types 2/3) crops in domain coordinates, after the shift, as the
`topotools` "Region terminology" docstring says. netCDF (type 4) crops in file
coordinates, before the shift, in both Python and Fortran. That docstring also
says `crop_bounds` subtracts `x_shift`, which the code does not do. Checked
2026-10-15 against `0ba4db7b`:

- The test file is global in longitude at 1°, with `z = 1000 cos(lon)`, run on
  the domain `[0, 20] × [-5, 5]` with `crop_extent = [0, 20, -5, 5]`.
- netCDF with no shift matches the exact cell averages to `2.5e-2`. That is
  inside the bilinear `O(dx²)` bound of about `3.8e-2`.
- ASCII with `x_shift = 10` matches the domain-frame field to `2.5e-2` and
  misses the file-frame field by 44.
- netCDF with `x_shift = 10` stops with "topo arrays do not cover domain"
  (overlap 100 of 200), because Fortran loaded file `x ∈ [0, 20]`.
- In Python, type 4 with `x_shift = 360` and crop `[180, 190]` returns a
  single column at `x = 540`. A type-4 crop that misses the file entirely
  returns the whole file. #763 is expected to make that second case an error;
  check once it merges.

#763's coverage check followed the same frames (ASCII shifted, type 4
unshifted), so PR A stacks on #763 and changes both sides.

**PR A, as implemented (2026-10-21).** Type 4 crops in domain coordinates in
Python (`gridded_input` picks the window on shifted axes; `TopographyData`
checks coverage the same way as ASCII) and in Fortran (the netCDF crop
always goes through the `topo_crop_extent` branch). `topo.data` has no
`crop_bounds` for topo, dtopo or met. A hand-built
`[4, path, TopoMetadata]` entry whose crop is only in its metadata has it
written as `crop_extent`. Fortran stops with "Regenerate topo.data" on an
old file that still has a `crop_bounds` line. Verified at the PR's commit:
`tests/test_topo_crop_frame.py` and the full Python suite, the `topo_crop`
regressions (including the netCDF + `x_shift` case and the old-file stop),
`bowl-slosh-netcdf` (including coarsened), `dtopo`, `topo0save`,
`force_dry`, and `chile2010` (including `lon_360`, the only end-to-end user of
hand-built entries with a −360 offset). The windows are documented, with a
diagram, in `topo_crop_windows.md`. Not yet done from the PR A row: the user
docs.

**ASCII and netCDF (Kyle, 2026-10-22).** The windowing logic and the names
apply to both formats wherever it makes sense, and any difference has to be
justified. `topo_crop_windows.md` lists the differences that remain. Three are
not justified:

- The antimeridian wrap is netCDF-only. This moves into PR B: a per-entry
  `lon_wrap_offset` line for every type.
- Fortran computes the window in two separate places.
- `fort.geo` reports the ASCII header but the netCDF window.

The last two can go in one small Fortran refactor PR, separate from B.

**Crop round-off tolerance (folded into PR A, Kyle 2026-10-21).** Python's
crop window compared exactly, while Fortran's had a `1e-6·dx` allowance. So a
crop written on a grid line could keep that line in Fortran and drop it in
Python, without a warning. Example: crop `[-99.5, -99.0]` on a 0.1° grid
written by `Topography.write` gave Python five points, x = -99.4 … -99.0. This
is because the writer stores the spacing as `9.9999999999994e-02` and the
point read back is `-99.50000000000003`. Over 1,905 grid-line crop edges
(spacings 0.1° to 1/240°, shifts including 0 and 360), Python and the Fortran
formula disagreed on 1,017. There is now one constant on each side, used
everywhere a crop edge meets grid points: `coordinate_tools.CROP_TOL` and the
Fortran `CROP_TOL` parameter in `topo_module.f90`, both `1e-6` of the spacing.
`tests/test_crop_tolerance.py` checks that the two values agree, that the
Python window matches the Fortran index formula on that sweep, and that the
-99.5 case keeps six points for ASCII and netCDF. Setting the constant to 0
makes five of those tests fail.

It guards round-off only. A user-facing tolerance, or Randy's earlier ~0.5·dx
(which snaps crops to cell edges and so changes which points are selected),
is a separate decision. It would need a `topo.data` line, so it belongs with
PR B and needs Randy's input. Follow-ups not in PR A:

- `Topography.write` takes `cellsize` from `x[1] - x[0]`, through the `delta`
  property. `(x[-1] - x[0]) / (n - 1)` would be more accurate, but `delta`
  feeds crop, align and every ASCII header written, so changing it needs its
  own check against golden files and regression outputs.
- The met reader's Fortran `crop_indices` still compares exactly. It should
  use the same constant in PR C.

**Decided (Kyle, 2026-10-15):**

- Drop `crop_bounds` from the topo, dtopo and met descriptors. Fortran derives
  each entry's file window as `crop_extent − lon_wrap_offset − shift`, which is
  the existing `topo_crop_extent` branch. Masking against the file coordinates
  clips the window to the file.
- Use `value =: key` everywhere, with an optional trailing `# comment`. This is
  what clawutil's `data_write(..., description=...)` already writes.
- Detect old files by checking the key. Fortran still reads values in a fixed
  order, but it also checks the key after `=:` and stops with "regenerate with
  `make data`" if it is wrong. Without this check, an old `key = value` line
  would be misread silently.
- Python wraps on read and returns one `Topography` across the seam. It must
  match what Fortran loads: the same points and values, and the same
  coarsen/align/buffer grid.
- Durable docs and a diagram of the windows go here in `dev/design/`. The
  user-facing pages go in the clawpack doc repository.

**PRs, in order:**

| PR | Scope |
|---|---|
| A | Document and diagram the windows. Drop `crop_bounds`. Make type 4 crop in domain coordinates in Python and Fortran, including #763's check. Failing test first: netCDF + `x_shift` must equal ASCII + `x_shift` (the setup above). Release note for netCDF users who crop with a shift. |
| A2 | Python wrap-on-read for geographic longitude axes, gated like `resolve_wrap`. Reuse `_compute_lon_entries`, and keep the seam column once where Fortran stores it twice. Verify with a periodic manufactured field: `Z` equals the exact field, `x` is strictly ascending with uniform `dx` across the seam, and the result equals the union of the Fortran entry windows. |
| B | `topo.data` and `dtopo.data` to `value =: key`. The descriptor gets a fixed set of lines (`fill_value` always written, plus a flag). The 1D dtopo `read(iunit,*) dtopofname, dtopotype` splits into two reads. Add `TopographyData.read` and a write → read round-trip test. |
| C | The met storm data file and its `&file_info` block move to the same convention, using B's Fortran helper. |

The windows to document and diagram:

1. The file extent.
2. `crop_extent`, in domain coordinates.
3. Each entry's file window (today's `crop_bounds`).
4. The aligned/buffered/coarsened index window that is actually loaded.
5. The extent of the loaded arrays.
6. The computational domain plus ghost cells, which is the fallback when no
   crop is given.

**Interactions:**

- **§2.1 (seam gap).** A2 runs into the same question. Joining a cell-centred
  global grid across the seam has no hole: `-179.9979 + 360 = 180.0021` sits
  one `dx` from `179.9979`. The two Fortran entries still leave the hole. A2
  is required to match Fortran, so §2.1 must be decided before or alongside
  A2.
- **§2.2 / mandli/geoclaw#15.** Its step 5 ("parse `crop_bounds`") becomes
  "use `crop_extent − lon_wrap_offset`" after A. `DTopoInspector.dtopo_entries`
  also folds the buffer into `crop_bounds`, which contradicts #15's own note
  that the buffer is applied in Fortran. Update #15 when it is picked up.
- **Met (separate PR, reproduction first).** These are suspected from reading
  the code and have not been reproduced. Fortran shifts a whole 0–360 file
  to ±180 (`gridded_met_forcing_module.f90`, `lon_wrap == 360`), and
  `crop_indices` then assumes longitude increases steadily. This would
  mis-select crops that straddle 0° on 0–360 files, and domains that cross
  180°.

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
- `fix-geotiff-reader`, checked 2026-10-11 against `master` at `6df4c00a`:
  `tests/ -m "not remote"` is 621 passed / 5 skipped / 1 xfailed with rasterio
  1.5.2 (bundled GDAL 3.12.2). The new tests cover a non-square raster with
  nodata, S→N row order, header-only reads, and type 3 ↔ type 5 equivalence
  where the GeoTIFF pixel centers sit on the type 3 grid points. **Not
  verified:** the `osgeo` backend. No GDAL Python bindings were available, so
  its parametrized case skipped. Independently, rjleveque (geoclaw#756 thread,
  2026-10-13) found `topo.x`/`topo.y` agree with a direct rasterio read, and
  NCEI CUDEM tiles carry `AREA_OR_POINT=Area` (cell-edge bounding box, points
  `delta/2` inside), matching the pixel-center convention used here.

**Update 2026-10-13.** Resolved by the merges: CI on a machine other than the
author's (`Test GeoClaw` on `master` at `1687a021` and `58d0290c`), and the
Fortran regression suite on the combined tree (42 passed under gcc debug and
optimized flags, locally, before merging). Still not verified here: the
remote-fetch tests (marked `remote`), and the five other `examples/1d_classic`
cases.

**Not verified (as of 2026-09-07):** that `topo-input-parity` passes CI on a machine other than the
author's, and the remote-fetch tests throughout (marked `remote`, skipped).
Also not re-run after #745's rebase: the other five `examples/1d_classic` cases,
and any Fortran regression suite (`-m regression`) on the combined tree.

**§2.8 (2026-10-15, against `0ba4db7b`).** The crop-frame results were
checked by running scratch GeoClaw problems: one level, spherical, built from
`examples/tsunami/bowl-slosh-netcdf`'s Makefile. The topography each run
built was read from `fort.a0000` and compared with exact cell averages. The
Python results came from direct `Topography.read` calls. PR A turns the
setup into `test_shifted_crop_is_in_domain_coordinates` in
`tests/regression/topo_crop`; its expected values use the cos(latitude)
area weighting of the spherical topo integral (`bilinearintegral.f`).

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
