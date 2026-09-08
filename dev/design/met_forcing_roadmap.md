# Meteorological Forcing — Post-Refactor Roadmap

> **Provenance.** Reconstructed 2026-08-23 from the Obsidian vault mirror
> *"Met Forcing Post-Refactor Roadmap"* (last synced 2026-07-30) + git history +
> merged code. The original `.plans/met_forcing_roadmap.md` was never tracked in
> git (unlike `met_forcing_refactor.md` / `met_forcing_docs_plan.md`, which were
> recovered verbatim from a blob), so the vault mirror is the only source; this is
> a reconstruction, not a recovered original. Every merged-status / PR claim below
> was verified against the geoclaw merge history: S1 = PR #715, S2 = PR #716,
> OWI/ASCII off-by-one = PR #717, Decision 1→B (met package) = PR #724 /
> commit `f39c0937` — all merged to `clawpack/master` by 2026-07-30. Obsidian
> `[[wikilinks]]` in the mirror have been rendered as plain text and `.plans/`
> paths repointed to `dev/design/`. Lines marked `VERIFY:` are not confirmed
> against code.

*Status: living document. Companion to the completed `met_forcing_refactor.md`
(object-model + Fortran rename) and `met_forcing_docs_plan.md` (documentation).*

**Progress (2026-07-30):** S1 (Fortran `met_*` de-prefix, PR #715), S2 (Python
OWI WIN/PRE reader/writer, PR #716) + the OWI/ASCII off-by-one fix (PR #717),
**Decision 1 → B** (the `clawpack.geoclaw.met` package + `surge` deprecation shim +
`MetData`/`rundata.met_data` aliases), and the **naming/docs conformance**
follow-up (geoclaw examples/tests + clawutil `met_data` + Sphinx docs) are all
merged to `clawpack/master`. **Next up: S3** (track-reader harvest), then the
**M1 parametric generator** — now set up as a tracked project, currently pinned.

**Progress (2026-09-01):** **S3a** (IBTrACS reader type/API hygiene) is done — see
below. It is the first of a short series driven by a downstream need to run a
45-year North Atlantic historical ensemble (IBTrACS v4 + HURDAT2, 1980–2025) through
GeoClaw, which surfaced that **the historical readers cannot currently produce a
runnable storm file at all**: `read_hurdat`/`read_ibtracs` mark missing radii `-1`
while `write_geoclaw` tests only `np.isnan`, so `fill_dict` never fires and `-1` is
written verbatim (a regression against v5.9.0, whose writer checked `== -1`).
**That is now also done** (S4-prerequisite, below): the missing-value contract is
unified on NaN and `test_storm_io[ibtracs]` is no longer `xfail`. **Multi-storm
ingestion is also done** (**S5**, new — `iter_hurdat` / `iter_ibtracs` /
`iter_atcf` plus `catalog_*` indices), and so is **S3** (quadrant radii,
`Basin`/`StormStatus` enums, the remaining `<U1` truncations, and HURDAT2's
post-2021 RMW column), and so is **S4** (`met/reconstruction.py`, PR #736), whose
coefficients were transcribed from the papers and verified against the 21,683
IBTrACS points that report an RMW. Note the scale of the gap S4 fills: RMW is
reported at only **55%** of North Atlantic 1980–2025 track points and just
**0.4%** in the 1980s–90s, so most of a 45-year ensemble has no reported storm
geometry.

**Progress (2026-09-04):** S1–S5 are closed. Three short-term items remain, all
opened by the S4 work rather than planned:

* **S6** — the Willoughby `X1` regression-family error. Patch and full
  verification prepared; the coefficient edit itself is the maintainer's, being
  a wind model.
* **S7** — `adjust_max_wind` over-corrects the motion asymmetry. Left alone
  deliberately: it moves every forward-model golden, *and* a downstream
  calibration was derived under the current value, so it is not a self-contained
  change.
* **S8** — seven parametric models have no Fortran-level regression coverage.
  This is the gap that let S6 survive, which is the argument for closing it.

After those, the **M1 parametric generator** (kernel territory, pinned).

## Context

The met-forcing refactor is merged to `clawpack/master`. The object model is now
met-generic (`Track` / `StormTrack` / `ParametricMetForcing` /
`GriddedMetForcing`), the Fortran modules are renamed
(`met_forcing_module` / `parametric_met_forcing_module` /
`gridded_met_forcing_module`), and the Sphinx/tutorial documentation is done.

This roadmap triages the *next* wave of work: finishing the "surge → met"
naming transition, trimming verbose Fortran names, adding meteorological-format
support (OWI, richer track readers), and a Python wind-field generator
(the parametric family now; the USACE/Chow solver later).

**Guiding principle (from the 2026-07 audit):** GeoClaw's recurring failure mode
is *forking a parallel implementation and letting it drift*. Every item below
favors **extending the one merged object model** over introducing a second one.

**Kernel-territory note:** the wind/pressure models (M1's parametric family —
Holland, CLE, … — and the L2 Chow PBL PDE solver) are numerical kernel code. Any
wind model must be verified with a real check (Python↔Fortran field equivalence,
gradient-wind far-field balance, stationary-storm symmetry, or a match to
published Chow/SLOSH output), not a plausible-looking run.

---

## Decision 1 — Python "surge" → "met" naming  ✅ DONE (Option B, merged 2026-07-30)

**Option B shipped:** `met/` is now the real package (implementation moved from
`surge/`); `surge/` is a thin `DeprecationWarning` shim; `MetData = SurgeData` in
`clawpack.geoclaw.data`; `rundata.met_data` aliases `surge_data` (clawutil);
examples/tests and the Sphinx docs conform to the `met` name. The `surge.data`
wire filename and the `Storm` class name are unchanged. The full rename (Option C
— `SurgeData`→`MetData`, `surge.data`→`met` wire, dir/label renames) remains the
eventual major-version end state, tracked as **L1**. Option table kept below for
history.

Contract surface a rename touches: package path `clawpack.geoclaw.surge`
(imported by community `setrun.py`); the `Storm` class; `SurgeData` +
`rundata.surge_data` (the attr is registered in **clawutil**, a separate repo);
the `surge.data` filename (**cross-language contract** with `met_forcing_module`);
and Sphinx labels `surgedata` / `quick_surge` / `setrun_surge` (doc repo).

| Option | End state | Consequences | Effort |
|---|---|---|---|
| **A. Docstrings/comments only** | `surge/` stays; wording reframed | Zero API/wire/build/downstream churn; does **not** fix the core inconsistency (met objects inside a `surge` package). | S |
| **B. Add `met` aliases, deprecate `surge` (RECOMMENDED)** | Canonical `clawpack.geoclaw.met` re-exports the model; `surge` becomes a `DeprecationWarning` shim; optional `MetData = SurgeData` + accept `rundata.met_data` | New users get clean naming; **no forced downstream break**; `surge.data` wire + Fortran untouched; clawutil change only if `met_data` attr wanted. Two names coexist for a deprecation cycle. | M |
| **C. Full rename to `met`** | `surge`→`met`, `SurgeData`→`MetData`, `surge.data`→`met` wire, `surge_data`→`met_data` | Cleanest, fully coherent; **breaks every community `setrun.py`**; needs coordinated geoclaw + clawutil + doc release + Fortran wire change. | L (multi-repo) |

**Recommendation: B now, C as the eventual major-version end state.** B resolves
the packaging inconsistency without breaking the storm-surge user base, and sets
up C as the deprecation-cycle endpoint. Keep the class name `Storm` (it *is* a
storm). Keep the `surge.data` filename for now (cross-language + cross-repo cost
not yet justified). Internal docstring/comment reframing (A) is absorbed into B.

---

## Decision 2 — Fortran `met_*` de-prefix  ✅ DONE

The only redundant prefixing was four public names + one subroutine inside
`gridded_met_forcing_module.f90` (a module whose name already says `met_forcing`),
none referenced by any external `use` site: `has_met_crop`→`has_crop`,
`met_crop_extent`→`crop_extent`, `met_x_shift`/`met_y_shift`→`x_shift`/`y_shift`
(now an exact match to the topo/dtopo names they parallel), subroutine
`met_crop_indices`→`crop_indices`. Stale pre-rename `.o` files removed. The module
*names* themselves are unchanged.

*Verified behavior-neutral:* the gridded netCDF + Holland aux regression tests
pass against existing goldens (`rtol=1e-14/atol=1e-8`, no `GEOCLAW_REGEN`).
Delivered on branch `met-deprefix-gridded`.

---

## Work items

### SHORT TERM

- **S1. Fortran `met_*` de-prefix** — ✅ done (Decision 2).

- **S2. OWI (Oceanweather WIN/PRE) reader/writer — ✅ DONE (PR #716).**
  Rewrote the orphaned `surge/data_storms.py` into a clean field reader/writer
  (`read_owi` / `write_owi` / `read_owi_start_time` + the `OWIData` dataclass; SI
  in memory, mbar/f10.4 on disk), registered it in `meson.build`, and added thin
  `GriddedMetForcing.from_owi` / `to_owi` helpers. Tests: object round-trip +
  real-`isaac` sample reads in `test_met_objects.py`, plus a format-1 gridded
  regression case in `tests/regression/met_forcing/`. Documented the f10.4
  precision limit.
  - **OWI/ASCII coordinate off-by-one — ✅ FIXED (PR #717).** The Fortran OWI
    reader placed grid point *i* at `sw + i*d` instead of `sw + (i-1)*d`, shifting
    the whole grid one cell NE and disagreeing with the NetCDF path. Proven with a
    new OWI-vs-NetCDF equivalence regression test (failed pre-fix by ~1535 Pa),
    fixed, and the `aux_owi` + isaac `owi_ascii` gauge goldens regenerated. Note:
    `test_isaac.py`'s OWI coordinate reconstruction was updated to match.

- **S3. Harvest low-risk track-reader improvements — ✅ DONE (PR #NNN, verified
  against `8f0984d8`).** Quadrant wind radii, `Basin`/`StormStatus` enums, the
  remaining `<U1` truncations, and the plot-swath guard. No parallel
  `CycloneDataset`: everything hangs off `StormTrack`.
  - **`wind_radii`** — shape `(n_times, n_thresholds, 4)` in metres, last axis
    NE/SE/SW/NW, middle axis matching `wind_radii_thresholds` (34/50/64 kt as
    m/s). Read for ATCF, HURDAT2 **and** IBTrACS. **Verification is
    cross-archive agreement, not shape:** all three archives report *identical*
    radii at Hurricane Ike's peak record — 34 kt `[213.0, 185.2, 166.7, 222.2]`
    km, 50 kt `[138.9, 92.6, 74.1, 111.1]`, 64 kt `[55.6, 46.3, 37.0, 46.3]` —
    through three independent parsers. A plumbing test would pass equally well
    with the quadrants transposed or the thresholds swapped; this would not.
    Plus invariants: radii shrink as the threshold rises, and a reported RMW sits
    inside the 34-kt envelope.
  - **Zero is data in HURDAT2 and missing in ATCF, and that asymmetry is real.**
    HURDAT2 has an explicit `-999`, so a zero quadrant radius there is a genuine
    zero (a weak system has no 34-kt winds); ATCF leaves an unavailable quadrant
    as `0`, indistinguishable from a real zero, so those become `NaN` — matching
    how this reader already treats ATCF's other radius fields. `wind_radii` is
    therefore self-consistent under "NaN means not reported"; documented on
    `StormTrack` and locked by a test.
  - **HURDAT2 carries a radius of maximum wind after all.** Field 21, added with
    the **2021 season** (100% populated 2021–2025, absent before — a handful of
    2018–2020 records), now read into `max_wind_radius`. Its meaning and units
    were established empirically rather than from the format PDF (which is not
    machine-readable): **2,676 of 2,676** coincident 2021–2025 North Atlantic
    records match IBTrACS `usa_rmw` *exactly* when read as nautical miles, max
    difference 0.000 m. This corrects the standing claim that HURDAT2 never
    supplies RMW — true for 1851–2020, false since. Pre-2021 revisions end the
    line at 20 fields, so `_sentinel_to_nan` now also treats an empty field as
    missing.
  - **`Basin` / `StormStatus` enums** plus `standardize_basin(code, source)` /
    `standardize_status(code)`, folding in `ATCF_basins`, `TCVitals_Basins`,
    `TC_designations` and IBTrACS's own codes. Added **alongside** the existing
    strings, reached via `basin_code`, `basin_standard()` and `status()`:
    `_STATE_FIELDS` is an explicit tuple, so additive attributes cause zero
    golden churn, whereas turning `basin` into an enum would have changed what
    every reader reports and how it serializes. Unknown codes degrade to
    `UNKNOWN` rather than raising.
  - **`<U1` truncation** fixed for `read_jma` and `read_tcvitals` (HURDAT2's was
    unavoidably fixed in S5). Asserted on content rather than dtype, because
    ATCF's classification comes back from pandas as object dtype — and, noted in
    passing, carries the source's leading whitespace (`' TD'`); `standardize_status`
    strips, so the enums are unaffected.
  - **`StormTrack.plot(plot_swath=True)`** masks non-finite and non-positive
    radii instead of feeding them to a matplotlib patch, and warns once when
    nothing is drawable. HURDAT2 supplies no outer radius at all, so this was
    reachable from a plain `read_hurdat` + `plot`.
  - **Zero existing goldens regenerated.** The ATCF quadrant pivot is
    deliberately a *separate* pass over the pre-`groupby` frame
    (`_atcf_quadrant_radii`), leaving the legacy pipeline untouched, because
    `atcf_geoclaw.txt` is the golden a mistake there would have moved silently.
    One new fixture, `tests/data/storm/hurdat_ike_mature.txt` (12 real
    mature-stage Ike records, NOAA/AOML public domain) — the bundled
    `hurdat.txt` is all genesis-stage zeros, so without it the offline
    cross-archive check would have had no content.
  - *The original "…and an IBTrACS netCDF path if the current reader is CSV-only"
    clause was **stale** and has been struck: `read_ibtracs` has always been an
    xarray/netCDF reader.* The bundled `ibtracs.nc` fixture lacks the `quadrant`
    dimension, so IBTrACS quadrant coverage is `remote`-only
    (`test_full_basin_sweep`), which also carries the RMW cross-check above.

- **S3a. IBTrACS reader type/API hygiene — ✅ DONE (PR #NNN, verified against
  `7ef7a816`).** Prerequisite for S4: the fill path could not run at all, because
  `read_ibtracs` left `self.t` as an `xarray.DataArray`, so `storm.t[n]` was a 0-d
  DataArray and `fill_rad_w_other_source`'s `interp` raised
  `ValueError: Could not convert object to NumPy datetime`.
  - `self.t` is now a `datetime64[s]` ndarray (matching every other reader),
    `classification` is decoded from bytes (`b'TD'` → `'TD'`), `ds.dims.keys()` →
    `ds.sizes` (xarray `FutureWarning`; `Dataset.dims` becomes a set of names), and
    the storm-dimension squeeze is explicit so a single-observation storm no longer
    has `date_time` collapsed along with it.
  - Truncating to seconds also drops the sub-second roundoff IBTrACS decodes
    (`…T06:00:00.000039936`). **This moved two rows of the written storm file back
    into exact agreement with the committed `ibtracs_geoclaw.txt` baseline**
    (`-3.60000003e+03` → `-3.60000000e+03` and `7.19999997e+03` → `7.20000000e+03`
    at data rows 103 and 105); the v5.9.0-era golden holds the exact values, so the
    roundoff was a post-v5.9.0 regression rather than the baseline being wrong.
    Verification is that agreement plus four contract tests in
    `tests/test_met_objects.py` — types, whole-second times, a read under
    `FutureWarning`-as-error, and `fill_rad_w_other_source` accepting a scalar
    *and* a 0-d DataArray. All four fail on `7ef7a816`.
  - Only golden regenerated: `characterization/read_ibtracs.json`, which was
    storing the *xarray repr string* for `t` (including `Size: 952B`, coords and
    attrs) and `"b'TD'"` for `classification` — i.e. exactly the type leak being
    fixed, and a snapshot that would have broken on any xarray repr change.
  - `test_storm.py::test_storm_io[ibtracs]` stays `xfail`; the `-1`/NaN sentinel
    bug that causes it is untouched here and is the subject of the next item.

- **S4-prerequisite. Unify the missing-data sentinel on NaN — ✅ DONE (PR #NNN,
  verified against `7ef7a816`).** The `fill_dict` machinery S4 is built on has
  been dead for HURDAT/JMA/IBTrACS since the write-path rewrite: those readers
  marked missing values `-1` while `write_geoclaw` tested only `np.isnan`
  (`met/parametric.py:172`), so no fill ever fired and `-1` was written verbatim
  into the storm file. v5.9.0's writer checked `== -1`, so this was a regression,
  not a design choice.
  - **Contract, now in three parts.** (1) *In memory, `np.nan` is the only
    missing marker*, in every reader — `read_atcf` already did this; hurdat, jma,
    ibtracs and tcvitals now agree, and archive sentinels (HURDAT2 `-99`/`-999`,
    tcvitals `-99`/`-999`) are normalized **before** unit conversion, so a
    sentinel is never scaled into a physical-looking value (`-999 mbar` was
    becoming `-99900.0 Pa`). (2) *At the write boundary a negative is also
    accepted as missing*, via `ParametricMetForcing._is_missing`, so hand-built
    objects and community fill functions still on the `-1` convention keep
    working — with a `DeprecationWarning`. (3) **Zero is never missing.** HURDAT2
    reports genuinely-zero quadrant wind radii for weak systems, so a "treat
    `<= 0` as missing" rule — the obvious alternative design — would corrupt real
    data. A test locks this against a future refactor.
  - Three further writer bugs found and fixed while in here: the built-in 500 km
    `storm_radius` default was applied with `dict.update()` *after* copying the
    caller's `fill_dict`, silently discarding a caller-supplied `storm_radius`
    fill (now `setdefault`); fills were written back onto the storm object, so
    writing the same storm twice with different `fill_dict`s reused the first
    fill (now applied to local copies, and the storm is not mutated); and a
    fill's return value was never re-checked, so a fill that legitimately fails
    wrote its sentinel into the file (now re-checked, and the row is skipped).
    `write_geoclaw` additionally warns when it emits a non-positive radius,
    naming the offending casts — the check that would have caught
    `hurdat_geoclaw.txt` being unrunnable.
  - **Verification: `test_storm_io[ibtracs]` is no longer `xfail`.** The
    committed `ibtracs_geoclaw.txt` baseline — an independent v5.9.0-era artifact
    — is reproduced **without regenerating it**. Seven new contract tests in
    `tests/test_met_objects.py` cover NaN/`-1` write equivalence, failed fills
    being skipped, zero-is-data, non-mutation, the `setdefault` fix, HURDAT
    sentinel normalization, and a positivity invariant over the bundled inputs;
    nine of them fail on `7ef7a816`. The Fortran side is untouched: the
    `met_forcing` regression storm input is **byte-identical**, so the aux
    goldens need no `GEOCLAW_REGEN`.
  - Goldens regenerated: `characterization/read_hurdat.json`,
    `read_jma.json`, `read_ibtracs.json` — the diffs are purely `-1.0` → `NaN`,
    which `test_storm_characterization.py:100` already anticipated in a comment.
    `hurdat_geoclaw.txt` / `jma_geoclaw.txt` are deliberately **not** regenerated
    here: neither format carries RMW/ROCI, so `test_storm.py`'s silent
    `radius[:] = 0.0` normalization was converted into an explicit, commented
    placeholder `fill_dict` rather than deleted (deleting it would skip every row
    and emit a zero-cast file). Those two baselines are still unrunnable and are
    S4's to fix.
  - *Release-note item:* a custom fill function that returns `-1` used to have
    that `-1` written into the storm file (a silently broken run) and will now
    have the row skipped instead. A correctness restoration, but it changes
    output cast counts for such users.

- **S4. Parametric geometry reconstruction (missing-data filling).** Restore and
  modernize the old code's "basic reconstruction from statistical and physical
  data" for the recurring pain of **missing storm geometry** in commonly-available
  best tracks (ATCF/HURDAT/IBTrACS): no radius of maximum winds (RMW), no
  outer/last-closed-isobar radius (ROCI), or only one of max-wind-speed /
  central-pressure. The readers mark these `NaN` and warn, but refuse to fill
  physically — the only automatic default is a flat constant
  (`storm_radius = 500 km`, `met/parametric.py:150`). *(The `fill_dict` path they
  plug into is live again as of the S4-prerequisite item above; before that it
  was dead for HURDAT/JMA/IBTrACS.)* This item adds **opt-in,
  explicit** estimator functions with the **existing `fill_dict` signature**
  `fn(t, storm) -> value (SI)`, so they drop straight into
  `ParametricMetForcing.write_geoclaw(fill_dict=...)` and compose with the existing
  `fill_rad_w_other_source` time-interpolator (`met/track.py:1035`). **No new
  object model, no Fortran, no wire change** (extend, don't fork). The default
  write path stays warn/skip — nothing is silently reconstructed. New module
  `met/reconstruction.py`; optional `default_fill_dict(track, models={...})`
  convenience assembler. Estimator families (each cites its paper; **maintainer
  confirms coefficients & provenance**):
  - **Wind–pressure relationship** — fill `max_wind_speed` ↔ `central_pressure`:
    Knaff & Zehr (2007) (recommended; operational standard), alt Atkinson &
    Holliday (1977). *WPR is the essential companion not in the reference note —
    flagged for sign-off.*
  - **RMW from intensity + latitude** — fill `max_wind_radius`: Willoughby et al.
    (2006) (recommended), alts Kimball & Mulekar (2004) / Carrasco et al. (2014),
    optional Vickery & Wadhera (2008) for Holland-B-consistent RMW.
  - **Outer/storm radius (ROCI) from environment** — replace the 500 km constant:
    Chavas et al. (2016), with a **climatological fallback** so it degrades to a
    physically-motivated default when environment (SST/outflow) is unknown.

  *Not numerical-kernel code* — empirical regressions from cited papers, not
  solver/Riemann/limiter math, so it does not hit the kernel STOP rule; coefficient
  choices remain the maintainer's. **Verification = reproduce the published
  relationship**, not a manufactured PDE: match published reference values
  (Willoughby RMW at documented `Vmax`/lat; Knaff–Zehr WPR table) to tolerance;
  physical invariants (RMW↓ with `Vmax`, deeper `dp`⇒higher `Vmax`, positivity/SI,
  latitude direction); distributional check vs observed IBTrACS/EBTRK; and a
  sparse-ATCF (RMW/ROCI stripped) → `fill_dict` → `write_geoclaw` round-trip
  yielding a runnable 7-column storm file with all fields finite
  (`tests/test_met_objects.py`). Opt-in ⇒ existing aux + bowl-slosh goldens stay
  byte-identical (no `GEOCLAW_REGEN`). Sequence **after S3** (can also draw on the
  harvested quadrant radii / basin) and as a **companion feeding M1** (the
  generator consumes complete parameter sets). Reference note (TC-radius section):
  *Rework TC Geometry Estimates* (vault project note).

  - **Status: ✅ DONE (PR #NNN, verified against `8f0984d8`).**
    `met/reconstruction.py` ships the estimator interface, the reporting, and
    coefficients transcribed from the papers (maintainer-verified).
    - **Willoughby Eq. (7a)**: `R_max = 46.4·exp(−0.0155·V_max + 0.0169·φ)`,
      km / m/s / degrees, Atlantic. Recorded with its range of validity
      (`valid_above_m_s = 26`) and a warn-once below it.
    - **Vickery & Wadhera (2008)**: `ln R_max = 3.015 − 6.291e−5·Δp² +
      0.0337·ψ`, plus its three-branch heteroscedastic `σ_lnRMW`. Included
      because it is pressure-driven and so usable where a wind is absent —
      **not** as the weak-storm estimator; see below.
    - **`rmw_sampled(estimator, sigma, seed)`** for propagating geometry
      uncertainty. Opt-in and seed-required: the default path stays
      deterministic so a storm is still a pure function of its track.
  - **The input basis was the subtle part.** Eq. (7a) is fit to *maximum
    azimuthally averaged flight-level* winds (850/700 hPa), which the paper
    contrasts explicitly with HURDAT's 10 m 1-minute point maximum.
    `to_willoughby_wind` performs both corrections. This also **confirms the
    Fortran Willoughby path is on the right basis** — `adjust_max_wind`'s
    boundary-layer division is a convention, not a bug.
  - **But the Fortran over-corrects the asymmetry.** It subtracts the *full*
    translation magnitude; the canonical reduction is ~0.5F at the RMW (Phadke
    et al. 2003). Since (7a) is monotone decreasing in wind, over-subtraction
    biases the reconstructed radius high. Python uses
    `ASYMMETRY_FRACTION = 0.5`; over the 18,028 points needing a fill this
    raises the median converted wind 16.3 → 20.1 m/s and drops the count driven
    to zero from 267 to 9. **The Fortran is deliberately left alone** — changing
    it moves every forward-model path (Holland included) and all their goldens.
    The divergence is pinned by
    `test_wind_conventions_versus_the_fortran` so it cannot become accidental.
  - **Verification is skill against observation, not the reference file.** The
    transcription cases in `reconstruction_reference.json` are computed *from*
    the coefficients, so they check the functional form, not the constants; the
    file says so. The real check is `test_skill_against_observed_rmw` (remote),
    which holds out the 21,683 IBTrACS-NA points that report an RMW:

    | | Willoughby (7a) | Vickery-Wadhera |
    |---|---|---|
    | ≥26 m/s (n≈8,600) | bias −3.2 km, MAE 12.4 km | +1.7 km, 15.2 km |
    | <26 m/s (n≈13,000) | bias **−35.5 km** | **−42.0 km** |

    A units or wind-convention error fails this immediately where it passes
    every invariant.
  - **Finding that overturned the plan: VW08 is the *worse* weak-storm
    estimator**, despite being the obvious candidate. At the median 13 hPa
    deficit its Δp² term contributes ~0.01, so it degenerates to a latitude-only
    constant and cannot track weak-storm spread. Neither raw fit is usable below
    ~26 m/s; the test asserts that limitation rather than describing it.
    Downstream (`atlantic-rp`) applies an archive-calibrated variant there — a
    calibration to one basin and one wind convention is a *study's* business, not
    geoclaw's, so it deliberately does not live here.
  - **`roci_climatology`** keeps 500 km, re-documented: `storm_radius` is the
    centre of a 100 km `tanh` taper on wind *and* pressure, i.e. an
    extent-of-forcing cutoff, not an ROCI. With `X1` ~250 km a
    (`R_max` 40 km, `V_max` 50 m/s) storm still has ~8 m/s azimuthal-mean wind
    at 500 km, so tapering at a physical NA ROCI (~350–400 km) would truncate
    surge-relevant outer winds. A size-varying default should be
    `ROCI + O(X1)`, not `ROCI`.
  - **Scope narrowed on evidence.** Measured over IBTrACS v4 North Atlantic
    1980–2025 (684 storms, 39,711 track points):

    | Field | Reported | 1980s | 2020s |
    |---|---|---|---|
    | `max_wind_speed` | 100% | — | — |
    | `central_pressure` | 95.9% | 78.1% | 100% |
    | `max_wind_radius` | 54.6% | **0.4%** | 100% |
    | `storm_radius` | 54.3% | 0.0% | 99.5% |

    Two items were **dropped from the MVP** as a result:
    - **`roci_from_r34` — dropped.** r34 is reported at 44.5% of points, but it
      fills only **1,831 points (4.6%)** that ROCI itself does not — and
      **zero** in the 1980s–90s, where r34 coverage is also 0.0%. A citation, a
      code path and a test surface to cover 4.6% of the record, none of it in
      the era that needs help.
    - **Wind–pressure relationship — deferred to its own PR** (`wind_pressure`
      is a declared `NotImplementedError` with the reasoning). It would recover
      **1,682 of 41,427 points (4.1%)** that are dropped for want of a
      pressure — 1,519 of them in the 1980s (21.9% of that decade), none after
      1999. The literature review added a decisive practical objection: KZ07 and
      CK09 both need *size* inputs (R34, or V₅₀₀ which CK09 estimates as
      R34/9 − 3) and an environmental pressure, and for those 1980s points the
      size inputs are themselves almost entirely absent. The implementation would
      collapse to a climatological-V₅₀₀ fallback, i.e. a latitude- and
      P_env-adjusted Dvorak-style relation. Revisit only if that 4.1% proves to
      matter.
  - **`roci_climatology` replaces the anonymous `500e3`** at
    `parametric.py`'s `setdefault`. Same value, so nothing moves — what changes
    is that the assumption is named, cited and overridable. Whether 500 km is
    the right climatological median for the basin in use is flagged for
    maintainer confirmation against Chavas et al. (2016); it is a science change,
    not a refactor.
  - **`rmw_constant(radius)`** is provided so a track can be made runnable
    *today*, explicitly labelling its geometry an assumption rather than an
    estimate, without waiting on the paper transcription. Downstream ensembles
    are therefore unblocked; they just have to say what they assumed.
  - **CI stays green; the guard is at runtime, not in a red test.** The earlier
    plan for this item was to ship `reconstruction_reference.json` empty with a
    `pytest.fail`, so "the PR cannot go green until a human opens the paper."
    Changed deliberately: a permanently failing test on master is worse than
    useless, and the runtime `CoefficientsNotTranscribed` is a *stronger*
    guarantee — it prevents an unverified coefficient being **used**, not merely
    merged. The reference test and the three Willoughby invariant tests skip
    with an instruction naming the paper, and switch themselves on when the
    coefficients land.
  - **Verified:** `default_fill_dict()` with no arguments reproduces
    `write_geoclaw`'s prior behaviour exactly (374 passed, **zero goldens
    regenerated**); a sparse-ATCF round trip (RMW/ROCI stripped) yields a
    7-column file with every row present, all values finite and both radii
    positive, without mutating the caller's track; `coverage_report` on real data
    takes a 1985 storm from **0/31 to 31/31 writable times**.
  - **What the maintainer must confirm.** The RMW regression is **Eq. (7)** of
    Willoughby et al. (2006) (maintainer-identified; the bibliographic record is
    Crossref-verified, the equation attribution is not). Alongside the constants,
    the reference JSON requires the basin, the published units, and **two
    independent wind conventions** that no invariant test can check:
    - *Reference height.* GeoClaw's Fortran `set_willoughby_fields` already
      evaluates the paper's **Section 4** profile coefficients (`n`, `X1`, `A`)
      on a **gradient-level** wind — `adjust_max_wind` subtracts translation
      speed then divides by `atmos_boundary_layer = 0.9`. If Eq. (7) shares that
      fitting basis the Python estimator must match it; a ~11% `Vmax`
      difference otherwise.
    - *Averaging period.* The Fortran treats storm-file winds as **1-minute**
      (HURDAT2 / ATCF / IBTrACS `usa_wind`) and converts to 10-minute on output
      via `sampling_time = 0.88`; SLOSH is the one model that also applies it
      upstream, because its own fit assumes 10-min. Willoughby & Rahn Part I
      (2004) evaluates against HURDAT 1-minute winds, which *suggests* but does
      not state the Part II convention. Wrong choice ⇒ ~12% `Vmax` error.

    **The same reading also settles a question about the existing Fortran.** If
    Eq. (7) and the Section 4 coefficients are both fit to 10 m 1-minute winds,
    then the boundary-layer division before applying `n`/`X1`/`A` is an
    unintended bias rather than a convention — worth checking while the paper is
    open, and reporting separately if so. Nothing in the Python can determine
    this. *(Noted also for Chavas et al. 2016: it works in 10 m winds but does
    not state an averaging period; that matters less here, since ROCI is used
    only as a climatological constant, not a wind-dependent regression.)*
  - **Provenance constraint restated:** these estimators exist in CLIMADA
    (GPL-3) and other copyleft TC toolkits. Nothing may be copied or adapted
    into this BSD tree — transcribe from the papers only.

- **S5. Multi-storm archive ingestion — ✅ DONE (PR #NNN, verified against
  `25ec8abc`).** HURDAT2 and IBTrACS both distribute *one file per basin holding
  every storm on record* (2,004 Atlantic storms / 55,605 records; 2,298 IBTrACS
  NA storms), so driving an ensemble means walking a file, not opening one per
  storm. Previously `read_hurdat` assumed a single storm per file and **crashed**
  on a real archive — it read the next storm's header as a data line and fed
  `"AL092008"` to `np.datetime64` — while `read_ibtracs` re-scanned the file for
  every storm.
  - **API: module-level generators, no collection class.** `iter_hurdat`,
    `iter_ibtracs` and `iter_atcf` yield `StormTrack`s; `catalog_hurdat` /
    `catalog_ibtracs` return a `pandas.DataFrame` index (id, name, season, basin,
    record count, time span, peak wind) so a driver can decide what to run before
    running it. A DataFrame index is what an ensemble driver actually needs; a
    `TrackCollection` would have been the parallel object model this roadmap
    warns against. **S5/M2 boundary:** M2 (`to_xarray()` / CF track I/O) is
    *serialization of our tracks*; S5 is *ingestion of foreign archives*. They
    compose — once M2 lands the multi-storm container is
    `xr.concat(t.to_xarray() for t in iter_*(...))` — and M2 need not re-solve
    ingestion.
  - **One parser per format.** `read_hurdat` gained `storm_id` / `name` / `year`
    selectors and is now built on the same `_hurdat_blocks` framing as
    `iter_hurdat`; `read_ibtracs` keeps its signature and is built on the same
    `_select_ibtracs_storm` + `_track_from_ibtracs` as `iter_ibtracs`;
    `make_multi_structure` is a thin wrapper over `iter_atcf` (which stages
    through a temp dir instead of littering the cwd) and its test passes
    unmodified. The verification is therefore *path equivalence* — bulk and
    single-storm reads produce identical objects — rather than two goldens that
    could drift.
  - **Record-count conservation** is the primary HURDAT2 integrity check: each
    header declares how many records follow, so declared totals must equal what
    parsed *and* the number of data lines in the file. Verified over the whole
    Atlantic archive (55,605 both ways). A header claiming more records than
    follow now raises with the line number instead of absorbing the next storm.
  - `iter_ibtracs(skip_invalid=True)` is load-bearing, not a convenience:
    **57 of 741** North Atlantic 1980-2025 storms have no time at which both a
    wind and a pressure were reported (`NoDataError`), and one of them must not
    abort a basin sweep. Measured throughput: 684 tracks in 8.8 s (13 ms/storm);
    `catalog_ibtracs` over the same selection takes 0.2 s.
  - Reader-correctness fixes that were **inseparable from rewriting the HURDAT2
    parser**, and so landed here rather than in S3: `ID` held the header's
    *record count* (`"06"`) and now holds the ATCF-style id (`"AL082008"`) — the
    count is what frames the block; `basin` was the raw `"AL"` where `read_atcf`
    gives `"Atlantic"`; and `classification`/`event` were `<U1`-truncated
    (`"LO"` → `"L"`, colliding with the landfall event code) because the arrays
    are allocated in the rewritten parser. One golden regenerated:
    `characterization/read_hurdat.json`, a three-field diff. **S3 retains** the
    quadrant radii, the `Basin`/`StormStatus` enums, the same `<U1` fix for
    `read_jma`/`read_tcvitals`, and the `plot` swath masking.
  - Also shared rather than duplicated: `_IBTRACS_AGENCY_PREF`, previously a
    mutable default argument on `read_ibtracs`, so the reader and the iterator
    cannot drift.
  - **Coverage gap stated plainly:** the bundled fixtures are a single-storm
    HURDAT2 excerpt and a stripped single-storm IBTrACS slice (no `quadrant`
    dim, no `season`, only `wmo_*`/`usa_*` agencies). Offline tests synthesize
    multi-storm files by relabeling those, which cannot exercise differing
    agencies, basin crossings or real invalid storms. The only test of real
    heterogeneity is `test_full_basin_sweep`
    (`@pytest.mark.remote @pytest.mark.slow`, driven by
    `GEOCLAW_TRACK_ARCHIVES`), which also cross-checks HURDAT2 against IBTrACS
    for Ike — the same storm from two independent archives through two
    independent readers agrees on eye position to 0.11 deg. Either commit a small
    real multi-storm NA slice, or accept that offline coverage is structural
    only.
  - *Note for whoever writes NaN-bearing test comparisons:* `_storm_state`
    dictionaries cannot be compared with `==`, because `nan != nan`. Compare
    `_dumps(...)` text, as the characterization tests do.

- **S6. Fortran Willoughby `X1` used the wrong regression family — ✅ FIXED
  (maintainer-applied, 2026-09-04).** Willoughby et al. (2006) gives two
  dual-exponential regression families: Eqs. (10a–c) predict from `V_max` and
  latitude, Eqs. (11a–c) add `ln(R_max)` as a predictor. They are not
  term-by-term interchangeable — in (10) the variance that would load on
  `ln(R_max)` instead loads on `V_max` and latitude through their correlations —
  so a routine must use one family throughout.

  `set_willoughby_fields` mixes them:

  | param | in `parametric_met_forcing_module.f90` | matches | verdict |
  |---|---|---|---|
  | `n` | `2.134 + 0.0077 V − 0.4522 ln(R) − 0.0038 φ` | Eq. (11b) | ✓ |
  | `A` | `0.5913 + 0.0029 V − 0.1361 ln(R) − 0.0042 φ` | Eq. (11c) | ✓ |
  | `X1` | `317.1 − 2.026 V + 1.915 φ` | **Eq. (10a)** | ✗ |

  ```fortran
  X1 = (287.6d0 - 1.942d0*mod_mws + 7.799d0*log(mwr/1.d3) &
        + 1.819d0*abs(sloc(2))) * 1.d3
  ```

  Kernel territory, so the coefficient edit was the maintainer's; the
  verification around it was prepared separately. All three coefficients now
  come from the Eq. (11) family, with the superseded Eq. (10) forms kept as
  commented reference.

  **Measured impact, so the cost is known before the edit.** Over the 21,683
  IBTrACS-NA points with a reported `R_max`: median change **+2.8 km** on an
  `X1` of ~330 km, 5th–95th percentile −5 to +10 km, maximum 20 km. Direction is
  as the mechanism predicts (compact storms −3.8 km, large storms +8.6 km) but
  the magnitude is low-single-digit percent, because `7.799·ln(R_max)` spans
  only ~11 km across observed sizes while the intercept difference partly
  offsets the slope changes. **Real and cheap to fix; not urgent.**

  **Why it survived: the model had no Fortran-level coverage at all** (see S8).
  This item therefore also adds the first Willoughby regression case, and the
  verification is built around not needing M1's Python generator, which is
  pinned:
  - `test_willoughby_n_and_A_use_the_ln_rmax_family` and
    `test_willoughby_x1_uses_the_same_family_as_n_and_A`
    (`tests/test_reconstruction.py`) assert one family throughout — catching
    exactly this bug class rather than this instance of it. The X1 one is
    `xfail(strict=True)` until the patch lands, at which point it XPASSes and
    forces removal of the marker.
  - `test_willoughby_field_shape` checks the produced aux field directly:
    finite and non-negative wind, a real pressure deficit, monotone decay
    outside the eyewall, and no jump across the sectionally-continuous
    0.9/1.1·`R_max` band. Statements about *shape*, since
    `post_process_wind_estimate` rescales and adds translation speed so the peak
    is not `V_max`.
  - `aux_willoughby.txt`, the new regression golden, generated **after** the
    fix. (Worth stating because the trap is live: adding the test first and
    running it auto-creates a golden from the unfixed code and enshrines the
    bug. That happened twice during this work and was caught both times by
    checking the file's mtime against the source's.)

  **Measured effect on the field**, from regenerating the golden either side of
  the change (409 cells carrying more than 0.5 m/s of wind): median **+1.37%**,
  5th-95th percentile +0.27% to +2.00%, maximum +2.08%; largest absolute change
  0.150 m/s on a 41.0 m/s peak, and pressure unchanged to the byte (`X1` enters
  only the wind profile). So the earlier estimate of "low single-digit percent
  in the outer wind" is confirmed by measurement rather than inference.
  - `X2 = 25 km` and the 0.9/1.1 transition band are cleared as non-bugs; the
    latter is a documented simplification of solving Eq. (3) for ξ.

  *Non-issue checked and dismissed:* the new `log(mwr/1.d3)` term does not add a
  failure mode for a non-positive `R_max` — `n` and `A` already call
  `log(mwr/1.d3)` unguarded two lines above, so that exposure is pre-existing
  and unchanged.

- **S7. Fortran `adjust_max_wind` over-corrects the motion asymmetry — OPEN.**
  It subtracts the *full* translation magnitude:

  ```fortran
  trans_speed = sqrt(tv(1)**2 + tv(2)**2)
  mod_mws = mws - trans_speed          ! full |t|
  ```

  The canonical reduction is about half the forward speed at the radius of
  maximum winds (Phadke et al. 2003, the same pipeline as Knaff et al. 2011);
  Lin & Chavas (2012) use a gentler 0.55–0.7 on a rotated motion vector. Full
  subtraction over-corrects, and since Willoughby's Eq. (7a) is monotone
  decreasing in `V_max`, it biases a reconstructed `R_max` high.

  **Deliberately not changed, for two reasons.** First, `adjust_max_wind` is
  shared by every model passing `convert_height=.true.` — Holland included — so
  the edit moves *every* forward-model golden, not just Willoughby's. Second,
  and less obvious: `clawpack.geoclaw.met.reconstruction.ASYMMETRY_FRACTION` is
  already 0.5, and the downstream `atlantic-rp` study's RMW calibration was
  **derived under that value** and raises `CalibrationConventionMismatch` if it
  changes. Aligning the Fortran is therefore not self-contained — it would
  require re-deriving that calibration. The divergence is pinned by
  `test_wind_conventions_versus_the_fortran` so it cannot become accidental.

  If it is changed, also revisit the bound-at-zero branch
  (`trans_mod = mws / trans_speed`), which fires far less often at 0.5F: over
  the North Atlantic 1980–2025 record, points driven to zero fall from 267 to 9.

- **S8. Seven parametric wind models have no Fortran-level regression coverage —
  OPEN.** `met_forcing_module.f90` dispatches nine
  (`holland80`, `holland2008`, `holland2010`, `cle`, `slosh`, `rankine`,
  `modified_rankine`, `demaria`, `willoughby`); the goldens on disk are
  `aux_holland80`, `aux_data`, `aux_owi`. **This gap is how S6's wrong
  regression family survived**, and it is the strongest argument for closing it.

  S6 adds the Willoughby case, leaving seven. The harness is already generic —
  `tests/regression/met_forcing/setrun.py` passes `forcing` straight to
  `storm_specification_type` for parametric models, and `_check_aux` creates a
  golden when absent — so each addition is close to mechanical, and S6's case is
  the worked example. Expect surprises: any of the seven could carry its own
  coefficient error, CLE most likely, being the heaviest. Note also that a
  golden characterizes *current* behaviour, so each new case should be paired
  with at least the field-shape invariants S6 introduces, or it will simply pin
  whatever is there.

### MEDIUM TERM

- **M1. Parametric-family offline wind-field generator (Python).** Turn a
  `StormTrack` + parameterization into gridded wind/pressure fields in Python and
  write netCDF (and/or OWI) that `GriddedMetForcing` already consumes — no
  live-solver changes, no new Fortran family. Fills the `met/storm.py`
  `holland_1980`/`holland_2010`/`cle_2015`/`construct_fields` stubs; adds
  `ParametricMetForcing.to_gridded(...)`. Phasing (G1–G5): **Holland-1980**
  (matches the existing `aux_holland80.txt` golden) → Holland-2010/2008 →
  **CLE15** (ODE solve; heaviest kernel; cross-check vs the Fortran CLE *and*
  `EmilioEchevarria/tcwindfields` CLE15, MIT — reference only, not vendored).
  **Verification is the point** (see the discipline note below). *Kernel territory
  — the model numerics are implemented/verified by the maintainer.* The **Chow**
  USACE PBL wind solver is the harder sub-case, deferred to L2 (nested
  telescoping-grid iterative solver: advection + diffusion + surface drag +
  pressure-gradient + implicit Coriolis; verify vs gradient-wind far-field,
  stationary-storm symmetry, published Chow output).
  - **Status: designed, tracked, pinned (not started).** Durable design +
    Fortran-fidelity spec + verification plan: `dev/design/met_wind_field_generator.md`.
    Tracked as an Obsidian project (currently pinned):
    *Implement parametric TC wind-field generator* — this roadmap remains the
    technical index; the project note tracks resumption.
  - **Sibling: see M3** (analytic/idealized forcing) — both add new parametric
    subtypes via the same registration path and share the Python↔Fortran
    field-equivalence verification harness.

- **M2. `StormTrack.to_xarray()` / CF netCDF track I/O** — a single canonical
  serialization for tracks (idea harvested from the artifact's CF schema), giving
  multi-storm datasets a home *inside* `StormTrack` rather than a rival class.
  Delivered alongside M1's gridding layer (recommend the same CF schema as
  `create_era5_storm_file`/`create_nws13_storm_file`).

- **M3. Analytic / idealized meteorological forcing (new parametric subtypes).**
  Closed-form pressure + wind evaluated **per cell, live in Fortran** — new
  `ParametricMetForcing` subtypes on a **generic `Track`** (a moving front), not
  gridded and not a new family. Primary motivation: meteotsunami modeling
  (issue #694; the hydrocoast/Miyashita `meteotsunami` fork) generalized to a
  reusable idealized-forcing capability. Fields: closed-form pressure with
  along-track shape Φ (sine/hat/morlet/gaussian/jump) × cross-track apodization
  Ψ⊥ × temporal ramp R(t); **wind first-class** (`none`/`geostrophic`/`gradient`).
  Integration = the **same subtype path as M1's models**: register in
  `data.py` `forcing_subtype_registry` (`("parametric", <new int>)`) +
  `forcing_spec_type`/`set_storm` case → new `set_analytic_fields` in
  `parametric_met_forcing_module.f90`; mimics the `test_topography` analytic
  bathymetry pattern (`topo_module.f90` ~L375). Source coupling (src2
  pressure-gradient + wind-drag) is **untouched** — kernel-territory note applies.
  **Verification is a strength:** analytic ocean response — **Proudman resonance**
  (moving pressure jump over constant depth) and **Greenspan edge waves** — gives
  end-to-end *manufactured-solution* tests, and the Python↔Fortran
  field-equivalence harness (shared with M1) covers the field evaluator itself.
  Deliverable: a Proudman moving-jump canonical example + regression/CI. Open
  design decisions (see the spec): subtypes-vs-sibling-family, wind-coupling
  formulation, the moving-front AMR criterion (**seam 2**), and subtype naming +
  legacy signed-int alias on the `surge.data` wire. Full spec:
  *Analytic Met Forcing Plan* (vault note); tracked project (status **active**,
  GeoClaw-Users'-Workshop-2026 target):
  *Implement analytical meteorological forcing for GeoClaw*.
  **Sibling to M1** — both are new parametric subtypes sharing the
  subtype-registration path and the verification harness.

### LONG TERM

- **L1. Full `met` rename** (Decision 1, Option C) at a GeoClaw major version once
  B's deprecation cycle has run — coordinated with clawutil + doc repos.
- **L2. Chow PBL wind solver** — first as the Python generator's hardest sub-case
  (offline, per M1), then as a native Fortran forcing family *only if* verification
  shows the live-solver cost is warranted. Explicit kernel work.
- **L3. De-duplicate track readers** — retire scattered prototype readers once
  S3/M2 land.
- **L4. TC precipitation estimation (placeholder — separate scoping).** The second
  half of the *Rework TC Geometry Estimates* reference note: estimate rainfall
  from a parametric storm — from vertical velocity (Emanuel et al. 2006, Zhu et al.
  2013, Lu et al. 2018, Gori et al. 2022) or from the empirical TC-precip↔radii
  relationship (Tuleya et al. 2007, Mudd et al. 2017). A **distinct, larger
  capability** than S4's geometry filling: it needs a rainfall forcing field plus
  overland-flow coupling that does not exist yet. Listed here only so the reference
  note's precip thread is not lost; scope separately before committing.

---

## Suggested sequencing

1. **S1** — ✅ done (PR #715).
2. **S2** (OWI reader/writer) + off-by-one fix — ✅ done (PR #716 / #717).
3. **Decision 1 → B** (met package + `surge` shim + `MetData`/`met_data`) +
   naming/docs conformance (examples/tests + clawutil + Sphinx) — ✅ done.
4. **S3** (track-reader harvest) — **NEXT**.
5. **S4 — parametric geometry reconstruction** (opt-in missing-data filling).
   After S3 (can use harvested quadrant radii / basin) and a **companion to M1**
   (feeds it complete parameter sets); pure Python, opt-in, low-risk.
6. **M1 — parametric-family Python generator** (the substantive new capability;
   designed, tracked as a project, currently pinned). Recommended **after S3** so
   it is authored against the richest `StormTrack`, but it can start independently
   if appetite comes first. Chow is the harder sub-case → **L2**.
7. **M3 — analytic / idealized forcing** — the currently **active** analytic thread
   (Workshop-2026 target), **independent of S3 and M1** (it uses a generic `Track`,
   not the enriched `StormTrack`), so it can proceed now / in parallel. Sibling to
   M1; shares the subtype-registration path + verification harness.
8. **M2 / L1 / L2 / L3 / L4** as evidence and appetite allow.

## Verification discipline (all items)

Behavior-neutral changes must keep the met-forcing aux goldens and bowl-slosh
canary byte-identical (no `GEOCLAW_REGEN`). New capabilities add their own
regression cases under `tests/regression/met_forcing/` and/or object tests in
`tests/test_met_objects.py`. Cross-repo changes (clawutil / riemann / doc) are
called out explicitly and merged in lockstep.

**Parametric generator = verification engine (M1).** A faithful Python parametric
generator must reproduce the Fortran forcing-aux goldens (e.g. `aux_holland80.txt`)
on the same grid/frames/storm — a Python↔Fortran equivalence test that makes the
two implementations *mutually verifying* and is the path to retiring the hand-rolled
analytic-vortex fixtures (`_netcdf_forcing`/`_owi_forcing`). **Anti-circularity
rule:** the equivalence golden always comes from the independent Fortran path, never
regenerated from the generator, so it never grades its own output.

**Analytic forcing = end-to-end manufactured solutions (M3).** The analytic
subtypes add *closed-form ocean-response* checks — Proudman resonance (known
amplification as `c → √(gh)`) and Greenspan edge-wave amplification — that verify
the full forcing→SWE coupling, not just the field. Their per-cell field evaluator
reuses M1's Python↔Fortran field-equivalence harness (same anti-circularity rule).
Together M1 and M3 give the met-forcing subsystem both field-level and
system-level manufactured verification.
