# Refactor of GeoClaw Meteorological Forcing Support

**Status:** Complete — all phases merged (2026-07-25). Phase 0 (verification
baseline, PR #706), Phase 1 (Python object model, #707), Phase F1 (center-assumption
guard, #708), Phase 2 (config family/subtype split, #710), and Phase F3 (Fortran
module/type rename, #711), plus the lockstep `clawpack/riemann` rename PR. The
Fortran modules/types are now `met_forcing_module` / `parametric_met_forcing_module`
/ `gridded_met_forcing_module` and `parametric_met_forcing_type` /
`gridded_met_forcing_type`. Remaining work is documentation (see
`met_forcing_docs_plan.md`).
**Scope:** `clawpack/geoclaw` storm forcing (Python `surge/` package and the Fortran modules listed above; formerly `storm_module` / `model_storm_module` / `data_storm_module`).
**Prerequisite:** netCDF met forcing capabilities are merged. This refactor branches off that state.

This note consolidates the prior design notes (Design Note Draft, Python API Design, Fortran Interface Evolution, Fortran Capability Map, Storm Module Inventory, Target Python capability map) into a single repository-facing document, and adds a development sequence, a verification strategy, and an explicit accounting of the seams where this work can affect general (non-storm) GeoClaw functionality.

---

## 1. Motivation

GeoClaw's storm support has grown to serve several distinct roles under one abstraction: parameterized tropical cyclone forcing, gridded meteorological forcing from external datasets, ingestion of many track and event formats, and newer applications that need atmospheric forcing but are not storm-centered (extratropical cyclones, meteo-tsunamis, asteroid-impact atmospheric forcing). As a result the term **storm** now overloads at least four ideas: a physical event with track metadata, a parameterized forcing model, a generic meteorological forcing source, and a descriptor for gridded external data.

The goal is to reorganize the code around the broader concept of **meteorological forcing**, with storms retained as one important subtype rather than the umbrella. This clarifies which capability a user is invoking, removes latent bugs that come from assuming every forcing source behaves like a storm, and gives a clean extension point for non-TC forcing.

---

## 2. Terminology (locked)

- **Meteorological forcing** ("met forcing"): the umbrella concept. Any source of wind and/or pressure forcing, parameterized or gridded. Concise code identifiers (`MetForcing`, `set_met_forcing`) are preferred over verbose ones.
- **Track**: an evolving spatial reference (center over time) for a feature.
- **StormTrack**: a track carrying storm-specific metadata. Supports both tropical and extratropical systems; it does not assume tropical-only semantics.
- **Parameterized forcing**: forcing generated from a compact set of evolving parameters.
- **Gridded forcing**: forcing read from external field datasets.
- **storm**: retained only as a subtype/specialization within met forcing, never as the umbrella.

---

## 3. Conceptual model

Four layers.

1. **Meteorological forcing** (umbrella). Owns forcing activation/timing, scaling, optional windowing/ramping, dispatch/orchestration, and the extension point for non-storm applications.
2. **Event/track data**. Center location over time, intensity-like parameters, classification where relevant. Does not own gridded dataset descriptors, file mappings, or variable mappings.
3. **Parameterized meteorological forcing**. Forcing from a compact evolving description: Holland family, CLE, SLOSH, Rankine variants, deMaria, Willoughby, and future ETC or idealized parameterizations. References a track/event object; does not require all parameterized forcing to reduce to one storm representation.
4. **Gridded meteorological forcing**. File-backed wind/pressure fields: OWI/ASCII, NetCDF, future reanalysis or model products. Owns file paths, variable/dimension mapping, interpolation metadata, and dataset access details.

---

## 4. Locked design decisions

- Umbrella concept is **meteorological forcing**; `storm` is a subtype only.
- **Track / StormTrack** split; `StormTrack` inherits `Track`. `center` is the spatial field name.
- **ParametricMetForcing / GriddedMetForcing** as the two forcing families.
- **`Storm` is a compatibility wrapper over `ParametricMetForcing` only.** It does not wrap `GriddedMetForcing`; the gridded path is newer and has no legacy user base to protect.
- Core objects are **NumPy structure-of-arrays**; no pandas/xarray dependency in the core model (gridded may use xarray internally).
- **Time representation:** ~~Python `datetime` is canonical on the `Track` layer (sparse, event-like)~~ **[amended in §6 — the track axis uses `datetime64`, uniform with the gridded axis]**; `GriddedMetForcing` holds its time axis as a native `datetime64` array (dense, CF-backed). Conversion happens only at the Fortran write boundary. See Section 6.
- **`time_offset`** is Python to Fortran conversion metadata only. Python remains datetime-aware; Fortran does not. It lives on the forcing/writer boundary, not on the bare track.
- **Asymmetric I/O is a design principle:** ingest broadly (many track/data formats), emit narrowly (only GeoClaw-consumable forcing formats). The apparent read/write symmetry in the current API is not a real requirement.
- **Python parametric models are not a prerequisite.** The current Python Holland/CLE/`construct_fields` methods are stubs; the real implementations live in Fortran. Python parametric implementations are future work, not part of this refactor.
- **Fortran: semantic clarification first, physical renames last.** Preserve the current logical split; improve naming/docs; rename modules and types only after the Python concepts and compatibility strategy are stable.
- **window / ramp_width** are implemented gridded-first, with the interface left flexible enough to support them on the parameterized path later.

---

## 5. Python target object model

New modules alongside the existing `storm.py`, introduced additively.

### `Track` (`track.py`)
Generic evolving feature with a spatial center over time.
- `t`: sequence of Python `datetime`
- `center`: `(n, 2)` array of `(lon, lat)`; for a TC this is the eye, for other systems the central point (e.g. pressure minimum)
- optional metadata: `name`, `id`, `event`
- SoA layout, NumPy-backed.

### `StormTrack(Track)`
Storm-specific track metadata.
- `max_wind_speed`, `max_wind_radius`, `central_pressure`, `storm_radius`
- optional: `classification`, `basin`, `wind_speeds` (format-specific; optional)
- Supports TC and ETC; not tropical-only.

### `ParametricMetForcing` (`parametric.py`)
Forcing from a parameterized model.
- `track`: a `Track` or `StormTrack` (does not duplicate track data; multiple models may share one track)
- `model`: identifier, e.g. `"holland80"`, `"holland2010"`
- model-specific configuration
- `write_geoclaw(...)` produces the compact GeoClaw parametric forcing file.

### `GriddedMetForcing` (`gridded.py`)
Forcing from external field datasets.
- `file_paths`, format metadata, variable/dimension mappings, interpolation metadata
- forcing controls: `window`, `ramp_width`, `window_type` (optional)
- `write_data(...)` produces the gridded forcing descriptor.
- Variable/dimension mapping is reader-specific (primarily netCDF), must be communicated to Fortran, is user-customizable in Python, and is ideally resolved internally by discovery. Discovery currently lives in `netcdf_utils.py`.
- May use xarray internally.

### `Storm` (compatibility wrapper)
- Wraps a `ParametricMetForcing`. Provides `read(...)`, `write(...)`, and compatibility accessors.
- Not the canonical internal model. Transitional only. Does not wrap gridded forcing.

### API relationships
```
Track
  StormTrack(Track)

ParametricMetForcing  -> holds Track / StormTrack
GriddedMetForcing     -> holds dataset metadata
Storm                 -> compatibility wrapper over ParametricMetForcing
```

### Reader / writer layer
Readers and writers are separated from the domain objects.
- Track readers return `StormTrack`: `read_atcf`, `read_hurdat`, `read_ibtracs`, `read_jma`, `read_tcvitals` (and `read_imd` if ever implemented).
- Parameterized: `read_geoclaw(...) -> ParametricMetForcing` (with an associated `StormTrack` view where the contents support it); `write_geoclaw(...)`.
- Gridded: `read_data(...) -> GriddedMetForcing`; `write_data(...)`.
- Primary write surface is only `write_geoclaw` and `write_data`. Track-format writers (ATCF/HURDAT/JMA/IMD/TCVITALS) are not primary supported outputs.

### Utilities
- `category(...)` becomes a `StormTrack` helper (categorization is storm-specific).
- `plot(...)` moves to a `Track`/`StormTrack` helper.
- `fill_rad_w_other_source(...)` is a track/data-cleaning utility (important for IBTrACS + ATCF workflows).
- `make_multi_structure(...)` moves to a separate `tools/` (workflow) namespace, with an import shim for compatibility.
- `construct_fields` and the Python Holland/CLE methods are removed from the object surface. A future `parametric` evaluator module can take a `StormTrack` and output wind/pressure on a space-time grid, mirroring the Fortran, but that is future work.
- `available_models()` splits into Fortran-supported vs Python-supported, both clearly labeled.

---

## 6. Data representation and the datetime question

> **Amendment (Phase 1):** the track axis uses `numpy.datetime64`, uniform with
> the gridded field axis. The track axis is already `datetime64` today, so
> keeping it is behavior-neutral and avoids conversion glue; `Storm.t` and
> `Track.t` remain `datetime64` arrays. This supersedes the "Python `datetime`
> canonical on the track" decision described in the rest of this section, which
> is retained below for historical context.

The instinct toward `numpy.datetime64` (because most netCDF forcing carries it) and the locked decision to use Python `datetime` are not in conflict once the two time axes are separated.

- The **track axis** is sparse and event-like (hundreds of points), where calendar semantics and readability matter. Keep scalar Python `datetime` canonical here.
- The **gridded field axis** is dense and CF-backed (thousands of steps), where `datetime64` is the efficient, lossless representation given CF units. Keep `datetime64` arrays native on `GriddedMetForcing`.
- Conversion happens only at the write boundary into Fortran, carried by `time_offset`. The Fortran side never needs datetime awareness.

The "locked" datetime decision governs the track layer; it was never intended to bind the gridded field axis.

---

## 7. Fortran target

The Fortran side already has the structural split we want (top-level dispatcher, parameterized implementation, gridded implementation). The work is renaming, interface clarification, and reducing overloaded terminology, not a new architecture. Semantic clarification comes first; physical renames come last.

### Module and type names (preferred, applied in the rename phase)

| Current | Target | Role |
|---|---|---|
| `storm_module` | `met_forcing_module` | top-level interface and dispatch |
| `model_storm_module` | `parametric_met_forcing_module` | parameterized forcing implementation |
| `data_storm_module` | `gridded_met_forcing_module` | gridded forcing implementation |
| `model_storm_type` | `parametric_met_forcing_type` | parameterized forcing state and model context |
| `data_storm_type` | `gridded_met_forcing_type` | gridded forcing dataset container |

No heavyweight top-level `met_forcing_type` unless a real need emerges; the dispatcher module is sufficient.

### Top-level procedures
`set_met_forcing`, `set_met_forcing_fields`. The top-level module owns shared forcing controls: wind/pressure toggles, drag law, AMR refinement thresholds. Helpers: `set_parametric_forcing`, `set_gridded_forcing`, `set_parametric_fields`, `set_gridded_fields`.

### Forcing selection semantics
Replace the overloaded `storm_specification_type` with an explicit split:
- **family:** `parametric`, `gridded`
- **parametric subtype:** `holland80`, `holland2008`, `holland2010`, `cle`, `slosh`, `rankine`, `modified_rankine`, `demaria`, `willoughby`
- **gridded subtype:** `owi_ascii`, `netcdf`

The legacy signed-integer encoding stays on the `surge.data` wire as an alias during transition so current Fortran keeps reading it (see Section 9, seam 4).

### Center / location abstraction
Do not require a universal storm-location abstraction at the top level. Parameterized forcing may provide `track_center` / `track_direction`; gridded forcing has no required center API. This assumption is a real source of current bugs on the gridded path and is addressed early (see Section 9, seam 2, and the Phase F1 note in Section 10).

### Fortran capability alignment (summary)

| Module | Alignment | Action |
|---|---|---|
| `storm_module` | strong structure, weak naming | keep structure; rename later; rethink `output_storm_location` |
| `model_storm_module` | strong; richer model support than Python exposes | keep; rename later; later split track/state from model identity |
| `data_storm_module` | strong; already clearly gridded | keep; rename later; do not force a center abstraction here |

---

## 8. Python capability triage (from the inventory)

**Keep and migrate early (Phase 1):** `Storm` wrapper; `read()`/`write()` dispatch; core track fields (`t`, `eye_location`/`center`, `max_wind_speed`, `max_wind_radius`, `central_pressure`, `storm_radius`); `time_offset`; `read_geoclaw`, `read_atcf`, `read_hurdat`, `read_ibtracs`, `read_data`; `write_geoclaw`, `write_data`; `fill_rad_w_other_source`; `make_multi_structure` (to tools/). `window` and `ramp_width` stay visible.

**Keep, de-emphasize (Phase 2):** `name`, `basin`, `ID`, `classification`, `event`; `wind_speeds` (format-specific, optional); `scaling`, `window_type`; `plot`, `category`; `available_formats()`, `available_models()` (split registries).

**Deprecate or stop advertising:** track-format writers (`write_atcf/hurdat/jma/imd/tcvitals`); Python `construct_fields` and Holland/CLE stubs as if functional. `read_imd` is a stub today; decide deprecate-now vs documented future-work placeholder.

---

## 9. General-functionality seams (must not break)

The maintainer owns storm forcing, so the automated verification gate carries the entire "nothing breaks" burden that a reviewer would otherwise share. That makes Phase 0 load-bearing. Characterization tests on `storm.py` I/O only prove the storm module reproduces its own outputs; they say nothing about the integration seams where forcing touches the rest of GeoClaw. Those seams are the general-functionality risk and each decision that would change them must be flagged, not made locally.

1. **Aux array contract.** Storm forcing owns aux slots for wind `(u, v)` and pressure, populated each step and consumed by `setaux`, `b4step2` (`set_storm_fields` / `set_met_forcing_fields`), `src2` (pressure-gradient and wind-drag source terms), output (`fort.a`), aux interpolation on regrid, and restart. **Invariant: aux count, component ordering, and indices unchanged.** Any unavoidable change is a flagged decision, not a local one.

2. **AMR refinement via storm center.** `flag2refine2` in the surge apps refines on distance-to-center (`R_refine`) alongside `wind_refine`. Optionalizing `storm_location` for the gridded case is therefore not free: it forces a decision about the gridded-forcing refinement criterion when there is no center (wind-speed threshold only, forcing-window bounding box, or other). This changes AMR behavior, which is general functionality. Pick it deliberately and flag it; do not let it fall out of the code.

   > **Decision (Phase F1):** gridded/no-center forcing refines on the
   > `wind_refine` wind-speed thresholds only; the distance-to-center `R_refine`
   > loop requires an eye and does not apply. This formalizes the pre-F1
   > behavior (`R_refine` was already inert for data storms because
   > `storm_location` returned `rinfinity`), so it is behavior-neutral. F1 makes
   > it intentional via a `storm_location_available` predicate (true iff
   > `storm_specification_type > 0`) that guards `flag2refine2` (R_refine),
   > `src2` (sector wind-drag angle → `phi = 0` when centerless), and `valout`
   > (`fort.track` is written for parameterized storms only).

3. **`src2` coupling.** Wind drag law selection and pressure forcing enter the momentum source terms alongside Manning friction. The drag pointer is a shared forcing control; if its selection or the source-term ordering moves, `src2` output moves. Every surge run is sensitive to this.

4. **`surge.data` format and `SurgeData` writer, in lockstep.** The `storm_specification_type` integer encoding lives on this wire. The Python writer and the Fortran reader change together or not at all. Keeping the legacy encoding as an alias is what protects this during the config split.

5. **`SurgeData` public surface and import paths.** Every surge `setrun.py` sets `SurgeData` attributes, and import paths (`from clawpack.geoclaw.surge.storm import Storm`, the `SurgeData` location) are themselves a compatibility contract. Module splits must keep those paths importable via shims. Treat the `SurgeData` surface and import paths as first-class must-keep, alongside the `Storm` class.

6. **Rename blast radius.** Every `use storm_module` / `use model_storm_module` / `use data_storm_module` site across the whole tree breaks on the rename, not just files under `surge/`. Enumerate them up front so the radius is known rather than discovered at link time.

7. **Non-shallow variants.** Verify whether multilayer SWE surge (different aux layout) and the Boussinesq path consume the storm aux contract. If either reads storm aux, the aux-contract invariant in seam 1 must hold for them too.

---

## 10. Development plan

### Phase 0 — verification harness (before any refactor)
Convert "best-effort compatibility" from a stance into a gate.
- **Contract audit** first: grep every `use` of the three Fortran modules tree-wide and every import of `SurgeData` and `clawpack.geoclaw.surge.storm`; report the current aux index/count contract and every routine that reads it (`setaux`, `b4step2`, `src2`, `flag2refine2`, output, restart); confirm whether multilayer and bouss consume it.
- **Python characterization tests:** snapshot object state from `read_geoclaw`, `read_atcf`, `read_hurdat`, `read_ibtracs`, `read_data`; golden-file the bytes from `write_geoclaw` and `write_data`, on a small fixture set. Assert identity through the entire refactor.
- **Fortran end-to-end regression:** one parametric (Holland 1980) and one gridded (netCDF) case that dumps the forcing aux arrays on a fixed grid at a few times, asserted bit- or tolerance-identical pre/post. The *Implement analytical storm surge specs* work is the natural verification oracle here and doubles as the guard on the gridded path.
- **Two-sided "nothing broke" check:**
  - a non-surge GeoClaw example (a tsunami/dtopo case) whose `fort.q`, gauges, and fgmax are asserted **byte-identical across every refactor commit** (storm work must have zero effect on it; the moment this canary moves, shared code has been touched);
  - the full surge example suite, with gauge/fgmax diffs held identical for the cases that should be preserved.

**Gate:** audit reported, baseline green, before any production code changes.

### Phase 1 — Python conceptual split (one PR, additive then swap)
Introduce `Track`/`StormTrack`/`ParametricMetForcing`/`GriddedMetForcing` alongside `Storm`; port readers to produce them; move `read_geoclaw -> ParametricMetForcing`, `read_data -> GriddedMetForcing`; put `write_geoclaw`/`write_data` on the objects with `time_offset` at the write boundary; reimplement `Storm` as a thin wrapper over `ParametricMetForcing`; preserve the must-keep surface and import paths via shims; demote stub writers and Python parametric stubs; split `available_models()`. Gate on the full Phase 0 suite.

### Phase 2 — config semantics (Python-first)
Replace the overloaded `storm_specification_type` at the API level with explicit `family` + `subtype`. Keep the legacy signed-integer encoding on the `surge.data` wire as an alias so current Fortran keeps reading it and `setrun.py` workflows do not break.

### Phase F1 — Fortran center-assumption fix (independent, early)
Guard/optionalize `storm_location` / `output_storm_location` so the gridded path no longer assumes a center. This is a real latent-bug fix, independent of the renames, and it de-risks gridded (ETC) runs. It carries the general AMR decision in seam 2, which must be surfaced and chosen here rather than at implementation time. Ships as its own small PR; changes no other logic and renames nothing.

### Phase F2 / F3 — Fortran interface split, then physical rename (last)
Split `storm_specification_type` into family+subtype at the Fortran config read (F2), then rename modules and types to the met-forcing vocabulary (F3). The rename is deliberately last and is required to be provably behavior-neutral: the end-to-end regressions and the non-surge canary must be identical across it, so any post-rename regression is definitionally not the rename.

### Phase 3+ — future work (not prerequisites)
Python-side parametric model implementations (a `StormTrack` to (wind, pressure) evaluator mirroring Fortran); ETC-specific parametric forcing on the non-tropical `StormTrack` path; Fortran track/state type separation. *Rework TC Geometry Estimates* (Kimball/Mulekar, Chavas radius; Emanuel/Tuleya precip) plugs in here as new parametric subtypes and radius-fill logic, riding on the clean parametric interface.

---

## 11. PR structure and rationale

Four separate PRs, in order: (1) Python split, (2) Fortran center-assumption fix, (3) config family/subtype split, (4) Fortran rename.

The reason is blast-radius isolation and bisectability, not review process. If a surge example regresses, `git bisect` should land on one conceptual change (the split, the config, or the rename) rather than a single commit that did all three. The rename PR in particular must be behavior-neutral, so isolating it means any post-rename regression is definitionally not the rename. Because the automated gate carries the full "nothing breaks" burden here, keeping the changes separable is what makes a regression attributable and cheap to revert.

---

## 12. Project gating

- **Model Alaskan ETCs:** the motivating case for the generalization, but not on the refactor's critical path. Its unblock was the netCDF merge (done) plus Derecho/`batch`. The refactor proceeds on its own branch without gating Alaska, provided Phase 0 guards the gridded path and Phase F1 lands the center-assumption fix (which helps it).
- **Model ETC Storm Surge:** gated on Derecho/`batch`; orthogonal to this refactor. Nothing here touches that path.
- **Run simulations for Atlantic risk assessment:** cares about `write_geoclaw` and Holland stability at scale; protected by Phase 0 characterization and the Phase 1 wrapper. Benefits later from Phase 3 Python parametric work and TC geometry.
- **Rework TC Geometry Estimates:** future parametric subtypes, not refactor structure; does not block this work.

---

## 13. Open questions (mostly settled; implementation-level)

- Exact target class names where still ambiguous (`Track` vs `TrackData`, `MetForcing` vs `BaseMetForcing`).
- New forcing file extensions/naming to distinguish parameterized forcing files from gridded descriptor files (a user-facing clarification opportunity).
- Degree of shared concrete base implementation for common forcing controls.
- `read_imd`: deprecate now or keep as a documented future-work placeholder.
- The gridded-forcing AMR refinement criterion (seam 2): the one genuinely general decision buried in Phase F1; decide it explicitly and early.
- Timing of Python object deprecations and file-naming transitions, weighed against test and documentation impact.
