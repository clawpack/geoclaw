# Met-Forcing Documentation, Transition Guide & Tutorial — Plan

**Status (2026-07-25):** Approved, execution deferred until the code work is fully
wrapped up. Companion to `met_forcing_refactor.md`. The refactor code is merged
(Phase 0/1/F1, Phase 2 = PR #710, F3 = PR #711, plus the lockstep `clawpack/riemann`
rename PR).

## Context

The met-forcing refactor reorganized GeoClaw storm forcing under a
**meteorological forcing** umbrella and added gridded netCDF met forcing. The
user-facing docs did not keep up. This effort (a) updates the sparse surge docs,
(b) adds a transition / "what's new" guide advertising the new capabilities and
easing migration, and (c) adds a capability tutorial notebook.

**Assessment of current docs.** GeoClaw's Sphinx docs live in the separate
`clawpack/doc` repo (`doc/doc/*.rst`, reST-only — no nbsphinx/MyST). The surge
cluster is thin:
- `quick_surge.rst` — a labeled stub (steps 2–7 empty, a dead link).
- New object model (`surge/track.py`, `parametric.py`, `gridded.py`, `tools.py`)
  and `SurgeData`'s new `storm_family`/`storm_subtype`, `t_ramp_on/off`,
  `rotation_override` — good docstrings/comments in code, **zero Sphinx surfacing**.
- `storm_module.rst` (lone surge autodoc page) is **orphaned** from the toctree;
  broken `:ref:` targets (`_storm_module`, `_surge_module`).
- `setrun_geoclaw.rst` surge section documents only the legacy integer
  `storm_specification_type` (−1..3); omits the newer attributes.
- Strong gridded-forcing content is **buried** in the general `netcdf.rst`.
- `examples/storm-surge/isaac/README.rst` is 3 stale lines despite its OWI/netCDF role.

**Decisions locked:** transition/"what's new" guide → new Sphinx surge page in
`clawpack/doc`; tutorial → notebook in `clawpack/apps` (`apps/notebooks/geoclaw/`),
covering **both** parametric and gridded netCDF forcing.

## Work item A — `clawpack/geoclaw` (docstrings & examples)
- `src/python/geoclaw/data.py`: expand the one-line `SurgeData` docstring into a
  full class docstring (attributes; `storm_family`/`storm_subtype` selection with
  `storm_specification_type` as the retained legacy path) so autodoc renders well.
- `examples/storm-surge/isaac/README.rst`: rewrite to describe the Holland /
  OWI-ASCII / netCDF (ERA5, NWS13) variants it exercises.
- Lightly annotate the `ike`/`isaac` `setrun.py` surge blocks to show the new
  `storm_family`/`storm_subtype` form (legacy string still works — do not break
  isaac's `== 'holland80'` comparisons).
- Reuse `clawpack.geoclaw.surge.plot` (`surgeplot`) helpers; do not rewrite plotting.

## Work item B — `clawpack/doc` (Sphinx surge cluster overhaul) — centerpiece
All under `doc/doc/`:
- Rewrite `quick_surge.rst` into a real quick start (end-to-end parametric Holland
  via the `ike` example, then a gridded pointer); drop the dead link and WIP warning.
- Update `setrun_geoclaw.rst` "Storm Specification Data": add `storm_family`/
  `storm_subtype`, `t_ramp_on`/`t_ramp_off`, `rotation_override`; reframe
  `storm_specification_type` as the supported legacy selector mapping to family/subtype.
- Add autodoc pages for `surge.track`, `surge.parametric`, `surge.gridded`,
  `surge.tools`, and `SurgeData`; wire them and the orphaned `storm_module.rst`
  into the `geoclaw.rst` toctree (model on dtopotools/fgmax_tools pages).
- Fix the broken `:ref:` targets (`_storm_module`, `_surge_module`).
- Surface gridded met forcing into the surge cluster (dedicated
  `gridded_met_forcing.rst` or promote the `netcdf.rst` section) so it is discoverable.
- New "What's new in met forcing" + migration page (`met_forcing_whatsnew.rst`) in
  the surge cluster + toctree. Advertises the object model
  (`Track`/`StormTrack`/`ParametricMetForcing`/`GriddedMetForcing`), explicit
  `family`+`subtype` selection, gridded netCDF forcing, temporal ramps. Migration
  notes: legacy `storm_specification_type` still works; the Fortran
  `storm_module` → `met_forcing_module` rename affects only user code that `use`s it
  (Python `setrun.py` unaffected). Source prose from `met_forcing_refactor.md` §7/§9.

## Work item C — `clawpack/apps` tutorial notebook (stretch)
- New `apps/notebooks/geoclaw/met_forcing_tutorial.ipynb`, complementing the existing
  `katrina/katrina.ipynb` and `data_derived_storms.ipynb`.
- Follow the `nbtools` convention (`clawpack.clawutil.nbtools.make_exe` /
  `make_output_and_plots`), as in `examples/tsunami/eta_init_force_dry/run_geoclaw.ipynb`:
  markdown intro (+ gallery link, version note), then
  - **Parametric:** read an ATCF track (`Storm.read(file_format="ATCF")`), write the
    GeoClaw file, build+run, plot wind/pressure/surface + gauges via `surgeplot`.
  - **Gridded netCDF:** build a small CF `u10`/`v10`/`msl` dataset (mirroring
    `tests/regression/met_forcing/`), write the `.storm` descriptor
    (`write(file_format="data")`), build with `-DNETCDF`, run, plot.
- The apps Makefile renders via `jupyter nbconvert --to html --execute` (1200 s
  timeout), so the notebook is executed at render time — it must run top-to-bottom,
  including a netCDF-enabled `xgeoclaw` for the gridded part.
- Add a gallery link in `doc/gallery/notebooks.rst` to the rendered HTML.

## Sequencing & PRs
Three coordinated PRs (one per repo):
1. `geoclaw` docstring/README/setrun polish (A).
2. `doc` surge-cluster overhaul + what's-new page (B) — centerpiece.
3. `apps` tutorial notebook (C) + `doc/gallery` link.

## Verification
- **doc build:** `cd $CLAW/doc/doc && make html` — new autodoc pages import against
  the installed `clawpack.geoclaw`; no new Sphinx warnings; broken `:ref:`s resolve;
  new pages appear in the nav.
- **geoclaw examples:** `ike`/`isaac` still `make all` and pass `test_ike.py` /
  `test_isaac.py` after the setrun annotations (legacy selection unchanged).
- **notebook:** `make met_forcing_tutorial.html -f $CLAW/apps/notebooks/Makefile`
  executes end-to-end (needs a `-DNETCDF` build); review the rendered figures.
  Not CI-executed; the tested surrogate remains `tests/regression/met_forcing/`.
