# Topo and dtopo input recipes (clawpack/geoclaw#768)

*Status (2026-11-25): implemented on branch `topo-recipe-768`; Python suite
and byte-for-byte checks pass (see Verification). Target: v5.15, cut before AGU (week of
2026-12-07), and before PR B (`topo.data` format, roadmap §2.8). Branch
`topo-recipe-768`, off `master` at `5fc5e7cc` (#765 merged).*

## Problem

`Topography` and `DTopography` do two jobs:

- **Recipe.** A file path and type, plus the settings `setrun.py` uses to turn
  the file into what Fortran loads: `crop_extent`, `coarsen`, `buffer`,
  `align`, `x_shift`, `y_shift`, `z_shift`, `negate_z`.
- **Result.** The arrays `x`, `y`, `Z` that come out of the read.

Mixing them causes the confusions raised in #768 and its review threads:

- After `topo.crop(coarsen=2)` the returned object holds coarsened data while
  its `coarsen` attribute says 1. The attributes say nothing about how the data
  was made.
- After `read()` the shift has already been applied to `x`/`y`, but
  `x_shift` is still set. Writing the object out and listing the written file
  with the same shift would apply it twice.
- `extent` means the file's extent after `read_header()` and the loaded
  extent after `read()` (#769 had to pin this down).
- A misspelled setting (`topo.crop_exent = [...]`, Randy's #765 repro) is
  silently accepted as a new attribute. `crop_extent` stays `None`, which is
  written as `0 0 0 0` and read by Fortran as "the whole file".

None of these attributes are in v5.14.0 (checked with `git show v5.14.0:`),
and in v5.14.0 `topofiles`/`dtopofiles` held only `[topo_type, path]` lists.
So the attributes and the object form of `topofiles` can be moved now without
a deprecation cycle. After v5.15 they can't.

## Design

### One recipe object per file

| Class | Module | Fields |
|---|---|---|
| `TopoFile` | `topotools` | `path`, `topo_type`, `crop_extent`, `coarsen`, `buffer`, `align`, `x_shift`, `y_shift`, `z_shift`, `negate_z` |
| `DTopoFile` | `dtopotools` | `path`, `dtopo_type`, `crop_extent`, `coarsen`, `buffer`, `align`, `x_shift`, `y_shift`, `z_shift`, `negate_z` |

The two classes have the same fields (the type field is named as on
`Topography`/`DTopography`), the same defaults and the same validation.

- Two ways to build one, as with `Topography` (Randy, 2026-11-24): keywords at
  construction, `TopoFile(path, 4, crop_extent=[...])`, or attributes after
  `TopoFile()`. Every field has a default (`path` and the type default to
  `None`, and are required when the entry is written). The two give equal
  objects.
- A slotted dataclass: setting an unknown name such as `crop_exent` raises
  `AttributeError`, and an unknown keyword raises `TypeError`.
- Fields are validated on assignment (crop is 4 numbers with `x1 < x2`,
  `y1 < y2`, or `None` for the whole file; `coarsen >= 1`; `buffer >= 0`
  integer; `align` 2 numbers or `None`).
- `datum` (both classes, default `None`) is informational: the file's
  vertical datum. It is not passed to `read()` and not written to the `.data`
  file; `TopographyData.write`/`DTopoData.write` warn when listed files
  declare different datums, as they did for `Topography`/`DTopography`
  entries. `DTopoData.write` still reads a netCDF file's datum when the entry
  declares none, without changing the entry.
- `path` and the type may be `None` while an entry is being built; writing
  an entry with no `path` raises `ValueError`. A missing type is inferred from
  the file extension at write time, as before.
- The shared base class is `coordinate_tools.InputRecipe` (the format-neutral
  module that both `topotools` and `dtopotools` already import).
- `DTopoFile` has the crop/coarsen/buffer/align fields even though Fortran
  does not implement them for dtopo yet (mandli/geoclaw#15; today
  `read_dtopo_settings` stops if they are set). Python already applies them
  (`DTopography.read`, types 2/3/4), so `DTopoFile.read()` does too.
  `DTopoData.write` keeps raising `NotImplementedError` for them, naming #15,
  so `make data` fails before a run could. When #15 lands, only that guard and
  the Fortran stop are removed; the API and the `dtopo.data` layout stay as
  they are. Keeping the classes the same follows the rule that window names
  and logic apply to every input wherever they make sense.
- `TopoFile.read()` returns a `Topography` with the recipe applied: domain
  coordinates, the same points Fortran loads. `TopoFile.read_header()`
  returns the file's own header (file coordinates). `DTopoFile.read()` and
  `DTopoFile.read_header()` do the same with `DTopography`.
- Private, not user-facing: the netCDF metadata that `topo_entries()` attaches
  today through `Topography._netcdf_meta` and `[type, path, TopoMetadata]`
  entries moves to a private `TopoFile` field.

### `Topography` and `DTopography` become data + I/O only

- The eight recipe attributes are removed from both classes.
- `Topography.read()` keeps `crop_extent=`, `coarsen=`, `buffer=`, `align=`
  as one-off keyword arguments, and gains `x_shift=`, `y_shift=`, `z_shift=`,
  `negate_z=`. They act on that read and are not stored. `TopoFile.read()` is
  `Topography.read(path, topo_type, **recipe)`.
- `DTopography.read()` gains the same eight keywords. Today it takes none of
  them and reads the settings from the object's attributes.
  `DTopoFile.read()` is `DTopography.read(path, dtopo_type, **recipe)`.
- `crop()`, `write()`, `extent`, `delta` are unchanged.
- Setting one of the eight names on a `Topography`/`DTopography` raises
  `AttributeError` naming `read()` and the recipe class. Code written against
  `master` before this change would otherwise store the value as a plain
  attribute that `read()` silently ignores.
- `TopographyData` orders files by each entry's header spacing
  (`TopoFile.read_header().delta`), as it did for an unread `Topography`.
- `entry.read().write(path)` writes a file that already has the shift applied;
  a `TopoFile` for that new file has no shift. This settles the double-shift
  case in #768.

### What `topofiles` and `dtopofiles` accept

| Entry | Behavior |
|---|---|
| `TopoFile` / `DTopoFile` | Canonical. |
| `[topo_type, path]` | Kept, with no deprecation warning (decided 2026-11-24): it is the v5.14 form and every example `setrun.py` uses it. Converted to a `TopoFile` with defaults. The `DeprecationWarning` currently on `master` is removed. |
| `[topo_type, minlev, maxlev, t1, t2, path]` | Deprecated since v5.8, unchanged. |
| dict | Unreleased form on `master`; removed (`TypeError` naming `TopoFile`). |
| `Topography` / `DTopography` | Unreleased form on `master`; `TypeError` naming `TopoFile`. |

### Out of scope

- Met forcing: `GriddedMetForcing` (`met/gridded.py`) is already a recipe
  (path, `crop_extent`, shifts) and holds no arrays. Unchanged here. It gets
  the same checks in a follow-up that should land before v5.15, by sharing
  the checking part of `InputRecipe` as a mixin rather than becoming an entry
  class (`topo_input_roadmap.md` §2.9).
- Tile selection (#770). The 2026-11-24 meeting deferred it: the PR is to be
  converted to an issue and redesigned, because it loops over every topo file,
  including the coarse ones that should always be kept. A redesign should work
  on `TopoFile` entries.
- `read_header()` vs. `read()` crop extent for dtopo: Randy is opening a
  separate issue (2026-11-24 meeting).
- The `topo.data` format: PR B, stacked on this. The block writer
  (`data._write_preprocessing_block`) reads the same field names, so it
  works on a `TopoFile` unchanged here.

## Migration (users of `master` since v5.14)

```python
# before (master, unreleased)
topo = topotools.Topography()
topo.path = path
topo.topo_type = 4
topo.crop_extent = [-124.1, -123.9, 46.8, 46.9]
topo_data.topofiles.append(topo)

# after
topo_data.topofiles.append(
    topotools.TopoFile(path, 4, crop_extent=[-124.1, -123.9, 46.8, 46.9]))
```

## Verification

The refactor holds the data path fixed, so the evidence is equality with the
tree before it (`5fc5e7cc`):

- **Baseline comparison.** A script captured, on `5fc5e7cc`, 328 cases: topo
  types 2/3/4 crossed with crop (none, on-grid, off-grid), coarsen (1, 3),
  buffer (0, 2), align (none, set), shifts (none, set) and `z_shift`/
  `negate_z` (none, set); dtopo type 3 reads across the same settings; and
  `dtopo.data` for dtopo types 1/3 with shifts. For each it kept the loaded
  arrays and the `topo.data`/`dtopo.data` text. On this branch, through
  `TopoFile.read()`/`DTopoFile.read()` and the new writers: **0 of 328
  differ**, arrays compared with `np.array_equal` and text exactly. The
  comparison was checked to fail on injected changes (an `x` offset of
  1e-13, a relative `dZ` change of 1e-15, one extra space in the data file).
  The script is not committed.
- `tests/test_topo_data_golden.py` (byte-exact `topo.data`) passes with the
  goldens unchanged.
- `tests/test_topo_recipe.py` (new): both construction styles give equal
  entries, misspelled names and invalid values raise, `read()` equals
  `Topography.read`/`DTopography.read` with the same keywords, `read_header()`
  is the file, the list form does not warn, `Topography`/`DTopography`/dict
  entries raise, and the data classes refuse the eight names.
- Full Python suite: 795 passed, 1 xfailed, on `443da050` (no tests
  regenerated or deleted; tests that set the settings as attributes now pass
  them as `read()` keywords or build `TopoFile`/`DTopoFile` entries).
- Fortran examples and regressions: see Progress.

## Progress

- [x] Design agreed (this doc; names `TopoFile`/`DTopoFile` and keeping the
      list form confirmed by Kyle 2026-11-24)
- [x] Failing tests (`tests/test_topo_recipe.py`, `f919a755`)
- [x] `TopoFile`, `DTopoFile`; `Topography`/`DTopography` attributes removed
- [x] `data.py`, `gridded_input`, tests and examples moved to the recipe
      objects (`topo_entries()` lists are still accepted)
- [x] Byte-identical `topo.data`/`dtopo.data` check (328 cases, above)
- [x] Roadmap §1 rows and §2.7/§2.8 notes updated
- [x] Fortran example/regression run (macOS, gfortran): bowl-slosh
      preprocessing, dtopo, bowl-slosh-netcdf, isaac, chile2010 and
      `tests/regression`: 67 passed, 1 xfailed, 3 deselected (marker-gated),
      no regression data regenerated
- [ ] User docs: the `topochanges` doc branch (rjleveque/doc#5) still shows
      the attribute form in `setrun_geoclaw.rst`, `topotools_module.rst` and
      `topodata_format.rst`
- [ ] clawpack/geoclaw#764 (kmltools) reads `crop_extent`/shifts from
      `Topography` entries; it needs rebasing onto `TopoFile`
