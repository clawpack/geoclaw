# GeoClaw units policy

**Status:** stated and enforced as of PR B1. The table below is generated from
`UNITS_POLICY` in `src/python/geoclaw/units.py` and verified against the real
readers by `tests/test_units_policy.py`.

## Why this document exists

GeoClaw already had a coherent units policy. It was implemented in the NetCDF
readers, it was well reasoned, and it was written down nowhere so enforcement
drifted, the override argument acquired four different spellings, and two
docstrings ended up claiming the opposite of what the code does. Users could not
discover any of it, and neither could reviewers.

Writing the rules down is only half of it. Prose drifts from code, and a table
in a document is exactly the kind of thing that quietly stops being true. So the
table here is *generated* from a registry, and every row is *executed* against
the real reader. A row cannot claim behavior the code does not have, and the
document cannot disagree with the registry.

## The rules

1. **Units must be declared in the file.** They are never silently assumed. A
   file with no declaration is an error, not a guess.
2. **Overriding is explicit.** `assume_units` (NetCDF) and `input_units`
   (subfault files) mean "treat the file as if it had declared this". They are
   the only way to supply units GeoClaw cannot read from the file.
3. **A recognized non-contract unit is converted, and the conversion is
   announced.** Reading a file in km is fine; doing it silently is not.
4. **An unrecognized unit raises.** GeoClaw does not guess at unit strings it
   does not know.
5. **After conversion, magnitude is sanity-checked.** Units can be declared
   *wrongly*, and rules 1-4 cannot catch that. `_check_magnitude` does, within
   limits: it auto-corrects only the one unambiguous case (pressure ~1000x low
   is unmistakably hPa/mbar) and raises on everything else, because feet vs
   meters and knots vs m/s cannot be told apart by magnitude alone.

Rule 5 is the reason the policy is not simply "trust the declaration". Rules 1-4
protect against *missing* information; rule 5 protects against *wrong*
information, which is the more common failure in practice.

## ASCII formats cannot declare units at all

This is the question the policy is least obvious about, so it is worth stating
plainly. A `topo_type=2/3` header is:

```
ncols / nrows / xlower / ylower / cellsize / nodata_value
```

There is no units field, and none is proposed here. The same is true of ASCII
dtopo. So for these formats:

- **Rule 1 cannot apply.** There is no declaration to require, and refusing to
  read an undeclared file would mean refusing every ASCII file GeoClaw has ever
  read.
- **The contract unit is assumed** -- meters for elevation and deformation.
- **What makes that safe is rule 5, not rule 1.** The magnitude check is the
  only defense available, and it is doing the whole job on its own here.
- **Rule 2 is how you say otherwise**: `assume_units` states the file's real
  unit, and the value is converted.

"Assumed, then checked" is the behavior, and the two ASCII rows below now
conform. `Topography.read` and `DTopography.read` both take `assume_units`, and
unit resolution happens once in the product rather than per format, so ASCII and
NetCDF go down the same path.

**The two ASCII rows differ from each other, and the table shows it.** Topo
ASCII is `assume+check`: meters is assumed and rule 5 catches a gross error.
Dtopo ASCII is `warn+assume`, because `_check_magnitude` defines **no bound for
deformation** -- there is no rule-5 backstop there, so the assumption is
announced instead. Pass `assume_units='m'` to state it deliberately and silence
the warning.

Defining a deformation bound would let dtopo ASCII reach `assume+check` too. It
needs a defensible physical limit on coseismic displacement, which is a
judgement about earthquake physics rather than about code, so it is left open
rather than guessed at.

## Contract units

GeoClaw works internally in SI: elevation and deformation in meters, wind in
m/s, pressure in Pa, time in seconds. See `GEOCLAW_NETCDF_UNITS` in
`src/python/geoclaw/units.py`.

## What each input path actually does

Regenerate this table with `GEOCLAW_REGEN=1 pytest tests/test_units_policy.py`.
Do not edit it by hand -- edit `UNITS_POLICY` instead.

<!-- BEGIN GENERATED TABLE -->
| Path | Contract | Declared in file | Override | Missing | Non-contract | Unrecognized | Magnitude | Conforms |
|---|---|---|---|---|---|---|---|---|
| `Topography.read (topo_type=4)` | m | yes | assume_units (str) | raise | convert+warn | raise | yes | yes |
| `Topography.read (topo_type=1,2,3)` | m | no | assume_units (str) | assume+check | convert | raise | yes | yes |
| `DTopoInspector (deformation)` | m | yes | assume_units (str) | raise | convert+warn | raise | no | yes |
| `DTopoInspector (time axis)` | s | yes | none | warn+assume | convert | raise | no | yes |
| `DTopography.read (dtopo_type=1,2,3)` | m | no | assume_units (str) | warn+assume | convert | raise | no | yes |
| `MetInspector (wind, pressure)` | m/s, Pa | yes | assume_units (dict; bool deprecated), format_units (dict) | raise | convert+warn | raise | yes | yes |
| `Fault.read (columnar subfault files)` | m, Pa | no | input_units (dict) | warn+assume | convert | raise | no | yes |
| `CSVFault.read (units in column headings)` | m, Pa | yes | input_units (dict), overrides the heading | warn+assume | convert | raise | no | yes |
<!-- END GENERATED TABLE -->

A row marked **no** in the *Conforms* column is a known hole: the policy says
one thing and that reader does another. Its conformance test is marked
`xfail(strict=True)`, so it is visible in every test run and the marker cannot
be left behind once the hole is closed.

## Known deliberate looseness

- **Unit-string matching is a hand-maintained table**, not dimensional
  analysis (see `_CONTRACT_UNIT_CF_ALIASES` and `_CF_TO_UNITS_PY` in
  `netcdf_utils.py`). It is exact-string and case-sensitive, so `meters` is
  accepted and `Meters` is not, and `feet` has no conversion at all. Using
  `cf-units`/udunits instead would be correct but adds a dependency.
- **`MetInspector.assume_units` is a `bool` where the other inspectors take a
  `str`.** This is deliberate, not an oversight: Met covers several variables
  with *different* contract units, so a single string is meaningless. `True`
  means "each variable is already in its own contract unit"; the per-role string
  form is `format_units={role: unit}`, which takes precedence.

## Definition of done

Any change to unit handling must update `UNITS_POLICY` in the same PR. The
generated table and the conformance test will fail otherwise, which is the
point.
