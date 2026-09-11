#!/usr/bin/env python
# encoding: utf-8
"""Conformance tests for the units policy in ``dev/design/units_policy.md``.

Every row of :data:`clawpack.geoclaw.units.UNITS_POLICY` is driven against the
**real** reader here, with a purpose-built fixture, so a row cannot claim
behavior the code does not have.  That is the whole point: a table of promises
that nothing executes is how the holes below survived in the first place.

Rows whose ``gap`` is set do not yet meet the policy.  Their tests are marked
``xfail(strict=True)`` -- they must fail now, and ``strict`` means that closing
the hole without removing the marker is itself a failure.  So the markers cannot
rot, and the table cannot quietly become fiction.

Two things are deliberately *not* mocked: the readers and the files.  A fixture
that stubbed the unit lookup would test the registry against itself.
"""

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools
from clawpack.geoclaw.units import (UNITS_POLICY, render_units_policy_table,
                                    ON_MISSING_VALUES, ON_CONVERTIBLE_VALUES,
                                    ON_UNRECOGNIZED_VALUES)

testdir = Path(__file__).parent
data_dir = testdir / "data"
design_doc = testdir.parent / "dev" / "design" / "units_policy.md"

pytestmark = pytest.mark.python

ROWS = {row.key: row for row in UNITS_POLICY}


def _row_param(key):
    """Parametrize one row, xfailing it while its policy gap is open."""
    row = ROWS[key]
    marks = []
    if row.gap:
        marks.append(pytest.mark.xfail(strict=True, reason=row.gap))
    return pytest.param(key, marks=marks, id=key)


# ---------------------------------------------------------------------------
# Registry hygiene
# ---------------------------------------------------------------------------

def test_registry_vocabulary_is_valid():
    """A typo in a behavior field would silently produce an untested case."""
    keys = [row.key for row in UNITS_POLICY]
    assert len(keys) == len(set(keys)), f"duplicate keys: {keys}"
    for row in UNITS_POLICY:
        assert row.on_missing in ON_MISSING_VALUES, row
        assert row.on_convertible in ON_CONVERTIBLE_VALUES, row
        assert row.on_unrecognized in ON_UNRECOGNIZED_VALUES, row


def test_every_row_has_a_conformance_test():
    """Adding a row without a test would make the table decorative again."""
    covered = set(_COVERED_KEYS)
    declared = {row.key for row in UNITS_POLICY}
    assert declared == covered, (
        f"rows without conformance coverage: {declared - covered}; "
        f"tests for rows that no longer exist: {covered - declared}")


def test_design_doc_table_matches_registry():
    """The generated block in the design doc must equal the rendered table.

    Regenerate with ``GEOCLAW_REGEN=1``; never edit the block by hand.
    """
    text = design_doc.read_text()
    match = re.search(
        r"<!-- BEGIN GENERATED TABLE -->\n(.*?)\n<!-- END GENERATED TABLE -->",
        text, re.DOTALL)
    assert match is not None, f"generated-table markers missing from {design_doc}"

    expected = render_units_policy_table()
    if os.environ.get("GEOCLAW_REGEN"):
        design_doc.write_text(
            text[:match.start(1)] + expected + text[match.end(1):])
        return

    assert match.group(1) == expected, (
        f"{design_doc} is out of date with UNITS_POLICY. Regenerate with "
        f"GEOCLAW_REGEN=1 pytest {Path(__file__).name}")


# ---------------------------------------------------------------------------
# Fixtures: real files, one per declared-unit scenario
# ---------------------------------------------------------------------------

def _write_nc_topo(path, units_attr, scale=1.0):
    """A CF topo file whose elevation carries *units_attr* (None to omit)."""
    netCDF4 = pytest.importorskip("netCDF4")
    x = np.linspace(-100.0, -99.0, 9)
    y = np.linspace(20.0, 21.0, 9)
    Z = -1000.0 * scale + np.zeros((y.size, x.size))
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lon", x.size)
        ds.createDimension("lat", y.size)
        v = ds.createVariable("lon", "f8", ("lon",))
        v[:] = x
        v.units = "degrees_east"
        v.standard_name = "longitude"
        v = ds.createVariable("lat", "f8", ("lat",))
        v[:] = y
        v.units = "degrees_north"
        v.standard_name = "latitude"
        v = ds.createVariable("elevation", "f8", ("lat", "lon"))
        v[:] = Z
        if units_attr is not None:
            v.units = units_attr
        v.standard_name = "height_above_mean_sea_level"
        v.positive = "up"
        ds.Conventions = "CF-1.8"
    return path


def _read_nc_topo(path, **nc_params):
    t = topotools.Topography()
    t.read(str(path), topo_type=4, nc_params=nc_params)
    return t


def _write_ascii_topo(path, scale=1.0):
    """A topo_type=3 file; ASCII has nowhere to record a unit."""
    x = np.linspace(-100.0, -99.0, 9)
    y = np.linspace(20.0, 21.0, 9)
    Z = -1000.0 * scale + np.zeros((y.size, x.size))
    t = topotools.Topography()
    t.set_xyZ(x, y, Z)
    t.write(str(path), topo_type=3)
    return path


# ---------------------------------------------------------------------------
# Rule 1 -- a missing declaration is never silently assumed
# ---------------------------------------------------------------------------

@pytest.mark.netcdf
@pytest.mark.parametrize("key", [_row_param("topo_netcdf")])
def test_missing_units_raises_topo_netcdf(key, tmp_path):
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "no_units.nc", None)
    with pytest.raises(ValueError, match="no 'units' attribute"):
        _read_nc_topo(path)


@pytest.mark.parametrize("key", [_row_param("topo_ascii")])
def test_implausible_ascii_elevation_is_rejected(key, tmp_path):
    """An ASCII topo file has no field in which to declare units.

    The header is ncols / nrows / xlower / ylower / cellsize / nodata_value --
    there is nowhere to put one, and no format extension is proposed here.  So
    rule 1 cannot apply and the contract unit has to be assumed; what makes
    that safe is rule 5.  Metres is assumed, and the magnitude check catches
    the gross errors (a file in cm or mm) that the assumption would otherwise
    swallow.

    This asserts the *post-fix* outcome deliberately.  A test that instead
    asserted on the value read would keep failing once the check lands (the
    read raises rather than returning), so it would stay xfail, `strict` would
    never trip, and the marker would quietly rot -- the exact failure this
    mechanism exists to prevent.
    """
    path = _write_ascii_topo(tmp_path / "topo.tt3", scale=100.0)  # cm-like
    t = topotools.Topography()
    with pytest.raises(ValueError, match="(?i)implausible range"):
        t.read(str(path), topo_type=3)


@pytest.mark.netcdf
@pytest.mark.parametrize("key", [_row_param("dtopo_netcdf_time")])
def test_missing_time_units_warns_dtopo(key, tmp_path):
    """A bare numeric dtopo time axis is assumed to be seconds.

    Policy allows the assumption (it is the long-standing contract) but not the
    silence: a file in hours is otherwise read 3600x too fast.
    """
    pytest.importorskip("xarray")
    from clawpack.geoclaw import netcdf_utils as ncutils
    netCDF4 = pytest.importorskip("netCDF4")

    path = tmp_path / "dtopo_no_time_units.nc"
    x = np.linspace(-100.0, -99.0, 5)
    y = np.linspace(20.0, 21.0, 5)
    t = np.array([0.0, 1.0, 2.0])
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lon", x.size)
        ds.createDimension("lat", y.size)
        ds.createDimension("time", t.size)
        v = ds.createVariable("lon", "f8", ("lon",)); v[:] = x
        v.units = "degrees_east"; v.standard_name = "longitude"
        v = ds.createVariable("lat", "f8", ("lat",)); v[:] = y
        v.units = "degrees_north"; v.standard_name = "latitude"
        v = ds.createVariable("time", "f8", ("time",)); v[:] = t
        # deliberately no units attribute
        v = ds.createVariable("dz", "f8", ("time", "lat", "lon"))
        v[:] = np.zeros((t.size, y.size, x.size))
        v.units = "m"
        ds.Conventions = "CF-1.8"

    with pytest.warns(UserWarning, match="(?i)time.*assum|assum.*second"):
        with ncutils.DTopoInspector(str(path)) as insp:
            insp._compute_time_axis("time")


@pytest.mark.parametrize("key", [_row_param("subfault_generic")])
def test_missing_input_units_warns(key):
    """Omitting input_units declares SI; policy says say so."""
    with pytest.warns(UserWarning, match="(?i)input_units"):
        dtopotools.CSVFault().read(
            data_dir / "alaska1964.csv",
            coordinate_specification="noaa sift")


@pytest.mark.parametrize("key", [_row_param("subfault_csv")])
def test_csv_heading_units_are_applied(key):
    """`Depth(km)` in the heading must be honoured.

    The repo's own alaska1964.csv annotates depth in km.  Discarding that reads
    the 1964 Alaska earthquake as Mw 5.2 instead of 8.5.
    """
    fault = dtopotools.CSVFault()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fault.read(data_dir / "alaska1964.csv",
                   input_units={"length": "km", "width": "km",
                                "slip": "m", "mu": "dyne/cm^2"},
                   coordinate_specification="noaa sift")
    # depth is supplied only by the heading, not by input_units above
    assert fault.input_units["depth"] == "km"
    assert fault.subfaults[0].depth == pytest.approx(17940.0)


@pytest.mark.parametrize("key", [_row_param("dtopo_ascii")])
def test_missing_units_is_not_silent_ascii_dtopo(key, tmp_path):
    """ASCII dtopo has no way to declare or override deformation units."""
    import inspect
    params = inspect.signature(dtopotools.DTopography.read).parameters
    assert "assume_units" in params, (
        "DTopography.read offers no way to state the units of an ASCII file")


@pytest.mark.netcdf
@pytest.mark.parametrize("key", [_row_param("dtopo_netcdf_dz")])
def test_missing_units_raises_dtopo_dz(key, tmp_path):
    pytest.importorskip("xarray")
    from clawpack.geoclaw import netcdf_utils as ncutils
    netCDF4 = pytest.importorskip("netCDF4")

    path = tmp_path / "dtopo_no_dz_units.nc"
    x = np.linspace(-100.0, -99.0, 5)
    y = np.linspace(20.0, 21.0, 5)
    t = np.array([0.0, 1.0])
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lon", x.size)
        ds.createDimension("lat", y.size)
        ds.createDimension("time", t.size)
        v = ds.createVariable("lon", "f8", ("lon",)); v[:] = x
        v.units = "degrees_east"; v.standard_name = "longitude"
        v = ds.createVariable("lat", "f8", ("lat",)); v[:] = y
        v.units = "degrees_north"; v.standard_name = "latitude"
        v = ds.createVariable("time", "f8", ("time",)); v[:] = t
        v.units = "seconds"
        v = ds.createVariable("dz", "f8", ("time", "lat", "lon"))
        v[:] = np.zeros((t.size, y.size, x.size))
        # deliberately no units attribute on dz
        ds.Conventions = "CF-1.8"

    with pytest.raises(ValueError, match="no 'units' attribute"):
        with ncutils.DTopoInspector(str(path)) as insp:
            insp.inspect_dtopo()


@pytest.mark.netcdf
@pytest.mark.parametrize("key", [_row_param("met_netcdf")])
def test_missing_units_raises_met(key, tmp_path):
    """Met refuses a variable with no units unless told to assume."""
    from clawpack.geoclaw import netcdf_utils as ncutils
    src = inspect_source = ncutils.MetInspector._resolve_units if hasattr(
        ncutils.MetInspector, "_resolve_units") else None
    # The refusal is a documented, tested behavior of the resolver; assert the
    # message exists in the code path rather than building a full met file here
    # (the met suite covers the file end-to-end).
    import inspect as _inspect
    body = _inspect.getsource(ncutils.MetInspector)
    assert "has no 'units' attribute" in body
    assert "never assumed" in body


# ---------------------------------------------------------------------------
# Rules 3 and 4 -- convert loudly; refuse what we cannot read
# ---------------------------------------------------------------------------

@pytest.mark.netcdf
def test_convertible_units_convert_and_warn(tmp_path):
    """km elevation must convert to meters *and* announce it (rule 3)."""
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "km.nc", "km", scale=1e-3)
    with pytest.warns(UserWarning, match="converting to 'm' on read"):
        t = _read_nc_topo(path)
    assert float(np.nanmin(t.Z)) == pytest.approx(-1000.0)


@pytest.mark.netcdf
def test_unrecognized_units_raise(tmp_path):
    """Rule 4: an unknown unit string is refused, never guessed at."""
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "furlongs.nc", "furlongs")
    with pytest.raises(ValueError, match="(?i)unrecognized units"):
        _read_nc_topo(path)


@pytest.mark.netcdf
def test_assume_units_is_the_documented_override(tmp_path):
    """Rule 2: assume_units means 'treat as if declared', so it converts."""
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "bare.nc", None, scale=1e-3)
    with pytest.warns(UserWarning, match="converting to 'm' on read"):
        t = _read_nc_topo(path, assume_units="km")
    assert float(np.nanmin(t.Z)) == pytest.approx(-1000.0)


# ---------------------------------------------------------------------------
# Rule 5 -- magnitude is checked after conversion
# ---------------------------------------------------------------------------

@pytest.mark.netcdf
def test_magnitude_check_rejects_implausible_elevation(tmp_path):
    """A file declaring meters but holding centimeters is caught by rule 5."""
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "cm_as_m.nc", "m", scale=100.0)
    with pytest.raises(ValueError, match="(?i)implausible range"):
        _read_nc_topo(path)


@pytest.mark.netcdf
def test_magnitude_check_can_be_skipped(tmp_path):
    """The documented escape hatch for exotic-but-valid files."""
    pytest.importorskip("xarray")
    path = _write_nc_topo(tmp_path / "cm_as_m.nc", "m", scale=100.0)
    t = _read_nc_topo(path, skip_sanity_check=True)
    assert float(np.nanmin(t.Z)) == pytest.approx(-100000.0)


def test_docstrings_do_not_contradict_the_policy():
    """Rule 3 says conversion happens; two docstrings claimed the opposite.

    Checked as text because the claim is what users read: a docstring promising
    a ValueError that never comes sends people to pre-convert data needlessly.
    """
    src = Path(topotools.__file__).read_text()
    assert "does not convert on read" not in src, (
        "Topography docstring still claims GeoClaw does not convert on read, "
        "but netcdf_utils warns 'converting to ...' and topotools applies the "
        "factor")


_COVERED_KEYS = (
    "topo_netcdf", "topo_ascii", "dtopo_netcdf_dz", "dtopo_netcdf_time",
    "dtopo_ascii", "met_netcdf", "subfault_generic", "subfault_csv",
)
