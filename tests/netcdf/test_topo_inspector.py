"""
Tests for TopoInspector.

Covers: fill-in-crop detection, fill-outside-crop pass, unit verification,
unit-conversion warn path, missing units warn path, descriptor key presence,
all coordinate variant combinations, and topo_entries() lon wrapping.
"""
import warnings

import pytest

from clawpack.geoclaw.netcdf_utils import TopoInspector, TopoMetadata

from ._helpers import (
    make_topo_dataset,
    COORD_VARIANTS,
    TOPO_CF_STANDARD_NAME_VARIANTS,
    TOPO_DIM_ORDER_VARIANTS,
    TOPO_FALLBACK_NAME_VARIANTS,
)

pytestmark = [pytest.mark.python, pytest.mark.netcdf]


# ============================================================
# Fill value in / outside crop region
# ============================================================

def test_fill_in_crop_region_raises(topo_file_factory):
    """
    TopoInspector raises ValueError when NaN is present within the
    crop region.  Silent NaN in bathymetry would corrupt simulations.
    """
    path = topo_file_factory(
        lon_min=-10.0, lon_max=10.0,
        lat_min=-10.0, lat_max=10.0,
        has_fill_in_data=True,  # NaN placed at data.flat[0]
    )
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-10.0, 10.0, -10.0, 10.0))
    with pytest.raises(ValueError, match="[Ff]ill"):
        insp.inspect_topo()
    insp.close()


def test_fill_outside_crop_region_passes(topo_file_factory):
    """
    TopoInspector passes when NaN exists only outside the crop region.
    """
    # The dataset has NaN at data.flat[0] which corresponds to the corner
    # of the full grid.  We crop to the interior so the NaN is outside.
    ds = make_topo_dataset(
        lon_min=-10.0, lon_max=10.0,
        lat_min=-10.0, lat_max=10.0,
        nlat=10, nlon=10,
        has_fill_in_data=True,  # NaN is at lat[-10], lon[-10] (S-W corner)
        lat_direction="S_to_N",
    )
    # The NaN is at the south-west corner; crop to the north-east interior.
    path = topo_file_factory(ds=ds)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(5.0, 10.0, 5.0, 10.0))
    meta = insp.inspect_topo()   # must not raise
    insp.close()
    assert isinstance(meta, TopoMetadata)


def test_no_fill_no_crop_passes(topo_file_factory):
    """TopoInspector passes when there are no NaN values at all."""
    path = topo_file_factory(has_fill_in_data=False)
    with TopoInspector(path, var_name="z") as insp:
        meta = insp.inspect_topo()
    assert isinstance(meta, TopoMetadata)


# ============================================================
# Unit verification
# ============================================================

def test_units_meters_passes_without_warning(topo_file_factory):
    """
    A variable whose units attribute is 'm' (or any meters alias) passes
    unit verification without emitting a warning.
    """
    path = topo_file_factory(units="m")
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # any warning is a test failure
        with TopoInspector(path, var_name="z") as insp:
            meta = insp.inspect_topo()
    assert meta.source_units in {"m", "meter", "meters", "metre", "metres"}


def test_units_convertible_raises(topo_file_factory):
    """
    A recognised but non-meter unit (e.g. 'cm') is rejected, not silently
    misread: the read paths do not convert, so the file must be pre-converted
    to meters.  (Unit mixing must never be silently assumed.)
    """
    path = topo_file_factory(units="cm")
    with TopoInspector(path, var_name="z") as insp:
        with pytest.raises(ValueError, match="cm"):
            insp.inspect_topo()


def test_units_unrecognised_raises(topo_file_factory):
    """
    A variable with units not in the known conversion table raises ValueError.
    """
    path = topo_file_factory(units="furlongs")
    with TopoInspector(path, var_name="z") as insp:
        with pytest.raises(ValueError, match="[Uu]nrecogni"):
            insp.inspect_topo()


def test_units_missing_raises(topo_file_factory):
    """
    A variable with no units attribute raises ValueError -- units are required
    and never silently assumed to be meters.
    """
    path = topo_file_factory(units="")   # empty string → no 'units' attr written
    with TopoInspector(path, var_name="z") as insp:
        with pytest.raises(ValueError, match="no 'units'"):
            insp.inspect_topo()


def test_units_missing_assume_units_opt_in(topo_file_factory):
    """
    The explicit assume_units escape hatch lets a caller declare the unit for
    a file that has no units attribute, without raising.
    """
    path = topo_file_factory(units="")
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # opt-in must be silent, not warn
        with TopoInspector(path, var_name="z", assume_units="m") as insp:
            meta = insp.inspect_topo()
    assert meta.source_units == "m"


# ============================================================
# Descriptor output keys
# ============================================================

def test_inspect_topo_returns_topo_metadata(topo_file_factory):
    """inspect_topo() returns a TopoMetadata instance."""
    path = topo_file_factory()
    with TopoInspector(path, var_name="z") as insp:
        meta = insp.inspect_topo()
    assert isinstance(meta, TopoMetadata)


def test_topo_metadata_required_fields_present(topo_file_factory):
    """
    TopoMetadata contains all required fields with valid values after
    inspection of a well-formed file.
    """
    path = topo_file_factory(
        lon_min=-10.0, lon_max=10.0,
        lat_min=-10.0, lat_max=10.0,
        units="m",
    )
    # crop_bounds is passed to TopoInspector, not to the dataset factory.
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-5.0, 5.0, -5.0, 5.0))
    meta = insp.inspect_topo()
    insp.close()

    assert meta.var_name == "z"
    assert meta.x_name is not None and len(meta.x_name) > 0
    assert meta.y_name is not None and len(meta.y_name) > 0
    assert meta.lon_wrap in (180, 360)
    assert meta.y_increasing in (True, False)
    assert len(meta.dim_order) == 2
    assert meta.fill_action == "abort"
    assert meta.source_units == "m"
    assert meta.crop_bounds == (-5.0, 5.0, -5.0, 5.0)


def test_topo_metadata_crop_bounds_none_when_not_specified(topo_file_factory):
    """When no crop_bounds are given, TopoMetadata.crop_bounds is None."""
    path = topo_file_factory()
    with TopoInspector(path, var_name="z") as insp:
        meta = insp.inspect_topo()
    assert meta.crop_bounds is None


# ============================================================
# All coordinate variant combinations
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants(topo_file_factory, coord_kwargs):
    """
    inspect_topo() completes without error for every coordinate variant
    and reports consistent lon_wrap / y_increasing values.
    """
    path = topo_file_factory(**coord_kwargs)
    with TopoInspector(path, var_name="z") as insp:
        meta = insp.inspect_topo()

    expected_conv = 360 if coord_kwargs["lon_max"] > 180 else 180
    assert meta.lon_wrap == expected_conv
    assert meta.y_increasing == (coord_kwargs["lat_direction"] == "S_to_N")


@pytest.mark.parametrize("dim_kwargs", TOPO_DIM_ORDER_VARIANTS)
def test_dim_order_variants(topo_file_factory, dim_kwargs):
    """
    inspect_topo() reports the correct dim_order regardless of how
    the variable axes are arranged in the file.
    """
    path = topo_file_factory(**dim_kwargs)
    with TopoInspector(path, var_name="z") as insp:
        meta = insp.inspect_topo()
    _role = {"lat": "y", "lon": "x", "time": "time"}
    assert meta.dim_order == [_role[d] for d in dim_kwargs["dim_order"]]


# ============================================================
# Variable auto-detection (no var_name supplied)
# ============================================================

@pytest.mark.parametrize("std_name", TOPO_CF_STANDARD_NAME_VARIANTS)
def test_autodetect_by_cf_standard_name(topo_file_factory, std_name):
    """
    When var_name is omitted, TopoInspector finds the elevation variable
    by matching the CF standard_name attribute.  The matched name is written
    back into the returned metadata.
    """
    ds = make_topo_dataset(var_name="topo_data", var_standard_name=std_name)
    path = topo_file_factory(ds=ds)
    with TopoInspector(path) as insp:
        meta = insp.inspect_topo()
    assert isinstance(meta, TopoMetadata)
    assert meta.var_name == "topo_data"


@pytest.mark.parametrize("var_name", TOPO_FALLBACK_NAME_VARIANTS)
def test_autodetect_by_fallback_name(topo_file_factory, var_name):
    """
    When var_name is omitted, TopoInspector falls back to matching the
    variable name against the built-in list of common elevation names.
    """
    path = topo_file_factory(var_name=var_name)
    with TopoInspector(path) as insp:
        meta = insp.inspect_topo()
    assert isinstance(meta, TopoMetadata)
    assert meta.var_name == var_name


def test_autodetect_no_match_raises(topo_file_factory):
    """
    When var_name is omitted and no variable name or CF standard_name matches
    any known elevation identifier, inspect_topo() raises ValueError with
    a message that mentions auto-detect and lists available variables.
    """
    path = topo_file_factory(var_name="completely_unknown_field")
    with TopoInspector(path) as insp:
        with pytest.raises(ValueError, match="auto-detect"):
            insp.inspect_topo()


# ============================================================
# topo_entries() — longitude offset and wrap logic
# ============================================================

def test_topo_entries_no_wrap(topo_file_factory):
    """
    File [-180, 180], crop (-90, 0, -45, 45): single entry, offset 0,
    crop_bounds in metadata equal the requested file-coord crop.
    """
    path = topo_file_factory(lon_min=-180.0, lon_max=180.0,
                              lat_min=-45.0, lat_max=45.0)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-90.0, 0.0, -45.0, 45.0))
    entries = insp.topo_entries()
    insp.close()

    assert len(entries) == 1
    ttype, fpath, meta = entries[0]
    assert ttype == 4
    assert meta.lon_wrap_offset == pytest.approx(0.0)
    assert meta.crop_bounds is not None
    lon0, lon1, lat0, lat1 = meta.crop_bounds
    assert lon0 == pytest.approx(-90.0)
    assert lon1 == pytest.approx(0.0)
    assert lat0 == pytest.approx(-45.0)
    assert lat1 == pytest.approx(45.0)


def test_topo_entries_simple_shift(topo_file_factory):
    """
    File [0, 360], crop (0, 180, -45, 45): single entry with lon_offset 0.0.
    """
    path = topo_file_factory(lon_min=0.0, lon_max=360.0,
                              lat_min=-45.0, lat_max=45.0)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(0.0, 180.0, -45.0, 45.0))
    entries = insp.topo_entries()
    insp.close()

    assert len(entries) == 1
    assert entries[0][2].lon_wrap_offset == pytest.approx(0.0)


def test_topo_entries_simple_shift_negative(topo_file_factory):
    """
    File [0, 360], crop (-180, 0, -45, 45): 1 or 2 entries are both valid;
    combined domain coverage must equal 180 degrees.
    """
    path = topo_file_factory(lon_min=0.0, lon_max=360.0,
                              lat_min=-45.0, lat_max=45.0)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-180.0, 0.0, -45.0, 45.0))
    entries = insp.topo_entries()
    insp.close()

    assert 1 <= len(entries) <= 2
    # Total lon coverage across all entries must equal the requested domain width.
    total = sum(
        meta.crop_bounds[1] - meta.crop_bounds[0]
        for _, _, meta in entries
        if meta.crop_bounds is not None
    )
    assert total == pytest.approx(180.0)


def test_topo_entries_wrap_required(topo_file_factory):
    """
    File [-180, 180], domain (-360, 0): requires two entries.
    Entry 1: file_crop [-180, 0], lon_offset 0.0.
    Entry 2: file_crop [0, 180], lon_offset -360.0.
    Combined domain coverage == 360 degrees.
    """
    path = topo_file_factory(lon_min=-180.0, lon_max=180.0,
                              lat_min=-90.0, lat_max=90.0)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-360.0, 0.0, -90.0, 90.0))
    entries = insp.topo_entries()
    insp.close()

    assert len(entries) == 2

    # Collect (file_crop_min, file_crop_max, lon_offset) pairs (order not mandated).
    pairs = {
        (meta.crop_bounds[0], meta.crop_bounds[1], meta.lon_wrap_offset)
        for _, _, meta in entries
        if meta.crop_bounds is not None
    }
    assert (-180.0, 0.0, 0.0) in pairs
    assert (0.0, 180.0, -360.0) in pairs


def test_topo_entries_no_coverage(topo_file_factory):
    """
    File [0, 90], domain [-180, -90]: no candidate offset can cover the
    domain; _compute_lon_entries raises ValueError.
    """
    path = topo_file_factory(lon_min=0.0, lon_max=90.0,
                              lat_min=-45.0, lat_max=45.0)
    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-180.0, -90.0, -45.0, 45.0))
    with pytest.raises(ValueError):
        insp.topo_entries()
    insp.close()


def test_topo_entries_near_global_gap(topo_file_factory):
    """
    Near-global file (e.g. GEBCO style) with lons [-179.99, 179.99] has a
    gap of ~0.0042 degrees at the dateline — one grid cell.  A domain that
    straddles the dateline (crop [-211, -99, 40, 75]) must still produce two
    entries without raising ValueError, since the gap is within one grid
    spacing.
    """
    import numpy as np
    import xarray as xr

    # 121 points from -179.99 to 179.99 gives spacing ~3.0 degrees; the
    # "missing" wrap column sits between 179.99 and -179.99 + 360 = 180.01.
    lons = np.linspace(-179.99, 179.99, 121)
    lats = np.linspace(40.0, 75.0, 10)
    data = np.full((len(lats), len(lons)), -500.0, dtype=np.float32)
    coords = {
        "lon": xr.DataArray(lons, dims=["lon"]),
        "lat": xr.DataArray(lats, dims=["lat"]),
    }
    ds = xr.Dataset(
        {"z": xr.DataArray(data, dims=["lat", "lon"], coords=coords,
                            attrs={"units": "m"})}
    )
    path = topo_file_factory(ds=ds)

    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-211.0, -99.0, 40.0, 75.0))
    entries = insp.topo_entries()
    insp.close()

    assert len(entries) == 2, (
        f"Expected 2 entries for dateline-straddling crop, got {len(entries)}"
    )
    # The seam gap (domain width minus combined file-coord coverage) must be
    # within one grid spacing.  crop_bounds widths equal domain-coord widths
    # because lon_offset is a constant shift.
    lon_spacing = float(lons[1] - lons[0])
    total_covered = sum(
        meta.crop_bounds[1] - meta.crop_bounds[0]
        for _, _, meta in entries
        if meta.crop_bounds is not None
    )
    domain_width = -99.0 - (-211.0)
    gap = domain_width - total_covered
    assert gap <= lon_spacing, (
        f"Seam gap {gap:.6f} exceeds grid spacing {lon_spacing:.6f}"
    )


def test_topo_entries_genuine_gap_still_errors(topo_file_factory):
    """
    File [-90, 90] genuinely cannot cover domain [-200, -50]: the uncovered
    gap (~60 degrees) is far larger than the grid spacing, so ValueError
    must still be raised even with the near-global-gap tolerance.
    """
    import numpy as np
    import xarray as xr

    lons = np.linspace(-90.0, 90.0, 91)
    lats = np.linspace(40.0, 75.0, 10)
    data = np.full((len(lats), len(lons)), -500.0, dtype=np.float32)
    coords = {
        "lon": xr.DataArray(lons, dims=["lon"]),
        "lat": xr.DataArray(lats, dims=["lat"]),
    }
    ds = xr.Dataset(
        {"z": xr.DataArray(data, dims=["lat", "lon"], coords=coords,
                            attrs={"units": "m"})}
    )
    path = topo_file_factory(ds=ds)

    insp = TopoInspector(path, var_name="z",
                             crop_bounds=(-200.0, -50.0, 40.0, 75.0))
    with pytest.raises(ValueError, match="[Gg]ap"):
        insp.topo_entries()
    insp.close()


def test_autodetect_cf_standard_name_takes_priority(topo_file_factory):
    """
    CF standard_name detection takes priority over the fallback name list.
    A dataset with both a CF-tagged variable and an 'elevation' variable must
    return the CF-tagged one.
    """
    import numpy as np
    import xarray as xr

    lons = np.linspace(-10.0, 10.0, 5)
    lats = np.linspace(-10.0, 10.0, 5)
    coords = {
        "lon": xr.DataArray(lons, dims=["lon"]),
        "lat": xr.DataArray(lats, dims=["lat"]),
    }
    data = np.full((5, 5), -100.0, dtype=np.float32)
    ds = xr.Dataset(
        {
            # CF-tagged variable: should be selected first
            "cf_topo": xr.DataArray(
                data, dims=["lat", "lon"], coords=coords,
                attrs={"units": "m", "standard_name": "surface_altitude"},
            ),
            # Fallback-name variable: present but must not be selected
            "elevation": xr.DataArray(
                data, dims=["lat", "lon"], coords=coords,
                attrs={"units": "m"},
            ),
        }
    )
    path = topo_file_factory(ds=ds)
    with TopoInspector(path) as insp:
        meta = insp.inspect_topo()
    assert meta.var_name == "cf_topo"


# ============================================================
# Geographic-axis gating of the 0-360 longitude wrap (Phase 3)
# ============================================================

def test_projected_axis_disables_lon_wrap(topo_file_factory):
    """A non-geographic x axis (projected meters) reports lon_wrap=None and
    never triggers the +/-360 wrap, even when cropped."""
    ds = make_topo_dataset(
        lon_name="x", lat_name="y",
        lon_min=0.0, lon_max=500000.0,
        lat_min=0.0, lat_max=400000.0,
        nlon=6, nlat=5,
        lon_axis_attr=True, lat_axis_attr=True,
    )
    ds["x"].attrs["units"] = "m"
    ds["y"].attrs["units"] = "m"
    path = topo_file_factory(ds=ds)

    with TopoInspector(path, var_name="z") as insp:
        assert insp.inspect_topo().lon_wrap is None

    # A crop must yield a single identity-offset entry (no 360 wrap attempt).
    with TopoInspector(path, var_name="z",
                       crop_bounds=(100000.0, 300000.0, 50000.0, 350000.0)) as insp:
        entries = insp.topo_entries()
    assert len(entries) == 1
    _, _, cropped_meta = entries[0]
    assert cropped_meta.lon_wrap_offset == 0.0


def test_length_units_mark_projected_even_with_small_values(topo_file_factory):
    """Length units flag a projected axis even when the coordinate values
    happen to fall within the degree range."""
    ds = make_topo_dataset(
        lon_name="x", lat_name="y",
        lon_min=0.0, lon_max=300.0,
        lat_min=0.0, lat_max=200.0,
        lon_axis_attr=True, lat_axis_attr=True,
    )
    ds["x"].attrs["units"] = "km"
    path = topo_file_factory(ds=ds)
    with TopoInspector(path, var_name="z") as insp:
        assert insp.inspect_topo().lon_wrap is None


def test_geographic_axis_without_units_keeps_wrap(topo_file_factory):
    """An unit-less lon/lat axis in [0, 360] is still treated as geographic,
    so the wrap convention is detected as before."""
    path = topo_file_factory(lon_min=170.0, lon_max=190.0,
                             lat_min=-10.0, lat_max=10.0)
    with TopoInspector(path, var_name="z") as insp:
        assert insp.inspect_topo().lon_wrap == 360
