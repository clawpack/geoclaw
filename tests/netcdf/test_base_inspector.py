"""
Tests for NetCDFInspector base class.

Covers: coordinate discovery, lon convention, lat order, dim order,
fill value resolution, crop bound validation, and laziness guarantee.
"""
import warnings

import numpy as np
import pytest

from clawpack.geoclaw.netcdf_utils import NetCDFInspector, FileMetadata

from ._helpers import (
    make_topo_dataset,
    COORD_VARIANTS,
    TOPO_DIM_ORDER_VARIANTS,
    FILL_VALUE_VARIANTS,
)

pytestmark = [pytest.mark.python, pytest.mark.netcdf]


# ============================================================
# Helpers
# ============================================================

def _open(path):
    """Return an inspector; caller is responsible for .close()."""
    return NetCDFInspector(path)


# ============================================================
# Lon convention detection
# ============================================================

@pytest.mark.parametrize("lon_min,lon_max,expected_convention", [
    (-10.0,  10.0, 180),   # [-180, 180] convention
    (170.0, 190.0, 360),   # [0, 360] convention
])
def test_lon_convention_detected(topo_file_factory, lon_min, lon_max,
                                 expected_convention):
    """Inspector reports the correct longitude convention."""
    path = topo_file_factory(lon_min=lon_min, lon_max=lon_max)
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    assert meta.lon_wrap == expected_convention, (
        f"Expected lon_wrap={expected_convention} for lon range "
        f"[{lon_min}, {lon_max}], got {meta.lon_wrap}"
    )


# ============================================================
# Lat order detection
# ============================================================

@pytest.mark.parametrize("lat_direction,expected_increasing", [
    ("S_to_N", True),
    ("N_to_S", False),
])
def test_lat_order_detected(topo_file_factory, lat_direction,
                            expected_increasing):
    """Inspector reports the correct y_increasing flag."""
    path = topo_file_factory(lat_direction=lat_direction)
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    assert meta.y_increasing == expected_increasing, (
        f"Expected y_increasing={expected_increasing!r} for lat_direction="
        f"{lat_direction!r}, got {meta.y_increasing!r}"
    )


# ============================================================
# Dim order detection
# ============================================================

@pytest.mark.parametrize("dim_kwargs", TOPO_DIM_ORDER_VARIANTS)
def test_dim_order_detected_topo(topo_file_factory, dim_kwargs):
    """
    Inspector reports canonical dim order using role names
    ('y', 'x') regardless of how the variable is stored.
    """
    path = topo_file_factory(**dim_kwargs)
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    _role = {"lat": "y", "lon": "x", "time": "time"}
    expected = [_role[d] for d in dim_kwargs["dim_order"]]
    assert meta.dim_order == expected, (
        f"Expected dim_order={expected}, got {meta.dim_order}"
    )


# ============================================================
# Coordinate discovery via axis attribute
# ============================================================

def test_coord_discovery_via_axis_attr(topo_file_factory):
    """
    Coordinate discovery succeeds when only the CF axis attribute is present
    (no standard_name, non-standard variable name).
    """
    path = topo_file_factory(
        lon_name="x_coord", lat_name="y_coord",
        lon_axis_attr=True, lat_axis_attr=True,
    )
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    assert meta.x_name == "x_coord"
    assert meta.y_name == "y_coord"


def test_coord_discovery_via_fallback_name(topo_file_factory):
    """
    Coordinate discovery falls back to common names when no axis/standard_name
    is present ('lon', 'lat' are in the fallback list).
    """
    path = topo_file_factory(lon_name="lon", lat_name="lat")
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    assert meta.x_name == "lon"
    assert meta.y_name == "lat"


def test_missing_lon_coord_raises(topo_file_factory):
    """
    inspect() raises ValueError when no longitude coordinate can be found.
    """
    import xarray as xr
    # Build a dataset whose coordinate names cannot be recognised
    ds = make_topo_dataset(lon_name="xpos", lat_name="lat")
    # Drop the axis/standard_name hints so it can't be found
    ds["xpos"].attrs.clear()
    path = topo_file_factory(ds=ds)
    with NetCDFInspector(path) as insp:
        with pytest.raises(ValueError, match="longitude"):
            insp.inspect("z")


# ============================================================
# Fill value resolution
# ============================================================

@pytest.mark.parametrize("fill_cfg", FILL_VALUE_VARIANTS)
def test_fill_value_resolution(nc4_topo_file_factory, fill_cfg):
    """
    _resolve_fill_value returns the correct value across all attribute
    placement combinations.  Conflicting _FillValue / missing_value must not
    raise or warn from the inspector (that is the normalizer's job).

    Note: xarray itself emits SerializationWarning when it opens a file that
    has conflicting fill values; that warning is suppressed here because it
    comes from xarray internals, not from any inspector code.
    """
    path = nc4_topo_file_factory(
        fill_value=fill_cfg["fill_value"],
        missing_value=fill_cfg["missing_value"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # any inspector warning is a failure
        # xarray emits SerializationWarning for conflicting fill values during
        # open_dataset; suppress it so it does not mask real test failures.
        warnings.filterwarnings("ignore", message=".*multiple fill values.*")
        with NetCDFInspector(path) as insp:
            meta = insp.inspect("z")

    expected = fill_cfg["expected"]
    if expected is None:
        assert meta.fill_value is None, (
            f"Expected fill_value=None, got {meta.fill_value}"
        )
    else:
        assert meta.fill_value == pytest.approx(expected), (
            f"Expected fill_value≈{expected}, got {meta.fill_value}"
        )


def test_conflicting_fill_values_no_inspector_warning(nc4_topo_file_factory):
    """
    Conflicting _FillValue and missing_value must not produce a warning from
    NetCDFInspector.  (Warnings are the CFNormalizer's responsibility.)

    xarray itself emits SerializationWarning when opening such a file; that
    warning is excluded from the assertion because it originates in xarray
    internals, not in any inspector code.
    """
    path = nc4_topo_file_factory(fill_value=-9999.0, missing_value=-8888.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with NetCDFInspector(path) as insp:
            insp.inspect("z")

    # Filter out xarray's own SerializationWarning about multiple fill values.
    inspector_warnings = [
        w for w in caught
        if "multiple fill values" not in str(w.message)
    ]
    assert not inspector_warnings, (
        f"Inspector emitted unexpected warning(s): "
        f"{[str(w.message) for w in inspector_warnings]}"
    )


# ============================================================
# Crop bound validation
# ============================================================

def test_crop_bounds_inside_extent_passes(topo_file_factory):
    """
    crop_bounds fully inside the file extent does not raise.
    """
    path = topo_file_factory(lon_min=-10.0, lon_max=10.0,
                             lat_min=-10.0, lat_max=10.0)
    with NetCDFInspector(path, crop_bounds=(-5.0, 5.0, -5.0, 5.0)) as insp:
        meta = insp.inspect("z")
    assert meta.crop_bounds == (-5.0, 5.0, -5.0, 5.0)


@pytest.mark.parametrize("crop_bounds,description", [
    ((-20.0, 5.0, -5.0, 5.0), "lon_min too small"),
    ((-5.0, 20.0, -5.0, 5.0), "lon_max too large"),
    ((-5.0, 5.0, -20.0, 5.0), "lat_min too small"),
    ((-5.0, 5.0, -5.0, 20.0), "lat_max too large"),
])
def test_crop_bounds_outside_extent_raises(topo_file_factory, crop_bounds,
                                           description):
    """
    crop_bounds that exceed the file spatial extent raise ValueError with a
    message that identifies which bound is out of range.
    """
    path = topo_file_factory(lon_min=-10.0, lon_max=10.0,
                             lat_min=-10.0, lat_max=10.0)
    with NetCDFInspector(path, crop_bounds=crop_bounds) as insp:
        with pytest.raises(ValueError):
            insp.inspect("z")


# ============================================================
# Laziness: no data arrays loaded during inspection
# ============================================================

def test_variables_remain_lazy_after_inspection(topo_file_factory):
    """
    Opening the inspector with chunks={} keeps variables as Dask arrays.
    inspect() must not trigger a full data load.
    """
    pytest.importorskip("dask")
    import dask.array as da

    path = topo_file_factory()
    with NetCDFInspector(path) as insp:
        insp.inspect("z")
        # The raw data backing the variable should still be a dask array
        assert isinstance(insp.ds["z"].data, da.Array), (
            "Variable 'z' should remain a Dask array after inspection "
            "(no full .compute() should have been triggered)"
        )


# ============================================================
# Return type
# ============================================================

def test_inspect_returns_file_metadata(topo_file_factory):
    """inspect() returns a FileMetadata instance."""
    path = topo_file_factory()
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")
    assert isinstance(meta, FileMetadata)


def test_inspect_missing_variable_raises(topo_file_factory):
    """inspect() raises KeyError when the requested variable is absent."""
    path = topo_file_factory()
    with NetCDFInspector(path) as insp:
        with pytest.raises(KeyError, match="not_a_var"):
            insp.inspect("not_a_var")


# ============================================================
# All coordinate variant combinations (smoke tests)
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants_smoke(topo_file_factory, coord_kwargs):
    """
    inspect() completes without error for every coord variant and
    reports the expected convention and order.
    """
    path = topo_file_factory(**coord_kwargs)
    with NetCDFInspector(path) as insp:
        meta = insp.inspect("z")

    lon_max = coord_kwargs["lon_max"]
    lat_dir = coord_kwargs["lat_direction"]
    expected_conv = 360 if lon_max > 180 else 180
    expected_increasing = (lat_dir == "S_to_N")

    assert meta.lon_wrap == expected_conv
    assert meta.y_increasing == expected_increasing
