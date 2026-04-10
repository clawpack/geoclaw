"""
Tests for NetCDFInterrogator base class.

Covers: coordinate discovery, lon convention, lat order, dim order,
fill value resolution, crop bound validation, and laziness guarantee.
"""
import warnings

import numpy as np
import pytest

from clawpack.geoclaw.netcdf_utils import NetCDFInterrogator, FileMetadata

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
    """Return an interrogator; caller is responsible for .close()."""
    return NetCDFInterrogator(path)


# ============================================================
# Lon convention detection
# ============================================================

@pytest.mark.parametrize("lon_min,lon_max,expected_convention", [
    (-10.0,  10.0, 180),   # [-180, 180] convention
    (170.0, 190.0, 360),   # [0, 360] convention
])
def test_lon_convention_detected(topo_file_factory, lon_min, lon_max,
                                 expected_convention):
    """Interrogator reports the correct longitude convention."""
    path = topo_file_factory(lon_min=lon_min, lon_max=lon_max)
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    assert meta.lon_convention == expected_convention, (
        f"Expected lon_convention={expected_convention} for lon range "
        f"[{lon_min}, {lon_max}], got {meta.lon_convention}"
    )


# ============================================================
# Lat order detection
# ============================================================

@pytest.mark.parametrize("lat_direction,expected_order", [
    ("S_to_N", "S_to_N"),
    ("N_to_S", "N_to_S"),
])
def test_lat_order_detected(topo_file_factory, lat_direction, expected_order):
    """Interrogator reports the correct latitude order."""
    path = topo_file_factory(lat_direction=lat_direction)
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    assert meta.lat_order == expected_order, (
        f"Expected lat_order={expected_order!r} for lat_direction="
        f"{lat_direction!r}, got {meta.lat_order!r}"
    )


# ============================================================
# Dim order detection
# ============================================================

@pytest.mark.parametrize("dim_kwargs", TOPO_DIM_ORDER_VARIANTS)
def test_dim_order_detected_topo(topo_file_factory, dim_kwargs):
    """
    Interrogator reports canonical dim order using role names
    ('lat', 'lon') regardless of how the variable is stored.
    """
    path = topo_file_factory(**dim_kwargs)
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    expected = list(dim_kwargs["dim_order"])
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
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    assert meta.lon_name == "x_coord"
    assert meta.lat_name == "y_coord"


def test_coord_discovery_via_fallback_name(topo_file_factory):
    """
    Coordinate discovery falls back to common names when no axis/standard_name
    is present ('lon', 'lat' are in the fallback list).
    """
    path = topo_file_factory(lon_name="lon", lat_name="lat")
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    assert meta.lon_name == "lon"
    assert meta.lat_name == "lat"


def test_missing_lon_coord_raises(topo_file_factory):
    """
    interrogate() raises ValueError when no longitude coordinate can be found.
    """
    import xarray as xr
    # Build a dataset whose coordinate names cannot be recognised
    ds = make_topo_dataset(lon_name="xpos", lat_name="lat")
    # Drop the axis/standard_name hints so it can't be found
    ds["xpos"].attrs.clear()
    path = topo_file_factory(ds=ds)
    with NetCDFInterrogator(path) as intr:
        with pytest.raises(ValueError, match="longitude"):
            intr.interrogate("z")


# ============================================================
# Fill value resolution
# ============================================================

@pytest.mark.parametrize("fill_cfg", FILL_VALUE_VARIANTS)
def test_fill_value_resolution(nc4_topo_file_factory, fill_cfg):
    """
    _resolve_fill_value returns the correct value across all attribute
    placement combinations.  Conflicting _FillValue / missing_value must not
    raise or warn from the interrogator (that is the normalizer's job).

    Note: xarray itself emits SerializationWarning when it opens a file that
    has conflicting fill values; that warning is suppressed here because it
    comes from xarray internals, not from any interrogator code.
    """
    path = nc4_topo_file_factory(
        fill_value=fill_cfg["fill_value"],
        missing_value=fill_cfg["missing_value"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # any interrogator warning is a failure
        # xarray emits SerializationWarning for conflicting fill values during
        # open_dataset; suppress it so it does not mask real test failures.
        warnings.filterwarnings("ignore", message=".*multiple fill values.*")
        with NetCDFInterrogator(path) as intr:
            meta = intr.interrogate("z")

    expected = fill_cfg["expected"]
    if expected is None:
        assert meta.fill_value is None, (
            f"Expected fill_value=None, got {meta.fill_value}"
        )
    else:
        assert meta.fill_value == pytest.approx(expected), (
            f"Expected fill_value≈{expected}, got {meta.fill_value}"
        )


def test_conflicting_fill_values_no_interrogator_warning(nc4_topo_file_factory):
    """
    Conflicting _FillValue and missing_value must not produce a warning from
    NetCDFInterrogator.  (Warnings are the CFNormalizer's responsibility.)

    xarray itself emits SerializationWarning when opening such a file; that
    warning is excluded from the assertion because it originates in xarray
    internals, not in any interrogator code.
    """
    path = nc4_topo_file_factory(fill_value=-9999.0, missing_value=-8888.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with NetCDFInterrogator(path) as intr:
            intr.interrogate("z")

    # Filter out xarray's own SerializationWarning about multiple fill values.
    interrogator_warnings = [
        w for w in caught
        if "multiple fill values" not in str(w.message)
    ]
    assert not interrogator_warnings, (
        f"Interrogator emitted unexpected warning(s): "
        f"{[str(w.message) for w in interrogator_warnings]}"
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
    with NetCDFInterrogator(path, crop_bounds=(-5.0, 5.0, -5.0, 5.0)) as intr:
        meta = intr.interrogate("z")
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
    with NetCDFInterrogator(path, crop_bounds=crop_bounds) as intr:
        with pytest.raises(ValueError):
            intr.interrogate("z")


# ============================================================
# Laziness: no data arrays loaded during interrogation
# ============================================================

def test_variables_remain_lazy_after_interrogation(topo_file_factory):
    """
    Opening the interrogator with chunks={} keeps variables as Dask arrays.
    interrogate() must not trigger a full data load.
    """
    pytest.importorskip("dask")
    import dask.array as da

    path = topo_file_factory()
    with NetCDFInterrogator(path) as intr:
        intr.interrogate("z")
        # The raw data backing the variable should still be a dask array
        assert isinstance(intr.ds["z"].data, da.Array), (
            "Variable 'z' should remain a Dask array after interrogation "
            "(no full .compute() should have been triggered)"
        )


# ============================================================
# Return type
# ============================================================

def test_interrogate_returns_file_metadata(topo_file_factory):
    """interrogate() returns a FileMetadata instance."""
    path = topo_file_factory()
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")
    assert isinstance(meta, FileMetadata)


def test_interrogate_missing_variable_raises(topo_file_factory):
    """interrogate() raises KeyError when the requested variable is absent."""
    path = topo_file_factory()
    with NetCDFInterrogator(path) as intr:
        with pytest.raises(KeyError, match="not_a_var"):
            intr.interrogate("not_a_var")


# ============================================================
# All coordinate variant combinations (smoke tests)
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants_smoke(topo_file_factory, coord_kwargs):
    """
    interrogate() completes without error for every coord variant and
    reports the expected convention and order.
    """
    path = topo_file_factory(**coord_kwargs)
    with NetCDFInterrogator(path) as intr:
        meta = intr.interrogate("z")

    lon_max = coord_kwargs["lon_max"]
    lat_dir = coord_kwargs["lat_direction"]
    expected_conv = 360 if lon_max > 180 else 180
    expected_order = lat_dir

    assert meta.lon_convention == expected_conv
    assert meta.lat_order == expected_order
