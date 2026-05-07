"""
Tests for TopoInterrogator.

Covers: fill-in-crop detection, fill-outside-crop pass, unit verification,
unit-conversion warn path, missing units warn path, descriptor key presence,
and all coordinate variant combinations.
"""
import warnings

import pytest

from clawpack.geoclaw.netcdf_utils import TopoInterrogator, TopoMetadata

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
    TopoInterrogator raises ValueError when NaN is present within the
    crop region.  Silent NaN in bathymetry would corrupt simulations.
    """
    path = topo_file_factory(
        lon_min=-10.0, lon_max=10.0,
        lat_min=-10.0, lat_max=10.0,
        has_fill_in_data=True,  # NaN placed at data.flat[0]
    )
    intr = TopoInterrogator(path, var_name="z",
                             crop_bounds=(-10.0, 10.0, -10.0, 10.0))
    with pytest.raises(ValueError, match="[Ff]ill"):
        intr.interrogate_topo()
    intr.close()


def test_fill_outside_crop_region_passes(topo_file_factory):
    """
    TopoInterrogator passes when NaN exists only outside the crop region.
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
    intr = TopoInterrogator(path, var_name="z",
                             crop_bounds=(5.0, 10.0, 5.0, 10.0))
    meta = intr.interrogate_topo()   # must not raise
    intr.close()
    assert isinstance(meta, TopoMetadata)


def test_no_fill_no_crop_passes(topo_file_factory):
    """TopoInterrogator passes when there are no NaN values at all."""
    path = topo_file_factory(has_fill_in_data=False)
    with TopoInterrogator(path, var_name="z") as intr:
        meta = intr.interrogate_topo()
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
        with TopoInterrogator(path, var_name="z") as intr:
            meta = intr.interrogate_topo()
    assert meta.source_units in {"m", "meter", "meters", "metre", "metres"}


def test_units_convertible_warns_and_records_source(topo_file_factory):
    """
    A variable with units convertible to meters (e.g. 'cm') triggers a
    UserWarning and records the original unit string in source_units so
    that a downstream conversion factor can be applied.
    """
    path = topo_file_factory(units="cm")
    with pytest.warns(UserWarning, match="cm"):
        with TopoInterrogator(path, var_name="z") as intr:
            meta = intr.interrogate_topo()
    assert meta.source_units == "cm", (
        f"Expected source_units='cm', got {meta.source_units!r}"
    )


def test_units_unrecognised_raises(topo_file_factory):
    """
    A variable with units not in the known conversion table raises ValueError.
    """
    path = topo_file_factory(units="furlongs")
    with TopoInterrogator(path, var_name="z") as intr:
        with pytest.raises(ValueError, match="[Uu]nrecogni"):
            intr.interrogate_topo()


def test_units_missing_warns_and_assumes_meters(topo_file_factory):
    """
    A variable with no units attribute emits a UserWarning and assumes 'm'.
    The interrogation itself must not raise.
    """
    path = topo_file_factory(units="")   # empty string → no 'units' attr written
    with pytest.warns(UserWarning):
        with TopoInterrogator(path, var_name="z") as intr:
            meta = intr.interrogate_topo()
    assert meta.source_units == "m", (
        f"Expected assumed source_units='m', got {meta.source_units!r}"
    )


# ============================================================
# Descriptor output keys
# ============================================================

def test_interrogate_topo_returns_topo_metadata(topo_file_factory):
    """interrogate_topo() returns a TopoMetadata instance."""
    path = topo_file_factory()
    with TopoInterrogator(path, var_name="z") as intr:
        meta = intr.interrogate_topo()
    assert isinstance(meta, TopoMetadata)


def test_topo_metadata_required_fields_present(topo_file_factory):
    """
    TopoMetadata contains all required fields with valid values after
    interrogation of a well-formed file.
    """
    path = topo_file_factory(
        lon_min=-10.0, lon_max=10.0,
        lat_min=-10.0, lat_max=10.0,
        units="m",
    )
    # crop_bounds is passed to TopoInterrogator, not to the dataset factory.
    intr = TopoInterrogator(path, var_name="z",
                             crop_bounds=(-5.0, 5.0, -5.0, 5.0))
    meta = intr.interrogate_topo()
    intr.close()

    assert meta.var_name == "z"
    assert meta.lon_name is not None and len(meta.lon_name) > 0
    assert meta.lat_name is not None and len(meta.lat_name) > 0
    assert meta.lon_convention in (180, 360)
    assert meta.lat_order in ("N_to_S", "S_to_N")
    assert len(meta.dim_order) == 2
    assert meta.fill_action == "abort"
    assert meta.source_units == "m"
    assert meta.crop_bounds == (-5.0, 5.0, -5.0, 5.0)


def test_topo_metadata_crop_bounds_none_when_not_specified(topo_file_factory):
    """When no crop_bounds are given, TopoMetadata.crop_bounds is None."""
    path = topo_file_factory()
    with TopoInterrogator(path, var_name="z") as intr:
        meta = intr.interrogate_topo()
    assert meta.crop_bounds is None


# ============================================================
# All coordinate variant combinations
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants(topo_file_factory, coord_kwargs):
    """
    interrogate_topo() completes without error for every coordinate variant
    and reports consistent lon_convention / lat_order values.
    """
    path = topo_file_factory(**coord_kwargs)
    with TopoInterrogator(path, var_name="z") as intr:
        meta = intr.interrogate_topo()

    expected_conv = 360 if coord_kwargs["lon_max"] > 180 else 180
    assert meta.lon_convention == expected_conv
    assert meta.lat_order == coord_kwargs["lat_direction"]


@pytest.mark.parametrize("dim_kwargs", TOPO_DIM_ORDER_VARIANTS)
def test_dim_order_variants(topo_file_factory, dim_kwargs):
    """
    interrogate_topo() reports the correct dim_order regardless of how
    the variable axes are arranged in the file.
    """
    path = topo_file_factory(**dim_kwargs)
    with TopoInterrogator(path, var_name="z") as intr:
        meta = intr.interrogate_topo()
    assert meta.dim_order == list(dim_kwargs["dim_order"])


# ============================================================
# Variable auto-detection (no var_name supplied)
# ============================================================

@pytest.mark.parametrize("std_name", TOPO_CF_STANDARD_NAME_VARIANTS)
def test_autodetect_by_cf_standard_name(topo_file_factory, std_name):
    """
    When var_name is omitted, TopoInterrogator finds the elevation variable
    by matching the CF standard_name attribute.  The matched name is written
    back into the returned metadata.
    """
    ds = make_topo_dataset(var_name="topo_data", var_standard_name=std_name)
    path = topo_file_factory(ds=ds)
    with TopoInterrogator(path) as intr:
        meta = intr.interrogate_topo()
    assert isinstance(meta, TopoMetadata)
    assert meta.var_name == "topo_data"


@pytest.mark.parametrize("var_name", TOPO_FALLBACK_NAME_VARIANTS)
def test_autodetect_by_fallback_name(topo_file_factory, var_name):
    """
    When var_name is omitted, TopoInterrogator falls back to matching the
    variable name against the built-in list of common elevation names.
    """
    path = topo_file_factory(var_name=var_name)
    with TopoInterrogator(path) as intr:
        meta = intr.interrogate_topo()
    assert isinstance(meta, TopoMetadata)
    assert meta.var_name == var_name


def test_autodetect_no_match_raises(topo_file_factory):
    """
    When var_name is omitted and no variable name or CF standard_name matches
    any known elevation identifier, interrogate_topo() raises ValueError with
    a message that mentions auto-detect and lists available variables.
    """
    path = topo_file_factory(var_name="completely_unknown_field")
    with TopoInterrogator(path) as intr:
        with pytest.raises(ValueError, match="auto-detect"):
            intr.interrogate_topo()


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
    with TopoInterrogator(path) as intr:
        meta = intr.interrogate_topo()
    assert meta.var_name == "cf_topo"
