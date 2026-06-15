"""
Tests for CFNormalizer.

Covers: coordinate renaming, axis/standard_name/units attribute injection,
_FillValue / missing_value conflict resolution, unknown variable passthrough,
and idempotency.
"""
import warnings

import numpy as np
import pytest
import xarray as xr

from clawpack.geoclaw.netcdf_utils import CFNormalizer

from ._helpers import make_topo_dataset

pytestmark = [pytest.mark.python, pytest.mark.netcdf]


# ============================================================
# Helpers
# ============================================================

def _simple_ds(lon_name="lon", lat_name="lat", time_name=None) -> xr.Dataset:
    """
    Build a minimal Dataset whose coord names may be non-standard.
    No axis / standard_name attributes.

    When time_name is given, the data variable includes time as a dimension so
    that CFNormalizer.rename() (which renames dimensions too) does not raise a
    CoordinateValidationError about dimensions not belonging to any variable.
    """
    lons = np.linspace(-10.0, 10.0, 4)
    lats = np.linspace(-10.0, 10.0, 4)
    coords: dict = {
        lon_name: xr.DataArray(lons, dims=[lon_name]),
        lat_name: xr.DataArray(lats, dims=[lat_name]),
    }

    if time_name is not None:
        import pandas as pd
        times = pd.date_range("2020-01-01", periods=3, freq="6h")
        coords[time_name] = xr.DataArray(times, dims=[time_name])
        data = np.ones((3, 4, 4), dtype=np.float32)
        da = xr.DataArray(data, dims=[time_name, lat_name, lon_name],
                          coords=coords, attrs={"units": "m"})
    else:
        data = np.ones((4, 4), dtype=np.float32)
        da = xr.DataArray(data, dims=[lat_name, lon_name], coords=coords,
                          attrs={"units": "m"})

    return xr.Dataset({"z": da})


def _normalize(ds: xr.Dataset) -> xr.Dataset:
    return CFNormalizer(ds).normalize()


# ============================================================
# Coordinate renaming
# ============================================================

@pytest.mark.parametrize("alias", ["lon", "x", "nav_lon", "LON", "Longitude"])
def test_lon_alias_renamed_to_longitude(alias):
    """
    Common longitude aliases are renamed to 'longitude'.
    """
    ds = _simple_ds(lon_name=alias)
    ds_norm = _normalize(ds)
    assert "longitude" in ds_norm.coords, (
        f"Expected 'longitude' in coords after renaming '{alias}'"
    )


@pytest.mark.parametrize("alias", ["lat", "y", "nav_lat", "LAT", "Latitude"])
def test_lat_alias_renamed_to_latitude(alias):
    """
    Common latitude aliases are renamed to 'latitude'.
    """
    ds = _simple_ds(lat_name=alias)
    ds_norm = _normalize(ds)
    assert "latitude" in ds_norm.coords, (
        f"Expected 'latitude' in coords after renaming '{alias}'"
    )


@pytest.mark.parametrize("alias", ["t", "TIME", "valid_time", "Time"])
def test_time_alias_renamed_to_time(alias):
    """
    Common time aliases are renamed to 'time'.
    """
    ds = _simple_ds(time_name=alias)
    ds_norm = _normalize(ds)
    assert "time" in ds_norm.coords, (
        f"Expected 'time' in coords after renaming '{alias}'"
    )


def test_already_canonical_names_not_duplicated():
    """
    A dataset already using 'longitude', 'latitude', 'time' is left unchanged
    and no duplicate coordinates are introduced.
    """
    ds = _simple_ds(lon_name="longitude", lat_name="latitude",
                    time_name="time")
    ds_norm = _normalize(ds)
    assert "longitude" in ds_norm.coords
    assert "latitude" in ds_norm.coords
    assert len(ds_norm.coords) == len(ds.coords)


def test_unrecognised_coord_name_not_renamed():
    """
    A coordinate with a name that does not match any alias (e.g. 'x_pos')
    is not renamed.
    """
    ds = _simple_ds(lon_name="x_pos", lat_name="y_pos")
    ds_norm = _normalize(ds)
    assert "x_pos" in ds_norm.coords
    assert "longitude" not in ds_norm.coords


# ============================================================
# Attribute injection
# ============================================================

def test_standard_name_added_to_longitude():
    """standard_name='longitude' is added to a renamed longitude coord."""
    ds = _simple_ds(lon_name="lon")
    ds_norm = _normalize(ds)
    assert ds_norm["longitude"].attrs.get("standard_name") == "longitude"


def test_standard_name_added_to_latitude():
    """standard_name='latitude' is added to a renamed latitude coord."""
    ds = _simple_ds(lat_name="lat")
    ds_norm = _normalize(ds)
    assert ds_norm["latitude"].attrs.get("standard_name") == "latitude"


def test_axis_x_added_to_longitude():
    """axis='X' is added to the longitude coordinate."""
    ds = _simple_ds()
    ds_norm = _normalize(ds)
    assert ds_norm["longitude"].attrs.get("axis") == "X"


def test_axis_y_added_to_latitude():
    """axis='Y' is added to the latitude coordinate."""
    ds = _simple_ds()
    ds_norm = _normalize(ds)
    assert ds_norm["latitude"].attrs.get("axis") == "Y"


def test_axis_t_added_to_time():
    """axis='T' is added to the time coordinate."""
    ds = _simple_ds(time_name="time")
    ds_norm = _normalize(ds)
    assert ds_norm["time"].attrs.get("axis") == "T"


def test_units_degrees_east_added_to_longitude():
    """units='degrees_east' is added to the longitude coordinate."""
    ds = _simple_ds()
    ds_norm = _normalize(ds)
    assert ds_norm["longitude"].attrs.get("units") == "degrees_east"


def test_units_degrees_north_added_to_latitude():
    """units='degrees_north' is added to the latitude coordinate."""
    ds = _simple_ds()
    ds_norm = _normalize(ds)
    assert ds_norm["latitude"].attrs.get("units") == "degrees_north"


def test_existing_attrs_not_overwritten():
    """
    Attributes already present on a coordinate are not overwritten by the
    normalizer (setdefault semantics).
    """
    ds = _simple_ds(lon_name="lon")
    ds["lon"].attrs["standard_name"] = "custom_longitude"
    ds["lon"].attrs["axis"] = "custom_axis"
    ds_norm = _normalize(ds)
    assert ds_norm["longitude"].attrs["standard_name"] == "custom_longitude"
    assert ds_norm["longitude"].attrs["axis"] == "custom_axis"


# ============================================================
# Fill value conflict resolution
# ============================================================

def test_only_missing_value_promoted_to_fill_value():
    """
    When only missing_value is present it is promoted to _FillValue and
    missing_value is removed.
    """
    ds = _simple_ds()
    ds["z"].attrs["missing_value"] = -8888.0

    ds_norm = _normalize(ds)

    assert ds_norm["z"].attrs.get("_FillValue") == -8888.0, (
        "missing_value should have been promoted to _FillValue"
    )
    assert "missing_value" not in ds_norm["z"].attrs, (
        "missing_value should have been removed after promotion"
    )


def test_conflicting_fill_values_warns_and_uses_fill_value():
    """
    Conflicting _FillValue and missing_value emit a UserWarning and resolve
    to _FillValue (CF precedence); missing_value is removed.
    """
    ds = _simple_ds()
    ds["z"].attrs["_FillValue"] = -9999.0
    ds["z"].attrs["missing_value"] = -8888.0

    with pytest.warns(UserWarning, match="_FillValue|missing_value|conflict"):
        ds_norm = _normalize(ds)

    assert ds_norm["z"].attrs.get("_FillValue") == -9999.0
    assert "missing_value" not in ds_norm["z"].attrs


def test_agreeing_fill_values_no_warning():
    """
    When _FillValue and missing_value agree, no warning is emitted and the
    value is kept.
    """
    ds = _simple_ds()
    ds["z"].attrs["_FillValue"] = -9999.0
    ds["z"].attrs["missing_value"] = -9999.0

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ds_norm = _normalize(ds)

    assert not caught, (
        f"Unexpected warnings when fill values agree: "
        f"{[str(w.message) for w in caught]}"
    )
    assert "missing_value" not in ds_norm["z"].attrs


def test_neither_fill_value_no_change():
    """When no fill-value attributes are present the variable is unchanged."""
    ds = _simple_ds()
    assert "_FillValue" not in ds["z"].attrs
    assert "missing_value" not in ds["z"].attrs

    ds_norm = _normalize(ds)

    assert "_FillValue" not in ds_norm["z"].attrs
    assert "missing_value" not in ds_norm["z"].attrs


# ============================================================
# Unknown variables left untouched
# ============================================================

def test_unknown_data_variable_attrs_unchanged():
    """
    Data variables with names that are not part of CF coordinate standards
    are not modified (no silent attribute injection on arbitrary variables).
    """
    ds = _simple_ds()
    ds["z"].attrs["custom_attr"] = "do_not_touch"
    original_attrs = dict(ds["z"].attrs)

    ds_norm = _normalize(ds)

    # custom_attr must still be present and unchanged
    assert ds_norm["z"].attrs.get("custom_attr") == "do_not_touch"


# ============================================================
# Original dataset not modified (copy semantics)
# ============================================================

def test_original_dataset_not_modified():
    """CFNormalizer operates on a copy; the input dataset is not modified."""
    ds = _simple_ds(lon_name="lon")
    assert "longitude" not in ds.coords   # confirm pre-condition

    _ = _normalize(ds)

    assert "longitude" not in ds.coords, (
        "CFNormalizer should not modify the original dataset"
    )
    assert "lon" in ds.coords


# ============================================================
# Idempotency
# ============================================================

def test_normalize_is_idempotent():
    """
    Calling CFNormalizer.normalize() twice produces the same result as once.
    No attributes are doubled or changed on a second pass.
    """
    ds = _simple_ds(lon_name="lon", lat_name="lat", time_name="time")
    ds["z"].attrs["missing_value"] = -8888.0   # will be promoted on first pass

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds_once = _normalize(ds)
        ds_twice = _normalize(ds_once)

    # Coordinate attributes should be the same
    for coord in ("longitude", "latitude", "time"):
        if coord in ds_once.coords:
            assert ds_once[coord].attrs == ds_twice[coord].attrs, (
                f"Coord '{coord}' attrs differ between first and second normalize"
            )

    # Data variable attributes should be stable
    for var in ds_once.data_vars:
        assert ds_once[var].attrs == ds_twice[var].attrs, (
            f"Variable '{var}' attrs differ between first and second normalize"
        )
