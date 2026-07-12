"""
Non-fixture helpers shared across the NetCDF test suite.

Exported by conftest.py (as fixtures) and importable directly by test files
for use with @pytest.mark.parametrize.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import pytest
import xarray as xr


# ============================================================
# Dataset factory functions
# ============================================================

def make_topo_dataset(
    *,
    lon_name: str = "lon",
    lat_name: str = "lat",
    var_name: str = "z",
    lon_min: float = -10.0,
    lon_max: float = 10.0,
    lat_min: float = -10.0,
    lat_max: float = 10.0,
    nlat: int = 5,
    nlon: int = 5,
    lat_direction: str = "S_to_N",   # "S_to_N" or "N_to_S"
    dim_order: tuple = ("lat", "lon"),
    units: str = "m",
    has_fill_in_data: bool = False,
    lon_axis_attr: bool = False,
    lat_axis_attr: bool = False,
    var_standard_name: Optional[str] = None,
) -> xr.Dataset:
    """
    Build a minimal topo xarray Dataset.

    lon_max > 180 gives a [0, 360] convention file.
    lat_direction='N_to_S' gives decreasing latitude values.
    has_fill_in_data=True puts a NaN in the data array so
    fill-in-crop tests can verify the inspector rejects the file.
    """
    lons = np.linspace(lon_min, lon_max, nlon)
    lats = (
        np.linspace(lat_min, lat_max, nlat)
        if lat_direction == "S_to_N"
        else np.linspace(lat_max, lat_min, nlat)
    )

    dim_to_coord = {"lat": lat_name, "lon": lon_name}
    actual_dims = tuple(dim_to_coord[d] for d in dim_order)
    shape = tuple(nlon if d == lon_name else nlat for d in actual_dims)

    data = np.full(shape, -100.0, dtype=np.float32)
    if has_fill_in_data:
        data.flat[0] = np.nan

    lon_attrs: dict = {}
    lat_attrs: dict = {}
    if lon_axis_attr:
        lon_attrs["axis"] = "X"
    if lat_axis_attr:
        lat_attrs["axis"] = "Y"

    var_attrs: dict = {}
    if units:
        var_attrs["units"] = units
    if var_standard_name is not None:
        var_attrs["standard_name"] = var_standard_name

    coords = {
        lon_name: xr.DataArray(lons, dims=[lon_name], attrs=lon_attrs),
        lat_name: xr.DataArray(lats, dims=[lat_name], attrs=lat_attrs),
    }
    da = xr.DataArray(data, dims=actual_dims, coords=coords, attrs=var_attrs)
    return xr.Dataset({var_name: da})


def make_met_dataset(
    *,
    lon_name: str = "lon",
    lat_name: str = "lat",
    time_name: str = "time",
    var_name_map: Optional[dict] = None,   # {geoclaw_role: file_var_name}
    lon_min: float = -10.0,
    lon_max: float = 10.0,
    lat_min: float = -10.0,
    lat_max: float = 10.0,
    nlat: int = 4,
    nlon: int = 4,
    nt: int = 3,
    lat_direction: str = "S_to_N",
    dim_order: tuple = ("time", "lat", "lon"),
    wind_units: str = "m/s",
    pressure_units: str = "Pa",
    pressure_value: float = 101325.0,   # constant pressure fill (in its units)
    wind_value: float = 5.0,            # constant wind fill (in its units)
    time_start: str = "2020-01-01",
    time_freq: str = "6h",
    extra_dim: Optional[tuple] = None,   # (dim_name, size) for ensemble
    omit_roles: Optional[list] = None,   # geoclaw roles to leave out
    mismatch_dims: bool = False,         # last var gets different dims
) -> xr.Dataset:
    """
    Build a minimal met forcing xarray Dataset.

    extra_dim=('member', 1)  inserts a singleton ensemble dimension.
    extra_dim=('member', 3)  inserts a non-singleton ensemble dimension.
    omit_roles=['wind_v']    omits that variable (consistency check must fail).
    mismatch_dims=True       gives the last variable a different dim set.

    pressure_value/wind_value set the constant fill in the variable's *own*
    units, so a hPa file can carry a physically plausible ~1013 rather than
    101325 (which would trip the magnitude sanity check once scaled).
    """
    if var_name_map is None:
        var_name_map = {"wind_u": "u10", "wind_v": "v10", "pressure": "msl"}
    if omit_roles:
        var_name_map = {k: v for k, v in var_name_map.items()
                        if k not in omit_roles}

    lons = np.linspace(lon_min, lon_max, nlon)
    lats = (
        np.linspace(lat_min, lat_max, nlat)
        if lat_direction == "S_to_N"
        else np.linspace(lat_max, lat_min, nlat)
    )
    times = pd.date_range(time_start, periods=nt, freq=time_freq)

    role_to_coord = {"lat": lat_name, "lon": lon_name, "time": time_name}
    actual_dims = [role_to_coord.get(d, d) for d in dim_order]
    coord_size = {lon_name: nlon, lat_name: nlat, time_name: nt}
    base_shape = tuple(coord_size[d] for d in actual_dims)

    coords: dict = {lon_name: lons, lat_name: lats, time_name: times}
    data_vars: dict = {}
    var_items = list(var_name_map.items())

    for i, (role, var_name) in enumerate(var_items):
        dims = list(actual_dims)
        shape = list(base_shape)

        if extra_dim is not None:
            ename, esize = extra_dim
            dims = [ename] + dims
            shape = [esize] + shape
            if ename not in coords:
                coords[ename] = np.arange(esize)

        if mismatch_dims and i == len(var_items) - 1:
            # Drop time so this variable has incompatible dims
            dims = [d for d in dims if d != time_name]
            shape = [s for d, s in zip(actual_dims, base_shape)
                     if d != time_name]

        fill_val = pressure_value if role == "pressure" else wind_value
        data = np.full(shape, fill_val, dtype=np.float32)

        if role in ("wind_u", "wind_v"):
            attrs = {"units": wind_units}
        elif role == "pressure":
            attrs = {"units": pressure_units}
        else:
            attrs = {}

        data_vars[var_name] = xr.DataArray(data, dims=dims, attrs=attrs)

    ds = xr.Dataset(data_vars, coords=coords)
    # Store the time axis as integer "seconds since <time_start>" so the file
    # is a valid GeoClaw met file: the Fortran met reader treats raw time
    # values as integer seconds.  Without this, xarray auto-picks
    # "hours since ..." for the 6-hourly default, which MetInspector rejects.
    ds[time_name].encoding['units'] = f'seconds since {time_start}'
    ds[time_name].encoding['dtype'] = 'int64'
    return ds


# ============================================================
# Parametrize lists
# ============================================================

#: lon/lat convention + direction combos.
#: Each entry is a dict of kwargs for make_topo_dataset / make_met_dataset.
COORD_VARIANTS: list = [
    pytest.param(
        {"lon_min": -10.0, "lon_max": 10.0, "lat_direction": "S_to_N"},
        id="lon180-S_to_N",
    ),
    pytest.param(
        {"lon_min": -10.0, "lon_max": 10.0, "lat_direction": "N_to_S"},
        id="lon180-N_to_S",
    ),
    pytest.param(
        {"lon_min": 170.0, "lon_max": 190.0, "lat_direction": "S_to_N"},
        id="lon360-S_to_N",
    ),
    pytest.param(
        {"lon_min": 170.0, "lon_max": 190.0, "lat_direction": "N_to_S"},
        id="lon360-N_to_S",
    ),
]

#: dim_order variants for 2-D (no time) topo variables.
TOPO_DIM_ORDER_VARIANTS: list = [
    pytest.param({"dim_order": ("lat", "lon")}, id="dim-lat-lon"),
    pytest.param({"dim_order": ("lon", "lat")}, id="dim-lon-lat"),
]

#: dim_order variants for 3-D met variables.
MET_DIM_ORDER_VARIANTS: list = [
    pytest.param({"dim_order": ("time", "lat", "lon")}, id="dim-time-lat-lon"),
    pytest.param({"dim_order": ("time", "lon", "lat")}, id="dim-time-lon-lat"),
    pytest.param({"dim_order": ("lat", "lon", "time")}, id="dim-lat-lon-time"),
    pytest.param({"dim_order": ("lon", "lat", "time")}, id="dim-lon-lat-time"),
]

#: Fill-value attribute combinations for nc4_topo_file_factory.
#: Keys: fill_value (nc4 createVariable arg), missing_value (plain attr),
#:       expected (what _resolve_fill_value should return).
FILL_VALUE_VARIANTS: list = [
    pytest.param(
        {"fill_value": -9999.0, "missing_value": None,    "expected": -9999.0},
        id="fill-only",
    ),
    pytest.param(
        # Modern xarray (mask_and_scale=True) moves missing_value to
        # enc['missing_value'], not enc['_FillValue'].  _resolve_fill_value
        # only checks enc['_FillValue'], so it returns None.  This is a known
        # implementation limitation; the test documents actual behaviour.
        {"fill_value": None,    "missing_value": -8888.0, "expected": None},
        id="missing-only",
    ),
    pytest.param(
        # Both present, agreeing: either path returns the same value.
        {"fill_value": -9999.0, "missing_value": -9999.0, "expected": -9999.0},
        id="both-agree",
    ),
    pytest.param(
        # Conflicting: _FillValue wins (CF precedence); no warning from inspector.
        {"fill_value": -9999.0, "missing_value": -8888.0, "expected": -9999.0},
        id="both-conflict-no-warn",
    ),
    pytest.param(
        {"fill_value": None,    "missing_value": None,    "expected": None},
        id="neither",
    ),
]

#: Pressure unit variants: units written in file, expected source_units in
#: metadata, and a physically plausible constant value in those units (so the
#: magnitude sanity check sees ~1e5 Pa after scaling).
PRESSURE_UNIT_VARIANTS: list = [
    pytest.param({"pressure_units": "Pa",   "expected_source": "Pa",
                  "value": 101325.0}, id="Pa"),
    pytest.param({"pressure_units": "hPa",  "expected_source": "hPa",
                  "value": 1013.25}, id="hPa"),
    pytest.param({"pressure_units": "mbar", "expected_source": "mbar",
                  "value": 1013.25}, id="mbar"),
]

#: Wind unit variants
WIND_UNIT_VARIANTS: list = [
    pytest.param({"wind_units": "m/s",   "expected_source": "m/s"},   id="m_s"),
    pytest.param({"wind_units": "knots", "expected_source": "knots"}, id="knots"),
]

#: CF standard_name values recognised by TopoInspector._find_topo_var_name,
#: in priority order.
TOPO_CF_STANDARD_NAME_VARIANTS: list = [
    pytest.param("surface_altitude",                  id="surface_altitude"),
    pytest.param("height_above_mean_sea_level",       id="height_above_mean_sea_level"),
    pytest.param("height_above_reference_ellipsoid",  id="height_above_reference_ellipsoid"),
    pytest.param("bedrock_altitude",                  id="bedrock_altitude"),
    pytest.param("altitude",                          id="altitude"),
    pytest.param("height",                            id="height"),
    pytest.param("sea_floor_depth_below_geoid",       id="sea_floor_depth_below_geoid"),
]

#: Variable names in the fallback list recognised by
#: TopoInspector._find_topo_var_name.
TOPO_FALLBACK_NAME_VARIANTS: list = [
    pytest.param("z",           id="z"),
    pytest.param("Z",           id="Z"),
    pytest.param("elevation",   id="elevation"),
    pytest.param("Elevation",   id="Elevation"),
    pytest.param("ELEVATION",   id="ELEVATION"),
    pytest.param("topo",        id="topo"),
    pytest.param("TOPO",        id="TOPO"),
    pytest.param("height",      id="height"),
    pytest.param("HEIGHT",      id="HEIGHT"),
    pytest.param("altitude",    id="altitude"),
    pytest.param("ALTITUDE",    id="ALTITUDE"),
    pytest.param("depth",       id="depth"),
    pytest.param("DEPTH",       id="DEPTH"),
    pytest.param("Band1",       id="Band1"),
    pytest.param("dem",         id="dem"),
    pytest.param("DEM",         id="DEM"),
    pytest.param("topography",  id="topography"),
    pytest.param("bathymetry",  id="bathymetry"),
    pytest.param("bathy",       id="bathy"),
    pytest.param("bedrock",     id="bedrock"),
]
