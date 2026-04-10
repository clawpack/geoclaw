"""
Shared pytest fixtures for GeoClaw NetCDF interrogator tests.

Factory functions, parametrize lists, and helper constants live in
_helpers.py so they can be imported directly by test files without
going through the pytest fixture machinery.

Fixtures provided here
----------------------
topo_file_factory     write a topo Dataset to a tmp_path .nc file
met_file_factory      write a met Dataset to a tmp_path .nc file
nc4_topo_file_factory write a topo file via netCDF4 for precise
                      attribute placement (fill value tests)
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from ._helpers import make_topo_dataset, make_met_dataset

try:
    import netCDF4 as nc4  # type: ignore
    _HAS_NC4 = True
except ImportError:
    _HAS_NC4 = False


@pytest.fixture
def topo_file_factory(tmp_path):
    """
    Return a callable that writes a topo Dataset to a .nc file in tmp_path.

    Signature::

        path = topo_file_factory(**make_topo_dataset_kwargs)
        # or supply an already-built dataset:
        path = topo_file_factory(ds=my_ds)

    Each call writes a distinct file (topo_1.nc, topo_2.nc, …).
    """
    counter = [0]

    def _make(ds=None, name: Optional[str] = None, **kwargs) -> Path:
        if ds is None:
            ds = make_topo_dataset(**kwargs)
        if name is None:
            counter[0] += 1
            name = f"topo_{counter[0]}.nc"
        path = tmp_path / name
        ds.to_netcdf(path)
        return path

    return _make


@pytest.fixture
def met_file_factory(tmp_path):
    """
    Return a callable that writes a met Dataset to a .nc file in tmp_path.

    Signature::

        path = met_file_factory(**make_met_dataset_kwargs)
        # or supply an already-built dataset:
        path = met_file_factory(ds=my_ds)
    """
    counter = [0]

    def _make(ds=None, name: Optional[str] = None, **kwargs) -> Path:
        if ds is None:
            from ._helpers import make_met_dataset as _make_met
            ds = _make_met(**kwargs)
        if name is None:
            counter[0] += 1
            name = f"met_{counter[0]}.nc"
        path = tmp_path / name
        ds.to_netcdf(path)
        return path

    return _make


@pytest.fixture
def nc4_topo_file_factory(tmp_path):
    """
    Return a callable that writes a topo file via the netCDF4 library,
    giving precise control over which fill-related attributes are present.

    This is required for fill value tests because xarray (mask_and_scale=True)
    moves _FillValue to encoding on read, so the only reliable way to exercise
    specific attribute combinations is to write at the netCDF4 level.

    Signature::

        path = nc4_topo_file_factory(
            fill_value=-9999.0,   # set at createVariable time; None = omit
            missing_value=-8888.0,  # written as plain attr; None = omit
            **other_kwargs,
        )

    Skips automatically if netCDF4 is not importable.
    """
    if not _HAS_NC4:
        pytest.skip("netCDF4 not available")

    counter = [0]

    def _make(
        *,
        fill_value: Optional[float] = None,
        missing_value: Optional[float] = None,
        lon_min: float = -10.0,
        lon_max: float = 10.0,
        lat_min: float = -10.0,
        lat_max: float = 10.0,
        nlat: int = 5,
        nlon: int = 5,
        lat_direction: str = "S_to_N",
        var_name: str = "z",
        lon_name: str = "lon",
        lat_name: str = "lat",
        units: str = "m",
        has_fill_in_data: bool = False,
        name: Optional[str] = None,
    ) -> Path:
        counter[0] += 1
        if name is None:
            name = f"nc4topo_{counter[0]}.nc"
        path = tmp_path / name

        lons = np.linspace(lon_min, lon_max, nlon)
        lats = (
            np.linspace(lat_min, lat_max, nlat)
            if lat_direction == "S_to_N"
            else np.linspace(lat_max, lat_min, nlat)
        )
        data = np.full((nlat, nlon), -100.0, dtype=np.float32)
        if has_fill_in_data:
            # Write the fill sentinel value so xarray masks it as NaN.
            sentinel = fill_value if fill_value is not None else 9.969209968386869e+36
            data[0, 0] = sentinel

        with nc4.Dataset(str(path), "w") as ds:
            ds.createDimension(lon_name, nlon)
            ds.createDimension(lat_name, nlat)

            lon_var = ds.createVariable(lon_name, "f4", (lon_name,))
            lon_var[:] = lons

            lat_var = ds.createVariable(lat_name, "f4", (lat_name,))
            lat_var[:] = lats

            # _FillValue must be set at createVariable time in netCDF4
            create_kwargs: dict = {}
            if fill_value is not None:
                create_kwargs["fill_value"] = fill_value

            z_var = ds.createVariable(
                var_name, "f4", (lat_name, lon_name), **create_kwargs
            )
            z_var.units = units
            if missing_value is not None:
                z_var.missing_value = missing_value
            z_var[:] = data

        return path

    return _make
