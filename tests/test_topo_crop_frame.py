#!/usr/bin/env python
# encoding: utf-8
"""A topo ``crop_extent`` is in domain coordinates for every topo type.

``crop_extent`` is the region of the *domain* to keep, so it is compared
against the topography after ``x_shift``/``y_shift`` are applied
(dev/design/topo_input_roadmap.md §2.8).  ASCII (types 2/3) already does this.
NetCDF (type 4) applied the crop in the file's own coordinates and shifted
afterwards, in Python and in Fortran alike, so the same ``Topography`` settings
selected a different region depending on the file format.

The Fortran half is in tests/regression/topo_crop/test_topo_descriptor_crop.py.
"""

import warnings

import numpy as np
import pytest

from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import TopographyData

pytestmark = pytest.mark.python

# File grid: node-registered, DX a power of two so every coordinate is exact.
FILE_X = (-100.0, -80.0)
FILE_Y = (20.0, 30.0)
DX = 0.125

# Shifts are whole multiples of DX, so shifted coordinates stay on the grid.
X_SHIFT = 2.0
Y_SHIFT = 1.0

# Crop in domain coordinates, inside the shifted file [-98, -78] x [21, 31].
CROP = (-95.0, -90.0, 22.0, 25.0)

TOPO_TYPES = [2, 3, 4]


def _field(x_file, y_file):
    """Linear in both directions, so any frame error changes the values."""
    return -100.0 + 2.0 * (x_file - FILE_X[0]) + 3.0 * (y_file - FILE_Y[0])


def _write_file(tmp_path, topo_type):
    x = np.arange(FILE_X[0], FILE_X[1] + DX / 2, DX)
    y = np.arange(FILE_Y[0], FILE_Y[1] + DX / 2, DX)
    X, Y = np.meshgrid(x, y)
    Z = _field(X, Y)
    if topo_type == 4:
        xr = pytest.importorskip("xarray")
        pytest.importorskip("netCDF4")
        path = tmp_path / "topo.nc"
        xr.Dataset(
            {"z": (("lat", "lon"), Z, {"units": "m"})},
            coords={
                "lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"}),
            }).to_netcdf(path)
        return path
    t = topotools.Topography()
    t.x, t.y, t.Z = x, y, Z
    path = tmp_path / f"topo.tt{topo_type}"
    t.write(str(path), topo_type=topo_type)
    return path


def _shifted_topo(path, topo_type, crop):
    t = topotools.Topography()
    t.path = str(path)
    t.topo_type = topo_type
    t.crop_extent = list(crop)
    t.x_shift = X_SHIFT
    t.y_shift = Y_SHIFT
    return t


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_read_crops_in_domain_coordinates(tmp_path, topo_type):
    """read() with a shift returns exactly CROP, holding the shifted field."""
    path = _write_file(tmp_path, topo_type)
    t = _shifted_topo(path, topo_type, CROP)
    with warnings.catch_warnings():
        warnings.simplefilter("error")      # no clipping warning either
        t.read(path=str(path), topo_type=topo_type, crop_extent=list(CROP))

    assert (t.x[0], t.x[-1], t.y[0], t.y[-1]) == pytest.approx(CROP)
    X, Y = np.meshgrid(t.x, t.y)
    np.testing.assert_allclose(t.Z, _field(X - X_SHIFT, Y - Y_SHIFT),
                               rtol=0, atol=1e-10)


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_write_accepts_crop_inside_shifted_file(tmp_path, topo_type):
    """Inside the file once shifted, outside it unshifted: valid."""
    path = _write_file(tmp_path, topo_type)
    crop = (-79.0, -78.5, 22.0, 25.0)     # file x ends at -80 before the shift
    td = TopographyData()
    td.topofiles = [_shifted_topo(path, topo_type, crop)]
    td.write(out_file=str(tmp_path / "topo.data"))


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_write_rejects_crop_outside_shifted_file(tmp_path, topo_type):
    """Inside the file unshifted, outside it once shifted: an error."""
    path = _write_file(tmp_path, topo_type)
    crop = (-99.5, -99.0, 22.0, 25.0)     # shifted file x starts at -98
    td = TopographyData()
    td.topofiles = [_shifted_topo(path, topo_type, crop)]
    with pytest.raises(ValueError, match="not covered by"):
        td.write(out_file=str(tmp_path / "topo.data"))
