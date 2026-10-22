#!/usr/bin/env python
# encoding: utf-8
"""One round-off allowance where a crop edge meets a grid point (§2.8).

Python's crop window used exact comparisons, while Fortran's had a 1e-6*dx
allowance.  So a crop written on a grid line could keep that line in Fortran
and drop it in Python, without a warning.  The usual cause is a coordinate
rebuilt with round-off: ``Topography.write`` stores a 0.1 spacing as
``9.9999999999994e-02``, and the point read back for -99.5 is
``-99.50000000000003``.  Both sides now use ``CROP_TOL`` (as a fraction of the
spacing): ``coordinate_tools.CROP_TOL`` and ``CROP_TOL`` in
src/2d/shallow/topo_module.f90.
"""

import math
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from clawpack.geoclaw import coordinate_tools, topotools

pytestmark = pytest.mark.python

TOPO_MODULE = (Path(__file__).resolve().parents[1]
               / "src" / "2d" / "shallow" / "topo_module.f90")

SPACINGS = [0.1, 1 / 60, 1 / 120, 0.25, 0.05, 1 / 240]
SHIFTS = [0.0, 0.1, 1 / 3, 360.0, -0.05]
XLL = -100.0
SPAN = 20.0


def test_python_and_fortran_use_the_same_tolerance():
    text = TOPO_MODULE.read_text()
    m = re.search(r"parameter\s*::\s*CROP_TOL\s*=\s*([0-9.dDeE+-]+)", text)
    assert m is not None, "CROP_TOL parameter not found in topo_module.f90"
    fortran = float(m.group(1).lower().replace("d", "e"))
    assert fortran == coordinate_tools.CROP_TOL


def _read_back_x(tmp_path, dx):
    """x as Python rebuilds it from a .tt3 that Topography.write produced."""
    n = int(round(SPAN / dx)) + 1
    t = topotools.Topography()
    t.x = XLL + dx * np.arange(n)
    t.y = np.array([0.0, dx, 2 * dx])
    t.Z = np.zeros((3, n))
    path = tmp_path / "grid.tt3"
    t.write(str(path), topo_type=3)
    r = topotools.Topography()
    r.read(path=str(path), topo_type=3)
    return path, r.x


def _fortran_window(lo, hi, shift, xlow, dx, n):
    """Index window from the types 2/3 formula in topo_module.f90
    (read_topo_settings), with the header values Fortran reads."""
    tol = coordinate_tools.CROP_TOL
    i_start = max(0, math.ceil((lo - shift - xlow) / dx - tol))
    i_end = min(n - 1, math.floor((hi - shift - xlow) / dx + tol))
    return i_start, i_end


@pytest.mark.parametrize("dx", SPACINGS, ids=lambda d: f"dx={d:.6g}")
def test_python_window_matches_fortran_formula(tmp_path, dx):
    """For crops on grid lines, Python picks the same first and last points
    as Fortran's index formula, for every spacing and shift tried."""
    path, x = _read_back_x(tmp_path, dx)
    header = path.read_text().splitlines()
    xlow = float(header[2].split()[0])
    dx_file = float(header[4].split()[0])
    n = x.size

    disagree = []
    for shift in SHIFTS:
        for k in range(5, n - 15, max(1, n // 60)):
            lo = round(XLL + k * dx + shift, 10)
            hi = round(XLL + (k + 10) * dx + shift, 10)
            delta = abs(x[1] - x[0])
            window = coordinate_tools.crop_indices(x + shift, lo, hi, delta)
            py = (window[0], window[1] - 1)
            fo = _fortran_window(lo, hi, shift, xlow, dx_file, n)
            if py != fo:
                disagree.append((shift, lo, hi, py, fo))
    assert not disagree, (
        f"{len(disagree)} crops differ (shift, lo, hi, python, fortran): "
        f"{disagree[:5]}")


@pytest.mark.parametrize("topo_type", [3, 4])
def test_crop_on_a_grid_line_keeps_that_line(tmp_path, topo_type):
    """The reported case: crop [-99.5, -99.0] on a 0.1-degree grid keeps
    both edge columns, six points, for ASCII and NetCDF."""
    dx = 0.1
    path, _x = _read_back_x(tmp_path, dx)
    if topo_type == 4:
        xr = pytest.importorskip("xarray")
        pytest.importorskip("netCDF4")
        src = topotools.Topography()
        src.read(path=str(path), topo_type=3)
        nc = tmp_path / "grid.nc"
        xr.Dataset(
            {"z": (("lat", "lon"), src.Z, {"units": "m"})},
            coords={"lon": ("lon", src.x, {"units": "degrees_east",
                                           "standard_name": "longitude"}),
                    "lat": ("lat", src.y, {"units": "degrees_north",
                                           "standard_name": "latitude"})},
        ).to_netcdf(nc)
        path = nc

    t = topotools.Topography()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t.read(path=str(path), topo_type=topo_type,
               crop_extent=[-99.5, -99.0, 0.0, 0.2])
    assert t.x.size == 6, t.x
    assert t.x[0] == pytest.approx(-99.5, abs=1e-9)
    assert t.x[-1] == pytest.approx(-99.0, abs=1e-9)
