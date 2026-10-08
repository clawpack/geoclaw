#!/usr/bin/env python
# encoding: utf-8
"""Byte-exact goldens for the ``topo.data`` writer.

Every other test of :meth:`TopographyData.write` checks substrings or line
counts, so the two write paths could diverge in field width, whitespace,
quoting or float formatting and still pass.  Fortran reads this file with
fixed-format ``read`` statements, so those details are the contract.

These goldens were generated **before** ``write()`` was restructured into
resolve-then-write (PR A2).  That restructuring must be a pure refactor for
every case that does not wrap across the antimeridian, and this file is what
pins it: if a byte moves in any case below, the refactor changed behavior it
was not supposed to touch.

Regenerate deliberately with ``GEOCLAW_REGEN=1`` and review the diff -- never
because a test failed.

The absolute path GeoClaw writes into ``topo.data`` depends on ``tmp_path``, so
paths are replaced with ``<PATH>`` before comparison.  Everything else,
including the exact spacing around ``# topo_type`` and the ``repr`` of every
float, is compared verbatim.
"""

import os
import re
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools
from clawpack.geoclaw.data import TopographyData

testdir = Path(__file__).parent
golden_dir = testdir / "data" / "topo_data_golden"

pytestmark = pytest.mark.python


def _write_tt3(path, x0=-100.0, x1=-90.0, y0=20.0, y1=30.0, delta=0.5):
    """A small ASCII topo_type=3 file on an exactly-representable lattice."""
    x = np.arange(x0, x1 + 1e-9, delta)
    y = np.arange(y0, y1 + 1e-9, delta)
    Z = -100.0 + np.outer(np.linspace(0.0, 50.0, y.size), np.ones_like(x))
    topo = topotools.Topography()
    topo.set_xyZ(x, y, Z)
    topo.write(str(path), topo_type=3)
    return path


def _write_nc(path, x0=-100.0, x1=-90.0, y0=20.0, y1=30.0, delta=0.5):
    """A CF-compliant NetCDF topo file (topo_type=4)."""
    netCDF4 = pytest.importorskip("netCDF4")
    x = np.arange(x0, x1 + 1e-9, delta)
    y = np.arange(y0, y1 + 1e-9, delta)
    Z = -100.0 + np.outer(np.linspace(0.0, 50.0, y.size), np.ones_like(x))
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
        v.units = "m"
        v.standard_name = "height_above_mean_sea_level"
        v.positive = "up"
        ds.Conventions = "CF-1.8"
    return path


def _normalize(text, tmp_path):
    """Replace machine-specific absolute paths with a stable placeholder."""
    text = text.replace(str(Path(tmp_path).resolve()), "<TMP>")
    text = text.replace(str(tmp_path), "<TMP>")
    # The data_source header line carries a timestamp/path in some setups.
    return re.sub(r"<TMP>[^\s']*/", "<TMP>/", text)


def _check(name, tmp_path, out_file):
    # ".txt", not ".data": the repo's .gitignore excludes "*.data", so a
    # golden named topo.data would be silently left out of the commit and the
    # test would fail for everyone else.  Matches the met_forcing goldens.
    golden = golden_dir / f"{name}.txt"
    actual = _normalize(Path(out_file).read_text(), tmp_path)

    if os.environ.get("GEOCLAW_REGEN"):
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(actual)
        return

    assert golden.exists(), (
        f"Missing golden {golden}. Generate it with GEOCLAW_REGEN=1 from the "
        f"code you intend to pin -- not from code you are about to change.")
    expected = golden.read_text()
    assert actual == expected, (
        f"topo.data for '{name}' changed.\n--- expected ---\n{expected}\n"
        f"--- actual ---\n{actual}")


def test_golden_ascii_single(tmp_path):
    """One ASCII file, no preprocessing: the simplest possible topo.data."""
    path = _write_tt3(tmp_path / "a.tt3")
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = 3

    td = TopographyData()
    td.topofiles = [topo]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    _check("ascii_single", tmp_path, out)


def test_golden_ascii_all_preprocessing(tmp_path):
    """Every preprocessing attribute non-default at once.

    Pins the float formatting of the crop_extent and align lines: these are
    written with repr() so coordinates reach Fortran at full precision, and a
    change to %g would silently truncate them.
    """
    path = _write_tt3(tmp_path / "b.tt3")
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = 3
    topo.crop_extent = [-98.25, -92.75, 22.5, 27.5]
    topo.coarsen = 2
    topo.buffer = 3
    topo.align = [-100.0, 20.0]
    topo.x_shift = 0.125
    topo.y_shift = -0.25
    topo.z_shift = 1.5
    topo.negate_z = True

    td = TopographyData()
    td.topofiles = [topo]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    _check("ascii_all_preprocessing", tmp_path, out)


def test_golden_two_files_priority_order(tmp_path):
    """Two files of different resolution: pins the coarsest-first ordering."""
    coarse = _write_tt3(tmp_path / "coarse.tt3", delta=1.0)
    fine = _write_tt3(tmp_path / "fine.tt3", delta=0.25)

    topos = []
    for p in (fine, coarse):          # deliberately not in priority order
        t = topotools.Topography()
        t.path = str(p)
        t.topo_type = 3
        topos.append(t)

    td = TopographyData()
    td.topofiles = topos
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    _check("two_files_priority", tmp_path, out)


@pytest.mark.netcdf
def test_golden_netcdf_no_crop(tmp_path):
    """type-4 with no crop: pins the whole CF descriptor block."""
    pytest.importorskip("xarray")
    path = _write_nc(tmp_path / "c.nc")
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = 4

    td = TopographyData()
    td.topofiles = [topo]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    _check("netcdf_no_crop", tmp_path, out)


@pytest.mark.netcdf
def test_golden_netcdf_cropped(tmp_path):
    """type-4 with an ordinary in-extent crop and a buffer.

    This is the case the resolve-then-write refactor touches most, and the one
    most likely to change by accident: it must keep emitting exactly one entry,
    with crop_bounds in file coordinates and lon_wrap_offset 0.0.
    """
    pytest.importorskip("xarray")
    path = _write_nc(tmp_path / "d.nc")
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = 4
    topo.crop_extent = [-98.0, -92.0, 22.0, 28.0]
    topo.buffer = 2
    topo.coarsen = 2

    td = TopographyData()
    td.topofiles = [topo]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    _check("netcdf_cropped", tmp_path, out)
