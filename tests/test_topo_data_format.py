#!/usr/bin/env python
# encoding: utf-8
"""``topo.data`` and ``dtopo.data`` use one notation: ``value =: key``.

Before this, one file mixed three notations: ``value =: key`` for the header,
``value   # key`` for each file's settings, and a ``key = value`` netCDF
descriptor ended by a blank line (dev/design/topo_input_roadmap.md §2.8,
PR B).  Now every value line is ``value =: key`` with an optional
``# comment``, in a fixed order, and Fortran checks each key as it reads, so
a file from an older GeoClaw stops with "regenerate with make data" instead of
being misread.

Unset values are written as ``None``, the convention ``surge.data`` already
uses.  This matters for ``align``: the old all-zero sentinel could not tell
"no align" from ``align = [0, 0]``, so Fortran dropped an align through the
origin that Python applied.

The Fortran half is in tests/regression/topo_crop/test_topo_descriptor_crop.py.
"""

import re

import numpy as np
import pytest

from clawpack.geoclaw import dtopotools, topotools
from clawpack.geoclaw.data import DTopoData, TopographyData

pytestmark = pytest.mark.python

KEYED = re.compile(r"^(?P<value>.*?)\s+=:\s+(?P<key>\w+)\s*(#.*)?$")

TOPO_HEADER = ["topo_missing", "test_topography", "ntopofiles"]
ENTRY = ["crop_extent", "coarsen", "buffer", "align", "x_shift", "y_shift",
         "z_shift", "negate_z", "lon_wrap_offset"]
TOPO_ENTRY = ["topo_path", "topo_type"] + ENTRY
TOPO_NC = ["var_name", "x_name", "y_name", "y_increasing", "dim_order",
           "scale_factor", "fill_value", "fill_action"]
DTOPO_ENTRY = ["dtopo_path", "dtopo_type"] + ENTRY
DTOPO_NC = ["var_name", "x_name", "y_name", "time_name", "y_increasing",
            "dim_order", "scale_factor", "t0", "dt"]


def _keyed_lines(path):
    """[(key, value)] for every line that is not blank or a # comment."""
    out = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        m = KEYED.match(s)
        assert m, f"line is not 'value =: key': {line!r}"
        out.append((m.group("key"), m.group("value").strip()))
    return out


def _keys(path):
    return [k for k, _ in _keyed_lines(path)]


def _values(path):
    return dict(_keyed_lines(path))


# --------------------------------------------------------------------------
# Input files
# --------------------------------------------------------------------------

def _topo_files(tmp_path):
    x = np.arange(-100.0, -79.99, 0.25)
    y = np.arange(20.0, 30.01, 0.25)
    X, Y = np.meshgrid(x, y)
    Z = -100.0 + 0.0 * X
    t = topotools.Topography()
    t.x, t.y, t.Z = x, y, Z
    tt3 = tmp_path / "a.tt3"
    t.write(str(tt3), topo_type=3)
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    nc = tmp_path / "a.nc"
    xr.Dataset(
        {"z": (("lat", "lon"), Z, {"units": "m"})},
        coords={"lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"})},
    ).to_netcdf(nc)
    return tt3, nc


def _dtopo_files(tmp_path):
    x = np.linspace(-99.0, -81.0, 5)
    y = np.linspace(21.0, 29.0, 5)
    X, Y = np.meshgrid(x, y)
    d = dtopotools.DTopography()
    d.x, d.y, d.X, d.Y = x, y, X, Y
    d.times = np.array([0.0, 1.0])
    d.dZ = np.zeros((2, 5, 5))
    dtt3 = tmp_path / "d.dtt3"
    d.write(str(dtt3), dtopo_type=3)
    pytest.importorskip("netCDF4")
    dnc = tmp_path / "d.nc"
    d.write(str(dnc), dtopo_type=4)
    return dtt3, dnc


def _write_topo(tmp_path, entries, num_dim=2):
    td = TopographyData(num_dim=num_dim)
    td.topofiles = entries
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    return out


def _write_dtopo(tmp_path, entries, num_dim=2):
    dd = DTopoData(num_dim=num_dim)
    dd.dtopofiles = entries
    out = tmp_path / "dtopo.data"
    dd.write(out_file=str(out))
    return out


# --------------------------------------------------------------------------
# One notation, fixed order
# --------------------------------------------------------------------------

def test_topo_data_every_line_is_keyed_in_order(tmp_path):
    tt3, nc = _topo_files(tmp_path)
    out = _write_topo(tmp_path, [[3, str(tt3)], [4, str(nc)]])
    assert _keys(out) == (TOPO_HEADER + TOPO_ENTRY + TOPO_ENTRY + TOPO_NC)


def test_dtopo_data_every_line_is_keyed_in_order(tmp_path):
    dtt3, dnc = _dtopo_files(tmp_path)
    out = _write_dtopo(tmp_path, [[3, str(dtt3)], [4, str(dnc)]])
    assert _keys(out) == (["mdtopofiles"] + DTOPO_ENTRY + DTOPO_ENTRY
                          + DTOPO_NC + ["dt_max_dtopo"])


def test_topo_data_1d_is_keyed(tmp_path):
    """1D keeps its own layout (no per-file settings), in the same notation."""
    tt3, _ = _topo_files(tmp_path)
    out = _write_topo(tmp_path, [[3, str(tt3)]], num_dim=1)
    assert _keys(out) == ["topo_missing", "test_topography", "ntopofiles",
                          "override_order", "topo_path", "topo_type"]


def test_dtopo_data_1d_is_keyed(tmp_path):
    dtt3, _ = _dtopo_files(tmp_path)
    out = _write_dtopo(tmp_path, [[3, str(dtt3)]], num_dim=1)
    assert _keys(out) == ["mdtopofiles", "dtopo_path", "dtopo_type",
                          "dt_max_dtopo"]


# --------------------------------------------------------------------------
# Values
# --------------------------------------------------------------------------

def test_unset_crop_and_align_are_none(tmp_path):
    tt3, _ = _topo_files(tmp_path)
    v = _values(_write_topo(tmp_path, [[3, str(tt3)]]))
    assert v["crop_extent"] == "None"
    assert v["align"] == "None"


def test_align_through_the_origin_is_written(tmp_path):
    """[0, 0] is a valid align, not "unset" (the old sentinel lost it)."""
    tt3, _ = _topo_files(tmp_path)
    entry = topotools.TopoFile(str(tt3), 3, crop_extent=[-95.2, -89.7, 21.7, 25.3],
                               coarsen=2, align=[0.0, 0.0])
    v = _values(_write_topo(tmp_path, [entry]))
    assert [float(a) for a in v["align"].split()] == [0.0, 0.0]
    assert [float(a) for a in v["crop_extent"].split()] == [-95.2, -89.7, 21.7, 25.3]


def test_netcdf_strings_are_quoted_and_logicals_are_t_f(tmp_path):
    """dim_order 'y,x' must survive a list-directed Fortran read, where a bare
    comma ends the value."""
    _, nc = _topo_files(tmp_path)
    v = _values(_write_topo(tmp_path, [[4, str(nc)]]))
    assert v["var_name"] == "'z'"
    assert v["dim_order"] == "'y,x'"
    assert v["y_increasing"] in ("T", "F")
    assert v["negate_z"] == "F"
    assert float(v["lon_wrap_offset"]) == 0.0


def test_netcdf_fill_value_none_when_absent(tmp_path):
    """The descriptor has a fixed set of lines: an absent fill value is
    written as None rather than left out."""
    _, nc = _topo_files(tmp_path)
    import xarray as xr
    with xr.open_dataset(nc) as ds:
        ds = ds.load()
    ds["z"].encoding["_FillValue"] = None
    nofill = tmp_path / "nofill.nc"
    ds.to_netcdf(nofill)
    v = _values(_write_topo(tmp_path, [[4, str(nofill)]]))
    assert v["fill_value"] == "None"


def test_wrapped_entries_carry_lon_wrap_offset_in_the_entry(tmp_path):
    """lon_wrap_offset is a per-entry line for every type (0 for ASCII), not
    part of the netCDF descriptor."""
    x = np.arange(-180.0, 180.01, 1.0)
    y = np.arange(-5.0, 5.01, 1.0)
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    g = tmp_path / "g.nc"
    xr.Dataset(
        {"z": (("lat", "lon"), np.zeros((y.size, x.size)), {"units": "m"})},
        coords={"lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"})},
    ).to_netcdf(g)
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [topotools.TopoFile(str(g), 4,
                                       crop_extent=[175.0, 185.0, -3.0, 3.0])]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    offsets = sorted(float(v) for k, v in _keyed_lines(out)
                     if k == "lon_wrap_offset")
    assert offsets == [0.0, 360.0]


# --------------------------------------------------------------------------
# Reading back
# --------------------------------------------------------------------------

def test_topography_data_read_round_trip(tmp_path):
    tt3, nc = _topo_files(tmp_path)
    entries = [
        topotools.TopoFile(str(tt3), 3, crop_extent=[-95.0, -90.0, 22.0, 25.0],
                           coarsen=2, buffer=1, align=[0.0, 0.0],
                           x_shift=0.0, z_shift=1.5, negate_z=True),
        topotools.TopoFile(str(nc), 4),
    ]
    td = TopographyData()
    td.override_order = True
    td.topo_missing = -12345.0
    td.topofiles = entries
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    back = TopographyData()
    back.read(str(out))
    assert back.topo_missing == -12345.0
    assert back.test_topography == 0
    assert back.topofiles == entries


def test_dtopo_data_read_round_trip(tmp_path):
    dtt3, dnc = _dtopo_files(tmp_path)
    entries = [dtopotools.DTopoFile(str(dtt3), 3, x_shift=1.0, y_shift=-0.5,
                                    z_shift=0.25, negate_z=True),
               dtopotools.DTopoFile(str(dnc), 4)]
    out = _write_dtopo(tmp_path, entries)
    back = DTopoData()
    back.read(str(out))
    assert back.dtopofiles == entries


def test_path_with_quote_and_space_round_trips(tmp_path):
    """Paths are quoted Fortran strings: an embedded quote is doubled."""
    d = tmp_path / "it's a dir"
    d.mkdir()
    tt3, _ = _topo_files(d)
    out = _write_topo(tmp_path, [[3, str(tt3)]])
    assert "'it''s a dir/" in out.read_text()
    back = TopographyData()
    back.read(str(out))
    assert back.topofiles[0].path == str(tt3)
