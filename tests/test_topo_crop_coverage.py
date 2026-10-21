#!/usr/bin/env python
# encoding: utf-8
"""A topo crop_extent must be covered by its file, for every topo type (#760).

Rule (decided 2026-10-19, dev/design/topo_input_roadmap.md §2.7): at
``make data`` (``TopographyData.write``) a ``crop_extent`` that the file does
not cover is an error, whatever the topo type.  "Covered" allows half a cell
beyond the outermost grid points, so a crop written at a tile's cell edges is
accepted.  Longitude on a geographic file that can wrap is the exception: a
type-4 crop across the antimeridian is split into two entries and stays valid.

Before the rule, types 2/3 wrote any crop and left Fortran to clip it (or stop,
when disjoint), while type 4 rejected even a crop at the tile's own edges,
because it compared against the outermost cell *centers*.
"""

import warnings

import numpy as np
import pytest

from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import TopographyData

pytestmark = pytest.mark.python

# A 0.25 x 0.25 degree regional tile with cell-centered points, like the NCEI
# CUDEM tiles in #760.  DX is a power of two so every coordinate is exact.
DX = 1.0 / 64
X_EDGES = (-124.25, -124.0)
Y_EDGES = (46.75, 47.0)

# Shared by every rejection so the tests pin one message, not three.
NOT_COVERED = "not covered by"

TOPO_TYPES = [2, 3, 4]


def _write_tile(tmp_path, topo_type, x_edges=X_EDGES, y_edges=Y_EDGES, dx=DX):
    """Write a cell-centered tile spanning x_edges x y_edges; return its path."""
    x = np.arange(x_edges[0] + dx / 2, x_edges[1], dx)
    y = np.arange(y_edges[0] + dx / 2, y_edges[1], dx)
    Z = -10.0 + np.zeros((y.size, x.size))

    path = tmp_path / f"tile.tt{topo_type}"
    if topo_type == 4:
        xr = pytest.importorskip("xarray")
        pytest.importorskip("netCDF4")
        path = tmp_path / "tile.nc"
        xr.Dataset(
            {"z": (("lat", "lon"), Z, {"units": "m"})},
            coords={
                "lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"}),
            }).to_netcdf(path)
    else:
        t = topotools.Topography()
        t.x, t.y, t.Z = x, y, Z
        t.write(str(path), topo_type=topo_type)
    return path


def _write_topo_data(tmp_path, path, topo_type, crop_extent, **attrs):
    """Run TopographyData.write for one cropped file; return topo.data text."""
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = topo_type
    topo.crop_extent = list(crop_extent)
    for name, value in attrs.items():
        setattr(topo, name, value)

    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [topo]
    out = tmp_path / "topo.data"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            td.write(out_file=str(out))
        finally:
            if getattr(td, "_out_file", None) is not None:
                td._out_file.close()
    return out.read_text()


# ---------------------------------------------------------------------------
# Accepted
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_crop_inside_file_is_accepted(tmp_path, topo_type):
    """Control: a crop well inside the tile writes for every type."""
    path = _write_tile(tmp_path, topo_type)
    text = _write_topo_data(tmp_path, path, topo_type,
                            [-124.2, -124.1, 46.8, 46.9])
    assert path.name in text


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_crop_at_tile_cell_edges_is_accepted(tmp_path, topo_type):
    """A crop at the tile's own cell edges is the tile; it must be accepted.

    The outermost points are half a cell inside these edges.  Type 4 used to
    reject this, comparing against the points rather than the cells.
    """
    path = _write_tile(tmp_path, topo_type)
    text = _write_topo_data(tmp_path, path, topo_type,
                            [*X_EDGES, *Y_EDGES])
    assert path.name in text


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_crop_is_checked_in_domain_coordinates(tmp_path, topo_type):
    """crop_extent is in domain coordinates, after x/y_shift, for every type.

    The tile shifted by +10 in x covers [-114.25, -114.0]; a crop there is
    covered even though it lies outside the file's own coordinates.  (Type 4
    was checked in file coordinates until dev/design/topo_input_roadmap.md
    §2.8; tests/test_topo_crop_frame.py covers the read side.)
    """
    path = _write_tile(tmp_path, topo_type)
    text = _write_topo_data(tmp_path, path, topo_type,
                            [-114.2, -114.1, 46.8, 46.9], x_shift=10.0)
    assert path.name in text


# ---------------------------------------------------------------------------
# Rejected
# ---------------------------------------------------------------------------

REJECTED = {
    "partly_west": [-124.3, -124.1, 46.8, 46.9],
    "partly_north": [-124.2, -124.1, 46.9, 47.22],
    "randy_760": [-124.25, -124.1, 46.95, 47.22],
    "disjoint": [-124.25, -124.1, 47.05, 47.22],
    # Just past the half-cell allowance on one side only.
    "edge_plus_0.6_cell": [X_EDGES[0] - 0.6 * DX, X_EDGES[1], *Y_EDGES],
}


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
@pytest.mark.parametrize("case", sorted(REJECTED))
def test_crop_not_covered_is_rejected(tmp_path, topo_type, case):
    """Partly or wholly outside the file is one error, for every type.

    Types 2/3 used to write these and let Fortran clip silently (or stop, if
    disjoint); type 4 raised one of three different messages.
    """
    path = _write_tile(tmp_path, topo_type)
    with pytest.raises(ValueError, match=NOT_COVERED) as err:
        _write_topo_data(tmp_path, path, topo_type, REJECTED[case])
    assert path.name in str(err.value)


@pytest.mark.parametrize("topo_type", TOPO_TYPES)
def test_rejection_reports_file_extent_and_overlap(tmp_path, topo_type):
    """The message gives the cell-edge extent of the file and the overlap,
    which is what the user needs to fix the crop."""
    path = _write_tile(tmp_path, topo_type)
    with pytest.raises(ValueError, match=NOT_COVERED) as err:
        _write_topo_data(tmp_path, path, topo_type,
                         [-124.25, -124.1, 46.95, 47.22])
    msg = str(err.value)
    assert "-124.25" in msg and "47.0" in msg      # file cell-edge extent
    assert "overlap" in msg and "46.95" in msg     # the usable overlap


@pytest.mark.parametrize("topo_type", [2, 3])
def test_ascii_crop_across_the_antimeridian_is_rejected(tmp_path, topo_type):
    """ASCII types have no antimeridian split at make data, so a crop past
    the file's longitude range is not covered.  It used to be written and then
    clipped by Fortran, silently losing the part across the seam."""
    path = _write_tile(tmp_path, topo_type, x_edges=(-180.0, 180.0),
                       y_edges=(-40.0, 0.0), dx=1.0)
    with pytest.raises(ValueError, match=NOT_COVERED):
        _write_topo_data(tmp_path, path, topo_type,
                         [-190.0, -150.0, -30.0, -10.0])


# ---------------------------------------------------------------------------
# Antimeridian: type 4 still splits
# ---------------------------------------------------------------------------

def test_netcdf_crop_across_the_antimeridian_still_splits(tmp_path):
    """A geographic netCDF file that can wrap covers a cross-seam crop through
    two entries.  The coverage rule must not reject it."""
    path = _write_tile(tmp_path, 4, x_edges=(-180.0, 180.0),
                       y_edges=(-40.0, 0.0), dx=1.0)
    text = _write_topo_data(tmp_path, path, 4, [-190.0, -150.0, -30.0, -10.0])
    assert text.count("lon_wrap_offset") == 2


def test_netcdf_crop_across_the_antimeridian_checks_latitude(tmp_path):
    """Longitude wraps; latitude never does."""
    path = _write_tile(tmp_path, 4, x_edges=(-180.0, 180.0),
                       y_edges=(-40.0, 0.0), dx=1.0)
    with pytest.raises(ValueError, match=NOT_COVERED):
        _write_topo_data(tmp_path, path, 4, [-190.0, -150.0, -30.0, 5.0])


def test_netcdf_crop_at_the_points_keeps_the_outermost_points(tmp_path):
    """A crop equal to a point-registered file's point extent must keep the
    outermost points.  Clipping to a point rebuilt as ``(x0 - dx/2) + dx/2``
    once landed a round-off inside ``x0``, Fortran's ``xlocs >= crop_bounds``
    then dropped the first column, and a coarsened read flipped parity (caught
    by examples/tsunami/bowl-slosh-netcdf[coarsened]).

    topo.data no longer has crop_bounds (dev/design/topo_input_roadmap.md
    §2.8): the crop reaches Fortran as written, on the crop_extent line, so
    no clipped value can round inside the outermost points.  Fortran's own
    comparison has a 1e-6*dx tolerance, as for types 2/3; the example above
    still covers that end to end."""
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    x = np.linspace(-3.1, 3.1, 619)
    y = np.linspace(-3.5, 2.5, 599)
    path = tmp_path / "bowl.nc"
    xr.Dataset(
        {"z": (("lat", "lon"), np.zeros((y.size, x.size)), {"units": "m"})},
        coords={"lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"})},
    ).to_netcdf(path)

    text = _write_topo_data(tmp_path, path, 4, [-3.1, 3.1, -3.5, 2.5])
    assert "crop_bounds" not in text
    line = next(l for l in text.splitlines() if "# crop_extent" in l)
    bounds = [float(v) for v in line.split("#")[0].split()]
    assert bounds[0] <= x[0] and bounds[1] >= x[-1], bounds
    assert bounds[2] <= y[0] and bounds[3] >= y[-1], bounds
