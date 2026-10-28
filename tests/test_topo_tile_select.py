#!/usr/bin/env python
# encoding: utf-8
"""Pick the tiles a crop uses from a predefined tile list (#765 review, #760).

A setrun can list every tile that might be needed (all NCEI tiles along a
coast, say) and one crop region; ``TopographyData.select_tiles`` decides which
tiles the crop uses.  ``crop_extent`` then takes one of:

* ``None`` or ``'all'`` -- the whole file (``'all'`` is the explicit spelling);
* ``[x1, x2, y1, y2]`` -- this window, under the per-entry coverage rule of
  test_topo_crop_coverage.py (a crop the file does not cover is an error);
* ``'none'`` -- this tile is not used.

``select_tiles`` reads headers only, marks a tile that misses the crop
``'none'``, clips a tile that partly overlaps to the overlap (so the clipped
crop passes the per-entry rule), and raises when the tiles together leave a
gap in the crop.  Until the topo.data format change (PR B), ``'none'`` entries
are left out of topo.data -- the current Fortran reader cannot parse them --
and ``make data`` reports them.
"""

import warnings

import numpy as np
import pytest

from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import TopographyData

pytestmark = pytest.mark.python

# 0.25 degree tiles with cell-centered points, like the NCEI CUDEM tiles in
# #760.  DX is a power of two so every coordinate is exact.
DX = 1.0 / 64
SIZE = 0.25

# A 2 x 2 block of tiles with south-west corners at these points.
SW = {"sw": (-124.25, 46.75), "se": (-124.0, 46.75),
      "nw": (-124.25, 47.0), "ne": (-124.0, 47.0)}

# Inside the two southern tiles only, crossing the seam between them.
SOUTH_CROP = [-124.1, -123.9, 46.8, 46.9]


def _write_tile(tmp_path, name, sw, topo_type=3):
    """Write a SIZE x SIZE cell-centered tile with south-west corner *sw*."""
    x = np.arange(sw[0] + DX / 2, sw[0] + SIZE, DX)
    y = np.arange(sw[1] + DX / 2, sw[1] + SIZE, DX)
    Z = -10.0 + np.zeros((y.size, x.size))
    if topo_type == 4:
        xr = pytest.importorskip("xarray")
        pytest.importorskip("netCDF4")
        path = tmp_path / f"{name}.nc"
        xr.Dataset(
            {"z": (("lat", "lon"), Z, {"units": "m"})},
            coords={
                "lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"}),
            }).to_netcdf(path)
    else:
        path = tmp_path / f"{name}.tt{topo_type}"
        t = topotools.Topography()
        t.x, t.y, t.Z = x, y, Z
        t.write(str(path), topo_type=topo_type)
    return path


def _entry(path, topo_type, **attrs):
    topo = topotools.Topography()
    topo.path = str(path)
    topo.topo_type = topo_type
    for name, value in attrs.items():
        setattr(topo, name, value)
    return topo


def _tile_data(tmp_path, names=("sw", "se", "nw", "ne"), topo_type=3,
               **attrs):
    """A TopographyData listing the named tiles of the 2 x 2 block."""
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [
        _entry(_write_tile(tmp_path, n, SW[n], topo_type), topo_type,
               **attrs)
        for n in names]
    return td


def _write(tmp_path, td):
    """Run TopographyData.write; return the topo.data text."""
    out = tmp_path / "topo.data"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            td.write(out_file=str(out))
        finally:
            if getattr(td, "_out_file", None) is not None:
                td._out_file.close()
    return out.read_text()


def _ntopofiles(text):
    line = next(l for l in text.splitlines() if "ntopofiles" in l)
    return int(line.split()[0])


def _crop(td, name):
    return next(t.crop_extent for t in td.topofiles
                if t.path.endswith(name))


# ---------------------------------------------------------------------------
# crop_extent = 'all' / 'none' on a single entry
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("topo_type", [3, 4])
def test_all_is_the_whole_file_like_none(tmp_path, topo_type):
    """'all' writes the same topo.data as the default None."""
    path = _write_tile(tmp_path, "sw", SW["sw"], topo_type)

    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [_entry(path, topo_type)]
    expected = _write(tmp_path, td)

    td.topofiles = [_entry(path, topo_type, crop_extent="all")]
    assert _write(tmp_path, td) == expected


@pytest.mark.parametrize("topo_type", [3, 4])
def test_none_entry_is_left_out_of_topo_data(tmp_path, topo_type):
    """A 'none' tile is not written and not counted; the others are written
    exactly as without it."""
    used = _write_tile(tmp_path, "sw", SW["sw"], topo_type)
    skipped = _write_tile(tmp_path, "ne", SW["ne"], topo_type)

    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [_entry(used, topo_type)]
    expected = _write(tmp_path, td)

    td.topofiles = [_entry(used, topo_type),
                    _entry(skipped, topo_type, crop_extent="none")]
    text = _write(tmp_path, td)
    assert skipped.name not in text.split("ntopofiles")[1]
    assert _ntopofiles(text) == 1
    assert text == expected


def test_none_entry_file_need_not_exist(tmp_path):
    """A 'none' tile is never opened, so a run can be shared without it."""
    used = _write_tile(tmp_path, "sw", SW["sw"])
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [_entry(used, 3),
                    _entry(tmp_path / "not_here.tt3", 3, crop_extent="none")]
    assert _ntopofiles(_write(tmp_path, td)) == 1


def test_make_data_reports_none_entries(tmp_path, capsys):
    """Until topo.data can carry 'none' (PR B), make data names the skipped
    tiles so the run's output still records them."""
    used = _write_tile(tmp_path, "sw", SW["sw"])
    skipped = _write_tile(tmp_path, "ne", SW["ne"])
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [_entry(used, 3), _entry(skipped, 3, crop_extent="none")]
    _write(tmp_path, td)
    out = capsys.readouterr()
    assert str(skipped) in out.out + out.err


def test_only_none_entries_is_an_error(tmp_path):
    """Every tile skipped leaves no topography at all."""
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [_entry(_write_tile(tmp_path, "sw", SW["sw"]), 3,
                           crop_extent="none")]
    with pytest.raises(ValueError, match="none"):
        _write(tmp_path, td)


@pytest.mark.parametrize("bad", ["All", "nothing", ""])
def test_other_strings_are_rejected(tmp_path, bad):
    """Only 'all' and 'none' are accepted, spelled exactly."""
    td = _tile_data(tmp_path, names=("sw",), crop_extent=bad)
    with pytest.raises(ValueError, match="crop_extent"):
        _write(tmp_path, td)


# ---------------------------------------------------------------------------
# select_tiles
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("topo_type", [3, 4])
def test_select_tiles_marks_and_clips(tmp_path, topo_type):
    """Tiles the crop misses become 'none'; tiles it crosses get the crop
    clipped to their own cells, so make data then accepts every entry."""
    td = _tile_data(tmp_path, topo_type=topo_type)
    td.select_tiles(SOUTH_CROP)

    ext = "nc" if topo_type == 4 else f"tt{topo_type}"
    assert _crop(td, f"nw.{ext}") == "none"
    assert _crop(td, f"ne.{ext}") == "none"
    np.testing.assert_allclose(_crop(td, f"sw.{ext}"),
                               [-124.1, -124.0, 46.8, 46.9])
    np.testing.assert_allclose(_crop(td, f"se.{ext}"),
                               [-124.0, -123.9, 46.8, 46.9])

    text = _write(tmp_path, td)
    assert _ntopofiles(text) == 2


def test_select_tiles_keeps_the_list_and_other_settings(tmp_path):
    """The list keeps every tile, in order, and select_tiles changes only
    crop_extent: coarsen, buffer, shifts stay as the user set them."""
    td = _tile_data(tmp_path, coarsen=2, buffer=1, z_shift=0.5)
    before = [t.path for t in td.topofiles]
    td.select_tiles(SOUTH_CROP)
    assert [t.path for t in td.topofiles] == before
    for t in td.topofiles:
        assert (t.coarsen, t.buffer, t.z_shift) == (2, 1, 0.5)


def test_select_tiles_reads_only_headers(tmp_path, monkeypatch):
    """Tile lists can name many large files: only headers may be read."""
    td = _tile_data(tmp_path)

    def no_full_read(self, *args, **kwargs):
        raise AssertionError("Topography.read called by select_tiles")

    monkeypatch.setattr(topotools.Topography, "read", no_full_read)
    td.select_tiles(SOUTH_CROP)


def test_select_tiles_gap_is_an_error(tmp_path):
    """Without the south-east tile the crop has a hole the other tiles cannot
    fill: an error naming the uncovered part, not a silent hole."""
    td = _tile_data(tmp_path, names=("sw", "nw", "ne"))
    with pytest.raises(ValueError, match="not covered"):
        td.select_tiles(SOUTH_CROP)


def test_select_tiles_crop_outside_every_tile_is_an_error(tmp_path):
    td = _tile_data(tmp_path)
    with pytest.raises(ValueError, match="not covered"):
        td.select_tiles([-120.0, -119.9, 46.8, 46.9])


def test_select_tiles_edge_contact_is_not_a_use(tmp_path):
    """A crop ending exactly on a tile's edge does not use that tile: it is
    marked 'none', not given a zero-width crop."""
    td = _tile_data(tmp_path, names=("sw", "se"))
    td.select_tiles([-124.2, -124.0, 46.8, 46.9])
    assert _crop(td, "se.tt3") == "none"
    np.testing.assert_allclose(_crop(td, "sw.tt3"),
                               [-124.2, -124.0, 46.8, 46.9])


def test_select_tiles_uses_domain_coordinates(tmp_path):
    """The crop is in domain coordinates: a tile's x_shift moves it."""
    td = _tile_data(tmp_path, names=("sw", "se"), x_shift=360.0)
    crop = [SOUTH_CROP[0] + 360.0, SOUTH_CROP[1] + 360.0] + SOUTH_CROP[2:]
    td.select_tiles(crop)
    np.testing.assert_allclose(_crop(td, "sw.tt3"),
                               [235.9, 236.0, 46.8, 46.9])
    np.testing.assert_allclose(_crop(td, "se.tt3"),
                               [236.0, 236.1, 46.8, 46.9])


def test_select_tiles_mixed_types(tmp_path):
    """ASCII and netCDF tiles in one list are treated alike."""
    pytest.importorskip("xarray")
    td = TopographyData()
    td.coordinate_system = 2
    td.topofiles = [
        _entry(_write_tile(tmp_path, "sw", SW["sw"], 3), 3),
        _entry(_write_tile(tmp_path, "se", SW["se"], 4), 4),
        _entry(_write_tile(tmp_path, "ne", SW["ne"], 4), 4)]
    td.select_tiles(SOUTH_CROP)
    assert _crop(td, "ne.nc") == "none"
    np.testing.assert_allclose(_crop(td, "se.nc"),
                               [-124.0, -123.9, 46.8, 46.9])
    assert _ntopofiles(_write(tmp_path, td)) == 2


def test_select_tiles_leaves_hand_set_none(tmp_path):
    """A tile the user switched off stays off, even if the crop meets it --
    and then the crop has a gap only that tile would fill."""
    td = _tile_data(tmp_path, names=("sw", "se"))
    td.topofiles[1].crop_extent = "none"
    with pytest.raises(ValueError, match="not covered"):
        td.select_tiles(SOUTH_CROP)
    assert td.topofiles[1].crop_extent == "none"


# ---------------------------------------------------------------------------
# Python read agrees with make data
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("topo_type", [3, 4])
def test_python_read_crop_missing_file_is_an_error(tmp_path, topo_type):
    """Reading a crop the file does not touch used to return the whole tile
    (type 4) silently; make data rejects the same crop, so the read does
    too."""
    path = _write_tile(tmp_path, "sw", SW["sw"], topo_type)
    with pytest.raises(ValueError):
        topotools.Topography(str(path), topo_type,
                             crop_extent=[-120.0, -119.9, 46.8, 46.9])
