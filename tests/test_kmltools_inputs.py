#!/usr/bin/env python
# encoding: utf-8
"""kmltools.make_input_data_kmls with Topography / DTopography entries.

make_input_data_kmls read topofiles and dtopofiles only in the legacy list
form ([topo_type, path]).  A setrun that appends Topography or DTopography
objects -- the current style -- failed with "TypeError: 'Topography' object is
not subscriptable", outside any try, so no kml was written at all.  The topo
box should also show the region actually used: the crop_extent when one is set
(#760), the file extent shifted by x/y_shift otherwise.
"""

import re
import warnings

import numpy as np
import pytest

from clawpack.clawutil.data import ClawRunData
from clawpack.geoclaw import dtopotools, kmltools, topotools

pytestmark = pytest.mark.python

# A small tile, cell-centered points on [-124.25, -124.0] x [46.75, 47.0].
DX = 1.0 / 64
X = np.arange(-124.25 + DX / 2, -124.0, DX)
Y = np.arange(46.75 + DX / 2, 47.0, DX)


def _write_tt3(path):
    t = topotools.Topography()
    t.x, t.y, t.Z = X, Y, np.zeros((Y.size, X.size))
    t.write(str(path), topo_type=3)
    return path


def _write_dtt3(path):
    d = dtopotools.DTopography()
    d.x, d.y = X, Y
    d.X, d.Y = np.meshgrid(X, Y)
    d.times = np.array([0.0, 1.0])
    d.dZ = np.zeros((2, Y.size, X.size))
    d.write(str(path), dtopo_type=3)
    return path


def _rundata():
    rundata = ClawRunData('geoclaw', 2)
    rundata.topo_data.topofiles = []
    rundata.dtopo_data.dtopofiles = []
    return rundata


def _box(kml_path):
    """The (x1, x2, y1, y2) box2kml wrote into the description."""
    text = kml_path.read_text()
    m = re.search(r"x1 = (\S+), x2 = (\S+)\s+y1 = (\S+), y2 = (\S+?)<", text)
    if m is None:
        m = re.search(r"x1 = (\S+), x2 = (\S+)\s+y1 = (\S+), y2 = (\S+)",
                      text)
    assert m is not None, text
    return tuple(float(v) for v in m.groups())


def test_topography_object_without_crop_draws_file_extent(tmp_path,
                                                          monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = _write_tt3(tmp_path / "tile.tt3")
    rundata = _rundata()
    topo = topotools.Topography(path=str(path), topo_type=3)
    rundata.topo_data.topofiles.append(topo)

    kmltools.make_input_data_kmls(rundata)

    assert _box(tmp_path / "tile.kml") == pytest.approx(
        (X[0], X[-1], Y[0], Y[-1]), abs=1e-6)


def test_topography_object_with_crop_draws_the_crop(tmp_path, monkeypatch):
    """With crop_extent set, the box is the region used, not the file."""
    monkeypatch.chdir(tmp_path)
    path = _write_tt3(tmp_path / "tile.tt3")
    rundata = _rundata()
    topo = topotools.Topography(path=str(path), topo_type=3)
    topo.crop_extent = [-124.2, -124.1, 46.8, 46.9]
    rundata.topo_data.topofiles.append(topo)

    kmltools.make_input_data_kmls(rundata)

    assert _box(tmp_path / "tile.kml") == pytest.approx(
        (-124.2, -124.1, 46.8, 46.9), abs=1e-6)


def test_topography_object_with_shift_draws_shifted_extent(tmp_path,
                                                           monkeypatch):
    """x/y_shift move the file into domain coordinates; so does the box."""
    monkeypatch.chdir(tmp_path)
    path = _write_tt3(tmp_path / "tile.tt3")
    rundata = _rundata()
    topo = topotools.Topography(path=str(path), topo_type=3)
    topo.x_shift = 10.0
    topo.y_shift = -1.0
    rundata.topo_data.topofiles.append(topo)

    kmltools.make_input_data_kmls(rundata)

    assert _box(tmp_path / "tile.kml") == pytest.approx(
        (X[0] + 10.0, X[-1] + 10.0, Y[0] - 1.0, Y[-1] - 1.0), abs=1e-6)


def test_legacy_list_entry_still_works(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = _write_tt3(tmp_path / "tile.tt3")
    rundata = _rundata()
    rundata.topo_data.topofiles.append([3, str(path)])

    # The deprecation is TopographyData.write's to report at make data; the
    # kml pass must not repeat it.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        kmltools.make_input_data_kmls(rundata)

    assert _box(tmp_path / "tile.kml") == pytest.approx(
        (X[0], X[-1], Y[0], Y[-1]), abs=1e-6)


def test_dtopography_object_draws_its_extent(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = _write_dtt3(tmp_path / "source.dtt3")
    rundata = _rundata()
    dtopo = dtopotools.DTopography(path=str(path), dtopo_type=3)
    rundata.dtopo_data.dtopofiles.append(dtopo)

    kmltools.make_input_data_kmls(rundata)

    assert _box(tmp_path / "source.kml") == pytest.approx(
        (X[0], X[-1], Y[0], Y[-1]), abs=1e-6)
