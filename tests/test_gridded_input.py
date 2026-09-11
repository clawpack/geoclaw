#!/usr/bin/env python
# encoding: utf-8
"""Tests for the gridded_input format-reader registry (the ingest seam).

The central guarantee is *cross-reader equivalence*: the same field written in
different on-disk formats must read back into an identical in-memory grid, so
routing ``Topography.read`` through the registry did not change what any format
produces.
"""

import numpy as np
import pytest

from clawpack.geoclaw import gridded_input, topotools


@pytest.mark.python
def test_registry_dispatch_and_unknown_type():
    assert isinstance(gridded_input.get_reader(1), gridded_input.ASCIIReader)
    assert isinstance(gridded_input.get_reader(2), gridded_input.ASCIIReader)
    assert isinstance(gridded_input.get_reader(3), gridded_input.ASCIIReader)
    assert isinstance(gridded_input.get_reader(4), gridded_input.NetCDFReader)
    assert isinstance(gridded_input.get_reader(5), gridded_input.GeoTIFFReader)
    # Sign only flips the elevation convention; dispatch is by abs(topo_type).
    assert isinstance(gridded_input.get_reader(-3), gridded_input.ASCIIReader)
    with pytest.raises(IOError):
        gridded_input.get_reader(99)


def _make_topo():
    topo = topotools.Topography()
    topo.x = np.linspace(-3.0, 5.0, 9)
    topo.y = np.linspace(10.0, 16.0, 7)
    X, Y = np.meshgrid(topo.x, topo.y)
    topo.Z = (np.sin(0.3 * X) + np.cos(0.2 * Y)).astype(float)
    return topo


@pytest.mark.python
def test_cross_reader_equivalence_ascii_vs_netcdf(tmp_path):
    """ASCII type 3 and NetCDF type 4 ingest of the same field agree exactly."""
    pytest.importorskip("netCDF4")
    ref = _make_topo()

    p3 = tmp_path / "field.tt3"
    p4 = tmp_path / "field.nc"
    ref.write(p3, topo_type=3, Z_format="%22.15e")
    ref.write(p4, topo_type=4, z_dtype="float64")

    a = topotools.Topography(str(p3), topo_type=3)
    a.read()
    b = topotools.Topography(str(p4), topo_type=4)
    b.read()

    np.testing.assert_allclose(a.x, b.x, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.y, b.y, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.Z, b.Z, rtol=0, atol=1e-12)
    # And both reproduce the source field.
    np.testing.assert_allclose(a.Z, ref.Z, rtol=0, atol=1e-12)


@pytest.mark.python
def test_cross_reader_equivalence_under_crop(tmp_path):
    """A crop_extent read agrees between the ASCII and NetCDF adapters."""
    pytest.importorskip("netCDF4")
    ref = _make_topo()

    p3 = tmp_path / "field.tt3"
    p4 = tmp_path / "field.nc"
    ref.write(p3, topo_type=3, Z_format="%22.15e")
    ref.write(p4, topo_type=4, z_dtype="float64")

    crop = [-1.0, 3.0, 11.0, 15.0]
    a = topotools.Topography(str(p3), topo_type=3)
    a.crop_extent = list(crop)
    a.read()
    b = topotools.Topography(str(p4), topo_type=4)
    b.crop_extent = list(crop)
    b.read()

    np.testing.assert_allclose(a.x, b.x, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.y, b.y, rtol=0, atol=1e-9)
    assert a.Z.shape == b.Z.shape
    np.testing.assert_allclose(a.Z, b.Z, rtol=0, atol=1e-12)
