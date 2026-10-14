#!/usr/bin/env python
# encoding: utf-8
"""Tests for the gridded_input format-reader registry (the ingest seam).

The central guarantee is *cross-reader equivalence*: the same field written in
different on-disk formats must read back into an identical in-memory grid, so
routing ``Topography.read`` through the registry did not change what any format
produces.
"""

import sys

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


# ---------------------------------------------------------------------------
# GeoTIFF (topo_type 5)
# ---------------------------------------------------------------------------

def _write_geotiff(path, z, gt, nodata=None):
    """Write *z* (file row order, row 0 first) as a one-band GeoTIFF with GDAL
    geotransform *gt*, using whichever backend is installed."""
    try:
        import rasterio
        from affine import Affine
    except ImportError:
        rasterio = None
    if rasterio is not None:
        with rasterio.open(path, "w", driver="GTiff", height=z.shape[0],
                           width=z.shape[1], count=1, dtype=z.dtype,
                           transform=Affine.from_gdal(*gt),
                           nodata=nodata) as ds:
            ds.write(z, 1)
        return
    gdal = pytest.importorskip("osgeo.gdal")
    gdal_type = {np.dtype("float32"): gdal.GDT_Float32,
                 np.dtype("float64"): gdal.GDT_Float64}[z.dtype]
    ds = gdal.GetDriverByName("GTiff").Create(str(path), z.shape[1],
                                              z.shape[0], 1, gdal_type)
    ds.SetGeoTransform(gt)
    band = ds.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    band.WriteArray(z)
    ds = None


@pytest.mark.python
@pytest.mark.parametrize("backend", ["rasterio", "osgeo"])
def test_geotiff_non_square_geometry(tmp_path, backend):
    """Pixel centers, S->N row order and nodata -> NaN on a 3x5 raster.

    Non-square on purpose: the old reader swapped the row/column counts, which
    a square raster cannot detect.
    """
    pytest.importorskip(backend)
    nrows, ncols = 3, 5
    z = (10.0 * np.arange(nrows)[:, None]
         + np.arange(ncols)[None, :]).astype("float32")
    z[0, 4] = -9999.0                       # north-east pixel is missing
    gt = (100.0, 2.0, 0.0, 50.0, 0.0, -0.5)
    path = tmp_path / "dem.tif"
    _write_geotiff(path, z, gt, nodata=-9999.0)

    x, y, Z, delta = gridded_input.read_geotiff(path, backends=(backend,))

    np.testing.assert_allclose(x, [101.0, 103.0, 105.0, 107.0, 109.0])
    np.testing.assert_allclose(y, [48.75, 49.25, 49.75])
    assert delta == (2.0, 0.5)
    assert Z.shape == (nrows, ncols)
    # Row 0 of Z is the south edge, i.e. the last row of the file.
    np.testing.assert_array_equal(Z[0], [20.0, 21.0, 22.0, 23.0, 24.0])
    np.testing.assert_array_equal(Z[2, :4], [0.0, 1.0, 2.0, 3.0])
    assert np.isnan(Z[2, 4])

    # Header-only read agrees with the full read and does not need the data.
    xh, yh, Zh, _ = gridded_input.read_geotiff(path, read_data=False,
                                               backends=(backend,))
    assert Zh is None
    np.testing.assert_array_equal(xh, x)
    np.testing.assert_array_equal(yh, y)


@pytest.mark.python
def test_geotiff_topography_read_and_header(tmp_path):
    """Topography.read() and read_header() both work for topo_type 5."""
    z = np.arange(15, dtype="float32").reshape(3, 5)
    path = tmp_path / "dem.tif"
    if not (_has("rasterio") or _has("osgeo")):
        pytest.skip("needs rasterio or osgeo")
    _write_geotiff(path, z, (0.0, 1.0, 0.0, 3.0, 0.0, -1.0))

    topo = topotools.Topography(str(path), topo_type=5)
    topo.read()
    np.testing.assert_allclose(topo.x, [0.5, 1.5, 2.5, 3.5, 4.5])
    np.testing.assert_allclose(topo.y, [0.5, 1.5, 2.5])
    np.testing.assert_array_equal(topo.Z, np.flipud(z))

    hdr = topotools.Topography(str(path), topo_type=5)
    assert hdr.read_header() == (5, 3)
    assert hdr.extent == [0.5, 4.5, 0.5, 2.5]


@pytest.mark.python
def test_cross_reader_equivalence_ascii_vs_geotiff(tmp_path):
    """A GeoTIFF whose pixel centers sit on the type 3 grid points reads back
    to the same in-memory grid as that type 3 file."""
    if not (_has("rasterio") or _has("osgeo")):
        pytest.skip("needs rasterio or osgeo")
    ref = _make_topo()
    dx = ref.x[1] - ref.x[0]
    dy = ref.y[1] - ref.y[0]

    p3 = tmp_path / "field.tt3"
    p5 = tmp_path / "field.tif"
    ref.write(p3, topo_type=3, Z_format="%22.15e")
    gt = (ref.x[0] - dx / 2, dx, 0.0, ref.y[-1] + dy / 2, 0.0, -dy)
    _write_geotiff(p5, np.flipud(ref.Z).astype("float64"), gt)

    a = topotools.Topography(str(p3), topo_type=3)
    a.read()
    b = topotools.Topography(str(p5), topo_type=5)
    b.read()

    np.testing.assert_allclose(a.x, b.x, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.y, b.y, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.Z, b.Z, rtol=0, atol=1e-12)


@pytest.mark.python
def test_geotiff_rotated_transform_refused(tmp_path):
    if not (_has("rasterio") or _has("osgeo")):
        pytest.skip("needs rasterio or osgeo")
    path = tmp_path / "rot.tif"
    _write_geotiff(path, np.zeros((2, 2), dtype="float32"),
                   (0.0, 1.0, 0.1, 2.0, 0.0, -1.0))
    with pytest.raises(ValueError, match="rotated"):
        gridded_input.read_geotiff(path)


@pytest.mark.python
def test_geotiff_no_backend_message(monkeypatch, tmp_path):
    """With neither library importable the error names both install routes."""
    monkeypatch.setitem(sys.modules, "rasterio", None)
    monkeypatch.setitem(sys.modules, "osgeo", None)
    with pytest.raises(ImportError, match="pip install rasterio"):
        gridded_input.read_geotiff(tmp_path / "missing.tif")


def _has(module):
    try:
        __import__(module)
    except ImportError:
        return False
    return True
