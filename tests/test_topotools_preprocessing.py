"""
Tests for Topography preprocessing attributes and TopoData priority ordering.

Existing tests modified
-----------------------
tests/test_topotools.py
  test_read_write_topo_bowl: removed topo_type=1 from parametrize — the
    type-1 write and read now emit DeprecationWarning; covered by Group 7 here.
  test_read_write_topo_bowl_hill: same reason.

Behaviors already tested (not duplicated here)
----------------------------------------------
- Topography.__init__ basic attribute setting (topo_func, topo_type, path)
- read/write round-trip for topo_type 2 and 3
- crop() geometry: coordinate values, extent, Z values after crop
- plot() smoke tests
- Unstructured point interpolation
- NetCDF round-trip via etopo1 / kahului fixture (test_topotools.py)
- TopoInspector fill / unit / CF detection (tests/netcdf/)

New groups with no prior coverage
----------------------------------
Group 1:  preprocessing attribute defaults
Group 2:  individual preprocessing operations (negate_z, z_shift, x_shift, crop, coarsen, buffer)
Group 3:  operation order (negate before shift; shifts before crop)
Group 4:  negate_z vs topo_type < 0 interaction
Group 5:  stride × coarsen interaction (NetCDF type 4)
Group 6:  read_header() lazy loading for topo_type 4
Group 7:  topo_type=1 deprecation (read warn, preprocessing guard, write warn, read_header IOError)
Group 8:  TopoData._normalize_topofiles
Group 9:  _compute_priority_order
Group 10: TopoData.write() output format (sentinels, 3-line header, per-file block)
Group 11: backward-compat round-trip for legacy list/dict topofile entries
"""

from __future__ import annotations

import re
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools
from clawpack.geoclaw.topotools import Topography
from clawpack.geoclaw.data import TopographyData


# ===========================================================================
# Constants and helpers
# ===========================================================================

# Analytic Z for a 10×10 grid: Z[i, j] = i + 10*j (integer, exact arithmetic)
_NX = 10  # number of x (column) points
_NY = 10  # number of y (row) points
_ORIGIN_X = 0.0
_ORIGIN_Y = 0.0
_DELTA = 1.0
_TOPO_MISSING = 99999.0

# NetCDF coordinate name variants to parametrize over.
_COORD_VARIANTS = [
    ("lon", "lat"),
    ("longitude", "latitude"),
    ("x", "y"),
]
_COORD_IDS = ["lon/lat", "longitude/latitude", "x/y"]

# The 8 preprocessing attributes and a non-default value for each.
_PREPROCESSING_NONDEFAULTS = [
    ("crop_extent", [1.0, 8.0, 1.0, 8.0]),
    ("coarsen", 2),
    ("buffer", 1),
    ("align", (0.0, 0.0)),
    ("x_shift", 5.0),
    ("y_shift", 7.0),
    ("z_shift", 10.0),
    ("negate_z", True),
]


# ===========================================================================
# Shared fixtures
# ===========================================================================

def _analytic_Z() -> np.ndarray:
    """Return the 10×10 analytic Z array: Z[i, j] = i + 10*j."""
    i_idx = np.arange(_NY, dtype=np.float64).reshape(_NY, 1)
    j_idx = np.arange(_NX, dtype=np.float64).reshape(1, _NX)
    return i_idx + 10.0 * j_idx


def _write_tt2(path: Path, Z: np.ndarray | None = None,
               cellsize: float = _DELTA) -> Path:
    """Write a minimal topo_type=2 file.  Z defaults to _analytic_Z()."""
    if Z is None:
        Z = _analytic_Z()
    ny, nx = Z.shape
    with open(path, "w") as f:
        f.write(f"{nx}          ncols\n")
        f.write(f"{ny}          nrows\n")
        f.write(f"{_ORIGIN_X}  xllcenter\n")
        f.write(f"{_ORIGIN_Y}  yllcenter\n")
        f.write(f"{cellsize}     cellsize\n")
        f.write(f"{_TOPO_MISSING}   nodata_value\n")
        # type-2: one value per line, written NW→SE (top row first, flipped)
        for row in np.flipud(Z):
            for val in row:
                f.write(f"{val:.1f}\n")
    return path


@pytest.fixture
def tt2_path(tmp_path):
    """Path to a type-2 10×10 topo file with Z[i,j] = i + 10*j."""
    return _write_tt2(tmp_path / "analytic.tt2")


@pytest.fixture
def tt2_path_with_missing(tmp_path):
    """Path to a type-2 10×10 file with Z[5,5] = topo_missing sentinel."""
    Z = _analytic_Z()
    Z[5, 5] = _TOPO_MISSING
    return _write_tt2(tmp_path / "analytic_missing.tt2", Z)


@pytest.fixture
def topo_direct():
    """In-memory Topography with _x/_y/_Z set directly (no file I/O)."""
    t = Topography()
    t._x = np.arange(_NX, dtype=np.float64)
    t._y = np.arange(_NY, dtype=np.float64)
    t._Z = _analytic_Z()
    t.topo_type = 2
    return t


def _make_nc_topo(path: Path, lon_name: str, lat_name: str,
                  lat_south_to_north: bool = True) -> Path:
    """Write a CF-compliant NetCDF topo file using xarray."""
    xr = pytest.importorskip("xarray")
    np_mod = np

    lons = np_mod.arange(_NX, dtype=np.float64) + _ORIGIN_X
    lats = np_mod.arange(_NY, dtype=np.float64) + _ORIGIN_Y
    if not lat_south_to_north:
        lats = lats[::-1]

    # Z[i,j] = i + 10*j, where i=row(lat), j=col(lon)
    # In memory: Z shape = (nlat, nlon) following NetCDF convention
    lat_idx = np_mod.arange(_NY, dtype=np.float64).reshape(_NY, 1)
    lon_idx = np_mod.arange(_NX, dtype=np.float64).reshape(1, _NX)
    data = lat_idx + 10.0 * lon_idx  # float64, shape (10, 10)

    # units / standard_name attributes for CF detection
    lon_attrs: dict = {"units": "degrees_east"}
    lat_attrs: dict = {"units": "degrees_north"}
    if lon_name in ("x",):
        lon_attrs["standard_name"] = "longitude"
    if lat_name in ("y",):
        lat_attrs["standard_name"] = "latitude"

    coords = {
        lon_name: xr.DataArray(lons, dims=[lon_name], attrs=lon_attrs),
        lat_name: xr.DataArray(lats, dims=[lat_name], attrs=lat_attrs),
    }
    da = xr.DataArray(
        data, dims=[lat_name, lon_name], coords=coords,
        attrs={"units": "m", "positive": "up"},
    )
    ds = xr.Dataset({"elevation": da})
    ds.to_netcdf(path)
    return path


@pytest.fixture(params=_COORD_VARIANTS, ids=_COORD_IDS)
def nc_topo_path(request, tmp_path):
    """NetCDF topo file parametrized over 3 coordinate name variants."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    lon_name, lat_name = request.param
    fname = f"topo_{lon_name}_{lat_name}.nc"
    return _make_nc_topo(tmp_path / fname, lon_name, lat_name), lon_name, lat_name


# ===========================================================================
# Group 1 — Preprocessing attribute defaults
# ===========================================================================

def test_preprocessing_defaults_all_attributes():
    t = Topography()
    assert t.crop_extent is None
    assert t.coarsen == 1
    assert t.buffer == 0.0
    assert t.align is None
    assert t.x_shift == 0.0
    assert t.y_shift == 0.0
    assert t.z_shift == 0.0
    assert t.negate_z is False


def test_preprocessing_noop_with_defaults(tt2_path):
    t_ref = Topography()
    t_ref.read(tt2_path, topo_type=2)

    t_proc = Topography()
    t_proc.read(tt2_path, topo_type=2)  # all preprocessing attrs at defaults

    np.testing.assert_array_equal(t_proc._x, t_ref._x)
    np.testing.assert_array_equal(t_proc._y, t_ref._y)
    np.testing.assert_array_equal(t_proc.Z, t_ref.Z)


# ===========================================================================
# Group 2 — Individual preprocessing operations
# ===========================================================================

def test_preprocessing_negate_z_inverts_sign(tt2_path):
    t = Topography()
    t.negate_z = True
    t.read(tt2_path, topo_type=2)

    expected = -_analytic_Z()
    np.testing.assert_array_equal(t.Z, expected)
    np.testing.assert_array_equal(t._x, np.arange(_NX, dtype=np.float64))
    np.testing.assert_array_equal(t._y, np.arange(_NY, dtype=np.float64))


def test_preprocessing_z_shift_offsets_valid_cells(tt2_path_with_missing):
    """z_shift shifts only real cells; missing cells are NaN in memory and stay
    NaN under the offset (the numeric sentinel is only used on file)."""
    t = Topography()
    t.z_shift = 5.0
    t.no_data_value = _TOPO_MISSING
    t.read(tt2_path_with_missing, topo_type=2)

    Z_expected = _analytic_Z() + 5.0

    # Non-missing cells (everything but the injected missing cell at [5, 5])
    # are shifted by z_shift.
    real = np.ones(t.Z.shape, dtype=bool)
    real[5, 5] = False
    np.testing.assert_allclose(t.Z[real], Z_expected[real], rtol=1e-12)
    # The missing cell is NaN in memory: not shifted, and the numeric sentinel
    # does not leak into the array.
    assert np.isnan(t.Z[5, 5])
    assert not np.any(t.Z == _TOPO_MISSING)


def test_preprocessing_z_shift_negative(tt2_path):
    t = Topography()
    t.z_shift = -3.0
    t.read(tt2_path, topo_type=2)

    expected = _analytic_Z() - 3.0
    np.testing.assert_allclose(t.Z, expected, rtol=1e-12)


def test_preprocessing_x_shift_translates_coordinates(tt2_path):
    original_x = np.arange(_NX, dtype=np.float64)

    t = Topography()
    t.x_shift = 10.0
    t.read(tt2_path, topo_type=2)

    np.testing.assert_allclose(t._x, original_x + 10.0, rtol=1e-12)
    np.testing.assert_array_equal(t._y, np.arange(_NY, dtype=np.float64))
    np.testing.assert_array_equal(t.Z, _analytic_Z())
    # x_shift modifies _x permanently; verify it is not a view of a fresh array
    assert not np.shares_memory(t._x, original_x)


def test_preprocessing_y_shift_translates_coordinates(tt2_path):
    original_y = np.arange(_NY, dtype=np.float64)

    t = Topography()
    t.y_shift = 10.0
    t.read(tt2_path, topo_type=2)

    np.testing.assert_allclose(t._y, original_y + 10.0, rtol=1e-12)
    np.testing.assert_array_equal(t._x, np.arange(_NX, dtype=np.float64))
    np.testing.assert_array_equal(t.Z, _analytic_Z())
    # y_shift modifies _y permanently; verify it is not a view of a fresh array
    assert not np.shares_memory(t._y, original_y)


def test_preprocessing_crop_extent_clips_domain(tt2_path):
    t = Topography()
    t.crop_extent = [2.0, 7.0, 3.0, 8.0]
    t.read(tt2_path, topo_type=2)

    assert t.extent[0] >= 2.0
    assert t.extent[1] <= 7.0
    assert t.extent[2] >= 3.0
    assert t.extent[3] <= 8.0
    assert t.Z.shape[1] == len(t._x)
    assert t.Z.shape[0] == len(t._y)


def test_preprocessing_coarsen_subsamples_not_averages(tt2_path):
    t = Topography()
    t.coarsen = 2
    t.read(tt2_path, topo_type=2)

    assert t.Z.shape == (5, 5)
    # Coarsening uses stride-2 subsampling, NOT averaging
    Z_orig = _analytic_Z()
    for i in range(5):
        for j in range(5):
            np.testing.assert_array_equal(t.Z[i, j], Z_orig[2 * i, 2 * j])
    # Explicit spot checks from the prompt
    assert t.Z[0, 0] == 0.0   # Z_orig[0,0] = 0+10*0
    assert t.Z[1, 0] == 2.0   # Z_orig[2,0] = 2+10*0, NOT average of rows 0-1


def test_preprocessing_buffer_expands_crop_region(tt2_path):
    t_no_buf = Topography()
    t_no_buf.crop_extent = [3.0, 7.0, 3.0, 7.0]
    t_no_buf.read(tt2_path, topo_type=2)

    t_buf = Topography()
    t_buf.crop_extent = [3.0, 7.0, 3.0, 7.0]
    t_buf.buffer = 2
    t_buf.read(tt2_path, topo_type=2)

    # Buffered result must be strictly wider than unbuffered in both dimensions
    assert t_buf.extent[0] < t_no_buf.extent[0] or t_buf.extent[1] > t_no_buf.extent[1]
    assert t_buf.extent[2] < t_no_buf.extent[2] or t_buf.extent[3] > t_no_buf.extent[3]


def test_preprocessing_buffer_float_truncates_to_int(tt2_path):
    """buffer is integer grid points; float values are truncated via int().

    buffer=0.5 is truncated to int(0.5)=0, so no expansion occurs.
    """
    t_zero = Topography()
    t_zero.crop_extent = [3.0, 7.0, 3.0, 7.0]
    t_zero.buffer = 0
    t_zero.read(tt2_path, topo_type=2)

    t_float = Topography()
    t_float.crop_extent = [3.0, 7.0, 3.0, 7.0]
    t_float.buffer = 0.5  # int(0.5) == 0
    t_float.read(tt2_path, topo_type=2)

    np.testing.assert_array_equal(t_float._x, t_zero._x)
    np.testing.assert_array_equal(t_float._y, t_zero._y)
    np.testing.assert_array_equal(t_float.Z, t_zero.Z)


# ===========================================================================
# Group 3 — Operation order
# ===========================================================================

def test_preprocessing_order_negate_before_z_shift(tt2_path):
    """negate_z fires before z_shift: result is -(Z_orig) + z_shift."""
    t = Topography()
    t.negate_z = True
    t.z_shift = 1.0
    t.read(tt2_path, topo_type=2)

    # Any cell where Z != 0 distinguishes negate-then-shift from shift-then-negate
    Z_orig = _analytic_Z()
    expected = -Z_orig + 1.0
    np.testing.assert_allclose(t.Z, expected, rtol=1e-12)

    # Verify the two orderings differ where Z != 0
    wrong_order = -(Z_orig + 1.0)
    assert not np.allclose(expected[Z_orig != 0], wrong_order[Z_orig != 0])


def test_preprocessing_order_shifts_before_crop(tt2_path):
    """z_shift is applied before crop_extent; the cropped region sees shifted values."""
    t = Topography()
    t.z_shift = 100.0
    t.crop_extent = [2.0, 7.0, 3.0, 8.0]
    t.read(tt2_path, topo_type=2)

    # All Z values in the cropped result must equal original_Z + 100.0
    # Reconstruct expected values for the cropped x/y range
    xv = np.round(t._x).astype(int)
    yv = np.round(t._y).astype(int)
    for ri, row_y in enumerate(yv):
        for ci, col_x in enumerate(xv):
            expected_val = float(row_y + 10 * col_x) + 100.0
            np.testing.assert_allclose(t.Z[ri, ci], expected_val, rtol=1e-12)


# ===========================================================================
# Group 4 — negate_z vs topo_type < 0 interaction
# ===========================================================================

def test_preprocessing_negative_topotype_negates_z(tt2_path):
    """topo_type < 0 negates Z via the existing sign convention.

    This is the pre-existing behavior (Fortran topo_type sign convention).
    negate_z is not involved here.
    """
    t = Topography()
    t.read(tt2_path, topo_type=-2)

    np.testing.assert_array_equal(t.Z, -_analytic_Z())


def test_preprocessing_negate_z_independent_of_topotype(tt2_path):
    """negate_z=True negates Z independently of topo_type sign.

    When topo_type=2 (positive), negate_z flips sign.
    """
    t = Topography()
    t.negate_z = True
    t.read(tt2_path, topo_type=2)

    np.testing.assert_array_equal(t.Z, -_analytic_Z())


def test_preprocessing_double_negate_is_identity(tt2_path):
    """topo_type=-2 AND negate_z=True applies two sign flips — net identity.

    negate_z=True does NOT undo the topo_type sign convention.
    Both negations are applied independently. Use topo_type < 0 OR negate_z,
    not both, unless the double-negate identity is intentional.
    """
    t = Topography()
    t.negate_z = True
    t.read(tt2_path, topo_type=-2)

    np.testing.assert_array_equal(t.Z, _analytic_Z())


# ===========================================================================
# Group 5 — coarsen for NetCDF type 4 and the deprecated `stride` alias
# ===========================================================================

@pytest.mark.netcdf
def test_preprocessing_stride_deprecated_maps_to_coarsen(nc_topo_path, tmp_path):
    """`stride` is deprecated: it warns and maps onto the scalar `coarsen`, so
    the result is identical to reading with `coarsen=2`."""
    pytest.importorskip("xarray")
    path, lon_name, lat_name = nc_topo_path

    t = Topography()
    with pytest.warns(DeprecationWarning):
        t.read(path, topo_type=4, stride=[2, 2])

    assert t.Z.shape == (5, 5)
    # Values must match stride-2 subsampling of the analytic Z
    Z_orig = _analytic_Z()
    for i in range(5):
        for j in range(5):
            np.testing.assert_allclose(t.Z[i, j], Z_orig[2 * i, 2 * j], rtol=1e-6)


@pytest.mark.netcdf
def test_preprocessing_coarsen_only_netcdf(nc_topo_path, tmp_path):
    pytest.importorskip("xarray")
    path, lon_name, lat_name = nc_topo_path

    t = Topography()
    t.coarsen = 2
    t.read(path, topo_type=4)

    assert t.Z.shape == (5, 5)
    Z_orig = _analytic_Z()
    for i in range(5):
        for j in range(5):
            np.testing.assert_allclose(t.Z[i, j], Z_orig[2 * i, 2 * j], rtol=1e-6)


@pytest.mark.netcdf
def test_preprocessing_coarsen_param_matches_attribute(nc_topo_path):
    """Passing coarsen= to read() is equivalent to setting the attribute."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    t_attr = Topography()
    t_attr.coarsen = 2
    t_attr.read(path, topo_type=4)

    t_arg = Topography()
    t_arg.read(path, topo_type=4, coarsen=2)

    np.testing.assert_array_equal(t_arg.Z, t_attr.Z)
    assert t_arg.coarsen == 2


@pytest.mark.netcdf
def test_preprocessing_stride_conflicts_raise(nc_topo_path):
    """Per-axis stride is unsupported; conflicting stride+coarsen is an error."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    # Unequal per-axis stride cannot be expressed by scalar coarsen.
    t = Topography()
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ValueError):
            t.read(path, topo_type=4, stride=[2, 3])

    # stride and a conflicting coarsen at once.
    t2 = Topography()
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ValueError):
            t2.read(path, topo_type=4, coarsen=3, stride=[2, 2])


# ===========================================================================
# Group 6 — read_header() for topo_type=4
# ===========================================================================

@pytest.mark.netcdf
def test_read_header_netcdf_no_z_loaded(nc_topo_path, tmp_path):
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    t = Topography()
    t.path = str(path)
    t.topo_type = 4
    t.read_header()

    assert t._Z is None
    assert t.extent is not None
    assert t.delta is not None


@pytest.mark.netcdf
def test_read_header_netcdf_deferred_z_load(nc_topo_path, tmp_path):
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    t = Topography()
    t.path = str(path)
    t.topo_type = 4
    t.read_header()

    # Accessing Z must trigger deferred read
    Z = t.Z
    assert Z is not None
    np.testing.assert_allclose(Z, _analytic_Z(), rtol=1e-6)


@pytest.mark.netcdf
def test_read_header_netcdf_sn_normalization(tmp_path):
    """read_header() normalizes lat to S→N regardless of file storage order."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    path_sn = _make_nc_topo(tmp_path / "sn.nc", "lon", "lat",
                             lat_south_to_north=True)
    path_ns = _make_nc_topo(tmp_path / "ns.nc", "lon", "lat",
                             lat_south_to_north=False)

    t_sn = Topography()
    t_sn.path = str(path_sn)
    t_sn.topo_type = 4
    t_sn.read_header()

    t_ns = Topography()
    t_ns.path = str(path_ns)
    t_ns.topo_type = 4
    t_ns.read_header()

    np.testing.assert_array_equal(t_sn._y, t_ns._y)
    np.testing.assert_array_equal(t_sn.extent, t_ns.extent)
    assert t_sn.delta == t_ns.delta


@pytest.mark.netcdf
def test_read_header_netcdf_no_double_flip(tmp_path):
    """read_header() + deferred t.Z == direct read() with no read_header() first."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    path = _make_nc_topo(tmp_path / "flip_check.nc", "lon", "lat",
                         lat_south_to_north=True)

    # Path A: read_header() then deferred Z
    t_deferred = Topography()
    t_deferred.path = str(path)
    t_deferred.topo_type = 4
    t_deferred.read_header()
    Z_deferred = t_deferred.Z

    # Path B: direct read()
    t_direct = Topography()
    t_direct.read(str(path), topo_type=4)
    Z_direct = t_direct.Z

    np.testing.assert_allclose(Z_deferred, Z_direct, rtol=1e-12)


# ===========================================================================
# Group 7 — topo_type=1 deprecation
# ===========================================================================

def _write_tt1(path: Path) -> Path:
    """Write a minimal type-1 (x, y, z) ASCII file for deprecation tests."""
    Z = _analytic_Z()
    with open(path, "w") as f:
        for i in range(_NY - 1, -1, -1):  # top row first (y decreasing)
            for j in range(_NX):
                x = _ORIGIN_X + j * _DELTA
                y = _ORIGIN_Y + i * _DELTA
                f.write(f"{x:.1f} {y:.1f} {Z[i,j]:.1f}\n")
    return path


def test_deprecation_type1_read_warns(tmp_path):
    path = _write_tt1(tmp_path / "analytic.tt1")

    t = Topography()
    with pytest.warns(DeprecationWarning, match="deprecated"):
        t.read(str(path), topo_type=1)

    # Type-1 read is deprecated but still functional
    assert t.Z is not None


@pytest.mark.parametrize("attr,value", _PREPROCESSING_NONDEFAULTS)
def test_deprecation_type1_preprocessing_raises(tmp_path, attr, value):
    path = _write_tt1(tmp_path / "analytic.tt1")

    t = Topography()
    setattr(t, attr, value)
    with pytest.warns(DeprecationWarning):
        with pytest.raises(NotImplementedError, match="Convert"):
            t.read(str(path), topo_type=1)


def test_deprecation_type1_write_warns(tmp_path):
    # Read a type-1 file (suppress the read deprecation warning)
    path_in = _write_tt1(tmp_path / "in.tt1")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        t = Topography()
        t.read(str(path_in), topo_type=1)

    path_out = tmp_path / "out.tt1"
    with pytest.warns(DeprecationWarning):
        t.write(str(path_out), topo_type=1)

    assert path_out.exists()


def test_deprecation_type1_read_header_raises():
    """read_header() for topo_type=1 raises IOError (pre-existing behavior)."""
    t = Topography()
    t.path = "dummy.tt1"
    t.topo_type = 1
    with pytest.raises(IOError):
        t.read_header()


# ===========================================================================
# Group 8 — TopoData._normalize_topofiles
# ===========================================================================

def _td_with(entry) -> TopographyData:
    td = TopographyData()
    td.topofiles = [entry]
    return td


def test_normalize_topography_passthrough(tt2_path):
    t = Topography()
    t.path = str(tt2_path)
    t.topo_type = 2

    td = _td_with(t)
    result = td._normalize_topofiles()

    assert len(result) == 1
    assert result[0] is t  # same object, not a copy


def test_normalize_list_deprecated():
    with pytest.warns(DeprecationWarning):
        td = _td_with([2, "/path/to/file.tt2"])
        result = td._normalize_topofiles()

    assert len(result) == 1
    assert result[0].topo_type == 2
    assert result[0].path == "/path/to/file.tt2"


def test_normalize_tuple_deprecated():
    with pytest.warns(DeprecationWarning):
        td = _td_with((2, "/path/to/file.tt2"))
        result = td._normalize_topofiles()

    assert len(result) == 1
    assert result[0].topo_type == 2
    assert result[0].path == "/path/to/file.tt2"


def test_normalize_dict_deprecated_attr_mapping():
    """dict 'extent' key must map to Topography.crop_extent, not .extent."""
    entry = {
        "topo_type": 2,
        "topo_path": "/path/to/file.tt2",
        "extent": [0.0, 5.0, 0.0, 5.0],
        "coarsen": 3,
        "z_shift": 2.5,
    }
    with pytest.warns(DeprecationWarning):
        td = _td_with(entry)
        result = td._normalize_topofiles()

    t = result[0]
    assert t.topo_type == 2
    assert t.path == "/path/to/file.tt2"
    assert t.crop_extent == [0.0, 5.0, 0.0, 5.0]  # dict 'extent' → crop_extent
    assert t.coarsen == 3
    assert t.z_shift == 2.5


def test_normalize_invalid_type_raises():
    td = _td_with(42)
    with pytest.raises(ValueError):
        td._normalize_topofiles()


# ===========================================================================
# Group 9 — _compute_priority_order
# ===========================================================================

def _make_topo_with_delta(dx: float, tmp_path: Path, suffix: str = "") -> Topography:
    """Write a minimal type-2 file and return a Topography with path+topo_type set."""
    n = 10
    t = Topography()
    t._x = np.arange(n) * dx
    t._y = np.arange(n) * dx
    t._Z = np.zeros((n, n))
    t.topo_type = 2
    path = tmp_path / f"topo_dx{dx}{suffix}.tt2"
    t.path = str(path)
    _write_tt2(path, t._Z, cellsize=dx)
    return t


def test_priority_order_finest_last(tmp_path):
    """Finest resolution (smallest dx*dy) must end up last (highest priority).

    Convention: cell area descending = coarsest first, finest last = written
    last = highest Fortran priority (the last file maps to rank 1).
    """
    coarse = _make_topo_with_delta(1.0, tmp_path, "c")
    medium = _make_topo_with_delta(0.5, tmp_path, "m")
    fine = _make_topo_with_delta(0.25, tmp_path, "f")

    td = TopographyData()
    result = td._compute_priority_order([coarse, medium, fine])

    assert result[0].delta == pytest.approx((1.0, 1.0))
    assert result[1].delta == pytest.approx((0.5, 0.5))
    assert result[2].delta == pytest.approx((0.25, 0.25))


def test_priority_order_stable_sort_equal_resolution(tmp_path):
    """Equal-area files preserve their relative input order (stable sort)."""
    A = _make_topo_with_delta(0.5, tmp_path, "A")
    B = _make_topo_with_delta(0.5, tmp_path, "B")

    td = TopographyData()
    result = td._compute_priority_order([A, B])

    assert result[0] is A
    assert result[1] is B


def test_priority_order_override_preserves_user_order(tmp_path):
    """override_order=True returns the list unchanged.

    The user is responsible for placing the highest-priority (finest) file
    last.  override_order=True returns the user list unchanged.
    """
    coarse = _make_topo_with_delta(1.0, tmp_path, "c")
    fine = _make_topo_with_delta(0.25, tmp_path, "f")

    td = TopographyData()
    td.override_order = True
    result = td._compute_priority_order([coarse, fine])

    assert result[0] is coarse
    assert result[1] is fine


def test_priority_order_override_vs_sorted_differ(tmp_path):
    """override_order=True with wrong input gives the wrong winner.

    The winner is the LAST file listed (highest Fortran priority).  Correct
    order is finest-last, which sorting produces; override keeps user order.
    """
    fine = _make_topo_with_delta(0.25, tmp_path, "f")
    coarse = _make_topo_with_delta(1.0, tmp_path, "c")

    td = TopographyData()

    # Wrong user order (finest first) with override: finest is not last, so
    # the coarse file ends up as the winner.
    td.override_order = True
    r_override = td._compute_priority_order([fine, coarse])
    assert r_override[-1] is coarse

    # Sorting reorders to coarsest-first, so the finest file is last and wins.
    td.override_order = False
    r_sorted = td._compute_priority_order([fine, coarse])
    assert r_sorted[-1] is fine

    # The winners (last file listed) differ.
    assert r_override[-1] is not r_sorted[-1]


def test_priority_order_without_z_loaded(tmp_path):
    """_compute_priority_order must not load Z; uses read_header() only."""
    t_fine = _make_topo_with_delta(0.25, tmp_path, "f")
    t_coarse = _make_topo_with_delta(1.0, tmp_path, "c")

    # Clear all cached in-memory data to force read_header() call.
    # _Z is also cleared so the assertion below tests that _compute_priority_order
    # genuinely does not load Z (not merely that the fixture-set zeros survive).
    t_fine._x = None
    t_coarse._x = None
    t_fine._Z = None
    t_coarse._Z = None

    td = TopographyData()
    result = td._compute_priority_order([t_coarse, t_fine])

    # Finest must be last (highest priority)
    assert result[-1].path == t_fine.path
    # Z must NOT have been loaded
    assert t_fine._Z is None
    assert t_coarse._Z is None


# ===========================================================================
# Group 10 — TopoData.write() format
# ===========================================================================

def _parse_topo_data(path: Path) -> dict:
    """Parse a topo.data file into {header: list, files: list-of-str-blocks}."""
    lines = [l for l in path.read_text().splitlines()
             if l.strip() and not l.strip().startswith("#")]
    # First 3 lines: topo_missing, test_topography, ntopofiles
    header = lines[:3]
    # Remaining lines: per-file blocks
    file_blocks = lines[3:]
    return {"header": header, "file_blocks": file_blocks, "all_lines": lines,
            "raw": path.read_text()}


def test_write_global_header_three_lines(tmp_path, tt2_path):
    t = Topography()
    t.path = str(tt2_path)
    t.topo_type = 2

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    parsed = _parse_topo_data(out)
    # Exactly 3 header data lines
    assert len(parsed["header"]) == 3
    # override_order must not appear anywhere
    assert "override_order" not in parsed["raw"]
    assert "override" not in parsed["raw"].lower()


def test_write_per_file_block_10_lines_ascii(tmp_path, tt2_path):
    t = Topography()
    t.path = str(tt2_path)
    t.topo_type = 2

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    raw = out.read_text()
    # Collect lines from the first per-file block (after the global header)
    # Using the inline comment markers to find per-file lines
    per_file_lines = [l for l in raw.splitlines()
                      if "topo_path" in l or "topo_type" in l
                      or "crop_extent" in l or "coarsen" in l
                      or "buffer" in l or "align" in l
                      or "x_shift" in l or "y_shift" in l or "z_shift" in l
                      or "negate_z" in l]
    assert len(per_file_lines) == 10

    # Sentinel values
    assert any("0. 0. 0. 0." in l for l in per_file_lines), \
        "crop_extent sentinel should be '0. 0. 0. 0.'"
    assert any("0. 0." in l and "align" in l for l in per_file_lines), \
        "align sentinel should be '0. 0.'"
    assert any("F" in l and "negate_z" in l for l in per_file_lines), \
        "negate_z default should be 'F'"


def test_write_per_file_block_non_default_values(tmp_path, tt2_path):
    t = Topography()
    t.path = str(tt2_path)
    t.topo_type = 2
    t.crop_extent = [1.0, 5.0, 2.0, 6.0]
    t.coarsen = 3
    t.z_shift = 2.5
    t.negate_z = True

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    raw = out.read_text()
    assert "1" in raw and "5" in raw and "2" in raw and "6" in raw  # crop_extent values
    assert "3" in raw   # coarsen
    assert "2.5" in raw  # z_shift
    crop_line = next(l for l in raw.splitlines() if "crop_extent" in l)
    assert "0. 0. 0. 0." not in crop_line  # not the sentinel
    negate_line = next(l for l in raw.splitlines() if "negate_z" in l)
    assert "T" in negate_line


def test_write_netcdf_descriptor_after_preprocessing_lines(tmp_path):
    """NetCDF descriptor block must appear AFTER the 8 preprocessing lines."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    nc_path = _make_nc_topo(tmp_path / "nc.nc", "lon", "lat")

    t = Topography()
    t.path = str(nc_path)
    t.topo_type = 4

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    raw = out.read_text()
    # negate_z must appear before var_name (descriptor key)
    assert raw.index("negate_z") < raw.index("var_name")


def test_write_ntopofiles_matches_list_length(tmp_path, tt2_path):
    topos = []
    for k in range(3):
        t = Topography()
        t.path = str(tt2_path)
        t.topo_type = 2
        topos.append(t)

    td = TopographyData()
    td.topofiles = topos
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    raw = out.read_text()
    # Third header line is ntopofiles
    data_lines = [l.strip() for l in raw.splitlines()
                  if l.strip() and "#" not in l.split()[0]]
    ntopofiles_line = data_lines[2]  # 0=topo_missing, 1=test_topo, 2=ntopofiles
    assert int(ntopofiles_line.split()[0]) == 3


def test_write_priority_order_finest_last_in_file(tmp_path):
    """Finest topo written last in topo.data = highest Fortran priority."""
    coarse = _make_topo_with_delta(1.0, tmp_path, "c")
    fine = _make_topo_with_delta(0.5, tmp_path, "f")

    td = TopographyData()
    td.topofiles = [coarse, fine]
    td.override_order = False
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    raw = out.read_text()
    coarse_pos = raw.find(Path(coarse.path).name)
    fine_pos = raw.find(Path(fine.path).name)
    # Coarse file must appear before fine file in the output
    assert coarse_pos < fine_pos


def test_write_warns_on_mismatched_datums(tmp_path):
    """TopographyData.write warns when topo files carry different datums."""
    a = _make_topo_with_delta(1.0, tmp_path, "a")
    b = _make_topo_with_delta(1.0, tmp_path, "b")
    a.datum = "NAVD88"
    b.datum = "MSL"

    td = TopographyData()
    td.topofiles = [a, b]
    with pytest.warns(UserWarning, match="mismatched vertical datums"):
        td.write(out_file=str(tmp_path / "topo.data"))


def test_write_no_datum_warning_when_consistent(tmp_path, recwarn):
    """No datum warning when datums agree or only one is set."""
    a = _make_topo_with_delta(1.0, tmp_path, "a")
    b = _make_topo_with_delta(1.0, tmp_path, "b")
    a.datum = "NAVD88"
    b.datum = "NAVD88"          # same datum
    # (b could also be left as None -- a single distinct datum is fine.)

    td = TopographyData()
    td.topofiles = [a, b]
    td.write(out_file=str(tmp_path / "topo.data"))
    assert not any("mismatched vertical datums" in str(w.message)
                   for w in recwarn.list)


# ===========================================================================
# Group 11 — Backward compatibility (deprecation shim round-trip)
# ===========================================================================

def test_backward_compat_list_format_round_trip(tmp_path, tt2_path):
    """Legacy [topo_type, path] list normalizes to Topography and writes correctly."""
    td = TopographyData()
    td.topofiles.append([2, str(tt2_path)])

    out = tmp_path / "topo.data"
    with pytest.warns(DeprecationWarning):
        td.write(out_file=str(out))

    raw = out.read_text()
    # Path appears in output
    assert str(tt2_path) in raw or Path(tt2_path).name in raw
    # topo_type=2 appears
    topo_type_line = next(l for l in raw.splitlines() if "topo_type" in l)
    assert "2" in topo_type_line
    # New per-file block format markers present
    assert "crop_extent" in raw
    assert "negate_z" in raw


def test_backward_compat_mixed_formats(tmp_path, tt2_path):
    """Mix of Topography, list, and dict entries all normalize; dict extent → crop_extent."""
    t_obj = Topography()
    t_obj.path = str(tt2_path)
    t_obj.topo_type = 2

    t_list = [3, str(tt2_path)]

    t_dict = {
        "topo_type": 2,
        "topo_path": str(tt2_path),
        "extent": [0.0, 5.0, 0.0, 5.0],
    }

    td = TopographyData()
    td.topofiles = [t_obj, t_list, t_dict]

    out = tmp_path / "topo.data"
    with pytest.warns(DeprecationWarning):
        td.write(out_file=str(out))

    raw = out.read_text()
    # 3 per-file blocks → 3 topo_path lines
    assert raw.count("topo_path") == 3
    # The dict's crop_extent [0,5,0,5] appears (not all-zero sentinel)
    assert "0 5" in raw or "0. 5." in raw or "0.0 5.0" in raw


# ===========================================================================
# Group 8 — crop_extent / stride pushdown to the NetCDF read (type 4)
#
# A large global NetCDF must not be fully materialized when only a subset is
# requested.  read() pushes both `stride` and a `crop_extent`-derived index
# window down to xarray's lazy indexing so the backend reads only that
# hyperslab.  The window is a *superset* of what the post-read crop() selects,
# so the visible result must be byte-identical to reading the whole file and
# cropping in memory.  These tests pin that equivalence (the correctness
# contract) and guard that a bounded window is actually used (the perf/memory
# contract).  End-to-end memory/time on a real 7 GB global DEM was confirmed
# manually (a strided read that previously OOM'd now returns in ~0.1 s).
# ===========================================================================

@pytest.mark.netcdf
def test_crop_pushdown_matches_full_read(nc_topo_path):
    """crop_extent applied during read == full read then crop() in memory."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path
    crop = [2.0, 7.0, 3.0, 8.0]

    ref = Topography()
    ref.read(path, topo_type=4)
    ref = ref.crop(crop_extent=crop)

    t = Topography()
    t.crop_extent = crop
    t.read(path, topo_type=4)

    np.testing.assert_array_equal(t.x, ref.x)
    np.testing.assert_array_equal(t.y, ref.y)
    np.testing.assert_array_equal(t.Z, ref.Z)


@pytest.mark.netcdf
def test_crop_pushdown_with_buffer_matches_full_read(nc_topo_path):
    """buffer is preserved: the pushed-down window is a superset of crop()'s."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path
    crop = [3.0, 7.0, 3.0, 7.0]

    ref = Topography()
    ref.read(path, topo_type=4)
    ref = ref.crop(crop_extent=crop, buffer=1)

    t = Topography()
    t.crop_extent = crop
    t.buffer = 1
    t.read(path, topo_type=4)

    np.testing.assert_array_equal(t.x, ref.x)
    np.testing.assert_array_equal(t.y, ref.y)
    np.testing.assert_array_equal(t.Z, ref.Z)


@pytest.mark.netcdf
def test_crop_pushdown_with_coarsen_align_matches_full_read(nc_topo_path):
    """coarsen + align survive the pushdown unchanged."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path
    crop = [2.0, 7.0, 3.0, 8.0]

    ref = Topography()
    ref.read(path, topo_type=4)
    ref = ref.crop(crop_extent=crop, coarsen=2, align=(0.0, 0.0))

    t = Topography()
    t.crop_extent = crop
    t.coarsen = 2
    t.align = (0.0, 0.0)
    t.read(path, topo_type=4)

    np.testing.assert_array_equal(t.x, ref.x)
    np.testing.assert_array_equal(t.y, ref.y)
    np.testing.assert_array_equal(t.Z, ref.Z)


@pytest.mark.netcdf
@pytest.mark.parametrize("buffer", [0, 1])
def test_crop_pushdown_with_coarsen_buffer_matches_full_read(nc_topo_path, buffer):
    """coarsen + crop_extent + buffer survive the pushdown: the read hyperslab
    equals a full read then crop(coarsen, buffer).  buffer=1 exercises the
    buffer expansion of the read window."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path
    crop = [3.0, 8.0, 3.0, 8.0]

    ref = Topography()
    ref.read(path, topo_type=4)
    ref = ref.crop(crop_extent=crop, coarsen=2, buffer=buffer)

    t = Topography()
    t.read(path, topo_type=4, crop_extent=crop, coarsen=2, buffer=buffer)

    np.testing.assert_array_equal(t.x, ref.x)
    np.testing.assert_array_equal(t.y, ref.y)
    np.testing.assert_array_equal(t.Z, ref.Z)


@pytest.mark.netcdf
@pytest.mark.parametrize("s2n", [True, False], ids=["S→N", "N→S"])
def test_crop_pushdown_respects_storage_order(tmp_path, s2n):
    """crop_extent pushdown matches full-read-then-crop for either lat storage
    order, and always returns coordinates normalized S→N (y increasing).
    (The bundled fixture flips only the coordinate on N→S, not the data, so the
    invariant is pushdown-vs-full on the *same* file, not S→N-vs-N→S.)"""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    crop = [2.0, 7.0, 3.0, 8.0]

    path = _make_nc_topo(tmp_path / "ordered.nc", "lon", "lat",
                         lat_south_to_north=s2n)

    ref = Topography()
    ref.read(str(path), topo_type=4)
    ref = ref.crop(crop_extent=crop)

    t = Topography()
    t.crop_extent = crop
    t.read(str(path), topo_type=4)

    np.testing.assert_array_equal(t.x, ref.x)
    np.testing.assert_array_equal(t.y, ref.y)
    np.testing.assert_array_equal(t.Z, ref.Z)
    # Coordinates are always normalized to S→N regardless of file order.
    assert np.all(np.diff(t.y) > 0)


@pytest.mark.netcdf
def test_crop_no_overlap_keeps_full_grid(nc_topo_path):
    """A crop_extent that misses the file leaves the array uncropped, matching
    crop()'s own fall-back (does not crash or return an empty grid)."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    ref = Topography()
    ref.read(path, topo_type=4)

    t = Topography()
    t.crop_extent = [100.0, 200.0, 100.0, 200.0]
    # The fall-back is kept, but it must not be silent: the Fortran reader
    # treats the same condition as fatal, so a run that ignores this here
    # fails there instead.
    with pytest.warns(UserWarning, match="did not overlap"):
        t.read(path, topo_type=4)

    np.testing.assert_array_equal(t.Z, ref.Z)


@pytest.mark.netcdf
def test_crop_pushdown_reads_bounded_window(nc_topo_path, monkeypatch):
    """Regression tripwire: a crop_extent read must compute (and therefore read)
    a strict sub-window on each axis, not the whole axis.  If someone reverts to
    materializing the full variable and cropping afterward, this fails."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path

    calls = []
    orig = topotools._crop_indices

    def _spy(x, y, crop_extent, coarsen, buffer, align):
        result = orig(x, y, crop_extent, coarsen, buffer, align)
        calls.append((result, len(x), len(y)))
        return result

    monkeypatch.setattr(topotools, "_crop_indices", _spy)

    t = Topography()
    t.crop_extent = [3.0, 6.0, 3.0, 6.0]
    t.read(path, topo_type=4)

    assert calls, "crop_extent read did not use the window pushdown"
    for result, nx, ny in calls:
        assert result is not None, "crop_extent overlaps the file but got None"
        il, iu, jl, ju = result
        assert 0 <= il < iu <= nx
        assert 0 <= jl < ju <= ny
        assert (iu - il) < nx, "x window spans the whole axis (no pushdown)"
        assert (ju - jl) < ny, "y window spans the whole axis (no pushdown)"


# ===========================================================================
# Group 12 — ASCII/NetCDF read equivalence for coarsen+align (rjl report)
#
# The bug: reading the same DEM as NetCDF (type 4) vs ASCII (type 3) with the
# same coarsening produced grids misaligned by fractions of a coarse cell,
# because `stride` coarsened+aligned NetCDF only and used a different alignment
# convention than crop().  These tests pin the unified behavior: read() takes
# `coarsen`/`align` and both file types produce identical, lattice-aligned grids.
# ===========================================================================

@pytest.mark.netcdf
@pytest.mark.parametrize("align", [None, (0.0, 0.0)])
def test_ascii_netcdf_coarsen_align_identical(nc_topo_path, tmp_path, align):
    """Same data read as NetCDF and as ASCII, with the same coarsen+align, must
    yield identical grids (the core of rjl's report)."""
    pytest.importorskip("xarray")
    path, _, _ = nc_topo_path
    crop = [1.0, 8.0, 1.0, 8.0]
    coarsen = 2

    tn = Topography()
    tn.read(path, topo_type=4, crop_extent=crop, coarsen=coarsen, align=align)

    # Round-trip the same data through ASCII (type 3) and read it back the same way.
    full = Topography()
    full.read(path, topo_type=4)
    asc = tmp_path / "roundtrip.asc"
    full.write(str(asc), topo_type=3, Z_format="%.10f")

    ta = Topography()
    ta.read(str(asc), topo_type=3, crop_extent=crop, coarsen=coarsen, align=align)

    np.testing.assert_allclose(ta.x, tn.x)
    np.testing.assert_allclose(ta.y, tn.y)
    np.testing.assert_allclose(ta.Z, tn.Z)


def test_coarsen_align_lattice_invariant(tmp_path):
    """With a fixed align, shifting crop_extent by whole native cells keeps the
    coarsened grid on the align lattice -- (x0-align)/dx_new stays integer --
    rather than drifting by 1/coarsen of a coarse cell (the reported symptom)."""
    n = 30
    x = np.arange(n, dtype=float)
    y = np.arange(n, dtype=float)
    X, Y = np.meshgrid(x, y)
    Z = X + 10.0 * Y
    t = Topography()
    t.set_xyZ(X, Y, Z)
    asc = tmp_path / "invariant.asc"
    t.write(str(asc), topo_type=3, Z_format="%.6f")

    coarsen = 3
    align = (0.0, 0.0)
    for k in range(coarsen):
        crop = [6.0 + k, 20.0 + k, 6.0 + k, 20.0 + k]
        tk = Topography()
        tk.read(str(asc), topo_type=3, crop_extent=crop,
                coarsen=coarsen, align=align)
        # dx_new = dx*coarsen = 1*3; alignment => (x0-align)/dx_new integer
        xphase = (tk.x[0] - align[0]) / coarsen
        yphase = (tk.y[0] - align[1]) / coarsen
        assert abs(xphase - round(xphase)) < 1e-9, (k, tk.x[0])
        assert abs(yphase - round(yphase)) < 1e-9, (k, tk.y[0])


# ===========================================================================
# Group N — Antimeridian and degenerate crop windows
#
# Nothing exercised antimeridian cropping through Topography before these,
# which is how the failures below shipped silently.  The point of the group is
# that a Topography *never wraps*: a crop is either an ordinary ascending
# window, or it is not expressible on a single Topography at all and has to go
# through TopoInspector.topo_entries().
# ===========================================================================

_WRAPPED_CROP = [170.0, -170.0, -5.0, 5.0]      # crosses the seam, descending
_CONTINUOUS_CROP = [-211.0, -99.0, -5.0, 5.0]   # same region, continuous coords


def _write_global_tt3(path: Path) -> Path:
    """A 1-degree global file spanning the antimeridian, x in [-180, 180]."""
    x = np.linspace(-180.0, 180.0, 361)
    y = np.linspace(-10.0, 10.0, 21)
    Z = -1000.0 + 10.0 * np.cos(np.radians(x))[None, :] * np.ones((y.size, 1))
    t = Topography()
    t.set_xyZ(x, y, Z)
    t.write(str(path), topo_type=3)
    return path


@pytest.fixture
def global_tt3_path(tmp_path):
    return _write_global_tt3(tmp_path / "global.tt3")


def test_wrapped_crop_spelling_raises_not_empty_grid(global_tt3_path):
    """[170, -170] used to produce an *empty* Topography (Z.shape == (11, 0))
    whose .extent then raised an opaque "zero-size array to reduction" from
    numpy, far from the cause."""
    t = Topography()
    t.crop_extent = list(_WRAPPED_CROP)
    with pytest.raises(ValueError, match="must increase in both coordinates"):
        t.read(str(global_tt3_path), topo_type=3)


def test_wrapped_crop_error_names_the_supported_route(global_tt3_path):
    """The error has to say what to do instead, or it just moves the confusion."""
    t = Topography()
    t.crop_extent = list(_WRAPPED_CROP)
    with pytest.raises(ValueError) as excinfo:
        t.read(str(global_tt3_path), topo_type=3)
    assert "topo_entries" in str(excinfo.value)


def test_descending_latitude_crop_also_raises(global_tt3_path):
    """A descending *latitude* pair has no antimeridian excuse at all; it was
    equally silent."""
    t = Topography()
    t.crop_extent = [-100.0, -80.0, 5.0, -5.0]
    with pytest.raises(ValueError, match="must increase in both coordinates"):
        t.read(str(global_tt3_path), topo_type=3)


def test_continuous_crop_is_clipped_and_says_so(global_tt3_path):
    """The continuous spelling is accepted, but it is *not* wrapped -- it is
    reduced to the part of the file that exists (112 degrees requested, 81
    delivered).  That silent reduction is the whole antimeridian confusion, so
    it must warn."""
    t = Topography()
    t.crop_extent = list(_CONTINUOUS_CROP)
    with pytest.warns(UserWarning, match="extends past the data"):
        t.read(str(global_tt3_path), topo_type=3)

    # Clipped to the file's western edge, not wrapped around to +149.
    assert float(t.x[0]) == pytest.approx(-180.0)
    assert float(t.x[-1]) == pytest.approx(-99.0)


def test_ordinary_crop_does_not_warn(global_tt3_path):
    """The clipping warning must not fire for a crop wholly inside the file,
    or it becomes noise everyone filters."""
    t = Topography()
    t.crop_extent = [-100.0, -80.0, -5.0, 5.0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        t.read(str(global_tt3_path), topo_type=3)
    assert float(t.x[0]) == pytest.approx(-100.0)


def test_crop_between_grid_points_raises(global_tt3_path):
    """A window inside the extent but narrower than one cell contains no data.
    It used to return the *full file* -- the opposite of what was asked."""
    t = Topography()
    t.crop_extent = [10.2, 10.8, -5.0, 5.0]
    with pytest.raises(ValueError, match="lies between grid points"):
        t.read(str(global_tt3_path), topo_type=3)


def test_crop_no_overlap_ascii_warns_and_keeps_full_grid(global_tt3_path):
    """ASCII counterpart of test_crop_no_overlap_keeps_full_grid."""
    t = Topography()
    t.crop_extent = [300.0, 320.0, -5.0, 5.0]
    with pytest.warns(UserWarning, match="did not overlap"):
        t.read(str(global_tt3_path), topo_type=3)
    assert t.x.size == 361


def test_unstructured_with_crop_raises_not_typeerror(tmp_path):
    """This used to die with `TypeError: list indices must be integers` from
    indexing a Python list as an array, several frames from the cause.  crop()
    already refused unstructured input; read() now agrees.

    The fixture is a genuine 3-column xyz file (not a headed .tt3) so that the
    read gets far enough to hit that bug when the guard is removed -- a test
    whose "before" failure is an unrelated parse error would pin nothing.
    """
    xyz = tmp_path / "scattered.xyz"
    with open(xyz, "w") as f:
        for x in np.linspace(-110.0, -70.0, 9):
            for y in np.linspace(-8.0, 8.0, 5):
                f.write(f"{x} {y} {-1000.0 + x + y}\n")

    t = Topography()
    t.crop_extent = [-100.0, -80.0, -5.0, 5.0]
    with pytest.raises(NotImplementedError, match="unstructured"):
        t.read(str(xyz), topo_type=1, unstructured=True)


def test_cross_seam_crop_raises_at_topo_data_write(global_tt3_path, tmp_path):
    """Writing a descending crop_extent used to emit `crop_bounds = 170.0
    -170.0`, which Fortran resolves to mx=0, my=0: an empty topo, no error.

    The wrapped spelling stays an error even though the *continuous* spelling
    is now split across the seam automatically -- a descending pair has no
    unambiguous reading.  The message must offer the continuous equivalent
    rather than telling the caller to go build descriptors by hand.
    """
    t = Topography()
    t.path = str(global_tt3_path)
    t.topo_type = 3
    t.crop_extent = list(_WRAPPED_CROP)

    td = TopographyData()
    td.topofiles = [t]
    with pytest.raises(ValueError, match="descending in longitude") as excinfo:
        td.write(out_file=str(tmp_path / "topo.data"))

    # _WRAPPED_CROP is [170, -170]; the continuous equivalent is [-190, -170].
    assert "[-190.0, -170.0]" in str(excinfo.value)


def test_topo_type_none_inferred_from_suffix(global_tt3_path, tmp_path):
    """topo_type=None reached the `:3d` format and raised a TypeError naming
    neither the file nor the attribute.  A .tt3 suffix is unambiguous."""
    t = Topography()
    t.path = str(global_tt3_path)
    t.topo_type = None

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    assert t.topo_type == 3
    assert "  3   # topo_type" in out.read_text()


def test_topo_type_none_unknown_suffix_raises(tmp_path):
    """When the suffix carries no type either, say which attribute to set."""
    src = _write_global_tt3(tmp_path / "global.tt3")
    unknown = tmp_path / "global.dat"
    unknown.write_bytes(src.read_bytes())

    t = Topography()
    t.path = str(unknown)
    t.topo_type = None

    td = TopographyData()
    td.topofiles = [t]
    with pytest.raises(ValueError, match="topo_type is not set"):
        td.write(out_file=str(tmp_path / "topo.data"))


@pytest.mark.netcdf
def test_cropped_netcdf_read_after_read_header_updates_extent(tmp_path):
    """read_header() populates _extent from the file header; the topo_type=4
    read applies its crop while reading the hyperslab, bypassing the property
    setters that would invalidate it.  The object then reported the *full file*
    extent alongside cropped data -- and .extent is what _compute_priority_order
    and the plotting routines use."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    path, _, _ = tmp_path / "nc_extent.nc", None, None
    _make_nc_topo(path, "lon", "lat")

    t = Topography()
    t.path = str(path)
    t.topo_type = 4
    t.read_header()
    assert list(t.extent) == pytest.approx([_ORIGIN_X, _ORIGIN_X + _NX - 1,
                                            _ORIGIN_Y, _ORIGIN_Y + _NY - 1])

    t.crop_extent = [_ORIGIN_X + 2, _ORIGIN_X + 5,
                     _ORIGIN_Y + 2, _ORIGIN_Y + 5]
    t.read()

    assert list(t.extent) == pytest.approx(
        [float(t.x[0]), float(t.x[-1]), float(t.y[0]), float(t.y[-1])])
    assert float(t.extent[0]) == pytest.approx(_ORIGIN_X + 2)


def test_buffer_and_coarsen_give_absolute_output_shape(tt2_path):
    """The existing combined test is a *relative* netCDF-vs-ASCII equality, so
    nothing pinned whether buffer=1, coarsen=2 adds 1 or 2 output points per
    side.  buffer counts coarsened *output* points: the window is expanded by
    buffer*coarsen native points before the strided slice."""
    ref = Topography()
    ref.crop_extent = [_ORIGIN_X + 2, _ORIGIN_X + 7, _ORIGIN_Y + 2,
                       _ORIGIN_Y + 7]
    ref.coarsen = 2
    ref.read(str(tt2_path), topo_type=2)

    t = Topography()
    t.crop_extent = list(ref.crop_extent)
    t.coarsen = 2
    t.buffer = 1
    t.read(str(tt2_path), topo_type=2)

    # One extra coarsened point on each side of each axis.
    assert t.x.size == ref.x.size + 2
    assert t.y.size == ref.y.size + 2
    assert t.Z.shape == (ref.Z.shape[0] + 2, ref.Z.shape[1] + 2)
    # Coarsening is unchanged by the buffer: still every 2nd native point.
    assert float(t.x[1] - t.x[0]) == pytest.approx(2.0 * _DELTA)
    # And the buffered window still contains the unbuffered one.
    assert float(t.x[0]) == pytest.approx(float(ref.x[0]) - 2.0 * _DELTA)


# ===========================================================================
# Group N+1 — Remote sources, and cross-seam crops from ordinary setrun code
#
# Both of these were reported from the field: the most natural possible setrun
# (set .path, .crop_extent, .buffer; append to topofiles) failed, once with a
# mangled local path and once with "crop_bounds exceed file extent", while the
# machinery to handle each already existed and was tested one layer down.
# These pin the two paths being reachable, not just present.
# ===========================================================================

REMOTE_URL = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/"
              "30s_bed_elev_netcdf/ETOPO_2022_v1_30s_N90W180_bed.nc")


def _make_global_nc(path, delta=0.5, lon0=-180.0, lon1=180.0,
                    lat0=-70.0, lat1=10.0):
    """A CF-compliant near-global NetCDF file spanning the antimeridian."""
    netCDF4 = pytest.importorskip("netCDF4")
    x = np.arange(lon0, lon1 + 1e-9, delta)
    y = np.arange(lat0, lat1 + 1e-9, delta)
    Z = -1000.0 + np.outer(np.linspace(0.0, 200.0, y.size), np.ones_like(x))
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


def _entry_blocks(text):
    """Split a written topo.data into its per-file blocks."""
    return [b for b in text.split("# topo_path") if "topo_type" in b]


def _descriptor_values(text, key):
    """Every value written for descriptor *key*, in file order."""
    return re.findall(rf"^{re.escape(key)}\s*=\s*(.+)$", text, re.MULTILINE)


def test_is_remote_url_discriminates_urls_from_paths():
    """The regex must not mistake a Windows drive letter for a URL scheme."""
    from clawpack.geoclaw.netcdf_utils import is_remote_url

    assert is_remote_url("https://example.org/topo.nc")
    assert is_remote_url("http://example.org/topo.nc")
    assert not is_remote_url("/tmp/topo.nc")
    assert not is_remote_url("topo.nc")
    assert not is_remote_url(r"C:\data\topo.nc")
    assert not is_remote_url(Path("/tmp/topo.nc"))


def test_remote_url_in_topofiles_raises_with_the_recipe(tmp_path):
    """A URL used to be run through os.path.abspath, producing

        FileNotFoundError: /run/dir/https:/www.ngdc.noaa.gov/...

    naming a path the user never typed and giving no hint that the fix is to
    fetch it first.
    """
    t = Topography()
    t.path = REMOTE_URL
    t.topo_type = 4
    t.crop_extent = [-160.0, -120.0, -60.0, 0.0]

    td = TopographyData()
    td.topofiles = [t]
    with pytest.raises(ValueError) as excinfo:
        td.write(out_file=str(tmp_path / "topo.data"))

    msg = str(excinfo.value)
    assert "fetch_remote_topo" in msg          # the actionable part
    assert REMOTE_URL in msg                   # unmangled
    assert "https:/www" not in msg             # specifically not collapsed


def test_remote_url_in_dtopofiles_raises(tmp_path):
    """Same trap on the dtopo writer, which shares the abspath pattern."""
    import clawpack.geoclaw.dtopotools as dtopotools

    d = dtopotools.DTopography()
    d.path = REMOTE_URL
    d.dtopo_type = 4

    from clawpack.geoclaw.data import DTopoData

    dtd = DTopoData()
    dtd.dtopofiles = [d]
    with pytest.raises(ValueError, match="URL"):
        dtd.write(out_file=str(tmp_path / "dtopo.data"))


@pytest.mark.netcdf
def test_cross_seam_crop_writes_two_entries(tmp_path):
    """The reported case: a continuous crop spanning the date line.

    Previously raised `crop_bounds lon [...] exceed file extent` even though
    _compute_lon_entries could already cover it.  Must now produce one entry
    per side of the seam with complementary crop_bounds.
    """
    pytest.importorskip("xarray")
    nc = _make_global_nc(tmp_path / "gebco_like.nc")

    t = Topography()
    t.path = str(nc)
    t.topo_type = 4
    t.crop_extent = [-190.0, -120.0, -60.0, 0.0]

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    text = out.read_text()

    assert "2                    =: ntopofiles" in text
    assert len(_entry_blocks(text)) == 2

    offsets = [float(v) for v in _descriptor_values(text, "lon_wrap_offset")]
    assert offsets == [0.0, -360.0]

    bounds = _descriptor_values(text, "crop_bounds")
    # East side comes from the file as-is; west side is the +170..180 strip
    # read with a -360 shift so Fortran places it at -190..-180.
    assert bounds[0].split() == ["-180.0", "-120.0", "-60.0", "0.0"]
    assert bounds[1].split() == ["170.0", "180.0", "-60.0", "0.0"]


@pytest.mark.netcdf
def test_cross_seam_entries_keep_buffer_and_coarsen(tmp_path):
    """buffer and coarsen must reach *every* expanded entry.

    topo_entries() builds Topography objects carrying only _netcdf_meta, so
    routing through it naively writes `buffer = 0` -- which would silently
    undo the Fortran fix that made buffer work for descriptor crops at all.
    """
    pytest.importorskip("xarray")
    nc = _make_global_nc(tmp_path / "gebco_like.nc")

    t = Topography()
    t.path = str(nc)
    t.topo_type = 4
    t.crop_extent = [-190.0, -120.0, -60.0, 0.0]
    t.buffer = 1
    t.coarsen = 20

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))

    blocks = _entry_blocks(out.read_text())
    assert len(blocks) == 2
    for block in blocks:
        assert "1   # buffer" in block
        assert "20   # coarsen" in block


@pytest.mark.netcdf
def test_off_seam_crop_writes_one_entry_with_nonzero_offset(tmp_path):
    """A crop wholly on the far side of the cut needs a *single* entry with a
    non-zero offset.  lon_wrap_offset was hard-coded to 0.0, so this case was
    wrong even though it never needed splitting."""
    pytest.importorskip("xarray")
    nc = _make_global_nc(tmp_path / "g.nc")

    t = Topography()
    t.path = str(nc)
    t.topo_type = 4
    t.crop_extent = [185.0, 195.0, -60.0, 0.0]   # i.e. -175..-165

    td = TopographyData()
    td.topofiles = [t]
    out = tmp_path / "topo.data"
    td.write(out_file=str(out))
    text = out.read_text()

    assert "1                    =: ntopofiles" in text
    assert [float(v) for v in
            _descriptor_values(text, "lon_wrap_offset")] == [360.0]
    assert _descriptor_values(text, "crop_bounds")[0].split() == [
        "-175.0", "-165.0", "-60.0", "0.0"]


@pytest.mark.netcdf
def test_wrapping_write_does_not_scan_the_whole_file(tmp_path, monkeypatch):
    """topo_entries() inspects with crop_bounds unset, so its fill/magnitude
    checks would read the *entire* variable and reject NaN anywhere in it.

    On a global DEM read over OPeNDAP that turns `make data` into a full
    download.  The regression is invisible on a small local fixture, so it is
    pinned directly rather than by timing.
    """
    pytest.importorskip("xarray")
    from clawpack.geoclaw import netcdf_utils as ncutils

    nc = _make_global_nc(tmp_path / "g.nc")

    calls = []
    original = ncutils.TopoInspector._check_fill_in_crop

    def spy(self, *args, **kwargs):
        calls.append(args)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(ncutils.TopoInspector, "_check_fill_in_crop", spy)

    t = Topography()
    t.path = str(nc)
    t.topo_type = 4
    t.crop_extent = [-190.0, -120.0, -60.0, 0.0]

    td = TopographyData()
    td.topofiles = [t]
    td.write(out_file=str(tmp_path / "topo.data"))

    assert calls == [], (
        "write() triggered a fill scan; on a remote global DEM this "
        "downloads the whole file during `make data`.")


@pytest.mark.netcdf
def test_wrapping_crop_still_checks_latitude(tmp_path):
    """Longitude wraps; latitude does not.  Dropping crop_bounds validation to
    allow the wrap must not also drop the latitude check."""
    pytest.importorskip("xarray")
    nc = _make_global_nc(tmp_path / "g.nc")   # lat spans -70..10

    t = Topography()
    t.path = str(nc)
    t.topo_type = 4
    t.crop_extent = [-190.0, -120.0, -60.0, 45.0]   # 45N is off the file

    td = TopographyData()
    td.topofiles = [t]
    with pytest.raises(ValueError, match="latitude extent"):
        td.write(out_file=str(tmp_path / "topo.data"))


# ===========================================================================
# _crop_indices: two-axis behavior preservation
#
# _crop_indices used to inline the whole index computation; the per-axis math
# now lives in coordinate_tools.crop_indices (shared with dtopo).  These pin
# that the move changed nothing observable, including the precedence between
# its two failure modes.
# ===========================================================================

def _ref_crop_indices_two_axis(x, y, crop_extent, coarsen, buffer, align):
    """Pre-refactor ``topotools._crop_indices``, inlined verbatim as an oracle.

    Returns the four-tuple, ``None`` for a non-overlap, or raises for a
    sub-cell window -- exactly as before the per-axis math moved out.  The
    clipping warning is omitted; only the index result is compared.
    """
    dx = np.round(abs(x[1] - x[0]), 15)
    dy = np.round(abs(y[1] - y[0]), 15)
    dx_new = dx * coarsen
    dy_new = dy * coarsen

    try:
        ilower = (x >= crop_extent[0]).nonzero()[0][0]
        iupper = (x <= crop_extent[1]).nonzero()[0][-1]
        jlower = (y >= crop_extent[2]).nonzero()[0][0]
        jupper = (y <= crop_extent[3]).nonzero()[0][-1]
    except IndexError:
        return None
    if iupper < ilower or jupper < jlower:
        raise ValueError("lies between grid points")

    if (coarsen > 1) and (align is not None):
        xs = np.array([x[ilower + i] for i in range(coarsen)])
        offsets = (xs - align[0]) / dx_new
        offsets_frac = offsets - np.round(offsets)
        ioffset = np.argmin(abs(offsets_frac))
        ilower = ilower + ioffset
        iupper = iupper - np.remainder(iupper - ilower, coarsen)

        ys = np.array([y[jlower + j] for j in range(coarsen)])
        offsets = (ys - align[1]) / dy_new
        offsets_frac = offsets - np.round(offsets)
        joffset = np.argmin(abs(offsets_frac))
        jlower = jlower + joffset
        jupper = jupper - np.remainder(jupper - jlower, coarsen)

    ilower = np.maximum(0, ilower - buffer * coarsen)
    jlower = np.maximum(0, jlower - buffer * coarsen)
    iupper = np.minimum(len(x) - 1, iupper + buffer * coarsen) + 1
    jupper = np.minimum(len(y) - 1, jupper + buffer * coarsen) + 1
    return int(ilower), int(iupper), int(jlower), int(jupper)


@pytest.mark.python
def test_crop_indices_two_axis_matches_pre_refactor():
    """Sweep crop x coarsen x buffer x align; windows must be identical.

    Behavior-preservation Vet for delegating the per-axis math: the two grids
    differ in spacing (dx != dy) and origin so an x/y mix-up cannot pass.
    """
    import itertools

    dx, dy = 0.25, 0.5
    x = np.arange(-5.0, 15.0 + dx / 2, dx)
    y = np.arange(2.0, 22.0 + dy / 2, dy)

    x_pairs = [(-5.0, 15.0), (-2.1, 6.0), (0.0, 9.9), (3.3, 12.5)]
    y_pairs = [(2.0, 22.0), (4.4, 13.0), (7.0, 19.9), (10.5, 21.0)]
    coarsens = [1, 2, 3, 5]
    buffers = [0, 1, 4]
    aligns = [None, (0.0, 0.0), (0.5, 0.25), (-1.0, 3.0)]

    checked = 0
    for (x1, x2), (y1, y2), coarsen, buffer, align in itertools.product(
            x_pairs, y_pairs, coarsens, buffers, aligns):
        crop_extent = [x1, x2, y1, y2]
        ref = _ref_crop_indices_two_axis(x, y, crop_extent, coarsen, buffer,
                                         align)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")   # clipping warning is not the subject
            got = topotools._crop_indices(x, y, crop_extent, coarsen, buffer,
                                          align)
        assert got == ref, (crop_extent, coarsen, buffer, align, got, ref)
        checked += 1

    assert checked == 4 * 4 * 4 * 3 * 4   # 768 combinations actually compared


@pytest.mark.python
def test_crop_indices_non_overlap_outranks_subcell_on_the_other_axis():
    """A non-overlapping axis still returns None, even when the other is empty.

    The original wrapped all four index lookups in one try/except, so an
    IndexError from *either* axis escaped before the sub-cell check ran.
    Resolving the axes independently would have turned this None into a
    ValueError -- a silent change in which failure the caller sees.
    """
    dx, dy = 0.5, 0.5
    x = np.arange(0.0, 10.0, dx)
    y = np.arange(0.0, 10.0, dy)

    # x: sub-cell (strictly between 1.0 and 1.5).  y: no overlap at all.
    crop_extent = [1.1, 1.4, 50.0, 60.0]
    assert _ref_crop_indices_two_axis(x, y, crop_extent, 1, 0, None) is None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert topotools._crop_indices(x, y, crop_extent, 1, 0, None) is None
