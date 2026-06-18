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

# The 7 preprocessing attributes and a non-default value for each.
_PREPROCESSING_NONDEFAULTS = [
    ("crop_extent", [1.0, 8.0, 1.0, 8.0]),
    ("coarsen", 2),
    ("buffer", 1),
    ("align", (0.0, 0.0)),
    ("x_shift", 5.0),
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

    This is the pre-existing behaviour (Fortran topo_type sign convention).
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
# Group 5 — stride × coarsen interaction (NetCDF type 4)
# ===========================================================================

@pytest.mark.netcdf
def test_preprocessing_stride_only(nc_topo_path, tmp_path):
    pytest.importorskip("xarray")
    path, lon_name, lat_name = nc_topo_path

    t = Topography()
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
def test_preprocessing_stride_and_coarsen_compound(nc_topo_path, tmp_path):
    """stride=[2,2] then coarsen=2: effective subsampling is stride-4.

    The shape is floor(10/4) → check against actual output shape.
    """
    pytest.importorskip("xarray")
    path, lon_name, lat_name = nc_topo_path

    t = Topography()
    t.coarsen = 2
    t.read(path, topo_type=4, stride=[2, 2])

    # After stride=[2,2], array is 5×5. After coarsen=2, it's 3×3 (floor(5/2)+1).
    # The exact shape depends on how crop() slices; just check stride-4 values.
    Z_orig = _analytic_Z()
    for i in range(t.Z.shape[0]):
        for j in range(t.Z.shape[1]):
            np.testing.assert_allclose(t.Z[i, j], Z_orig[4 * i, 4 * j], rtol=1e-6)


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
    """read_header() normalises lat to S→N regardless of file storage order."""
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
    """read_header() for topo_type=1 raises IOError (pre-existing behaviour)."""
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


def test_write_per_file_block_9_lines_ascii(tmp_path, tt2_path):
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
                      or "x_shift" in l or "z_shift" in l
                      or "negate_z" in l]
    assert len(per_file_lines) == 9

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
    """NetCDF descriptor block must appear AFTER the 9 preprocessing lines."""
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
    """Legacy [topo_type, path] list normalises to Topography and writes correctly."""
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
    # New 9-line format markers present
    assert "crop_extent" in raw
    assert "negate_z" in raw


def test_backward_compat_mixed_formats(tmp_path, tt2_path):
    """Mix of Topography, list, and dict entries all normalise; dict extent → crop_extent."""
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
