"""
Tests for DTopoInspector.

Covers: deformation-variable discovery (by name, by unique 3-D fallback,
ambiguity error), time-axis handling (t0/dt extraction, uniform-spacing
validation, datetime conversion via time_reference, missing time axis),
the time-slowest-dimension requirement, the fill-value warning, and the
descriptor writer output.
"""
import io
import warnings

import numpy as np
import pytest
import xarray as xr

from clawpack.geoclaw.netcdf_utils import (
    DTopoInspector,
    DTopoMetadata,
    DescriptorWriter,
)

pytestmark = [pytest.mark.python, pytest.mark.netcdf]


def make_dtopo_dataset(times=None, var_name="dz", dims=("time", "lat", "lon"),
                       nlat=6, nlon=5, units="meters"):
    """Build a small in-memory dtopo dataset with deterministic values.

    The deformation variable declares ``units`` (meters by default); the
    dtopo inspector requires a units attribute on the NetCDF read path.
    """
    if times is None:
        times = np.array([0.0, 0.5, 1.0])
    coords = {
        "time": np.asarray(times),
        "lat": np.linspace(0.0, 1.0, nlat),
        "lon": np.linspace(0.0, 2.0, nlon),
    }
    sizes = {"time": len(coords["time"]), "lat": nlat, "lon": nlon}
    shape = [sizes[d] for d in dims]
    data = np.arange(np.prod(shape), dtype=float).reshape(shape)
    attrs = {"units": units} if units else {}
    return xr.Dataset({var_name: (list(dims), data, attrs)},
                      coords={d: coords[d] for d in dims})


def write_dataset(ds, path, var_name="dz"):
    ds.to_netcdf(path, encoding={var_name: {"_FillValue": None}})
    return path


def test_inspect_dtopo_basics(tmp_path):
    path = write_dataset(make_dtopo_dataset(), tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()

    assert isinstance(meta, DTopoMetadata)
    assert meta.var_name == "dz"
    assert meta.x_name == "lon"
    assert meta.y_name == "lat"
    assert meta.time_name == "time"
    assert meta.dim_order == ["time", "y", "x"]
    assert np.isclose(meta.t0, 0.0)
    assert np.isclose(meta.dt, 0.5)
    assert meta.mt == 3


def test_var_discovery_unique_3d_fallback(tmp_path):
    """A 3-D variable with a non-standard name is found when unique."""
    ds = make_dtopo_dataset(var_name="seafloor_motion")
    path = write_dataset(ds, tmp_path / "dt.nc", var_name="seafloor_motion")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()
    assert meta.var_name == "seafloor_motion"


def test_var_discovery_ambiguous_raises(tmp_path):
    """Two unrecognized 3-D variables cannot be disambiguated."""
    ds = make_dtopo_dataset(var_name="motion_a")
    ds["motion_b"] = ds["motion_a"].copy()
    path = write_dataset(ds, tmp_path / "dt.nc", var_name="motion_a")
    with DTopoInspector(path) as insp:
        with pytest.raises(ValueError, match="var_name"):
            insp.inspect_dtopo()


def test_nonuniform_time_raises(tmp_path):
    ds = make_dtopo_dataset(times=[0.0, 0.5, 2.0])
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        with pytest.raises(ValueError, match="uniform"):
            insp.inspect_dtopo()


def test_time_not_slowest_raises(tmp_path):
    ds = make_dtopo_dataset(dims=("lat", "time", "lon"))
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        with pytest.raises(ValueError, match="slowest"):
            insp.inspect_dtopo()


def test_missing_time_axis_raises(tmp_path):
    ds = make_dtopo_dataset()
    ds = ds.isel(time=0).drop_vars("time")   # 2-D variable, no time coord
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path, var_name="dz") as insp:
        with pytest.raises(ValueError, match="time"):
            insp.inspect_dtopo()


def test_datetime_time_requires_reference(tmp_path):
    times = np.array(["2011-03-11T05:46:00", "2011-03-11T05:46:30",
                      "2011-03-11T05:47:00"], dtype="datetime64[ns]")
    ds = make_dtopo_dataset(times=times)
    path = write_dataset(ds, tmp_path / "dt.nc")

    with DTopoInspector(path) as insp:
        with pytest.raises(ValueError, match="time_reference"):
            insp.inspect_dtopo()

    with DTopoInspector(path,
                        time_reference="2011-03-11T05:46:00") as insp:
        meta = insp.inspect_dtopo()
    assert np.isclose(meta.t0, 0.0)
    assert np.isclose(meta.dt, 30.0)
    assert meta.mt == 3


def test_timedelta_time_axis_converts_correctly(tmp_path):
    """Regression test: a time coordinate xarray decodes to timedelta64 must
    convert to seconds via true unit conversion, not have its raw tick count
    reinterpreted as a plain float.

    Some xarray versions decode a bare duration-style ``units`` attribute
    (e.g. ``"seconds"``, as ``DTopography._write_netcdf`` writes) into
    ``timedelta64[ns]`` even though the on-disk data is a plain float.  Prior
    to this fix, ``_compute_time_axis`` fell through to
    ``np.asarray(tvals, dtype=float)`` for anything that wasn't
    ``datetime64``, which reinterprets the underlying nanosecond tick count
    as if it were already in seconds -- e.g. t=1 second (1e9 ns) read back
    as t=1e9 "seconds".  That is exactly the billion-fold inflation reported
    against the chile2010 example, where ``dt_max_dtopo`` never triggered
    because ``tfdtopo`` looked a billion seconds in the future.
    """
    ds = write_dataset(make_dtopo_dataset(times=[0.0, 0.5, 1.0]),
                       tmp_path / "dt.nc")
    with DTopoInspector(ds) as insp:
        # Simulate the timedelta64-decoded scenario directly, since whether
        # this happens at xr.open_dataset() time is xarray-version-dependent.
        insp.ds = insp.ds.assign_coords(
            time=(np.array([0.0, 0.5, 1.0]) * 1e9).astype("timedelta64[ns]"))
        t0, dt, mt = insp._compute_time_axis("time")

    assert np.isclose(t0, 0.0)
    assert np.isclose(dt, 0.5)
    assert mt == 3


def test_single_time_gives_zero_dt(tmp_path):
    ds = make_dtopo_dataset(times=[7.0])
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()
    assert np.isclose(meta.t0, 7.0)
    assert meta.dt == 0.0
    assert meta.mt == 1


def test_fill_value_warns(tmp_path):
    ds = make_dtopo_dataset()
    path = tmp_path / "dt.nc"
    ds.to_netcdf(path)   # default encoding declares a NaN _FillValue
    with DTopoInspector(path) as insp:
        with pytest.warns(UserWarning, match="fill"):
            insp.inspect_dtopo()


@pytest.mark.parametrize("cf_unit,factor", [
    ("seconds", 1.0), ("minutes", 60.0), ("hours", 3600.0), ("days", 86400.0),
])
def test_numeric_time_units_scaled_to_seconds(tmp_path, cf_unit, factor):
    """A bare numeric time axis is scaled to seconds using its units attribute,
    not assumed to already be seconds (regression for the finding that a
    'hours'/'minutes' axis was read too fast)."""
    ds = make_dtopo_dataset(times=[0.0, 1.0, 2.0])
    ds["time"].attrs["units"] = cf_unit
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()
    assert np.isclose(meta.t0, 0.0)
    assert np.isclose(meta.dt, 1.0 * factor)


def test_numeric_time_unrecognised_units_raise(tmp_path):
    """An unrecognised time-units string is rejected, never assumed seconds."""
    ds = make_dtopo_dataset(times=[0.0, 1.0, 2.0])
    ds["time"].attrs["units"] = "fortnights"
    path = write_dataset(ds, tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        with pytest.raises(ValueError, match="[Uu]nrecognised time units"):
            insp.inspect_dtopo()


def test_open_nonstandard_extension(tmp_path):
    """A valid NetCDF file with a non-NetCDF extension (e.g. '.dtt3') still
    opens: the inspector selects an explicit engine rather than relying on
    xarray's extension-based backend guessing."""
    ds = make_dtopo_dataset()
    path = write_dataset(ds, tmp_path / "quake.dtt3")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()
    assert meta.mt == 3


def test_descriptor_writer_output(tmp_path):
    path = write_dataset(make_dtopo_dataset(), tmp_path / "dt.nc")
    with DTopoInspector(path) as insp:
        meta = insp.inspect_dtopo()

    buf = io.StringIO()
    DescriptorWriter.write_dtopo_descriptor(buf, meta)
    text = buf.getvalue()

    for key in ("var_name", "x_name", "y_name", "time_name",
                "lon_wrap_offset", "y_increasing", "dim_order", "t0", "dt"):
        assert key in text, f"missing descriptor key {key}"
    assert "dim_order      = time,y,x" in text
    assert "dt             = 0.5" in text
    # Blank line terminates the block for the Fortran parser.
    assert text.endswith("\n\n")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
