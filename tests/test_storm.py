#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data."""

from pathlib import Path
import sys

import numpy as np
import pytest

import clawpack.geoclaw.surge.storm as storm
import clawpack.geoclaw.util as util

# Local test directory and bundled storm test data
testdir = Path(__file__).parent
data_dir = testdir / "data" / "storm"

# Current file-format tests
FILE_FORMAT_TESTS = ["atcf", 
                    "hurdat", 
                    "jma", 
                    "tcvitals", 
                    pytest.param("ibtracs", 
                                 marks=pytest.mark.xfail(
                                 reason=("IBTrACS GeoClaw output does not ,"
                                         "currently match baseline")))
]

DATA_FILE_FORMAT_MAP = {
    1: ["ascii", 1, "nws12", "owi"],
    2: ["netcdf", 2, "nws13"],
}


def create_netcdf_storm_file(path):
    """Create a small deterministic NetCDF storm dataset for testing."""
    xr = pytest.importorskip("xarray")
    size = (3, 5, 3)
    coords = [
        ("longitude", np.linspace(-1.0, 1.0, size[0])),
        ("latitude", np.linspace(-2.0, 2.0, size[1])),
        ("valid_time", np.linspace(0.0, 2.0, size[2])),
    ]

    values = np.arange(np.prod(size), dtype=float).reshape(size)
    wind_x = xr.DataArray(values, coords=coords)
    wind_y = xr.DataArray(values + 100.0, coords=coords)
    pressure = xr.DataArray(values + 1000.0, coords=coords)
    ds = xr.Dataset({"u": wind_x, "v": wind_y, "pressure": pressure})
    ds.to_netcdf(path)


def create_era5_storm_file(path, pressure_fields=None, u_fields=None,
                           v_fields=None, times=None, lon=None, lat=None):
    """Create an ERA5-style CF-compliant NetCDF met-forcing file.

    ERA5 (ECMWF Reanalysis v5) is one example of a CF-compliant netCDF
    met-forcing source; GeoClaw supports any CF-compliant netCDF data.
    This function produces the variable/dimension schema used by ECMWF ERA5
    reanalysis: dims ``(valid_time, latitude, longitude)``, variables
    ``msl`` (Pa), ``u10`` (m/s), ``v10`` (m/s), lon in [0, 360], lat S-to-N,
    and a ``Conventions = "CF-1.7"`` global attribute.

    The dimension order ``(valid_time, latitude, longitude)`` = ``(nt, nlat,
    nlon)`` matches the C-row-major layout expected by the Fortran NetCDF reader
    ``nf90_get_var`` when it fills its column-major ``pressure(nlon, nlat, nt)``
    array: the library reverses the dimension order automatically.

    Parameters
    ----------
    path : Path or str
        Output NetCDF file path.
    pressure_fields : ndarray, shape (nt, nlat, nlon), optional
        Pressure in Pa.  Defaults to 101300 Pa everywhere.
    u_fields : ndarray, shape (nt, nlat, nlon), optional
        Eastward 10-m wind in m/s.  Defaults to zero.
    v_fields : ndarray, shape (nt, nlat, nlon), optional
        Northward 10-m wind in m/s.  Defaults to zero.
    times : array-like of int, optional
        Time values in **seconds since the time_offset epoch** used in the
        accompanying .storm descriptor.  Defaults to ``[0, 3600, 7200]``.
    lon : array-like of float, optional
        Longitude values in [0, 360].  Defaults to a 5-point sample grid.
    lat : array-like of float, optional
        Latitude values, S-to-N.  Defaults to a 4-point sample grid.
    """
    # netCDF4 is required as xarray's NetCDF backend.
    pytest.importorskip("netCDF4")
    xr = pytest.importorskip("xarray")

    # --- defaults (synthetic values for unit tests) ---
    if lon is None:
        lon = np.linspace(260.0, 300.0, 5)   # [0, 360] convention
    else:
        lon = np.asarray(lon, dtype=float)
    if lat is None:
        lat = np.linspace(10.0, 35.0, 4)     # S-to-N
    else:
        lat = np.asarray(lat, dtype=float)
    if times is None:
        times = np.array([0, 3600, 7200], dtype=np.int64)
    else:
        times = np.asarray(times, dtype=np.int64)

    nt, nlat, nlon = len(times), len(lat), len(lon)
    shape = (nt, nlat, nlon)

    if pressure_fields is None:
        pressure_fields = np.full(shape, 101300.0)
    if u_fields is None:
        u_fields = np.zeros(shape)
    if v_fields is None:
        v_fields = np.zeros(shape)

    # Use a fixed epoch string; callers that need a specific offset should
    # convert their datetime values to integer seconds before the call.
    time_epoch_str = "2012-08-29T00:00:00"

    ds = xr.Dataset(
        {
            "msl": (
                ["valid_time", "latitude", "longitude"],
                pressure_fields,
                {
                    "units": "Pa",
                    "standard_name": "air_pressure_at_mean_sea_level",
                    "long_name": "Mean sea level pressure",
                },
            ),
            "u10": (
                ["valid_time", "latitude", "longitude"],
                u_fields,
                {
                    "units": "m s-1",
                    "standard_name": "eastward_wind",
                    "long_name": "10 metre U wind component",
                },
            ),
            "v10": (
                ["valid_time", "latitude", "longitude"],
                v_fields,
                {
                    "units": "m s-1",
                    "standard_name": "northward_wind",
                    "long_name": "10 metre V wind component",
                },
            ),
        },
        coords={
            "longitude": (
                ["longitude"],
                lon,
                {"units": "degrees_east", "axis": "X",
                 "standard_name": "longitude"},
            ),
            "latitude": (
                ["latitude"],
                lat,
                {"units": "degrees_north", "axis": "Y",
                 "standard_name": "latitude"},
            ),
            "valid_time": (
                ["valid_time"],
                times,
                {
                    "units": f"seconds since {time_epoch_str}",
                    "axis": "T",
                    "standard_name": "time",
                    "long_name": "time",
                },
            ),
        },
        attrs={"Conventions": "CF-1.7"},
    )
    ds.to_netcdf(path)


def create_nws13_storm_file(path, pressure_fields=None, u_fields=None,
                            v_fields=None, times=None, lon=None, lat=None):
    """Create a NWS13/OWI-NetCDF style met-forcing file.

    NWS13 (NOAA Weather Service format 13) is the OWI NetCDF met-forcing
    schema used by ADCIRC and GeoClaw.  This function produces the canonical
    NWS13 variable and dimension names: dims ``(time, lat, lon)``, variables
    ``uwnd`` (m/s), ``vwnd`` (m/s), ``press`` (mb).

    Pressure is stored in millibar (mb) to exercise the unit-awareness path
    in tests.  Note that the Fortran NetCDF reader does **not** perform unit
    conversion; callers responsible for running Fortran should convert to Pa
    before calling this function.

    The dimension order ``(time, lat, lon)`` = ``(nt, nlat, nlon)`` is
    correct for the Fortran column-major allocation ``pressure(nlon, nlat,
    nt)`` via the library-managed Fortran–C transpose.

    Parameters
    ----------
    path : Path or str
        Output NetCDF file path.
    pressure_fields : ndarray, shape (nt, nlat, nlon), optional
        Pressure in **mb**.  Defaults to 1013.0 mb everywhere.
    u_fields : ndarray, shape (nt, nlat, nlon), optional
        Eastward wind in m/s.  Defaults to zero.
    v_fields : ndarray, shape (nt, nlat, nlon), optional
        Northward wind in m/s.  Defaults to zero.
    times : array-like of int, optional
        Time in seconds (absolute UNIX seconds or relative to an epoch).
        Defaults to ``[0, 3600, 7200]``.
    lon : array-like of float, optional
        Longitude values.  Defaults to a 5-point sample.
    lat : array-like of float, optional
        Latitude values.  Defaults to a 4-point sample.
    """
    pytest.importorskip("netCDF4")
    xr = pytest.importorskip("xarray")

    if lon is None:
        lon = np.linspace(-99.0, -70.0, 5)
    else:
        lon = np.asarray(lon, dtype=float)
    if lat is None:
        lat = np.linspace(8.0, 32.0, 4)
    else:
        lat = np.asarray(lat, dtype=float)
    if times is None:
        times = np.array([0, 3600, 7200], dtype=np.int64)
    else:
        times = np.asarray(times, dtype=np.int64)

    nt, nlat, nlon = len(times), len(lat), len(lon)
    shape = (nt, nlat, nlon)

    if pressure_fields is None:
        pressure_fields = np.full(shape, 1013.0)  # mb
    if u_fields is None:
        u_fields = np.zeros(shape)
    if v_fields is None:
        v_fields = np.zeros(shape)

    ds = xr.Dataset(
        {
            "press": (
                ["time", "lat", "lon"],
                pressure_fields,
                {"units": "mb", "long_name": "Sea-level pressure"},
            ),
            "uwnd": (
                ["time", "lat", "lon"],
                u_fields,
                {"units": "m/s", "long_name": "U-component of wind"},
            ),
            "vwnd": (
                ["time", "lat", "lon"],
                v_fields,
                {"units": "m/s", "long_name": "V-component of wind"},
            ),
        },
        coords={
            "lon": (["lon"], lon, {"units": "degrees_east"}),
            "lat": (["lat"], lat, {"units": "degrees_north"}),
            "time": (["time"], times, {"units": "seconds since 1970-01-01T00:00:00"}),
        },
    )
    ds.to_netcdf(path)


def check_geoclaw(paths, check_header=False):
    """
    Check that two GeoClaw-formatted storm files are numerically equivalent.
    """
    paths = [Path(path) for path in paths]

    if check_header:
        with paths[0].open("r") as first_file, paths[1].open("r") as second_file:
            # Check for number of lines
            assert int(first_file.readline()) == int(second_file.readline())

            # Check for time offset
            assert first_file.readline() == second_file.readline()

    data = [np.loadtxt(path, skiprows=3) for path in paths]
    np.testing.assert_almost_equal(data[0], data[1])


def _storm_input_path(file_format):
    """Return the bundled input path for a given storm file format."""
    suffix = "nc" if file_format == "ibtracs" else "txt"
    return data_dir / f"{file_format}.{suffix}"


def _storm_check_path(file_format):
    """Return the bundled GeoClaw-format reference file for a storm format."""
    return data_dir / f"{file_format}_geoclaw.txt"


def _make_storm_from_format(file_format):
    """Read one bundled storm file and return the storm plus fill callbacks."""
    input_path = _storm_input_path(file_format)

    if file_format == "ibtracs":
        kwargs = {"sid": "2008245N17323", "agency_pref": ["wmo", "usa"]}

        atcf_path = data_dir / "atcf.txt"
        storm_atcf = storm.Storm(atcf_path, file_format="ATCF")

        def fill_mwr(t, this_storm):
            return storm.fill_rad_w_other_source(
                t, this_storm, storm_atcf, "max_wind_radius"
            )

        def fill_rad(t, this_storm):
            return storm.fill_rad_w_other_source(
                t, this_storm, storm_atcf, "storm_radius"
            )
    else:
        kwargs = {}
        fill_mwr = None
        fill_rad = None

    test_storm = storm.Storm(input_path, file_format=file_format, **kwargs)

    # Temporary normalization for formats that do not provide these radii.
    if file_format in ["hurdat", "jma"]:
        test_storm.max_wind_radius[:] = 0.0
        test_storm.storm_radius[:] = 0.0

    return test_storm, fill_mwr, fill_rad


@pytest.mark.python
@pytest.mark.parametrize("file_format", FILE_FORMAT_TESTS)
def test_storm_io(tmp_path, file_format):
    r"""Test reading one storm format and writing it back in GeoClaw format."""
    if file_format == "ibtracs":
        pytest.importorskip("xarray")
    if file_format == "atcf":
        pytest.importorskip("pandas")

    test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
    out_path = tmp_path / f"{file_format}_geoclaw.txt"
    check_path = _storm_check_path(file_format)
    write_kwargs = {"file_format": "geoclaw"}
    if fill_mwr is not None:
        write_kwargs["max_wind_radius_fill"] = fill_mwr
    if fill_rad is not None:
        write_kwargs["storm_radius_fill"] = fill_rad

    test_storm.write(out_path, **write_kwargs)

    check_geoclaw([out_path, check_path])


def save_storm_test_data(output_dir):
    """Utility helper to regenerate bundled GeoClaw-format storm baselines."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for file_format in FILE_FORMAT_TESTS:
        if file_format == "ibtracs":
            try:
                import xarray  # noqa: F401
            except ImportError:
                print("Skipping IBTrACS save, missing xarray.")
                continue
        if file_format == "atcf":
            try:
                import pandas  # noqa: F401
            except ImportError:
                print("Skipping ATCF save, missing pandas.")
                continue

        test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
        check_path = output_dir / f"{file_format}_geoclaw.txt"
        write_kwargs = {"file_format": "geoclaw"}
        if fill_mwr is not None:
            write_kwargs["max_wind_radius_fill"] = fill_mwr
        if fill_rad is not None:
            write_kwargs["storm_radius_fill"] = fill_rad
        test_storm.write(check_path, **write_kwargs)


@pytest.mark.python
@pytest.mark.parametrize("file_format",
                         ["ascii", "netcdf", "netcdf_era5", "netcdf_nws13"])
def test_data_storm_roundtrip(file_format, tmp_path):
    """Test round-trip read/write of data-driven storm metadata.

    Exercises:
    - ``ascii``      : OWI/NWS12 ASCII pair (two PRE + WIN files)
    - ``netcdf``     : generic NetCDF format with existing helper
    - ``netcdf_era5``: ERA5-style CF NetCDF (valid_time/latitude/longitude,
                       u10/v10/msl) — verifies that write_data discovers and
                       records ERA5 variable and dimension names correctly
    - ``netcdf_nws13``: NWS13-style OWI NetCDF (time/lat/lon, uwnd/vwnd/press,
                        pressure in mb) — verifies user-mapping path for
                        non-default variable names
    """
    storm_path = tmp_path / "test.storm"
    data_storm = storm.Storm()
    data_storm.time_offset = np.datetime64("2012-08-29")
    data_storm.met_crop_extent = [-100.0, -60.0, 10.0, 40.0]
    data_storm.ramp_width = 3

    if file_format == "ascii":
        data_storm.file_format = "ascii"
        data_storm.file_paths = [
            Path("storm_1.PRE"),
            Path("storm_1.WIN"),
            Path("storm_2.PRE"),
            Path("storm_2.WIN"),
        ]
        data_storm.write(storm_path, file_format="data")
        read_storm = storm.Storm(storm_path, file_format="data")

    elif file_format == "netcdf":
        data_storm.file_format = "netcdf"
        data_storm.file_paths = [tmp_path / "storm.nc"]
        create_netcdf_storm_file(data_storm.file_paths[0])
        data_storm.write(
            storm_path,
            file_format="data",
            dim_mapping={"t": "valid_time"},
        )
        read_storm = storm.Storm(storm_path, file_format="data")

    elif file_format == "netcdf_era5":
        # ERA5-style: dims (valid_time, latitude, longitude),
        # vars (u10, v10, msl).  write_data uses MetInspector, which
        # discovers "valid_time" via CF axis='T' automatically; dim_mapping
        # is accepted for backwards compatibility but not required.
        pytest.importorskip("netCDF4")
        data_storm.file_format = "netcdf"
        data_storm.file_paths = [tmp_path / "era5.nc"]
        create_era5_storm_file(data_storm.file_paths[0])
        data_storm.write(
            storm_path,
            file_format="data",
            dim_mapping={"t": "valid_time"},
        )
        read_storm = storm.Storm(storm_path, file_format="data")

        # Verify the descriptor body contains the correct coordinate and
        # variable names in the new &file_info / &variable_info format.
        desc_text = storm_path.read_text()
        assert "lon_name       = longitude" in desc_text
        assert "lat_name       = latitude" in desc_text
        assert "time_name      = valid_time" in desc_text
        assert "var_name=u10" in desc_text
        assert "var_name=v10" in desc_text
        assert "var_name=msl" in desc_text

    elif file_format == "netcdf_nws13":
        # NWS13-style: dims (time, lat, lon), vars (uwnd, vwnd, press, mb).
        # write_data uses MetInspector; variable names require an explicit
        # var_mapping because "uwnd"/"vwnd"/"press" are not in the default
        # fallback lists used by get_netcdf_names.
        pytest.importorskip("netCDF4")
        data_storm.file_format = "nws13"
        data_storm.file_paths = [tmp_path / "nws13.nc"]
        create_nws13_storm_file(data_storm.file_paths[0])
        data_storm.write(
            storm_path,
            file_format="data",
            var_mapping={
                "wind_u": "uwnd",
                "wind_v": "vwnd",
                "pressure": "press",
            },
        )
        read_storm = storm.Storm(storm_path, file_format="data")

        # Verify the descriptor body contains the correct coordinate and
        # variable names in the new &file_info / &variable_info format.
        desc_text = storm_path.read_text()
        assert "lon_name       = lon" in desc_text
        assert "lat_name       = lat" in desc_text
        assert "time_name      = time" in desc_text
        assert "var_name=uwnd" in desc_text
        assert "var_name=vwnd" in desc_text
        assert "var_name=press" in desc_text

    assert data_storm.time_offset == read_storm.time_offset
    assert data_storm.file_format in DATA_FILE_FORMAT_MAP[read_storm.file_format]
    assert data_storm.met_crop_extent == read_storm.met_crop_extent
    assert data_storm.ramp_width == read_storm.ramp_width
    assert data_storm.storm_time_scale == read_storm.storm_time_scale
    assert len(data_storm.file_paths) == len(read_storm.file_paths)
    for i, path in enumerate(data_storm.file_paths):
        assert read_storm.file_paths[i] == path


def _met_roles(meta):
    """Return {geoclaw_role: var_name} from a MetMetadata instance."""
    return {v.geoclaw_role: v.var_name for v in meta.variables}


@pytest.mark.python
@pytest.mark.netcdf
def test_netcdf_var_mapping(tmp_path):
    """MetInspector discovers coords/vars by common names (no attrs)."""
    pytest.importorskip("netCDF4")
    from clawpack.geoclaw.netcdf_utils import MetInspector

    storm_data_file = tmp_path / "storm.nc"
    create_netcdf_storm_file(storm_data_file)

    with MetInspector(storm_data_file) as mi:
        meta = mi.inspect_met()

    assert meta.lon_name == "longitude"
    assert meta.lat_name == "latitude"
    assert meta.time_name == "valid_time"
    assert _met_roles(meta) == {"wind_u": "u", "wind_v": "v",
                                "pressure": "pressure"}


@pytest.mark.python
@pytest.mark.netcdf
def test_netcdf_var_mapping_era5(tmp_path):
    """MetInspector discovers ERA5 coords/vars via CF attributes.

    ERA5 (ECMWF Reanalysis v5) is one example of a CF-compliant netCDF
    met-forcing source supported by GeoClaw; the same discovery logic applies
    to any CF-compliant dataset.  Coordinates are found via ``axis`` attrs
    (``valid_time``/``latitude``/``longitude``) and the wind/pressure
    variables via ``standard_name`` (``eastward_wind``, ``northward_wind``,
    ``air_pressure_at_mean_sea_level``) — no explicit mapping required.
    """
    pytest.importorskip("netCDF4")
    from clawpack.geoclaw.netcdf_utils import MetInspector

    nc_path = tmp_path / "era5.nc"
    create_era5_storm_file(nc_path)

    with MetInspector(nc_path) as mi:
        meta = mi.inspect_met()

    assert meta.lon_name == "longitude"
    assert meta.lat_name == "latitude"
    assert meta.time_name == "valid_time"
    assert _met_roles(meta) == {"wind_u": "u10", "wind_v": "v10",
                                "pressure": "msl"}


@pytest.mark.python
@pytest.mark.netcdf
def test_netcdf_var_mapping_nws13(tmp_path):
    """MetInspector on NWS13/OWI-NetCDF: coords auto, vars need a map.

    NWS13 is the OWI NetCDF met-forcing format used by ADCIRC and GeoClaw.
    Its coordinate names (``lon``/``lat``/``time``) are found by common-name
    fallback, but the variable names (``uwnd``/``vwnd``/``press``) are neither
    standard_name-tagged nor in the fallback lists, so they require an
    explicit ``variable_map``; omitting it raises.
    """
    pytest.importorskip("netCDF4")
    from clawpack.geoclaw.netcdf_utils import MetInspector

    nc_path = tmp_path / "nws13.nc"
    create_nws13_storm_file(nc_path)

    # Without a mapping, the non-standard variable names cannot be resolved.
    with MetInspector(nc_path) as mi:
        with pytest.raises(ValueError, match="met role"):
            mi.inspect_met()

    with MetInspector(nc_path, variable_map={
            "wind_u": "uwnd", "wind_v": "vwnd", "pressure": "press"}) as mi:
        meta = mi.inspect_met()

    assert meta.lon_name == "lon"
    assert meta.lat_name == "lat"
    assert meta.time_name == "time"
    assert _met_roles(meta) == {"wind_u": "uwnd", "wind_v": "vwnd",
                                "pressure": "press"}


@pytest.mark.python
def test_get_netcdf_names_deprecated(tmp_path):
    """util.get_netcdf_names still works but emits a DeprecationWarning."""
    pytest.importorskip("netCDF4")
    storm_data_file = tmp_path / "storm.nc"
    create_netcdf_storm_file(storm_data_file)

    with pytest.warns(DeprecationWarning, match="get_netcdf_names"):
        util.get_netcdf_names(storm_data_file, lookup_type="var")


@pytest.mark.python
@pytest.mark.skip(
    reason=(
        "WRF string-time axis and curvilinear grid require MetPreprocessor, "
        "not yet implemented"
    )
)
def test_netcdf_wrf_stub(tmp_path):
    """Stub for future WRF NetCDF met-forcing test (currently skipped).

    WRF output files differ from ERA5 and NWS13 in two ways that require
    pre-processing before the GeoClaw Fortran reader can consume them:

    1. **String time axis**: WRF stores time as ``char Times(Time, DateStrLen)``
       (ISO 8601 strings) rather than a numeric variable.  GeoClaw's Fortran
       reader expects integer seconds; a ``MetPreprocessor`` step is needed to
       decode and convert the time axis.

    2. **Curvilinear grid**: WRF uses a map-projected (Lambert conformal,
       polar stereographic, or Mercator) curvilinear grid described by 2-D
       ``XLAT``/``XLONG`` arrays rather than 1-D coordinate variables.
       GeoClaw's reader expects regular (lon, lat) grids; regridding to a
       regular lat/lon grid is required before use.

    Once a ``MetPreprocessor`` class (or equivalent) is implemented to handle
    these two cases, this test slot should be filled in with a synthetic WRF
    file (minimal curvilinear grid, string time) and an end-to-end round-trip
    assertion.
    """
    pass  # placeholder — see docstring for implementation notes


@pytest.mark.python
@pytest.mark.parametrize(
    "speeds_knots, expected_categories",
    [
        (np.array([20, 50, 70, 90, 100, 120, 140]), 
         np.array([-1,  0,  1,  2,   3,   4,   5])),
    ],
)
def test_storm_category_nhc(speeds_knots, expected_categories):
    """Test NHC categorization from known wind speeds."""
    s = storm.Storm()
    s.max_wind_speed = speeds_knots * 0.514444  # convert knots to m/s

    categories = s.category(categorization="NHC")

    print(f"{categories}")
    print(f"{expected_categories}")

    assert np.array_equal(categories, expected_categories)


@pytest.mark.python
@pytest.mark.parametrize(
    "speed_knots, expected_category, expected_name",
    [
        (20, -1, "Tropical Depression"),
        (50, 0, "Tropical Storm"),
        (80, 1, "Category 1 Hurricane"),
    ],
)
def test_storm_category_names(speed_knots, expected_category, expected_name):
    """Test NHC category-name output."""
    s = storm.Storm()
    s.max_wind_speed = np.array([speed_knots]) * 0.514444  # convert knots to m/s

    categories, names = s.category(categorization="NHC", cat_names=True)

    assert categories[0] == expected_category
    assert expected_name in names[0]


@pytest.mark.python
def test_storm_plot_smoke(tmp_path):
    """Smoke test for Storm.plot (no assertions, just ensure it runs)."""
    plt = pytest.importorskip("matplotlib.pyplot")

    s = storm.Storm()

    # Minimal valid track
    s.t = np.array([
        np.datetime64("2020-01-01"),
        np.datetime64("2020-01-02"),
    ])
    s.eye_location = np.array([
        [-90.0, 25.0],
        [-89.5, 25.5],
    ])

    fig, ax = plt.subplots()
    s.plot(ax)

    # Save to tmp_path just to ensure backend works
    fig.savefig(tmp_path / "storm_plot.png")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "save":
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else data_dir
        save_storm_test_data(output_dir)
    else:
        raise SystemExit(pytest.main([__file__]))
