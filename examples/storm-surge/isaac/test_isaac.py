"""Pytest regression test for the GeoClaw Isaac storm-surge example."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from datetime import datetime
import pytest

import numpy as np

import clawpack.geoclaw.test as test
from clawpack.geoclaw.surge.storm import Storm

# ---------------------------------------------------------------------------
# Import shared NetCDF file generators from tests/test_storm.py.
# No Python package __init__.py exists in tests/, so we insert the directory
# into sys.path rather than using a standard package import.
# ---------------------------------------------------------------------------
_tests_dir = Path(__file__).parents[3] / "tests"
if str(_tests_dir) not in sys.path:
    sys.path.insert(0, str(_tests_dir))
try:
    from test_storm import create_era5_storm_file, create_nws13_storm_file
except ImportError:
    create_era5_storm_file = None
    create_nws13_storm_file = None

# Generic helper functions
def days2seconds(t):
    """Convert days to seconds."""
    return t * 60.0**2 * 24.0


# OWI data format helper functions
def _owi_master_header(kind: str, start: datetime, end: datetime) -> str:
    if kind == "win":
        label = "OWI WWS Wind Output Ucomp,Vcomp in m/s"
    elif kind == "pre":
        label = "OWI WWS Pressure Output in mb"
    else:
        raise ValueError(f"Unknown OWI kind={kind}")
    return f"{label:<56}Start:{start:%Y%m%d%H} End:{end:%Y%m%d%H}"


def _owi_snapshot_header(
    nlat: int,
    nlon: int,
    dx: float,
    dy: float,
    swlat: float,
    swlon: float,
    when: datetime,
) -> str:
    return (
        f"iLat={nlat:4d}"
        f"iLong={nlon:4d}"
        f"DX={dx:6.4f}"
        f"DY={dy:6.4f}"
        f"SWLat={swlat:8.4f}"
        f"SWLon={swlon:8.4f}"
        f"DT={when:%Y%m%d%H%M}"
    )


def _write_owi_values(fobj, values: np.ndarray) -> None:
    flat = np.asarray(values).T.flatten()  # match Fortran/OWI ordering
    for i, value in enumerate(flat, start=1):
        fobj.write(f"{value:10.4f}")
        if i % 8 == 0 or i == len(flat):
            fobj.write("\n")


def write_owi_pressure(
    path: Path,
    pressure_fields: list[np.ndarray],
    times: list[datetime],
    lon: np.ndarray,
    lat: np.ndarray,
) -> None:
    with Path(path).open("w") as fobj:
        fobj.write(_owi_master_header("pre", times[0], times[-1]) + "\n")
        for field, when in zip(pressure_fields, times):
            fobj.write(
                _owi_snapshot_header(
                    nlat=len(lat),
                    nlon=len(lon),
                    dx=float(lon[1] - lon[0]),
                    dy=float(lat[1] - lat[0]),
                    swlat=float(lat[0]),
                    swlon=float(lon[0]),
                    when=when,
                )
                + "\n"
            )
            _write_owi_values(fobj, field)


def write_owi_wind(
    path: Path,
    u_fields: list[np.ndarray],
    v_fields: list[np.ndarray],
    times: list[datetime],
    lon: np.ndarray,
    lat: np.ndarray,
) -> None:
    with Path(path).open("w") as fobj:
        fobj.write(_owi_master_header("win", times[0], times[-1]) + "\n")
        for u_field, v_field, when in zip(u_fields, v_fields, times):
            fobj.write(
                _owi_snapshot_header(
                    nlat=len(lat),
                    nlon=len(lon),
                    dx=float(lon[1] - lon[0]),
                    dy=float(lat[1] - lat[0]),
                    swlat=float(lat[0]),
                    swlon=float(lon[0]),
                    when=when,
                )
                + "\n"
            )
            _write_owi_values(fobj, u_field)
            _write_owi_values(fobj, v_field)


def _read_owi_all_timesteps(pre_path: Path, win_path: Path):
    """Read all time steps from committed OWI PRE and WIN files.

    There is no Python OWI reader in the GeoClaw codebase (the Fortran
    ``read_OWI_ASCII`` routine in ``data_storm_module.f90`` is the only
    reader).  This function uses minimal line-by-line parsing of the
    fixed-format OWI ASCII files.

    The OWI format stores data with longitude varying fastest (inner loop =
    longitude).  The flat data array for one time step is reshaped to
    ``(nlat, nlon)`` in standard numpy ``[lat, lon]`` indexing.

    Returns
    -------
    pressure_mb : ndarray, shape (nt, nlat, nlon)
        Sea-level pressure in mb (millibar).
    u_wind : ndarray, shape (nt, nlat, nlon)
        Eastward wind speed in m/s.
    v_wind : ndarray, shape (nt, nlat, nlon)
        Northward wind speed in m/s.
    lon : ndarray, shape (nlon,)
        Longitude grid in degrees.  First value = SWLon + DX (matching the
        1-based grid indexing used by the Fortran OWI reader).
    lat : ndarray, shape (nlat,)
        Latitude grid in degrees.  First value = SWLat + DY.
    times : list of datetime
        Datetime objects for each time step (parsed from DT= field).
    """

    def _parse_snapshot_header(line):
        """Extract grid parameters and datetime from an OWI snapshot header."""
        nlat = int(re.search(r'iLat=\s*(\d+)', line).group(1))
        nlon = int(re.search(r'iLong=\s*(\d+)', line).group(1))
        # OWI snapshot headers are fixed-width with NO spaces between fields,
        # e.g. "DX=0.2500DY=0.2500SWLat=...".  Use digit/decimal-only patterns
        # so the match stops at the next field keyword, not at whitespace.
        dx = float(re.search(r'DX=([\d.]+)', line).group(1))
        dy = float(re.search(r'DY=([\d.]+)', line).group(1))
        swlat = float(re.search(r'SWLat=\s*([-\d.]+)', line).group(1))
        swlon = float(re.search(r'SWLon=\s*([-\d.]+)', line).group(1))
        dt_str = re.search(r'DT=(\d{12})', line).group(1)
        dt = datetime.strptime(dt_str, '%Y%m%d%H%M')
        return nlat, nlon, dx, dy, swlat, swlon, dt

    def _read_n_values(fobj, n):
        """Read exactly n whitespace-separated floats from an OWI file."""
        vals = []
        while len(vals) < n:
            vals.extend(float(x) for x in fobj.readline().split())
        return np.array(vals[:n])

    # --- Read pressure (all time steps) ------------------------------------
    # First pass: determine grid dimensions from the first snapshot header.
    with open(pre_path) as f:
        f.readline()                              # master header
        header = f.readline()
    nlat, nlon, dx, dy, swlat, swlon, _ = _parse_snapshot_header(header)
    n_per_step = nlat * nlon

    # The Fortran OWI reader computes coordinates as 1-based:
    #   longitude(j) = swlon + j * dx,  j = 1..mx
    # so the first grid point is swlon+dx, not swlon.
    lon = swlon + np.arange(1, nlon + 1) * dx
    lat = swlat + np.arange(1, nlat + 1) * dy

    pre_timesteps = []
    times = []
    with open(pre_path) as f:
        f.readline()  # master header
        while True:
            header = f.readline()
            if not header.strip():
                break
            _, _, _, _, _, _, dt = _parse_snapshot_header(header)
            times.append(dt)
            flat = _read_n_values(f, n_per_step)
            # OWI: lon varies fastest -> reshape as (nlat, nlon)
            pre_timesteps.append(flat.reshape(nlat, nlon))

    # --- Read wind (all time steps) ----------------------------------------
    u_timesteps = []
    v_timesteps = []
    with open(win_path) as f:
        f.readline()  # master header
        for _ in times:
            f.readline()                              # snapshot header
            flat_u = _read_n_values(f, n_per_step)
            flat_v = _read_n_values(f, n_per_step)
            u_timesteps.append(flat_u.reshape(nlat, nlon))
            v_timesteps.append(flat_v.reshape(nlat, nlon))

    pressure_mb = np.stack(pre_timesteps, axis=0)   # (nt, nlat, nlon)
    u_wind = np.stack(u_timesteps, axis=0)
    v_wind = np.stack(v_timesteps, axis=0)

    return pressure_mb, u_wind, v_wind, lon, lat, times


def _make_isaac_netcdf(fmt: str, tmp_path: Path, test_path: Path,
                       met_crop_extent=None) -> Path:
    """Generate a NetCDF met-forcing file from the committed OWI Isaac files.

    Reads all time steps from ``isaac.PRE`` and ``isaac.WIN`` using minimal
    OWI line-parsing (no Python OWI reader exists in the codebase; flagged
    above in ``_read_owi_all_timesteps``), converts units, writes a NetCDF
    file to ``tmp_path``, and writes the corresponding ``.storm`` descriptor
    via ``Storm.write(file_format='data')``.

    Both ERA5 and NWS13 variants write pressure in **Pa** (as required by the
    Fortran NetCDF reader, which has no unit conversion).  The ERA5 variant
    also stores wind in m/s; the NWS13 file uses the same values but with
    NWS13 variable/dimension names (``uwnd``, ``vwnd``, ``press``).

    Parameters
    ----------
    fmt : str
        One of ``"netcdf_era5"`` or ``"netcdf_nws13"``.
    tmp_path : Path
        Directory for generated files.
    test_path : Path
        Path to the Isaac example directory (for finding committed OWI files).

    Returns
    -------
    storm_path : Path
        Path to the written ``.storm`` descriptor.
    """
    if create_era5_storm_file is None or create_nws13_storm_file is None:
        pytest.skip("NetCDF file generators not importable from test_storm.py")

    pytest.importorskip("netCDF4")

    pre_path = test_path / "isaac.PRE"
    win_path = test_path / "isaac.WIN"

    pressure_mb, u_wind, v_wind, lon, lat, times = _read_owi_all_timesteps(
        pre_path, win_path
    )

    # Convert OWI mb pressure to Pa (the Fortran OWI reader multiplies by
    # 100; the NetCDF reader reads raw values, so we must pre-convert).
    pressure_pa = pressure_mb * 100.0

    # Convert datetime list to integer UNIX seconds (Fortran uses UNIX epoch
    # 1970-01-01 via seconds_from_epoch in utility_module.f90).
    unix_epoch = np.datetime64("1970-01-01T00:00:00", "s")
    times_unix = np.array(
        [
            int((np.datetime64(dt) - unix_epoch) / np.timedelta64(1, "s"))
            for dt in times
        ],
        dtype=np.int64,
    )

    # lon/lat already adjusted for 1-based OWI grid offset in
    # _read_owi_all_timesteps.

    time_offset = np.datetime64("2012-08-29T00:00:00")

    isaac = Storm()
    isaac.time_offset = time_offset
    isaac.met_crop_extent = met_crop_extent

    if fmt == "netcdf_era5":
        nc_path = tmp_path / "isaac_era5.nc"
        # create_era5_storm_file uses "seconds since 2012-08-29" as the time
        # epoch (see time_epoch_str in test_storm.py).  Convert Unix seconds
        # to seconds relative to that epoch so xarray decodes them correctly
        # and MetInspector computes the right nc_time_offset.
        time_offset_s = int(
            (time_offset - unix_epoch) / np.timedelta64(1, "s")
        )
        create_era5_storm_file(
            nc_path,
            pressure_fields=pressure_pa,
            u_fields=u_wind,
            v_fields=v_wind,
            times=times_unix - time_offset_s,
            lon=lon,
            lat=lat,
        )
        isaac.file_format = "netcdf"
        isaac.file_paths = [nc_path]
        storm_path = tmp_path / "isaac_era5.storm"
        isaac.write(
            storm_path,
            file_format="data",
            dim_mapping={"t": "valid_time"},
        )

    elif fmt == "netcdf_nws13":
        nc_path = tmp_path / "isaac_nws13.nc"
        # NWS13 file: same Pa pressure (Fortran no-conversion path),
        # NWS13 variable names.  The "mb" units attribute in
        # create_nws13_storm_file is overridden by passing Pa values;
        # the Fortran uses only the variable name, not the units attribute.
        create_nws13_storm_file(
            nc_path,
            pressure_fields=pressure_pa,
            u_fields=u_wind,
            v_fields=v_wind,
            times=times_unix,
            lon=lon,
            lat=lat,
        )
        isaac.file_format = "nws13"
        isaac.file_paths = [nc_path]
        storm_path = tmp_path / "isaac_nws13.storm"
        isaac.write(
            storm_path,
            file_format="data",
            var_mapping={
                "wind_u": "uwnd",
                "wind_v": "vwnd",
                "pressure": "press",
            },
        )

    else:
        raise ValueError(f"Unknown fmt={fmt!r}")

    return storm_path


def _check_netcdf_storm_descriptor(generated_path: Path, fmt: str) -> None:
    """Validate structure of a NetCDF data-storm descriptor file.

    Checks that ``generated_path`` contains the correct format number (2),
    the expected coordinate names in the ``&file_info`` block, the expected
    variable/role pairs in ``&variable_info`` blocks, and exactly one file
    path entry.

    This is a structural validation against known-good expected values for
    each format, not a comparison against a committed regression file (which
    would contain machine-specific absolute paths).

    Parameters
    ----------
    generated_path : Path
        Path to the ``.storm`` descriptor written by ``Storm.write``.
    fmt : str
        One of ``"netcdf_era5"`` or ``"netcdf_nws13"``.
    """
    text = generated_path.read_text()
    lines = [ln.rstrip() for ln in text.splitlines()]

    # Locate the format line (contains the integer file format number).
    format_line = next(
        ln for ln in lines if ln.startswith("2 ") or ln.strip().startswith("2")
        and "# File format" in ln
    )
    assert "2" in format_line, (
        f"Expected file format 2 in descriptor, got: {format_line!r}"
    )

    # Locate the "# Format Data Information" section.
    try:
        fmt_info_idx = next(
            i for i, ln in enumerate(lines)
            if "# Format Data Information" in ln
        )
    except StopIteration:
        raise AssertionError(
            f"'# Format Data Information' section not found in {generated_path}"
        )

    # Parse &file_info block for coordinate names.
    file_info: dict = {}
    in_file_info = False
    for ln in lines[fmt_info_idx:]:
        stripped = ln.strip()
        if stripped == "&file_info":
            in_file_info = True
            continue
        if in_file_info:
            if stripped == "/":
                break
            if "=" in stripped:
                key, _, val = stripped.partition("=")
                file_info[key.strip()] = val.strip()

    # Parse &variable_info blocks for role→var_name mapping.
    var_info: dict = {}
    for ln in lines[fmt_info_idx:]:
        stripped = ln.strip()
        if stripped.startswith("&variable_info"):
            kv = dict(re.findall(r'(\w+)=(\S+)', stripped))
            if "geoclaw_role" in kv and "var_name" in kv:
                var_info[kv["geoclaw_role"]] = kv["var_name"].rstrip("/").strip()

    if fmt == "netcdf_era5":
        assert file_info.get("lon_name") == "longitude", (
            f"ERA5 lon_name wrong: {file_info.get('lon_name')!r}"
        )
        assert file_info.get("lat_name") == "latitude", (
            f"ERA5 lat_name wrong: {file_info.get('lat_name')!r}"
        )
        assert file_info.get("time_name") == "valid_time", (
            f"ERA5 time_name wrong: {file_info.get('time_name')!r}"
        )
        assert var_info.get("wind_u") == "u10", (
            f"ERA5 wind_u var wrong: {var_info.get('wind_u')!r}"
        )
        assert var_info.get("wind_v") == "v10", (
            f"ERA5 wind_v var wrong: {var_info.get('wind_v')!r}"
        )
        assert var_info.get("pressure") == "msl", (
            f"ERA5 pressure var wrong: {var_info.get('pressure')!r}"
        )
    elif fmt == "netcdf_nws13":
        assert file_info.get("lon_name") == "lon", (
            f"NWS13 lon_name wrong: {file_info.get('lon_name')!r}"
        )
        assert file_info.get("lat_name") == "lat", (
            f"NWS13 lat_name wrong: {file_info.get('lat_name')!r}"
        )
        assert file_info.get("time_name") == "time", (
            f"NWS13 time_name wrong: {file_info.get('time_name')!r}"
        )
        assert var_info.get("wind_u") == "uwnd", (
            f"NWS13 wind_u var wrong: {var_info.get('wind_u')!r}"
        )
        assert var_info.get("wind_v") == "vwnd", (
            f"NWS13 wind_v var wrong: {var_info.get('wind_v')!r}"
        )
        assert var_info.get("pressure") == "press", (
            f"NWS13 pressure var wrong: {var_info.get('pressure')!r}"
        )
    else:
        raise ValueError(f"Unknown fmt={fmt!r}")

    # There should be exactly one file path entry.
    try:
        paths_idx = next(
            i for i, ln in enumerate(lines) if "# File paths" in ln
        )
    except StopIteration:
        raise AssertionError(
            f"'# File paths' section not found in {generated_path}"
        )
    path_entries = [
        ln for ln in lines[paths_idx + 1:] if ln.strip()
    ]
    assert len(path_entries) == 1, (
        f"Expected 1 file path, got {len(path_entries)}: {path_entries}"
    )


def _check_geoclaw_storm_descriptor(generated_path: Path, regression_path: Path) -> None:
    """Compare two GeoClaw storm files semantically using Storm readers."""
    generated = Storm(path=generated_path, file_format="geoclaw")
    regression = Storm(path=regression_path, file_format="geoclaw")

    assert generated.time_offset == regression.time_offset
    assert generated.t.shape == regression.t.shape
    np.testing.assert_array_equal(generated.t, regression.t)
    np.testing.assert_allclose(generated.eye_location, regression.eye_location)
    np.testing.assert_allclose(generated.max_wind_speed, regression.max_wind_speed)
    np.testing.assert_allclose(generated.max_wind_radius, regression.max_wind_radius)
    np.testing.assert_allclose(generated.central_pressure, regression.central_pressure)
    np.testing.assert_allclose(generated.storm_radius, regression.storm_radius)


def _check_data_storm_descriptor(generated_path: Path, regression_path: Path) -> None:
    """Compare two data-derived storm descriptor files using Storm readers."""
    generated = Storm(path=generated_path, file_format="data")
    regression = Storm(path=regression_path, file_format="data")

    assert generated.time_offset == regression.time_offset
    assert generated.file_format == regression.file_format
    assert len(generated.file_paths) == len(regression.file_paths)
    np.testing.assert_allclose(generated.scaling, regression.scaling)
    assert generated.met_crop_extent == regression.met_crop_extent
    np.testing.assert_allclose(generated.ramp_width, regression.ramp_width)
    assert [path.name for path in generated.file_paths] == [
        path.name for path in regression.file_paths
    ]


@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.parametrize(
    "data_file_format",
    [
        "holland80",
        "owi_ascii",
        pytest.param(
            "netcdf_era5",
            marks=[pytest.mark.netcdf],
        ),
        pytest.param(
            "netcdf_nws13",
            marks=[pytest.mark.netcdf],
        ),
    ],
)
def test_isaac(tmp_path: Path, data_file_format: str, save: bool) -> None:
    """Regression test for Isaac storm surge using several storm-input modes.

    The ``netcdf_era5`` and ``netcdf_nws13`` variants generate NetCDF met-
    forcing files from the committed OWI ASCII files (``isaac.PRE`` /
    ``isaac.WIN``) and verify that GeoClaw produces gauge output identical to
    the ``owi_ascii`` baseline.  The two variants differ only in variable and
    dimension naming conventions:

    ``netcdf_era5``
        CF-1.7 ERA5-style file: dims ``valid_time / latitude / longitude``,
        vars ``u10 / v10 / msl``, lon in [0, 360], lat S-to-N.
        ERA5 (ECMWF Reanalysis v5) is used here as one example of a
        CF-compliant netCDF met-forcing source; GeoClaw supports any
        CF-compliant netCDF data with compatible structure.

    ``netcdf_nws13``
        OWI NWS13-style file: dims ``time / lat / lon``, vars
        ``uwnd / vwnd / press``.  Exercises the explicit ``user_mapping``
        path in ``util.get_netcdf_names`` for non-default variable names.

    Both NetCDF variants store pressure in Pa (after converting from the OWI
    mb values) because the Fortran NetCDF reader has no unit conversion.
    The gauge regression data from ``owi_ascii`` is reused since the
    underlying forcing is identical once unit-converted.
    """
    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    runner.set_data()
    # TODO: Decide on time range for Isaac test, currently is just a half day
    runner.rundata.clawdata.t0 = days2seconds(-1)
    runner.rundata.clawdata.tfinal = days2seconds(-0.5)
    runner.rundata.clawdata.num_output_times = 1

    # TODO: May also want to change number of levels used, currently matches the
    # example data but may be shortened for testing purposes.
    runner.rundata.amrdata.amr_levels_max = 2

    surge_data = runner.rundata.surge_data
    surge_data.storm_file = tmp_path / "isaac.storm"

    # Use a committed local ATCF file as the canonical Isaac source.
    atcf_path = runner.test_path / "bal092012.dat"

    isaac = Storm(path=atcf_path, file_format="ATCF")
    isaac.time_offset = np.datetime64("2012-08-29")

    if data_file_format == "holland80":
        surge_data.storm_specification_type = "holland80"
        isaac.write(surge_data.storm_file, file_format="geoclaw", verbose=True)

    elif data_file_format == "owi_ascii":
        surge_data.storm_specification_type = "data"
        isaac.file_format = "NWS12"
        isaac.file_paths = [runner.test_path / "isaac.PRE",
                            runner.test_path / "isaac.WIN"]
        isaac.write(surge_data.storm_file, file_format="data")
        # TODO: Generate pressure and wind fields from the ATCF data, requires
        # storm field generation, which is not yet supported.

    elif data_file_format in ("netcdf_era5", "netcdf_nws13"):
        # Generate NetCDF file + descriptor from committed OWI data.
        surge_data.storm_file = _make_isaac_netcdf(
            data_file_format, tmp_path, runner.test_path
        )
        surge_data.storm_specification_type = "data"
        # Validate descriptor structure before running.
        _check_netcdf_storm_descriptor(surge_data.storm_file, data_file_format)

    else:
        raise ValueError(f"Unsupported data_file_format={data_file_format}")

    runner.write_data()

    assert Path(surge_data.storm_file).exists()

    # Check storm descriptor files semantically.
    check_path = runner.test_path / "regression_data" / data_file_format
    regression_storm_file = check_path / "isaac.storm"
    if data_file_format == "holland80":
        _check_geoclaw_storm_descriptor(surge_data.storm_file, regression_storm_file)
    elif data_file_format == "owi_ascii":
        _check_data_storm_descriptor(surge_data.storm_file, regression_storm_file)
    # For netcdf variants, structural validation already done above via
    # _check_netcdf_storm_descriptor; no committed regression .storm file.

    # Run geoclaw
    runner.build_executable()
    runner.run_code()

    # Gauge checks.
    # The netcdf variants use the owi_ascii regression gauge data because the
    # underlying forcing (from the same committed OWI files, unit-converted)
    # is identical.  If the simulation results differ from the owi_ascii
    # baseline, stop and investigate before generating new reference files.
    if data_file_format in ("netcdf_era5", "netcdf_nws13"):
        gauge_regression_path = runner.test_path / "regression_data" / "owi_ascii"
    else:
        gauge_regression_path = check_path
    runner.check_gauge(gauge_id=1, regression_path=gauge_regression_path, save=save)
    runner.check_gauge(gauge_id=2, regression_path=gauge_regression_path, save=save)


@pytest.mark.storm
@pytest.mark.netcdf
def test_isaac_netcdf_crop(tmp_path: Path) -> None:
    """met_crop_extent restricts the forcing to the cropped sub-region.

    The met file is cropped to the southern Gulf (lat 8-15), well south of
    both Isaac gauges (lat ~29).  With the forcing cropped out at the gauges,
    no wind or pressure perturbation ever reaches them, so the sea stays at
    rest: gauge momentum and surface elevation remain ~0 for the whole run.
    The uncropped netcdf variants (test_isaac[netcdf_*]) show the same gauges
    develop a real surge, so this isolates the crop's spatial restriction.
    """
    pytest.importorskip("netCDF4")
    import clawpack.pyclaw.gauges as gauges

    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)
    runner.set_data()
    runner.rundata.amrdata.amr_levels_max = 2
    # Short run: the crop's effect at the gauges is time-independent (forcing
    # there is always ambient), so a brief integration suffices.
    runner.rundata.clawdata.num_output_times = 1

    surge_data = runner.rundata.surge_data
    # Crop to the southern Gulf, excluding both gauges (lat ~29).
    surge_data.storm_file = _make_isaac_netcdf(
        "netcdf_era5", tmp_path, example_dir,
        met_crop_extent=[-99.0, -70.0, 8.0, 15.0],
    )
    surge_data.storm_specification_type = "data"

    # Confirm the crop line reached the descriptor.
    assert "8.0 15.0" in surge_data.storm_file.read_text() or \
           "8.0" in surge_data.storm_file.read_text()

    runner.write_data()
    runner.build_executable()
    runner.run_code()

    for gauge_id in (1, 2):
        g = gauges.GaugeSolution(gauge_id, path=runner.temp_path)
        # q = (h, hu, hv, eta).  With forcing cropped out here, the sea is at
        # rest: momentum and surface perturbation stay negligible.
        assert np.max(np.abs(g.q[1, :])) < 1e-3, \
            f"gauge {gauge_id} hu not at rest: {np.max(np.abs(g.q[1, :]))}"
        assert np.max(np.abs(g.q[2, :])) < 1e-3, \
            f"gauge {gauge_id} hv not at rest: {np.max(np.abs(g.q[2, :]))}"
        assert np.max(np.abs(g.q[3, :])) < 1e-3, \
            f"gauge {gauge_id} eta not at rest: {np.max(np.abs(g.q[3, :]))}"


@pytest.mark.storm
def test_isaac_model_storm_time_scale(tmp_path: Path) -> None:
    """storm_time_scale changes the model-storm (Holland) forcing.

    Slowing the storm 4x (storm_time_scale=4.0) shifts the eye position at
    every sim time relative to the unscaled run, so the gauge series must
    differ from the committed holland80 (scale=1.0) regression.  Confirms
    the model-storm time-scale path in set_storm_fields is active.
    """
    import clawpack.pyclaw.gauges as gauges

    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)
    runner.set_data()
    runner.rundata.clawdata.t0 = days2seconds(-1)
    runner.rundata.clawdata.tfinal = days2seconds(-0.5)
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.amrdata.amr_levels_max = 2

    surge_data = runner.rundata.surge_data
    surge_data.storm_file = tmp_path / "isaac.storm"
    surge_data.storm_specification_type = "holland80"
    surge_data.storm_time_scale = 4.0   # 4x slower storm

    isaac = Storm(path=example_dir / "bal092012.dat", file_format="ATCF")
    isaac.time_offset = np.datetime64("2012-08-29")
    isaac.write(surge_data.storm_file, file_format="geoclaw")

    # Confirm the parameter reached surge.data.
    runner.write_data()
    assert "storm_time_scale" in (runner.temp_path / "surge.data").read_text()

    runner.build_executable()
    runner.run_code()

    scaled = gauges.GaugeSolution(1, path=runner.temp_path)
    regression = gauges.GaugeSolution(
        1, path=example_dir / "regression_data" / "holland80")

    # Adaptive time-stepping gives different sample counts when the forcing
    # differs, so interpolate the scaled surface elevation onto the
    # regression's gauge times before comparing.
    eta_scaled = np.interp(regression.t, scaled.t, scaled.q[3, :])
    assert not np.allclose(eta_scaled, regression.q[3, :],
                           rtol=1e-3, atol=1e-3), \
        "storm_time_scale=4.0 produced the same gauge as the scale=1.0 baseline"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
