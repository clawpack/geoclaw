"""Pytest regression test for the GeoClaw Isaac storm-surge example."""

from __future__ import annotations

from pathlib import Path
from datetime import datetime
import pytest

import numpy as np

import clawpack.geoclaw.test as test
from clawpack.geoclaw.surge.storm import Storm

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
    assert generated.window_type == regression.window_type
    np.testing.assert_allclose(generated.ramp_width, regression.ramp_width)
    assert [path.name for path in generated.file_paths] == [
        path.name for path in regression.file_paths
    ]


# @pytest.mark.slow

@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.parametrize(
    "data_file_format",
    ["holland80", "owi_ascii"],
)
def test_isaac(tmp_path: Path, data_file_format: str, save: bool) -> None:
    """Regression test for Isaac storm surge using several storm-input modes."""
    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    runner.set_data()
    # TODO: Decide on time range for Isaac test; this is currently set to match
    # the Isaac test, but may need to be adjusted.
    runner.rundata.clawdata.t0 = days2seconds(-1)
    runner.rundata.clawdata.tfinal = days2seconds(-0.5)
    runner.rundata.clawdata.num_output_times = 1

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
        # write_owi_pressure(isaac.file_paths[0],
        #                    pressure_fields=[isaac.pressure_field],
        #                    times=[isaac.time_offset.astype(datetime)],
        #                    lon=isaac.lon,
        #                    lat=isaac.lat,
        # )
        # write_owi_wind(isaac.file_paths[1],
        #                u_fields=[isaac.u_wind_field],
        #                v_fields=[isaac.v_wind_field],
        #                times=[isaac.time_offset.astype(datetime)],
        #                lon=isaac.lon,
        #                lat=isaac.lat,
        # )
    # elif data_file_format == "owi_netcdf":
    #     surge_data.storm_specification_type = "data"
    #     isaac.data_file_format = "NWS13"
    #     isaac.file_paths = [tmp_path / "isaac.nc"]
    #     isaac.write(surge_data.storm_file, file_format="data")
    # TODO: Add test for generic NetCDF formatted storms (should be a super set
    # of NWS13)
    # elif data_file_format == "netcdf": 
    #     surge_data.storm_specification_type = "data"
    #     isaac.data_file_format ="netcdf" 
    #     isaac.file_paths = [tmp_path / "isaac.nc"]
    #     isaac.write(surge_data.storm_file, file_format="data")
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
    else:
        raise ValueError(f"Unsupported data_file_format={data_file_format}")

    # Run geoclaw
    runner.build_executable()
    runner.run_code()

    # Gauge checks - may need to be format specific for interpolation
    # differences
    runner.check_gauge(gauge_id=1, regression_path=check_path, save=save)
    runner.check_gauge(gauge_id=2, regression_path=check_path, save=save)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
