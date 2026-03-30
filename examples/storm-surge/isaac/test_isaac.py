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


# @pytest.mark.slow

@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.parametrize(
    "data_file_format",
    ["holland80", "owi_ascii"],
)
def test_isaac(tmp_path: Path, data_file_format: str, save: bool) -> None:
    """Regression test for Isaac storm surge using several storm-input modes."""
    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    runner.set_data()
    # TODO: Decide on time range for Isaac test; this is currently set to match
    # the Isaac test, but may need to be adjusted.
    runner.rundata.clawdata.t0 = days2seconds(-1)
    runner.rundata.clawdata.tfinal = days2seconds(-0.5)

    surge_data = runner.rundata.surge_data
    surge_data.storm_file = tmp_path / "isaac.storm"

    # Use a committed local ATCF file as the canonical Isaac source.
    atcf_path = example_dir / "bal092012.dat"

    isaac = Storm(path=atcf_path, file_format="ATCF")
    isaac.time_offset = np.datetime64("2012-08-29")
    
    if data_file_format == "holland80":
        surge_data.storm_specification_type = "holland80"
        isaac.write(surge_data.storm_file, file_format="geoclaw")
    elif data_file_format == "owi_ascii":
        surge_data.storm_specification_type = "data"
        isaac.file_format = "NWS12"
        isaac.file_paths = [example_dir / "isaac.PRE", 
                            example_dir / "isaac.WIN"]
        
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

    # Check storm files are correct
    check_path = example_dir / "regression_data" / data_file_format
    print(save)
    runner.check_files_equal(surge_data.storm_file, check_path / "isaac.storm", 
                             save=save)

    # Run geoclaw
    runner.build_executable()
    runner.run_code()

    # Gauge checks - may need to be format specific for interpolation
    # differences
    runner.check_gauge(gauge_id=1, regression_path=check_path, save=save)
    runner.check_gauge(gauge_id=2, regression_path=check_path, save=save)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
