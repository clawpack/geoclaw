"""Pytest regression test for the GeoClaw Isaac storm-surge example."""

from pathlib import Path
import gzip
import shutil

import numpy as np
import pytest

import clawpack.geoclaw.test as test
from clawpack.geoclaw.surge.storm import Storm


def days2seconds(t):
    """Convert days to seconds."""
    return t * 60.0**2 * 24.0

# @pytest.mark.slow

@pytest.mark.regression
@pytest.mark.parametrize(
    "data_file_format",
    ["geoclaw", "owi_ascii", "owi_netcdf"],
)
def test_isaac_formats(tmp_path: Path, data_file_format: str) -> None:
    """Regression test for Isaac storm surge using several storm-input modes."""
    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    runner.set_data()
    # TODO: Decide on time range for Isaac test; this is currently set to match the Ike test, but may need to be adjusted.
    runner.rundata.clawdata.t0 = days2seconds(-2)
    runner.rundata.clawdata.tfinal = days2seconds(1.5)

    surge_data = runner.rundata.surge_data
    surge_data.storm_file = str(tmp_path / "isaac.storm")

    # Use a committed local ATCF file as the canonical Isaac source.
    atcf_path = example_dir / "bal092012.dat"

    isaac = Storm(path=atcf_path, file_format="ATCF")
    isaac.time_offset = np.datetime64("2012-08-29")

    if data_file_format == "geoclaw":
        surge_data.storm_specification_type = "holland80"
        isaac.write(surge_data.storm_file, file_format="geoclaw")
    elif data_file_format == "owi_ascii":
        surge_data.storm_specification_type = "OWI"
        isaac.data_file_format = "NWS12"
        isaac.file_paths = [tmp_path / "isaac.PRE", tmp_path / "isaac.WIN"]
        isaac.write(surge_data.storm_file, file_format="data")
    elif data_file_format == "owi_netcdf":
        surge_data.storm_specification_type = "OWI"
        isaac.data_file_format = "NWS13"
        isaac.file_paths = [tmp_path / "isaac.nc"]
        isaac.write(surge_data.storm_file, file_format="data")
    elif data_file_format == "netcdf":
        isaac.data_file_format = "netcdf"
        isaac.file_paths = [tmp_path / "isaac.nc"]
        isaac.write(surge_data.storm_file, file_format="data")
    else:
        raise ValueError(f"Unsupported data_file_format={data_file_format}")

    runner.write_data()
    runner.build_executable()
    runner.run_code()

    # Gauge checks are re-enabled here; tighten or broaden later as needed.
    runner.check_gauge(gauge_id=1)
    runner.check_gauge(gauge_id=2)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
