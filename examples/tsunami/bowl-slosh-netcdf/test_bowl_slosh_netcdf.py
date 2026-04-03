"""Pytest regression test for the GeoClaw bowl-slosh NetCDF example."""

from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import pytest

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools


def _make_bowl_netcdf_topography(output_dir: Path) -> None:
    """
    Create the NetCDF topography input used by the example.
    
    We are directly creating the NetCDF file here rather than relying on the
    example's maketopo.py to test out a slightly generic NetCDF writing approach
    that doesn't rely on the topotools write method.™
    """
    netCDF4 = pytest.importorskip("netCDF4")

    a = 1.0
    h0 = 0.1
    topo_func = lambda x, y: h0 * (x**2 + y**2) / a**2 - h0

    topo = topotools.Topography(topo_func=topo_func)
    topo.x = np.linspace(-3.1, 3.1, 310)
    topo.y = np.linspace(-3.5, 2.5, 300)

    # Intentionally create latitude first to exercise dimension discovery.
    with netCDF4.Dataset(output_dir / "bowl.nc", "w") as out:
        out.createDimension("lat", len(topo.y))
        out.createDimension("lon", len(topo.x))

        latitudes = out.createVariable("lat", "f8", ("lat",))
        longitudes = out.createVariable("lon", "f8", ("lon",))
        elevations = out.createVariable("elevation", "f8", ("lat", "lon"))

        latitudes[:] = topo.y
        longitudes[:] = topo.x
        elevations[:] = topo.Z


@pytest.mark.regression
@pytest.mark.netcdf
def test_bowl_slosh_netcdf(tmp_path: Path, save: bool) -> None:
    """Regression test for bowl-slosh using NetCDF topography input."""
    pytest.importorskip("netCDF4")

    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    # Generate topography
    _make_bowl_netcdf_topography(tmp_path)

    runner.set_data()
    runner.rundata.clawdata.lower[1] = -3.0
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.5
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.5, 0.5, 0.0, 1e10])
    runner.write_data()

    runner.build_executable()
    runner.run_code()

    runner.check_gauge(gauge_id=1, indices=(2, 3), rtol=1.0e-4, atol=1.0e-4, save=save)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
