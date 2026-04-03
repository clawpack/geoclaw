#!/usr/bin/env python
"""Pytest regression test for multilayer plane wave example."""

from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.test as test

def transform_p2c(x, y, x0, y0, theta):
    return (
        x * np.cos(theta) + y * np.sin(theta) - x0,
        -x * np.sin(theta) + y * np.cos(theta) - y0,
    )


def bathy_step(x, y, location=0.15, angle=0.0, left=-1.0, right=-0.2):
    x_c, y_c = transform_p2c(x, y, location, 0.0, angle)
    return (x_c <= 0.0) * left + (x_c > 0.0) * right


@pytest.mark.multilayer
@pytest.mark.regression
def test_multilayer_plane_wave(tmp_path: Path, save: bool):
    """Regression test for multilayer shallow water plane wave example."""

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    # Setup data for test
    runner.set_data()
    
    runner.rundata.clawdata.lower[0] = -1.0
    runner.rundata.clawdata.upper[0] = 2.0
    runner.rundata.clawdata.lower[1] = -1.0
    runner.rundata.clawdata.upper[1] = 2.0

    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 1.0

    runner.rundata.amrdata.refinement_ratios_x = [2, 6]
    runner.rundata.amrdata.refinement_ratios_y = [2, 6]
    runner.rundata.amrdata.refinement_ratios_t = [2, 6]

    runner.rundata.geo_data.rho = [0.9, 1.0]

    runner.rundata.topo_data.topofiles = []
    topo_path = tmp_path / 'jump_topo.tt2'
    runner.rundata.topo_data.topofiles.append([2, topo_path])
    
    runner.rundata.multilayer_data.wave_tolerance = [0.1, 0.2]

    runner.rundata.qinit_data.angle = np.pi / 4.0

    runner.write_data()

    # Create topography file
    topo_func = lambda x, y: bathy_step(
        x, y,
        location=0.15,
        angle=np.pi / 8.0,
        left=-1.0,
        right=-0.2,
    )

    topo = topotools.Topography(topo_func=topo_func)
    topo.x = np.linspace(-1.16, 2.16, 166)
    topo.y = np.linspace(-1.16, 2.16, 166)
    topo.write(topo_path, topo_type=2)

    # Build and run test
    runner.build_executable()
    runner.run_code()

    # Check gauge outputs (surface heights)
    for gauge_id in range(5):
        runner.check_gauge(gauge_id=gauge_id, indices=(6, 7), atol=1e-5, save=save)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
