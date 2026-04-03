#!/usr/bin/env python
"""Regression test for storm surge based on Hurricane Ike"""

from pathlib import Path
import pytest
import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

CASES = [
    pytest.param(
        {"num_cells": [29, 24], 
         "amr_levels_max": 2,
         "num_output_times": 1},
        id="coarse",
    ),
    pytest.param(
        {"num_cells": [29 * 4, 24 * 4], 
         "amr_levels_max": 6, 
         "num_output_times": 16},
        id="fine",
        marks=pytest.mark.slow,
    ),
]

@pytest.mark.storm
@pytest.mark.regression
@pytest.mark.parametrize("case", CASES)
def test_ike(case: dict, tmp_path: Path, save: bool):
    r"""Hurricane Ike regression test for GeoClaw"""

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    # Setup data for test
    runner.set_data()

    runner.rundata.clawdata.num_cells = case["num_cells"]
    runner.rundata.clawdata.num_output_times = case["num_output_times"]
    runner.rundata.amrdata.amr_levels_max = case["amr_levels_max"]
    runner.rundata.surge_data.storm_file = runner.test_path / 'ike.storm'

    runner.write_data()

    # Build Topography - called the same as in the default setrun.py
    topo = topotools.Topography()
    topo.x = np.linspace(-100, -69, 125)
    topo.y = np.linspace(7.0, 33.0, 105)
    topo.Z = 25.0 * ((topo.X + 84.5)**2 + (topo.Y - 20.0)**2) - 4000.0
    topo.write(runner.temp_path / 'gulf_caribbean.tt3', topo_type=2,
                Z_format="%22.15e")

    # Build and run test
    runner.build_executable()
    runner.run_code()

    # Check results
    check_path = runner.test_path / "regression_data" / f"{case['amr_levels_max']}_levels"
    if case["amr_levels_max"] == 2:
        runner.check_gauge(save=save, gauge_id=2, regression_path=check_path)
        runner.check_gauge(save=save, gauge_id=3, regression_path=check_path)
    else:
        # For the finer test, we need to figure out somewhere to put the
        # regression data still
        # runner.check_gauge(save=save, gauge_id=2, regression_path=check_path)
        # runner.check_gauge(save=save, gauge_id=3, regression_path=check_path)
        pass

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
