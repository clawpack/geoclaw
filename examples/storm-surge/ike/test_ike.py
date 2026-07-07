#!/usr/bin/env python
"""Regression test for storm surge based on Hurricane Ike"""

from pathlib import Path
import pytest

import clawpack.geoclaw.test as test

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
        marks=[
            pytest.mark.slow,
            pytest.mark.skip(
                reason="Appears to hang/blow up (unbounded AMR grid growth "
                "on level 5) in CI; disabled pending investigation."
            ),
        ],
    ),
]

@pytest.mark.storm
@pytest.mark.regression
@pytest.mark.remote
@pytest.mark.parametrize("case", CASES)
def test_ike(case: dict, tmp_path: Path, download_cache: Path, save: bool):
    r"""Hurricane Ike regression test for GeoClaw"""

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    # Setup data for test.  setrun.py downloads gulf_caribbean.tt3; direct it
    # into the shared download_cache instead of $CLAW/geoclaw/scratch.
    runner.set_data(download_dir=download_cache)

    runner.rundata.clawdata.num_cells = case["num_cells"]
    runner.rundata.clawdata.num_output_times = case["num_output_times"]
    runner.rundata.amrdata.amr_levels_max = case["amr_levels_max"]
    runner.rundata.surge_data.storm_file = runner.test_path / 'ike.storm'

    runner.write_data()

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
