#!/usr/bin/env python
r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from pathlib import Path
import pytest

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.fgmax_tools as fgmax_tools
import clawpack.geoclaw.topotools as topotools

@pytest.mark.regression
def test_bowl_slosh(tmp_path: Path, save: bool):
    """Bowl-Slosh regression test for GeoClaw"""

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    # Setup data for test
    runner.set_data()

    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.5

    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.5, 0.5, 0, 1e10])

    # == fgmax.data values ==
    runner.rundata.fgmax_data.num_fgmax_val = 2
    fgmax_grids = runner.rundata.fgmax_data.fgmax_grids
    
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2       # will specify a 2d grid of points
    fg.x1 = -2.0
    fg.x2 = 2.0
    fg.y1 = -2.0
    fg.y2 = 2.0
    fg.dx = 0.1
    fg.tstart_max = 0.0       # when to start monitoring max values
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 0.1         # target time (sec) increment between updating
                            # max values
    fg.min_level_check = 2    # which levels to monitor max on
    fg.arrival_tol = 1.e-2    # tolerance for flagging arrival

    fgmax_grids.append(fg)  # written to fgmax_grids.data

    runner.write_data()

    # Build Topography
    a = 1.
    h0 = 0.1
    topo_func = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

    topo = topotools.Topography(topo_func=topo_func)
    topo.topo_type = 2
    topo.x = np.linspace(-2.0, 2.0, 200)
    topo.y = np.linspace(-2.0, 2.0, 200)
    topo.write(Path(runner.temp_path) / "bowl.topotype2", topo_type=2, 
                                                          Z_format="%22.15e")

    # Build and run test
    runner.build_executable()
    runner.run_code()

    # Check results
    runner.check_gauge(save=save, gauge_id=1, indices=(2, 3))
    runner.check_fgmax(save=save)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
