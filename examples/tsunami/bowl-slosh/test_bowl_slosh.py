#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from pathlib import Path
import sys
import unittest

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.fgmax_tools as fgmax_tools

class BowlSloshTest(test.GeoClawRegressionTest):

    r"""Bowl-Slosh regression test for GeoClaw"""

    def setUp(self):

        super(BowlSloshTest, self).setUp()

        # fgmax_grids.data created by setrun.py now contains all info
        #from . import make_fgmax_grid 
        #make_fgmax_grid.make_fgmax_grid1(self.temp_path)


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test bowl-slosh example

        Note that this stub really only runs the code and performs no tests.

        """

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.5

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 0.5, 0.5, 0, 1e10])

        self.rundata.refinement_data.deep_depth = 1e2
        self.rundata.refinement_data.max_level_deep = 3

        # == fgmax.data values ==
        self.rundata.fgmax_data.num_fgmax_val = 2
        fgmax_grids = self.rundata.fgmax_data.fgmax_grids  # empty list to start
        # Now append to this list objects of class fgmax_tools.FGmaxGrid
        # specifying any fgmax grids.

        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2       # will specify a 2d grid of points
        fg.x1 = -2.
        fg.x2 = 2.
        fg.y1 = -2.
        fg.y2 = 2.
        fg.dx = 0.1
        fg.tstart_max = 0.        # when to start monitoring max values
        fg.tend_max = 1.e10       # when to stop monitoring max values
        fg.dt_check = 0.1         # target time (sec) increment between updating
                               # max values
        fg.min_level_check = 2    # which levels to monitor max on
        fg.arrival_tol = 1.e-2    # tolerance for flagging arrival

        fgmax_grids.append(fg)  # written to fgmax_grids.data

        self.write_rundata_objects()

        # Make topography
        a = 1.
        h0 = 0.1
        topo_func = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-2.0, 2.0, 200)
        topo.y = numpy.linspace(-2.0, 2.0, 200)
        topo.write(Path(self.temp_path) / "bowl.topotype2", topo_type=2, 
                                                            Z_format="%22.15e")

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3))
        self.check_fgmax(save=save)
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = BowlSloshTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
