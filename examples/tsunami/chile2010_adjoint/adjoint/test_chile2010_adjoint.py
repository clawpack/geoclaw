"""
Test for the adjoint problem for a tsunami
"""

import sys
import unittest

import clawpack.geoclaw.test as test

class Chile2010AdjointTest(test.GeoClawRegressionTest):

    def runTest(self, save=False):

        # Make topography and qinit
        import maketopo
        maketopo.get_topo(path=self.temp_path)
        maketopo.makeqinit(path=self.temp_path, center=(-76.0, -36.0))

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_cells[0] = 240
        self.rundata.clawdata.num_cells[1] = 240

        self.rundata.clawdata.num_output_times = 4
        self.rundata.clawdata.tfinal = 3600.0

        self.rundata.amrdata.amr_levels_max = 1
        self.rundata.amrdata.flag2refine = False

        self.rundata.gaugedata.gauges.append([1, -76, -36., 0., 1.e10])

        self.rundata.flagregiondata.flagregions = []

        self.write_rundata_objects()

        self.run_code()

        # Perform Tests
        self.check_gauges(save=save, gauge_id=1)

        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Chile2010AdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
