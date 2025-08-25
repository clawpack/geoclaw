#!/usr/bin/env python

r"""Shallow Water Travelling Vortex Test

To create new regression data use
    `python test_vortex.py True`
"""

from pathlib import Path
import sys
import unittest

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.fgmax_tools as fgmax_tools

class VortexTest(test.GeoClawRegressionTest):
    r"""Traveling vortex regression test"""

    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        # self.rundata.clawdata.num_output_times = 1
        # self.rundata.clawdata.tfinal = 0.5

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        for i in range(5):
            self.check_gauges(save=save, gauge_id=i, indices=(2, 3))
        
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = VortexTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()