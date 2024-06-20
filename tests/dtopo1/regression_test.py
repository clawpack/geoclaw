#!/usr/bin/env python

r"""Particles regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from __future__ import absolute_import
import os
import sys
import unittest
import shutil

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools


class DTopoTests(test.GeoClawRegressionTest):

    r"""dtopo1 regression test for GeoClaw"""

    def setUp(self):

        super(DTopoTests, self).setUp()
        start_dir = os.getcwd()

        # Make topography and dtopo

        shutil.copy(os.path.join(self.test_path, "maketopo.py"),
                                 self.temp_path)
        shutil.copy(os.path.join(self.test_path, "dtopo1.csv"),
                                 self.temp_path)               
        shutil.copy(os.path.join(self.test_path, "dtopo2.csv"),
                                 self.temp_path)
        shutil.copy(os.path.join(self.test_path, "dtopo3.tt1"),
                                 self.temp_path)
        os.chdir(self.temp_path)
        os.system('python maketopo.py')
        os.chdir(start_dir)


    def runTest(self, save=False, indices=(2, 3)):
        r"""DTopography basic regression test
        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()
        
        import clawpack.pyclaw.gauges as gauges
        gauge = gauges.GaugeSolution(1, path=self.temp_path)
        #print('+++ Gauge 1:\n', gauge.q)
        gauge = gauges.GaugeSolution(2, path=self.temp_path)
        #print('+++ Gauge 2:\n', gauge.q)
        
        # Perform tests
        self.check_gauges(save=save, gauge_id=1, 
                          indices=(2, 3))
        print('gauge 1 ascii agrees')
        self.check_gauges(save=save, gauge_id=2,
                          indices=(2, 3))
        print('gauge 2 binary agrees')

        # If we have gotten here then we do not need to copy the run results
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = DTopoTests()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
