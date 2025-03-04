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


class ParticlesTest(test.GeoClawRegressionTest):

    r"""Particles regression test for GeoClaw"""

    def setUp(self):

        super(ParticlesTest, self).setUp()
        start_dir = os.getcwd()

        # Make topography

        shutil.copy(os.path.join(self.test_path, "maketopo.py"),
                                 self.temp_path)
        os.chdir(self.temp_path)
        os.system('python maketopo.py')
        os.chdir(start_dir)


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test particles example

        Note that this stub really only runs the code and performs no tests.

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        import clawpack.pyclaw.gauges as gauges
        gauge = gauges.GaugeSolution(1, path=self.temp_path)
        print('+++ Gauge 1:\n', gauge.q)
        gauge = gauges.GaugeSolution(2, path=self.temp_path)
        print('+++ Gauge 2:\n', gauge.q)

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(1, 2))
        self.check_gauges(save=save, gauge_id=2, indices=(1, 2))
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = ParticlesTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
