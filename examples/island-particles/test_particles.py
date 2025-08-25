#!/usr/bin/env python

r"""Particles regression test for GeoClaw

To create new regression data use
    `python test_particles.py True`
"""

import sys
import unittest
import shutil

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools


class ParticlesTest(test.GeoClawRegressionTest):
    r"""Particles regression test for GeoClaw"""

    def runTest(self, save=False, indices=(2, 3)):
        r"""Test particles example

        Note that this stub really only runs the code and performs no tests.

        """

        # Create topo and qinit
        import maketopo
        maketopo.maketopo(self.temp_path)
        maketopo.makeqinit(self.temp_path)

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 6.0

        self.rundata.amrdata.refinement_ratios_x = [2, 2]
        self.rundata.amrdata.refinement_ratios_y = [2, 2]
        self.rundata.amrdata.refinement_ratios_t = [2, 2]


        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gtype = {}
        self.rundata.gaugedata.gauges.append([1, 15., 20., 0., 1e10])
        self.rundata.gaugedata.gtype[1] = 'stationary'
        self.rundata.gaugedata.gauges.append([2, 15., 30., 0., 1e10])
        self.rundata.gaugedata.gtype[2] = 'lagrangian'

        self.write_rundata_objects()

        # Run code
        self.run_code()

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
