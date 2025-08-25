#!/usr/bin/env python

r"""chile2010_adjoint regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from pathlib import Path
import os
import sys
import unittest
import shutil

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

from adjoint.test_chile2010_adjoint import Chile2010AdjointTest

try:
    CLAW = Path(os.environ['CLAW'])
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = CLAW / 'geoclaw' / 'scratch'

class Chile2010AdjointForwardTest(test.GeoClawRegressionTest):
    r"""Chile2010AdjointTest regression test for GeoClaw"""

    def runTest(self, save=False, indices=(2, 3)):
        r"""Test chile2010_adjoint example"""

        # Run adjoint problem
        try:
            adjoint_run = Chile2010AdjointTest()    
            adjoint_run.setUp()
            adjoint_run.runTest()
            
            # Copy output to local directory
            adjoint_output = Path(self.temp_path) / "_adjoint_output"

            if adjoint_output.exists():
                shutil.rmtree(adjoint_output)
            shutil.copytree(adjoint_run.temp_path, adjoint_output)
        finally:
            adjoint_run.tearDown()

        # Make topo and dtopo
        import maketopo
        maketopo.get_topo(path=Path(self.temp_path))
        maketopo.make_dtopo(path=Path(self.temp_path))

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 2
        self.rundata.clawdata.tfinal = 3600.

        self.rundata.amrdata.flag2refine_tol = 0.0005

        self.rundata.regiondata.regions = []
        # all 3 levels anywhere, based on flagging:
        self.rundata.regiondata.regions.append([1, 3, 0., 1e9, -220,0,-90,90])

        # earthquake source region - force refinement initially:
        # dtopo region, replacing minlevel in dtopofile specification:
        self.rundata.regiondata.regions.append([3, 3, 0., 10., -77,-67,-40,-30])
        # later times from original test:
        self.rundata.regiondata.regions.append([3, 3, 0., 200., -85,-70,-38,-25])

        self.rundata.flagregiondata.flagregions = []

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, -76, -36., 0., 1.e10])

        self.rundata.adjointdata.adjoint_outdir = adjoint_output.resolve()

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3))
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Chile2010AdjointForwardTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
