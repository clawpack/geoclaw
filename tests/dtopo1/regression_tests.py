#!/usr/bin/env python

r"""Regression tests.  Execute via:
    python regression_tests.py
to test, or
    python regression_tests.py True
to create new regression data for archiving.
"""

from __future__ import absolute_import
import os
import sys
import unittest
import shutil

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools

class DTopoTests(test.GeoClawRegressionTest):
    
    def setUp(self):

        super(DTopoTests, self).setUp()

        # Make topography
        h0 = 1000.0
        topo_func = lambda x,y: -h0*(1 + 0.5 * numpy.cos(x - y))
        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-10.0, 10.0, 201)
        topo.y = numpy.linspace(-10.0, 10.0, 201)
        topo.write(os.path.join(self.temp_path, "topo1.topotype2"), \
                topo_type=2, Z_format="%22.15e")

        h0 = 1000.0
        topo_func = lambda x,y: -h0*(1. + numpy.exp(x+y))
        topo = topotools.Topography(topo_func=topo_func)
        topo.topo_type = 2
        topo.x = numpy.linspace(-0.5, -0.3, 21)
        topo.y = numpy.linspace(-0.1, 0.4, 51)
        topo.write(os.path.join(self.temp_path, "topo2.topotype2"), \
                topo_type=2, Z_format="%22.15e")

        # Make dtopography
        subfault_path = os.path.join(self.test_path, "dtopo1.csv")
        input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                coordinate_specification="top center")
        fault.rupture_type = 'dynamic'
        times = numpy.linspace(0.0, 1.0, 25)
        x = numpy.linspace(-0.4,0.6,151)
        y = numpy.linspace(-0.4,0.4,121)
        dtopo = fault.create_dtopography(x,y,times=times)
        dtopo.write(os.path.join(self.temp_path, "dtopo1.tt3"), dtopo_type=3)

        subfault_path = os.path.join(self.test_path, "dtopo2.csv")
        input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                    coordinate_specification="top center")
        fault.rupture_type = 'dynamic'
        times = numpy.linspace(0.5, 1.2, 25)
        x = numpy.linspace(-0.9,0.1,201)
        y = numpy.linspace(-0.4,0.4,161)
        dtopo = fault.create_dtopography(x,y,times=times)
        dtopo.write(os.path.join(self.temp_path, "dtopo2.tt3"), dtopo_type=3)

        # copy existing file:
        shutil.copy(os.path.join(self.test_path, "dtopo3.tt1"),
                                 self.temp_path)


    def check_gauges(self, save=False, gauge_id=1, regression_gauge_id=1, 
                     indices=[0], rtol=1e-14, atol=1e-8, tolerance=None):

        r"""Basic test to assert gauge equality

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
         - *gauge_id* (int) - The gauge to test.
         - *regression_gauge_id* (int) - The gauge to test against.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison.  Defaults to *(0)*.
         - *rtol* (float) - Relative tolerance used in the comparison, default 
           is *1e-14*.  Note that the old *tolerance* input is now synonymous 
           with this parameter.
         - *atol* (float) - Absolute tolerance used in the comparison, default
           is *1e-08*.
        """

        import clawpack.pyclaw.gauges as gauges
        import clawpack.clawutil.claw_git_status as claw_git_status

        if isinstance(tolerance, float):
            rtol = tolerance

        if not(isinstance(indices, tuple) or isinstance(indices, list)):
            indices = tuple(indices)

        # Get gauge data
        gauge = gauges.GaugeSolution(gauge_id, path=self.temp_path)

        # Get regression comparison data
        regression_data_path = os.path.join(self.test_path, "regression_data")
        if save:
            gauge_file_name = "gauge%s.txt" % str(regression_gauge_id).zfill(5)
            shutil.copy(os.path.join(self.temp_path, gauge_file_name), 
                                                           regression_data_path)
            claw_git_status.make_git_status_file(outdir=regression_data_path)

        regression_gauge = gauges.GaugeSolution(regression_gauge_id,
                                                path=regression_data_path)

        # Compare data
        try:
            for n in indices:
                numpy.testing.assert_allclose(gauge.q[n, :],
                                              regression_gauge.q[n, :], 
                                              rtol=rtol, atol=atol, 
                                              verbose=False)
        except AssertionError as e:
            err_msg = "\n".join((e.args[0], 
                                "Gauge Match Failed for gauge = %s" % gauge_id))
            err_msg = "\n".join((err_msg, "  failures in fields:"))
            failure_indices = []
            for n in indices:
                if ~numpy.allclose(gauge.q[n, :], regression_gauge.q[n, :], 
                                                          rtol=rtol, atol=atol):
                    failure_indices.append(str(n))
            index_str = ", ".join(failure_indices)
            raise AssertionError(" ".join((err_msg, index_str)))



    def runTest(self, save=False, indices=(2, 3)):
        r"""DTopography basic regression test

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Passed to *check_gauges*.  Default is
           *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison and passed to *check_gauges*.  Defaults to *(2, 3)*.

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, regression_gauge_id=1,
                          indices=(2, 3))
        print('gauge 1 ascii agrees')
        self.check_gauges(save=save, gauge_id=2, regression_gauge_id=1,
                          indices=(2, 3), rtol=1e-6, atol=1e-6)
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

