#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import unittest
import subprocess
import tempfile
import time

import numpy
import xarray as xr
import nose

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

build_failure_str = ("NetCDF topography test skipped due to failure to build "
                     "test program.")
class NetCDFBowlSloshTest(test.GeoClawRegressionTest):

    r"""NetCDF regression test for GeoClaw based on the bowl-slosh example"""

    def __init__(self, methodName="runTest"):

        super(NetCDFBowlSloshTest, self).__init__(methodName=methodName)


    def setUp(self):

        self.temp_path = tempfile.mkdtemp()

        self.stdout = open(os.path.join(self.temp_path, "run_output.txt"), "w")
        self.stdout.write("Output from Test %s\n" % self.__class__.__name__)
        # TODO - Should change this to use the time module's formatting
        # apparatus
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        self.stdout.write(date)
        self.stdout.write(("="*80 + "\n"))

        self.stderr = open(os.path.join(self.temp_path, "error_output.txt"), "w")
        self.stderr.write("Errors from Test %s\n" % self.__class__.__name__)
        self.stderr.write(date)
        self.stderr.write(("="*80 + "\n"))

        self.stdout.flush()
        self.stderr.flush()

        self.stdout.write("Paths:")
        self.stdout.write("  %s" % self.temp_path)
        self.stdout.write("  %s" % self.test_path)
        self.stdout.flush()

        # Make topography
        a = 1.
        h0 = 0.1
        z = lambda x,y: h0 * (x**2 + y**2) / a**2 - h0

        topo = topotools.Topography(topo_func=z)
        topo.x = numpy.linspace(-3.1, 3.1, 310)
        topo.y = numpy.linspace(-3.5,2.5, 300)

        try:
            import netCDF4
            this_path = os.path.join(self.temp_path, 'bowl.nc')

            # now mess with the order of the dimension IDs (lat, then lon)
            with netCDF4.Dataset(this_path,'w') as out:
                lat = out.createDimension('lat',len(topo.y))
                lon = out.createDimension('lon',len(topo.x))

                # create latitude first
                latitudes = out.createVariable('lat','f8',("lat",))
                longitudes = out.createVariable('lon','f8',("lon",))
                elevations = out.createVariable('elevation','f8',("lat","lon"))
                latitudes[:] = topo.y
                longitudes[:] = topo.x
                elevations[:] = topo.Z


        except ImportError:
            # Assume that NetCDF is not installed and move on
            raise nose.SkipTest(build_failure_str)

        except RuntimeError as e:
            print(e.message)
            raise nose.SkipTest("NetCDF topography test skipped due to " +
                                "runtime failure.")
        else:
            self.build_executable()

    def build_executable(self):
        try:
            self.stdout.write("Teting NetCDF output:\n")
            subprocess.check_call("cd %s ; make netcdf_test " % self.test_path,
                                                                stdout=self.stdout,
                                                                stderr=self.stderr,
                                                                shell=True)
        except subprocess.CalledProcessError:

            self.stdout.write(build_failure_str)
            raise nose.SkipTest(build_failure_str)

        else:
            # Force recompilation of topo_module to add NetCDF flags
            mod_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                    "shallow", "topo_module.mod")
            obj_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                    "shallow", "topo_module.o")
            if os.path.exists(mod_path):
                os.remove(mod_path)
            if os.path.exists(obj_path):
                os.remove(obj_path)

            super(NetCDFBowlSloshTest, self).build_executable()


    def runTest(self, save=False, indices=(2, 3)):
        r"""Test NetCDF topography support

        """

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3),
                          tolerance=1e-4)


    def tearDown(self):
        r"""Tear down test infrastructure.

        This version removes source that may have been compiled with NetCDF
        turned on so that it does not conflict with subsequent tests.

        """

        # Force recompilation of topo_module to add NetCDF flags
        mod_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                "shallow", "topo_module.mod")
        obj_path = os.path.join(os.environ["CLAW"], "geoclaw", "src", "2d",
                                "shallow", "topo_module.o")
        if os.path.exists(mod_path):
            os.remove(mod_path)
        if os.path.exists(obj_path):
            os.remove(obj_path)

        super(NetCDFBowlSloshTest, self).tearDown()



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = NetCDFBowlSloshTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
