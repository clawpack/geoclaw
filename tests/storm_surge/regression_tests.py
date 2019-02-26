#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

from __future__ import absolute_import
import sys
import os
import unittest
import gzip
import nose
import datetime

try:
    # For Python 3.0 and later
    from urllib.error import URLError
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import URLError

import numpy

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools
import clawpack.geoclaw.surge.storm as storm

class IkeTest(test.GeoClawRegressionTest):

    r"""Hurricane Ike regression test"""

    def setUp(self):

        super(IkeTest, self).setUp()

        # Fetch storm data
        remote_url = "http://ftp.nhc.noaa.gov/atcf/archive/2008/bal092008.dat.gz"
        atcf_path = os.path.join(self.temp_path, "bal092008.dat")
        storm_path = os.path.join(os.path.dirname(self.temp_path), 'ike.storm')
        
        try:
            path = self.get_remote_file(remote_url, unpack=False)
        except URLError:
            raise nose.SkipTest("Could not fetch remote file, skipping test.")

        # Need to additionally deal with the fact the file is gzipped
        with gzip.GzipFile(path, 'r') as gzip_file:
            file_content = gzip_file.read()
        
        with open(atcf_path, 'wb') as out_file:
            out_file.write(file_content)
        
        ike = storm.Storm(path=atcf_path, file_format="ATCF")
        ike.time_offset = datetime.datetime(2008, 9, 13, 7)
        ike.write(storm_path, file_format="geoclaw")

        

        # # Convert ATCF data to GeoClaw format
        # clawutil.data.get_remote_file(
        #                "http://ftp.nhc.noaa.gov/atcf/archive/2008/bal092008.dat.gz")
        # atcf_path = os.path.join(scratch_dir, "bal092008.dat")
        # # Note that the get_remote_file function does not support gzip files which
        # # are not also tar files.  The following code handles this
        # with gzip.open(".".join((atcf_path, 'gz')), 'rb') as atcf_file,    \
        #         open(atcf_path, 'w') as atcf_unzipped_file:
        #     atcf_unzipped_file.write(atcf_file.read().decode('ascii'))

        # # Uncomment/comment out to use the old version of the Ike storm file
        # ike = Storm(path=atcf_path, file_format="ATCF")

        # # Calculate landfall time - Need to specify as the file above does not
        # # include this info (9/13/2008 ~ 7 UTC)
        # ike.time_offset = datetime.datetime(2008, 9, 13, 7)

        # ike.write(data.storm_file, file_format='geoclaw')

        # Download file
        #self.get_remote_file(
        #   "http://www.columbia.edu/~ktm2132/bathy/gulf_caribbean.tt3.tar.bz2")

        # Create synthetic bathymetry - needs more work
        topo = clawpack.geoclaw.topotools.Topography()
        topo.x = numpy.linspace(-100, -69, 125)
        topo.y = numpy.linspace(7.0, 33.0, 105)
        topo.Z = 25.0 * ((topo.X + 84.5)**2 + (topo.Y - 20.0)**2) - 4000.0
        topo.write(os.path.join(self.temp_path, 'gulf_caribbean.tt3'), \
                topo_type=2, Z_format="%22.15e")


    def runTest(self, save=False, indices=(2, 3)):
        r"""Storm Surge Regression Test

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Passed to *check_gauges*.  Default is
           *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison and passed to *check_gauges*.  Defaults to *(2, 3)*.

        """

        # Write out data files
        self.load_rundata()
        self.rundata.surge_data.storm_file = os.path.join(
                                  os.path.dirname(self.temp_path), 'ike.storm')

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=indices)
        self.check_gauges(save=save, gauge_id=2, indices=indices)
        self.check_gauges(save=save, gauge_id=3, indices=indices)
        self.check_gauges(save=save, gauge_id=4, indices=indices)

        # If we have gotten here then we do not need to copy the run results
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = IkeTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
    