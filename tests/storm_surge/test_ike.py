#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

import sys
import os
import unittest
import pytest
import gzip
import urllib.error
import pytest

import numpy

import clawpack.geoclaw.test
import clawpack.geoclaw.topotools
import clawpack.clawutil
from clawpack.geoclaw.surge import storm

class IkeTest(clawpack.geoclaw.test.GeoClawRegressionTest):

    r"""Hurricane Ike regression test"""

    def setUp(self):

        super(IkeTest, self).setUp()

        # Download storm data
        remote_url = "http://ftp.nhc.noaa.gov/atcf/archive/2008/bal092008.dat.gz"
        try:
            path = self.get_remote_file(remote_url, unpack=False)
        except urllib.error.URLError as e:
            pytest.skip("Could not fetch remote file, skipping test.")
            raise e
        
        storm_path = os.path.join(os.path.dirname(path), 'ike.storm')

        # Need to additionally deal with the fact the file is gzipped
        with gzip.GzipFile(path, 'r') as gzip_file:
            file_content = gzip_file.read()
        
        with open(storm_path+'.atcf', 'wb') as out_file:
            out_file.write(file_content)
            
        # now convert to geoclaw format
        ike_storm = storm.Storm(storm_path+'.atcf', file_format='ATCF', verbose=True)
        ike_storm.write(storm_path)
        
        # Download file
        # clawpack.clawutil.data.get_remote_file(
        #    "https://depts.washington.edu/clawpack/geoclaw/topo/gulf_caribbean.tt3.tar.bz2")

        # Create synthetic bathymetry - needs more work
        topo = clawpack.geoclaw.topotools.Topography()
        topo.x = numpy.linspace(-100, -69, 125)
        topo.y = numpy.linspace(7.0, 33.0, 105)
        topo.Z = 25.0 * ((topo.X + 84.5)**2 + (topo.Y - 20.0)**2) - 4000.0
        topo.write(os.path.join(self.temp_path, 'gulf_caribbean.tt3'), \
                topo_type=2, Z_format="%22.15e")


    def runTest(self, save=False, indices=range(4)):
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
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=indices)

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