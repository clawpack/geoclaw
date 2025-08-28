#!/usr/bin/env python

"""Regression test for GeoClaw's storm surge functionality"""

from pathlib import Path
import sys
import unittest
import urllib
import gzip

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.clawutil.data
import clawpack.geoclaw.topotools
from clawpack.geoclaw.surge import storm

class IkeStormSurgeTest(test.GeoClawRegressionTest):
    r"""Hurricane Ike regression test"""

    # TODO: Need to enable wind and pressure field comparisons (4, 5, 6)), needs
    #       a PR for fixing aux gauge recording
    def runTest(self, save=False, indices=(0, 1, 2, 3, 4, 5, 6)):
        r"""Storm Surge Regression Test

        :Input:
         - *save* (bool) - If *True* will save the output from this test to
           the file *regresion_data.txt*.  Passed to *check_gauges*.  Default is
           *False*.
         - *indices* (tuple) - Contains indices to compare in the gague
           comparison and passed to *check_gauges*.  Defaults to *(2, 3)*.

        """

        # Download and write out storm
        remote_url = "http://ftp.nhc.noaa.gov/atcf/archive/2008/bal092008.dat.gz"
        try:
            atcf_path = Path(clawpack.clawutil.data.get_remote_file(remote_url))
        except urllib.error.URLError as e:
            pytest.skip(f"Could not fetch remote file {remote_url}.")
        storm_path = Path(self.temp_path) / 'ike.storm'
        ike = storm.Storm(atcf_path, file_format='ATCF')
        ike.time_offset = np.datetime64("2008-09-13T07")
        ike.write(storm_path, file_format="geoclaw")


        # Create synthetic bathymetry
        topo = clawpack.geoclaw.topotools.Topography()
        topo.x = np.linspace(-100, -69, 125)
        topo.y = np.linspace(7.0, 33.0, 105)
        topo.Z = 25.0 * ((topo.X + 84.5)**2 + (topo.Y - 20.0)**2) - 4000.0
        topo.write(Path(self.temp_path) / 'gulf_caribbean.tt3', topo_type=2,
                    Z_format="%22.15e")

        # Write out rundata
        self.load_rundata()

        self.rundata.clawdata.num_cells = [29, 24]
        self.rundata.clawdata.num_output_times = 0

        self.rundata.surge_data.storm_file = storm_path

        self.rundata.gaugedata.gauges = [[1, -90., 25.,
                        self.rundata.clawdata.t0, self.rundata.clawdata.tfinal]]

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=indices)

        # If we have gotten here then we do not need to copy the run results
        self.success = True


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = IkeStormSurgeTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
