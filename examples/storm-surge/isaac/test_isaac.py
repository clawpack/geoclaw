#!/usr/bin/env python
"""
Regression tests for Isaac storm surge
Primarily tests all of the input formats GeoClaw handles
"""

from pathlib import Path
import os
import gzip
import unittest
import pytest

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.clawutil as clawutil
from clawpack.geoclaw.surge.storm import Storm

days2seconds = lambda t: t * 60.0**2 * 24.0

# Look in scratch directory for storm file
scratch_dir = Path(os.environ['CLAW']) / 'geoclaw' / 'scratch'

# os.path.join(os.environ["CLAW"], 'geoclaw', 'scratch')

@pytest.mark.parametrize("data_file_format", 
                         ['geoclaw', 'owi_ascii', 'owi_netcdf'])
def test_isaac_formats(data_file_format):
    test = IsaacStormSurgeRun(file_format=data_file_format)
    try:
        test.setUp()
        test.runTest(save=True)
    finally:
        test.tearDown()

# :TODO: Restructure this to use pytest.mark.parameterize


class IsaacStormSurgeRun(test.GeoClawRegressionTest):
    r"""Regression test for Hurrican Isaac storm surge"""

    def __init__(self, methodName="runTest", file_format="geoclaw"):

        super(IsaacStormSurgeRun, self).__init__(methodName=methodName)

        if file_format.lower() not in ["geoclaw", "owi_ascii", "owi_netcdf"]:
            raise ValueError(f"Invalid test file fromat {self.file_format}.")
        self.file_format = file_format

    def setUp(self):


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.t0 = days2seconds(-3)
        self.rundata.clawdata.tfinal = days2seconds(1)

        surge_data = self.rundata.surge_data

        # Reload all of storm specific stuff
        # The fetching of the remote file should have already happened but we 
        # will check anyway
        clawutil.data.get_remote_file(
                   "http://ftp.nhc.noaa.gov/atcf/archive/2012/bal092012.dat.gz")
        atcf_path = scratch_dir / "bal092012.dat"
        # with gzip.open(".".join((atcf_path, 'gz')), 'rb') as atcf_file,    \
        with gzip.open(atcf_path.with_suffix(".dat.gz"), 'rb') as atcf_file,    \
                open(atcf_path, 'w') as atcf_unzipped_file:
            atcf_unzipped_file.write(atcf_file.read().decode('ascii'))

        # Uncomment/comment out to use the old version of the Ike storm file
        isaac = Storm(path=atcf_path, file_format="ATCF")

        # Calculate landfall time - Need to specify as the file above does not
        # include this info (~2345 UTC - 6:45 p.m. CDT - on August 28)
        # isaac.time_offset = datetime.datetime(2012, 8, 29, 0)
        isaac.time_offset = np.datetime64("2012-08-29")

        # Temporarily use generated files for OWI storms
        if self.file_format.lower == "geoclaw":
            surge_data.storm_specification_type = "holland80"
        elif self.file_format.lower == "owi_ascii":
            surge_data.storm_specification_type = 'OWI'
            isaac.data_file_format = 'NWS12'
            isaac.file_paths = [Path.cwd() / "isaac.PRE",
                                Path.cwd() / "isaac.WIN"]
            isaac.write(data.storm_file, file_format='OWI')
        elif self.file_format.lower == "owi_netcdf":
            surge_data.storm_specification_type = 'OWI'
            isaac.data_file_format = "NWS13"
            isaac.file_paths = [self.test_path / "isaac.nc"]
            isaac.write(data.storm_file, file_format='OWI')

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        # self.check_gauges(save=save, gauge_id=1)
        # self.check_gauges(save=save, gauge_id=2)

        self.success = True


if __name__=="__main__":
    # :TODO: Add each format here to loop over
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = IsaacStormSurgeRun()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()