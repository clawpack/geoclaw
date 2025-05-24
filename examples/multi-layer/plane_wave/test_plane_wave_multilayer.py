#!/usr/bin/env python

r"""Multilayer Shallow Water Test Case

To create new regression data use
    `python test_plane_wave_multilayer.py True`
"""

from pathlib import Path
import sys
import unittest

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

class PlaneWaveMultilayerTest(test.GeoClawRegressionTest):
    r"""Multilayer plane-wave regression test for GeoClaw

    initial condition angle = np.pi / 4.0
    bathy_angle = np.pi / 8.0

    """

    def runTest(self, save=False):
        r"""Test multi-layer basic plane-waves."""

        # Load and write data, change init-condition's starting angle
        self.load_rundata()

        self.rundata.clawdata.lower[0] = -1.0
        self.rundata.clawdata.upper[0] = 2.0

        self.rundata.clawdata.lower[1] = -1.0
        self.rundata.clawdata.upper[1] = 2.0

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 1.0

        self.rundata.amrdata.refinement_ratios_x = [2, 6]
        self.rundata.amrdata.refinement_ratios_y = [2, 6]
        self.rundata.amrdata.refinement_ratios_t = [2, 6]

        self.rundata.geo_data.rho = [0.9, 1.0]
    
        self.rundata.topo_data.topofiles = []
        topo_path = Path(self.temp_path) / 'jump_topo.topotype2'
        self.rundata.topo_data.topofiles.append([2, topo_path])
        

        self.rundata.multilayer_data.wave_tolerance = [0.1, 0.2]

        self.rundata.qinit_data.angle = np.pi / 4.0

        self.write_rundata_objects()

        # Create topography
        import setrun
        setrun.write_topo_file(self.rundata, topo_path, location=0.15, 
                                                        angle=np.pi / 8.0, 
                                                        left=-1.0, 
                                                        right=-0.2)

        # Run code and check surface heights
        self.run_code()
        self.check_gauges(save=save, gauge_id=0, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=1, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=2, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=3, indices=(6, 7), atol=1e-5)
        self.check_gauges(save=save, gauge_id=4, indices=(6, 7), atol=1e-5)

        # If we have gotten here then we do not need to copy the run results
        self.success = True


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = PlaneWaveMultilayerTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
