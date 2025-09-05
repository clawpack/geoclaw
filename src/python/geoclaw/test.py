r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

from pathlib import Path
import os
import shutil
import glob

import numpy as np

import clawpack.clawutil.test
import clawpack.pyclaw.util

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in [CLAW / "amrclaw" / "src" / "2d",
                 CLAW / "geoclaw" / "src" / "2d" / "shallow",
                 CLAW / "geoclaw" / "src" / "2d" / "shallow" / "multilayer",
                 CLAW / "geoclaw" / "src" / "2d" / "shallow" / "surge"]:
    for path in lib_path.glob("*.o"):
        path.unlink()
    for path in lib_path.glob("*.mod"):
        path.unlink()


class GeoClawRegressionTest(clawpack.clawutil.test.ClawpackRegressionTest):

    r"""Base GeoClaw regression test setup derived from ClawpackRegressionTest

    """

    __doc__ += clawpack.pyclaw.util.add_parent_doc(
                                  clawpack.clawutil.test.ClawpackRegressionTest)


    def build_executable(self, executable_name="xgeoclaw"):
        r"""Build executable by running `make .exe` in test directory.

        Moves the resulting executable to the temporary directory.


        """

        super(GeoClawRegressionTest, self).build_executable(
                                                executable_name=executable_name)


    def check_fgmax(self, fgno=1, save=False, **kwargs):
        r"""Basic test to assert fgmax equality
        Currently just records sum of fg.h and of fg.s.

        :Input:
         - *fgno* (int) - Which fixed grid to compare
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
         - *kwargs* (dict) Dictionary of key-word arguments passed to 
           *numpy.assert_allclose* such as *atol* and *rtol*.
        """

        from clawpack.geoclaw import fgmax_tools

        fg = fgmax_tools.FGmaxGrid()
        fname = Path(self.temp_path) / 'fgmax_grids.data'
        fg.read_fgmax_grids_data(fgno, fname)
        fg.read_output(outdir=self.temp_path)

        data_sum = np.array([fg.h.sum(), fg.s.sum()])

        # Get (and save) regression comparison data
        # :TODO: allow more than one fgmax file?  Maybe compare entire file
        regression_data_file = (Path(self.test_path) / "regression_data"
                                                 / "regression_data_fgmax.txt")
        if save:
            np.savetxt(regression_data_file, data_sum)
        regression_sum = np.loadtxt(regression_data_file)

        # Compare data
        kwargs.setdefault('atol', 1e-14)
        np.testing.assert_allclose(data_sum, regression_sum, **kwargs)
