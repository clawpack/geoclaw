r"""
Defines the GeoClaw Clawpack Test Runner class for running PyTest based
regression tests in GeoClaw.

Refer to the documentation for PyTest to manage output and reporting.
"""

from pathlib import Path
from typing import Optional
import os
import numpy as np

import clawpack.clawutil.test as test
import clawpack.geoclaw.fgmax_tools as fgmax_tools

# Set environment variable to avoid warning about missing CLAW variable
if "CLAW" in os.environ:
    CLAW = os.environ["CLAW"]
else:
    raise ValueError("Need to set CLAW environment variable.")

class GeoClawTestRunner(test.ClawpackTestRunner):
    r"""Class for running GeoClaw regression tests.

    """

    def __init__(self, path: Path, test_path: Optional[Path]=None):
        super(GeoClawTestRunner, self).__init__(path, test_path=test_path)
        self.executable_name = "xgeoclaw"

    def check_fgmax(self, fgno: int=1, save: bool=False, tolerance: float=1e-14):
        r"""Basic test to assert fgmax equality
        Currently just records sum of fg.h and of fg.s.

        :TODO: Add documentation
        """
        fg = fgmax_tools.FGmaxGrid()
        fname = self.temp_path / 'fgmax_grids.data'
        fg.read_fgmax_grids_data(fgno, fname)
        fg.read_output(outdir=self.temp_path)

        data_sum = np.array([fg.h.sum(), fg.s.sum()])

        # Get (and save) regression comparison data
        regression_data_file = (self.test_path/ "regression_data" /
                "regression_data_fgmax.txt")
        if save:
            np.savetxt(regression_data_file, data_sum)
        regression_sum = np.loadtxt(regression_data_file)

        # Compare data
        assert np.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)


# Useful for running tests that need more than one example to be run, e.g. an
# adjoint test that needs to run the forward problem first.
def run_example_for_test(
    runner_cls,
    output_path: Path,
    test_path: Path,
    *,
    setrun_path: Path | None = None,
    configure_runner=None,
    build_kwargs: dict | None = None,
):
    """
    Build and run one example in a specified output directory.

    Parameters
    ----------
    runner_cls
        Runner class, e.g. AMRClawTestRunner.
    output_path
        Temporary directory for data files, executable, and output.
    test_path
        Path to the example directory.
    setrun_path
        Optional explicit path to setrun.py.
    configure
        Optional callback taking the runner after set_data() and before
        write_data().
    build_kwargs
        Optional keyword arguments passed to build_executable().
    """
    output_path.mkdir(parents=True, exist_ok=True)

    runner = runner_cls(output_path, test_path=test_path)
    runner.set_data(setrun_path=setrun_path)

    if configure_runner is not None:
        configure_runner(runner)

    runner.write_data()
    runner.build_executable(**(build_kwargs or {}))
    runner.run_code()
    return runner
