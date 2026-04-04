#!/usr/bin/env python
"""Pytest helper regression test for the Chile 2010 adjoint stage."""

from pathlib import Path

import pytest

from clawpack.clawutil.util import fullpath_import
import clawpack.geoclaw.test as test

def set_adjoint_data(runner):
    """Set adjoint data values in the provided runner.
    
    Need to do this in a function so the forward problem can run 
    this to generate the necessary files.
    """
    # Create topo and qinit inputs
    maketopo_module = fullpath_import(runner.test_path / "maketopo.py")
    maketopo_module.get_topo(runner.temp_path)
    maketopo_module.makeqinit(runner.temp_path, center=(-76., -36.))
    
    runner.rundata.clawdata.num_cells[0] = 240
    runner.rundata.clawdata.num_cells[1] = 240
    runner.rundata.clawdata.num_output_times = 4
    runner.rundata.clawdata.tfinal = 3600.0
    runner.rundata.amrdata.amr_levels_max = 1
    runner.rundata.amrdata.flag2refine = False
    runner.rundata.gaugedata.gauges = [[1, -76, -36., 0., 1.e10]]
    runner.rundata.flagregiondata.flagregions = []


def run_chile2010_adjoint_stage(output_dir: Path):
    """Run the Chile 2010 adjoint stage and return the configured runner."""
    return clawtest.run_example_for_test(
        test.GeoClawTestRunner,
        output_dir,
        Path(__file__).parent,
        configure_runner=set_adjoint_data,
    )


@pytest.mark.regression
@pytest.mark.tsunami
@pytest.mark.remote
@pytest.mark.adjoint
def test_chile2010_adjoint(tmp_path: Path, save: bool) -> None:
    """Smoke test for the Chile 2010 adjoint stage."""
    runner = run_chile2010_adjoint_stage(tmp_path)
    runner.check_gauge(gauge_id=1, save=save)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
