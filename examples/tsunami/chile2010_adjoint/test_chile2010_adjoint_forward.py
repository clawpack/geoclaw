#!/usr/bin/env python
"""Pytest regression test for the Chile 2010 forward-adjoint example."""

from pathlib import Path

import pytest

import clawpack.clawutil.test as clawtest
import clawpack.clawutil.util as util
import clawpack.geoclaw.test as test

@pytest.mark.regression
@pytest.mark.tsunami
@pytest.mark.remote
@pytest.mark.xfail(reason="Slight gauge mismatch.")
def test_chile2010_adjoint_forward(tmp_path: Path, save: bool) -> None:
    """Regression test for the Chile 2010 forward run using adjoint output."""
    example_path = Path(__file__).parent
    adjoint_path = example_path / "adjoint"
    adjoint_output = tmp_path / "_adjoint_output"

    runner = test.GeoClawTestRunner(tmp_path, test_path=example_path)

    # Run adjoint test to create output for forward test to use
    adjoint_setup = util.fullpath_import(adjoint_path / "test_chile2010_adjoint.py")
    clawtest.run_example_for_test(
        test.GeoClawTestRunner,
        adjoint_output,
        adjoint_path,
        configure_runner=adjoint_setup.set_adjoint_data,
    )

    # Create topo and qinit inputs
    maketopo_module = util.fullpath_import(runner.test_path / "maketopo.py")
    maketopo_module.get_topo(tmp_path, verbose=True)
    maketopo_module.make_dtopo(tmp_path)

    runner.set_data()

    runner.rundata.clawdata.num_output_times = 2
    runner.rundata.clawdata.tfinal = 3600.0

    runner.rundata.amrdata.flag2refine_tol = 0.0005

    runner.rundata.regiondata.regions = []
    
    runner.rundata.regiondata.regions.append([1, 3, 0., 1e9, -220,0,-90,90])
    runner.rundata.regiondata.regions.append([3, 3, 0., 10., -77,-67,-40,-30])
    runner.rundata.regiondata.regions.append([3, 3, 0., 200., -85,-70,-38,-25])

    runner.rundata.flagregiondata.flagregions = []

    runner.rundata.gaugedata.gauges = [[1, -76, -36., 0., 1.e10]]

    runner.rundata.adjointdata.adjoint_outdir = str(adjoint_output)
    runner.write_data()

    runner.build_executable()
    runner.run_code()

    runner.check_gauge(save=save, gauge_id=1)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
