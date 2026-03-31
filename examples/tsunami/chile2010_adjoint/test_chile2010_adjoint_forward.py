#!/usr/bin/env python
"""Pytest regression test for the Chile 2010 forward-adjoint example."""

from pathlib import Path
import importlib.util
import random
import string
import sys
import pytest

import clawpack.geoclaw.test as test

def _load_adjoint_test_module(adjoint_path: Path):
    """Load the example-local test module under a unique module name."""
    test_path = adjoint_path / "test_chile2010_adjoint.py"
    mod_name = "_".join(
        (
            "test_chile2010_adjoint",
            "".join(random.choices(string.ascii_letters + string.digits, k=32)),
        )
    )
    spec = importlib.util.spec_from_file_location(mod_name, test_path)
    test_module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = test_module
    spec.loader.exec_module(test_module)
    return test_module

def _load_maketopo_module(example_dir: Path):
    """Load the example-local maketopo module under a unique module name."""
    maketopo_path = example_dir / "maketopo.py"
    mod_name = "_".join(
        (
            "maketopo",
            "".join(random.choices(string.ascii_letters + string.digits, k=32)),
        )
    )
    spec = importlib.util.spec_from_file_location(mod_name, maketopo_path)
    maketopo_module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = maketopo_module
    spec.loader.exec_module(maketopo_module)
    return maketopo_module

@pytest.mark.regression
@pytest.mark.tsunami
@pytest.mark.remote
def test_chile2010_adjoint_forward(tmp_path: Path, save: bool) -> None:
    """Regression test for the Chile 2010 forward run using adjoint output."""
    example_path = Path(__file__).parent
    adjoint_path = example_path / "adjoint"
    adjoint_output = tmp_path / "_adjoint_output"

    runner = test.GeoClawTestRunner(tmp_path, test_path=example_path)

    # Run adjoint test to create output for forward test to use
    adjoint_setup = _load_adjoint_test_module(adjoint_path)
    test.run_example_for_test(test.GeoClawTestRunner, adjoint_output, adjoint_path, 
                              configure_runner=adjoint_setup.set_adjoint_data)

    # Create topo and qinit inputs
    maketopo_module = _load_maketopo_module(runner.test_path)
    maketopo_module.get_topo(tmp_path)
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

    runner.rundata.adjointdata.adjoint_outdir = adjoint_output
    runner.write_data()

    runner.build_executable()
    runner.run_code()

    runner.check_gauge(save=save, gauge_id=1)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
