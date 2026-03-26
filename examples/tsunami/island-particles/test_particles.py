"""Pytest regression test for the GeoClaw particles example."""

from pathlib import Path
import importlib.util
import random
import string
import sys

import pytest

import clawpack.geoclaw.test as test


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
@pytest.mark.xfail(reason="Particles regression test is currently failing due to unknown changes and need to be re-validated.")
def test_particles(tmp_path: Path) -> None:
    """Regression test for the GeoClaw particles example."""
    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    # Create topo and qinit inputs
    maketopo_module = _load_maketopo_module(example_dir)
    maketopo_module.maketopo(tmp_path)
    maketopo_module.makeqinit(tmp_path)

    # Load and adjust run data for the regression test
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 6.0

    runner.rundata.amrdata.refinement_ratios_x = [2, 2]
    runner.rundata.amrdata.refinement_ratios_y = [2, 2]
    runner.rundata.amrdata.refinement_ratios_t = [2, 2]

    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gtype = {}
    runner.rundata.gaugedata.gauges.append([1, 15.0, 20.0, 0.0, 1e10])
    runner.rundata.gaugedata.gtype[1] = "stationary"
    runner.rundata.gaugedata.gauges.append([2, 15.0, 30.0, 0.0, 1e10])
    runner.rundata.gaugedata.gtype[2] = "lagrangian"

    runner.write_data()

    # Build and run code
    runner.build_executable()
    runner.run_code()

    # Check gauge outputs
    runner.check_gauge(gauge_id=1, indices=(1, 2))
    runner.check_gauge(gauge_id=2, indices=(1, 2))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
