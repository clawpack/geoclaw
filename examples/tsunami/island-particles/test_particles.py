"""Pytest regression test for the GeoClaw particles example."""

from pathlib import Path
import pytest

import clawpack.geoclaw.test as test
from clawpack.clawutil.util import fullpath_import


@pytest.mark.regression
@pytest.mark.xfail(reason="Particles regression test is currently failing due to unknown changes and need to be re-validated.")
def test_particles(tmp_path: Path, save: bool) -> None:
    """Regression test for the GeoClaw particles example."""
    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    # Create topo and qinit inputs
    maketopo_module = fullpath_import(example_dir / "maketopo.py")
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
    runner.check_gauge(gauge_id=1, indices=(1, 2), save=save)
    runner.check_gauge(gauge_id=2, indices=(1, 2), save=save)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
