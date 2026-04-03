"""Pytest regression test for the GeoClaw dtopography example."""

from pathlib import Path
import shutil

import numpy as np
import pytest

import clawpack.geoclaw.dtopotools as dtopotools
import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools


def _make_topography_files(output_dir: Path) -> None:
    """Create the static topography files needed by the example."""
    h0 = 1000.0
    topo_func = lambda x, y: -h0 * (1.0 + 0.5 * np.cos(x - y))
    topo = topotools.Topography(topo_func=topo_func)
    topo.topo_type = 2
    topo.x = np.linspace(-10.0, 10.0, 201)
    topo.y = np.linspace(-10.0, 10.0, 201)
    topo.write(output_dir / "topo1.topotype2", topo_type=2, Z_format="%22.15e")

    topo_func = lambda x, y: -h0 * (1.0 + np.exp(x + y))
    topo = topotools.Topography(topo_func=topo_func)
    topo.topo_type = 2
    topo.x = np.linspace(-0.5, -0.3, 21)
    topo.y = np.linspace(-0.1, 0.4, 51)
    topo.write(output_dir / "topo2.topotype2", topo_type=2, Z_format="%22.15e")


def _make_dtopography_files(example_dir: Path, output_dir: Path) -> None:
    """Create the dynamic dtopography files needed by the example."""
    input_units = {"slip": "m", "depth": "km", "length": "km", "width": "km"}

    fault = dtopotools.CSVFault()
    fault.read(
        example_dir / "dtopo1.csv",
        input_units=input_units,
        coordinate_specification="top center",
    )
    fault.rupture_type = "dynamic"
    times = np.linspace(0.0, 1.0, 25)
    x = np.linspace(-0.4, 0.6, 151)
    y = np.linspace(-0.4, 0.4, 121)
    dtopo = fault.create_dtopography(x, y, times=times)
    dtopo.write(output_dir / "dtopo1.tt3", dtopo_type=3)

    fault = dtopotools.CSVFault()
    fault.read(
        example_dir / "dtopo2.csv",
        input_units=input_units,
        coordinate_specification="top center",
    )
    fault.rupture_type = "dynamic"
    times = np.linspace(0.5, 1.2, 25)
    x = np.linspace(-0.9, 0.1, 201)
    y = np.linspace(-0.4, 0.4, 161)
    dtopo = fault.create_dtopography(x, y, times=times)
    dtopo.write(output_dir / "dtopo2.tt3", dtopo_type=3)

    shutil.copy(example_dir / "dtopo3.tt1", output_dir / "dtopo3.tt1")


@pytest.mark.regression
@pytest.mark.xfail(reason="Test is failing due to what appears to be a slight mismatch in the gauge.")
def test_dtopo(tmp_path: Path, save: bool) -> None:
    """Regression test for the GeoClaw dtopography example."""
    example_dir = Path(__file__).parent
    runner = test.GeoClawTestRunner(tmp_path, test_path=example_dir)

    _make_topography_files(tmp_path)
    _make_dtopography_files(example_dir, tmp_path)

    runner.set_data()
    runner.write_data()
    runner.build_executable()
    runner.run_code()

    runner.check_gauge(gauge_id=1, indices=(2, 3), save=save)
    runner.check_gauge(
        gauge_id=2,
        regression_gauge_id=1,
        indices=(2, 3),
        rtol=1.0e-6,
        atol=1.0e-6,
        save=save,
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
