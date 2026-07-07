#!/usr/bin/env python
"""Pytest helper regression test for the Chile 2010 adjoint stage."""

from pathlib import Path

import pytest

from clawpack.clawutil.util import fullpath_import
import clawpack.clawutil.test as clawtest
import clawpack.geoclaw.test as test

def set_adjoint_data(runner, topo_dir: Path | None = None):
    """Set adjoint data values in the provided runner.

    Need to do this in a function so the forward problem can run
    this to generate the necessary files.

    ``topo_dir`` is where the (downloaded) etopo file lives; it defaults to
    ``runner.temp_path`` but callers should normally pass a shared
    ``download_cache`` directory so the file isn't re-downloaded per test.
    """
    if topo_dir is None:
        topo_dir = runner.temp_path
    topo_dir = Path(topo_dir)

    # Create topo and qinit inputs
    maketopo_module = fullpath_import(runner.test_path / "maketopo.py")
    maketopo_module.get_topo(topo_dir)
    maketopo_module.makeqinit(runner.temp_path, center=(-76., -36.))

    # setrun.py hardcodes topofiles to a shared $CLAW/geoclaw/scratch path;
    # point it at wherever the etopo file was actually downloaded instead.
    runner.rundata.topo_data.topofiles = [
        [2, topo_dir / "etopo10min120W60W60S0S.asc"],
    ]

    runner.rundata.clawdata.num_cells[0] = 240
    runner.rundata.clawdata.num_cells[1] = 240
    runner.rundata.clawdata.num_output_times = 4
    runner.rundata.clawdata.tfinal = 3600.0
    runner.rundata.amrdata.amr_levels_max = 1
    runner.rundata.amrdata.flag2refine = False
    runner.rundata.gaugedata.gauges = [[1, -76, -36., 0., 1.e10]]
    runner.rundata.flagregiondata.flagregions = []


def run_chile2010_adjoint_stage(output_dir: Path, topo_dir: Path | None = None):
    """Run the Chile 2010 adjoint stage and return the configured runner."""
    return clawtest.run_example_for_test(
        test.GeoClawTestRunner,
        output_dir,
        Path(__file__).parent,
        configure_runner=lambda runner: set_adjoint_data(runner, topo_dir=topo_dir),
    )


@pytest.mark.regression
@pytest.mark.tsunami
@pytest.mark.remote
@pytest.mark.adjoint
@pytest.mark.slow
def test_chile2010_adjoint(tmp_path: Path, download_cache: Path, save: bool) -> None:
    """Smoke test for the Chile 2010 adjoint stage."""
    runner = run_chile2010_adjoint_stage(tmp_path, topo_dir=download_cache)
    runner.check_gauge(gauge_id=1, save=save)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
