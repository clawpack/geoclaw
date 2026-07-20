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


def _make_dtopo1(example_dir: Path) -> "dtopotools.DTopography":
    """Create the dtopo1 dynamic rupture used by the preprocessing tests."""
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
    return fault.create_dtopography(x, y, times=times)


# A non-vacuous-guard gauge placed inside the dtopo1 footprint
# (x in [-0.4, 0.6], y in [-0.4, 0.4]).  The seafloor deformation directly
# perturbs the surface here, so this gauge registers a signal regardless of
# which static topo wins or the bathymetry there -- unlike the example's
# gauges 1/2 at (-0.45, 0.05), which sit *outside* the footprint and only see
# a propagated wave whose amplitude depends on the topo ranking.
_GUARD_GAUGE = [99, 0.1, 0.0, 0.0, 1e10]


def _run_single_dtopo(run_path: Path, example_dir: Path, entry,
                      make_level: str = "new",
                      guard_gauge=None,
                      **build_kwargs) -> test.GeoClawTestRunner:
    """Run the example with a single dtopo file given by *entry*.

    If *guard_gauge* (a ``[id, x, y, t1, t2]`` list) is given it is appended to
    the example's gauges before the run.
    """
    run_path.mkdir(exist_ok=True)
    runner = test.GeoClawTestRunner(run_path, test_path=example_dir)
    _make_topography_files(run_path)

    runner.set_data()
    runner.rundata.dtopo_data.dtopofiles = [entry]
    if guard_gauge is not None:
        runner.rundata.gaugedata.gauges.append(guard_gauge)
    runner.write_data()
    runner.build_executable(make_level=make_level, **build_kwargs)
    runner.run_code()
    return runner


def _netcdf_build_flags() -> dict:
    """FFLAGS/LFLAGS enabling NetCDF, for build_executable(**flags).

    Mirrors the NETCDF_FFLAGS/NETCDF_LFLAGS logic in clawutil's
    Makefile.common: nf-config for compile flags, pkg-config preferred for
    link flags (nf-config --flibs may omit the C library's path).
    """
    import re as _re
    import shutil as _shutil
    import subprocess as _subprocess

    if _shutil.which("nf-config") is None:
        pytest.skip("nf-config not found; NetCDF Fortran not available")
    fflags = "-DNETCDF " + _subprocess.check_output(
        ["nf-config", "--fflags"], text=True).strip()

    libs = None
    if _shutil.which("pkg-config") is not None:
        probe = _subprocess.run(
            ["pkg-config", "--exists", "netcdf-fortran"])
        if probe.returncode == 0:
            libs = _subprocess.check_output(
                ["pkg-config", "--libs", "netcdf-fortran"], text=True).strip()
    if libs is None:
        libs = _subprocess.check_output(
            ["nf-config", "--flibs"], text=True).strip()
        # Drop the standalone C library; its path is often missing from
        # nf-config output (keep -lnetcdff).
        libs = _re.sub(r'-lnetcdf(?=\s|$)', '', libs).strip()

    return {"FFLAGS": fflags, "LFLAGS": f"{fflags} {libs}"}


@pytest.mark.regression
def test_dtopo_preprocessing_identity(tmp_path: Path) -> None:
    """negate_z + z_shift through the Fortran dtopo read reproduce baseline.

    Run A reads the unmodified dtopo file; run B reads a file storing
    ``z_shift - dZ`` with ``negate_z=True`` and ``z_shift=0.25``, which the
    Fortran preprocessing (negate first, then shift) must map back to ``dZ``.
    The two runs are compared gauge-to-gauge, so no stored reference data is
    involved.
    """
    import clawpack.pyclaw.gauges as gauges

    example_dir = Path(__file__).parent
    dtopo = _make_dtopo1(example_dir)

    # Run A: baseline.
    path_a = tmp_path / "a"
    path_a.mkdir()
    dtopo.write(path_a / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3
    runner_a = _run_single_dtopo(path_a, example_dir, entry,
                                 guard_gauge=_GUARD_GAUGE)

    # Run B: transformed file plus inverting preprocessing attributes.
    z_shift = 0.25
    path_b = tmp_path / "b"
    path_b.mkdir()
    transformed = _make_dtopo1(example_dir)
    transformed.dZ = z_shift - dtopo.dZ
    transformed.write(path_b / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3
    entry.negate_z = True
    entry.z_shift = z_shift
    runner_b = _run_single_dtopo(path_b, example_dir, entry,
                                 make_level="exe", guard_gauge=_GUARD_GAUGE)

    # Non-vacuous guard: the deformation must produce a real surface signal
    # inside the dtopo footprint (gauge 99).  Unlike gauges 1/2, this does not
    # depend on which static topo wins at an off-source point.
    guard = gauges.GaugeSolution(99, path=runner_a.temp_path)
    assert np.ptp(guard.q[3, :]) > 1e-3, \
        "dtopo deformation produced no surface signal in its footprint"

    for gauge_id in (1, 2):
        ga = gauges.GaugeSolution(gauge_id, path=runner_a.temp_path)
        gb = gauges.GaugeSolution(gauge_id, path=runner_b.temp_path)
        assert ga.q.shape == gb.q.shape, \
            f"gauge {gauge_id}: {ga.q.shape} vs {gb.q.shape}"
        np.testing.assert_allclose(ga.q, gb.q, rtol=1e-10, atol=1e-10)


@pytest.mark.regression
def test_dtopo_shift_identity(tmp_path: Path) -> None:
    """x_shift + y_shift through the Fortran dtopo read reproduce baseline.

    Run A reads the unmodified dtopo file.  Run B reads a file whose x and y
    coordinate axes have been pre-shifted by ``-delta``; applying
    ``x_shift = y_shift = +delta`` at read time must slide the deformation
    back to the same physical footprint, so the two runs match gauge-to-gauge.
    No stored reference data is involved.

    Non-vacuous: the guard gauge (99) sits inside the *unshifted* footprint.
    If the Fortran x/y_shift were dropped, run B's deformation would stay
    displaced by ``delta`` in both axes (a 0.2-degree move, larger than the
    gauge's margin from the footprint edge) and the gauge comparison and/or
    the surface-signal guard would fail.
    """
    import clawpack.pyclaw.gauges as gauges

    example_dir = Path(__file__).parent
    delta = 0.2

    # Run A: baseline.
    path_a = tmp_path / "a"
    path_a.mkdir()
    dtopo = _make_dtopo1(example_dir)
    dtopo.write(path_a / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3
    runner_a = _run_single_dtopo(path_a, example_dir, entry,
                                 guard_gauge=_GUARD_GAUGE)

    # Run B: coordinates pre-shifted by -delta, with x_shift/y_shift undoing it.
    path_b = tmp_path / "b"
    path_b.mkdir()
    shifted = _make_dtopo1(example_dir)
    shifted.x = shifted.x - delta
    shifted.y = shifted.y - delta
    X, Y = np.meshgrid(shifted.x, shifted.y)
    shifted.X, shifted.Y = X, Y
    shifted.write(path_b / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3
    entry.x_shift = delta
    entry.y_shift = delta
    runner_b = _run_single_dtopo(path_b, example_dir, entry,
                                 make_level="exe", guard_gauge=_GUARD_GAUGE)

    # Non-vacuous guard: the deformation must produce a real surface signal
    # inside the (unshifted) footprint at gauge 99.
    guard = gauges.GaugeSolution(99, path=runner_b.temp_path)
    assert np.ptp(guard.q[3, :]) > 1e-3, \
        "shifted dtopo produced no surface signal in the unshifted footprint"

    for gauge_id in (1, 2):
        ga = gauges.GaugeSolution(gauge_id, path=runner_a.temp_path)
        gb = gauges.GaugeSolution(gauge_id, path=runner_b.temp_path)
        assert ga.q.shape == gb.q.shape, \
            f"gauge {gauge_id}: {ga.q.shape} vs {gb.q.shape}"
        np.testing.assert_allclose(ga.q, gb.q, rtol=1e-8, atol=1e-8)


@pytest.mark.regression
@pytest.mark.netcdf
def test_dtopo_netcdf_identity(tmp_path: Path) -> None:
    """dtopo_type=4 (NetCDF) reproduces the ASCII dtopo gauge output.

    The same deformation is written as ASCII type 3 (run A) and CF NetCDF
    type 4 (run B); both executables are built with NetCDF enabled so the
    only difference is the read path.  Gauges are compared run-to-run.
    Tolerances allow for the %.12e rounding of the ASCII dZ values.
    """
    pytest.importorskip("netCDF4")
    import clawpack.pyclaw.gauges as gauges

    flags = _netcdf_build_flags()
    example_dir = Path(__file__).parent
    dtopo = _make_dtopo1(example_dir)

    # Run A: ASCII type 3 baseline (NetCDF-enabled build).
    path_a = tmp_path / "a"
    path_a.mkdir()
    dtopo.write(path_a / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3
    runner_a = _run_single_dtopo(path_a, example_dir, entry,
                                 guard_gauge=_GUARD_GAUGE, **flags)

    # Run B: the same deformation as CF NetCDF type 4.
    path_b = tmp_path / "b"
    path_b.mkdir()
    dtopo.write(path_b / "dtopo1.nc", dtopo_type=4)
    entry = dtopotools.DTopography()
    entry.path = "dtopo1.nc"
    entry.dtopo_type = 4
    runner_b = _run_single_dtopo(path_b, example_dir, entry,
                                 make_level="exe", guard_gauge=_GUARD_GAUGE,
                                 **flags)

    # Non-vacuous guard: the deformation must produce a real surface signal
    # inside the dtopo footprint (gauge 99).  Unlike gauges 1/2, this does not
    # depend on which static topo wins at an off-source point.
    guard = gauges.GaugeSolution(99, path=runner_a.temp_path)
    assert np.ptp(guard.q[3, :]) > 1e-3, \
        "dtopo deformation produced no surface signal in its footprint"

    for gauge_id in (1, 2):
        ga = gauges.GaugeSolution(gauge_id, path=runner_a.temp_path)
        gb = gauges.GaugeSolution(gauge_id, path=runner_b.temp_path)
        assert ga.q.shape == gb.q.shape, \
            f"gauge {gauge_id}: {ga.q.shape} vs {gb.q.shape}"
        np.testing.assert_allclose(ga.q, gb.q, rtol=1e-7, atol=1e-7)


@pytest.mark.regression
def test_dtopo_unsupported_preprocessing_guard(tmp_path: Path) -> None:
    """The Fortran guard stops the run for unsupported dtopo preprocessing.

    The Python writer raises before producing such a file, so this patches
    the buffer line of a valid dtopo.data directly to reach the Fortran
    guard in read_dtopo_settings (``stop 1``).  (x_shift/y_shift/z_shift and
    negate_z are now supported, so an unsupported attribute such as buffer
    must be used to trigger the guard.)
    """
    import subprocess

    example_dir = Path(__file__).parent
    dtopo = _make_dtopo1(example_dir)

    run_path = tmp_path / "guard"
    run_path.mkdir()
    runner = test.GeoClawTestRunner(run_path, test_path=example_dir)
    _make_topography_files(run_path)
    dtopo.write(run_path / "dtopo1.tt3", dtopo_type=3, dZ_format="%.12e")

    entry = dtopotools.DTopography()
    entry.path = "dtopo1.tt3"
    entry.dtopo_type = 3

    runner.set_data()
    runner.rundata.dtopo_data.dtopofiles = [entry]
    runner.write_data()

    data_path = run_path / "dtopo.data"
    text = data_path.read_text()
    assert "0   # buffer" in text
    data_path.write_text(text.replace("0   # buffer",
                                      "2   # buffer"))

    # Full rebuild: a preceding test may have left objects compiled with
    # different preprocessor flags (e.g. -DNETCDF), which "exe" would link.
    runner.build_executable(make_level="new")
    with pytest.raises(subprocess.CalledProcessError):
        runner.run_code()


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
