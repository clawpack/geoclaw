#!/usr/bin/env python
# encoding: utf-8
"""Phase 0 forcing-aux regression: characterize the wind/pressure aux fields.

Two dedicated single-grid (no-AMR) cases dump the storm forcing aux arrays
(wind_u, wind_v, pressure) on a fixed grid at a few fixed times and assert they
are reproduced to a tight tolerance across the meteorological-forcing refactor:

* ``holland80`` — parametric Holland-1980 forcing from a runtime-generated
  storm file (this test *specifies* the Holland fields).
* ``data``      — gridded netCDF forcing from a runtime-generated ``.nc`` +
  descriptor.

Because ``set_storm_fields`` makes the aux fields a pure function of (x, y, t),
independent of the evolving solution, these isolate the forcing contract from
AMR interpolation and solver feedback.  The forcing aux indices are the
GeoClaw defaults: Fortran slots 5 (wind_u), 6 (wind_v), 7 (pressure), i.e.
0-based Solution.aux components 4, 5, 6.

Tolerance ``rtol=1e-14, atol=1e-8`` (the GeoClaw regression convention): the
fields are deterministic in (x, y, t) with no chaotic feedback, so this absorbs
only benign compiler reassociation while catching any real contract change.

Golden arrays live in ``regression_data/aux_<forcing>.txt``.  Regenerate after
an intentional change with ``GEOCLAW_REGEN=1``.
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as gtest
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.surge.storm as storm
from clawpack.pyclaw import solution

testdir = Path(__file__).parent
regression_dir = testdir / "regression_data"

# 0-based Solution.aux components for wind_u, wind_v, pressure (Fortran 5,6,7).
AUX_WIND_U = 4
AUX_WIND_V = 5
AUX_PRESSURE = 6

OUTPUT_FRAMES = [0, 1, 2]           # matches setrun output_times [0, 3h, 6h]
RTOL, ATOL = 1e-14, 1e-8

# Storm time reference; storm track seconds are measured from here.
TIME_OFFSET = np.datetime64("2020-08-01T00:00:00")


def _flat_topo(path):
    """Write a flat (constant-depth) topotype-3 file covering the domain."""
    topo = topotools.Topography(topo_func=lambda x, y: -200.0 + 0.0 * x)
    topo.topo_type = 3
    topo.x = np.linspace(-6.0, 6.0, 25)
    topo.y = np.linspace(14.0, 26.0, 25)
    topo.write(path, topo_type=3, Z_format="%22.15e")


def _holland_storm(path):
    """Build and write a deterministic Holland-1980 storm file."""
    hours = np.array([0.0, 3.0, 6.0, 9.0])
    t = TIME_OFFSET + (hours * 3600.0).astype("timedelta64[s]")
    n = len(t)

    s = storm.Storm()
    s.t = t
    s.time_offset = TIME_OFFSET
    # Eye drifts slowly across the domain center.
    s.eye_location = np.empty((n, 2))
    s.eye_location[:, 0] = np.linspace(0.0, 1.0, n)     # lon
    s.eye_location[:, 1] = np.linspace(20.0, 20.5, n)   # lat
    s.max_wind_speed = np.full(n, 50.0)                 # m/s
    s.max_wind_radius = np.full(n, 50.0e3)              # m
    s.central_pressure = np.full(n, 95000.0)            # Pa
    s.storm_radius = np.full(n, 300.0e3)                # m
    s.write(path, file_format="geoclaw")


def _netcdf_forcing(nc_path, descriptor_path):
    """Build a gridded netCDF forcing file (Gaussian vortex) + descriptor."""
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    lon = np.linspace(-5.0, 5.0, 21)
    lat = np.linspace(15.0, 25.0, 21)
    hours = np.array([0.0, 3.0, 6.0, 9.0])
    # CF datetime axis relative to the same reference as the parametric case.
    time = TIME_OFFSET + (hours * 3600.0).astype("timedelta64[s]")

    LON, LAT = np.meshgrid(lon, lat, indexing="xy")
    u = np.empty((len(hours), len(lat), len(lon)))
    v = np.empty_like(u)
    p = np.empty_like(u)
    for k, _ in enumerate(hours):
        # Vortex center drifts, mirroring the parametric eye track.
        cx = 0.0 + (1.0 / 3.0) * (hours[k] / 3.0)
        cy = 20.0 + (0.5 / 3.0) * (hours[k] / 3.0)
        r2 = (LON - cx) ** 2 + (LAT - cy) ** 2
        env = np.exp(-r2 / 4.0)
        u[k] = -20.0 * (LAT - cy) * env
        v[k] = 20.0 * (LON - cx) * env
        p[k] = 101300.0 - 6000.0 * env      # dip to ~95300 Pa at center

    ds = xr.Dataset(
        {
            "u10": (("valid_time", "latitude", "longitude"), u,
                    {"units": "m/s", "standard_name": "eastward_wind"}),
            "v10": (("valid_time", "latitude", "longitude"), v,
                    {"units": "m/s", "standard_name": "northward_wind"}),
            "msl": (("valid_time", "latitude", "longitude"), p,
                    {"units": "Pa",
                     "standard_name": "air_pressure_at_mean_sea_level"}),
        },
        coords={
            "longitude": ("longitude", lon,
                          {"units": "degrees_east", "axis": "X"}),
            "latitude": ("latitude", lat,
                         {"units": "degrees_north", "axis": "Y"}),
            "valid_time": ("valid_time", time, {"axis": "T"}),
        },
    )
    ds.to_netcdf(nc_path)

    s = storm.Storm()
    s.time_offset = TIME_OFFSET
    s.file_format = "netcdf"
    s.file_paths = [nc_path]
    s.write(descriptor_path, file_format="data")


def _collect_aux(temp_path):
    """Stack (wind_u, wind_v, pressure) over the fixed grid at each frame.

    Returns a 2D array of shape (n_frames * 3, n_cells) for a compact,
    diff-friendly golden.
    """
    rows = []
    for frame in OUTPUT_FRAMES:
        sol = solution.Solution(frame, path=temp_path, read_aux=True)
        aux = sol.states[0].aux
        for comp in (AUX_WIND_U, AUX_WIND_V, AUX_PRESSURE):
            rows.append(np.asarray(aux[comp]).ravel())
    return np.array(rows)


def _check_aux(actual, forcing):
    golden = regression_dir / f"aux_{forcing}.txt"
    regen = bool(os.environ.get("GEOCLAW_REGEN"))
    if regen or not golden.exists():
        golden.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(golden, actual, fmt="%.16e")
        if not regen:
            pytest.skip(f"Baseline created: {golden.name} (rerun to assert)")
        return
    expected = np.loadtxt(golden)
    np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=ATOL)


def _netcdf_build_vars():
    """Return make_vars enabling a NetCDF build, or None if unavailable.

    GeoClaw must be compiled with ``-DNETCDF`` and linked against
    netcdf-fortran to read gridded netCDF forcing; probe ``nf-config`` (falling
    back to ``nc-config``) for the flags.
    """
    nf = shutil.which("nf-config")
    nc = shutil.which("nc-config")
    if nf is None:
        return None
    try:
        fflags = subprocess.check_output([nf, "--fflags"], text=True).strip()
        flibs = subprocess.check_output([nf, "--flibs"], text=True).strip()
        # nf-config --flibs references -lnetcdf (the C library) but not its
        # -L path (it lives in a separate netcdf prefix from netcdf-fortran);
        # append nc-config --libs so the C lib is found at link time.
        if nc is not None:
            flibs += " " + subprocess.check_output(
                [nc, "--libs"], text=True).strip()
    except (subprocess.CalledProcessError, OSError):
        return None
    return {"USE_NETCDF": "1", "NETCDF_FFLAGS": fflags, "NETCDF_LFLAGS": flibs}


def _run_case(tmp_path, forcing):
    """Generate inputs, build, run, and return the collected aux arrays."""
    topo_path = tmp_path / "flat.tt3"
    _flat_topo(topo_path)

    make_vars = None
    if forcing == "data":
        make_vars = _netcdf_build_vars()
        if make_vars is None:
            pytest.skip("NetCDF (nf-config/nc-config) unavailable; the gridded "
                        "end-to-end path is covered by the isaac suite.")
        nc_path = tmp_path / "met.nc"
        storm_path = tmp_path / "met.storm"
        _netcdf_forcing(nc_path, storm_path)
    else:
        storm_path = tmp_path / "test.storm"
        _holland_storm(storm_path)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data(forcing=forcing, topo_path=str(topo_path),
                    storm_path=str(storm_path))
    runner.write_data()
    runner.build_executable(make_vars=make_vars)
    runner.run_code()
    return _collect_aux(runner.temp_path)


@pytest.mark.regression
@pytest.mark.storm
def test_holland80_forcing_aux(tmp_path):
    """Parametric Holland-1980 forcing aux fields are reproduced identically."""
    _check_aux(_run_case(tmp_path, "holland80"), "holland80")


@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.netcdf
def test_netcdf_forcing_aux(tmp_path):
    """Gridded netCDF forcing aux fields are reproduced identically."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    _check_aux(_run_case(tmp_path, "data"), "data")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
