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
import clawpack.geoclaw.met.storm as storm
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


def _owi_forcing(pre_path, win_path, descriptor_path):
    """Build a small OWI WIN/PRE pair (Gaussian vortex) + format-1 descriptor.

    Mirrors the netCDF vortex so the OWI (ASCII, format-1) branch of the Fortran
    ``read_OWI_ASCII`` reader is exercised end-to-end.  Needs no NetCDF build.
    """
    from clawpack.geoclaw.met.data_storms import OWIData
    from clawpack.geoclaw.met.gridded import GriddedMetForcing

    # Grid spans slightly beyond the [-5,5]x[15,25] domain so it fully covers
    # it (the Fortran OWI reader applies a one-cell coordinate shift; see the
    # note in the test below).
    lon = np.linspace(-6.0, 6.0, 25)
    lat = np.linspace(14.0, 26.0, 25)
    hours = np.array([0.0, 3.0, 6.0, 9.0])
    time = TIME_OFFSET + (hours * 3600.0).astype("timedelta64[s]")

    LON, LAT = np.meshgrid(lon, lat, indexing="xy")   # (ny, nx)
    u = np.empty((len(hours), len(lat), len(lon)))
    v = np.empty_like(u)
    p = np.empty_like(u)
    for k in range(len(hours)):
        cx = 0.0 + (1.0 / 3.0) * (hours[k] / 3.0)
        cy = 20.0 + (0.5 / 3.0) * (hours[k] / 3.0)
        env = np.exp(-((LON - cx) ** 2 + (LAT - cy) ** 2) / 4.0)
        u[k] = -20.0 * (LAT - cy) * env
        v[k] = 20.0 * (LON - cx) * env
        p[k] = 101300.0 - 6000.0 * env

    data = OWIData(time=time, longitude=lon, latitude=lat,
                   wind_u=u, wind_v=v, pressure=p)
    forcing = GriddedMetForcing().to_owi(data, pre_path, win_path)
    forcing.write_data(descriptor_path)


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


@pytest.fixture(scope="module")
def plain_xgeoclaw(tmp_path_factory):
    """Build the non-NetCDF ``xgeoclaw`` once for the whole module.

    Serves every parametric/OWI-ASCII case (holland80, owi): the forcing kind is
    a runtime data choice, so a single ``make new`` replaces the per-test
    rebuilds that recompiled ~140 Fortran sources each time.
    """
    build_dir = tmp_path_factory.mktemp("met_plain_build")
    builder = gtest.GeoClawTestRunner(build_dir, test_path=testdir)
    builder.build_executable()
    return build_dir / builder.executable_name


@pytest.fixture(scope="module")
def netcdf_xgeoclaw(tmp_path_factory):
    """Build the NetCDF-enabled ``xgeoclaw`` once, or skip if NetCDF is absent.

    Serves the gridded-netCDF cases (data, gridded_no_center_amr, and the
    NetCDF half of the OWI/NetCDF equivalence test).
    """
    make_vars = _netcdf_build_vars()
    if make_vars is None:
        pytest.skip("NetCDF (nf-config/nc-config) unavailable; the gridded "
                    "end-to-end path is covered by the isaac suite.")
    build_dir = tmp_path_factory.mktemp("met_netcdf_build")
    builder = gtest.GeoClawTestRunner(build_dir, test_path=testdir)
    builder.build_executable(make_vars=make_vars)
    return build_dir / builder.executable_name


def _install_executable(runner, prebuilt: Path) -> None:
    """Place a shared prebuilt executable where ``run_code`` expects it."""
    shutil.copy(prebuilt, runner.temp_path / runner.executable_name)


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


def _run_case(tmp_path, forcing, prebuilt, return_path=False):
    """Generate inputs, install the shared build, run, return the aux arrays."""
    topo_path = tmp_path / "flat.tt3"
    _flat_topo(topo_path)

    if forcing == "data":
        nc_path = tmp_path / "met.nc"
        storm_path = tmp_path / "met.storm"
        _netcdf_forcing(nc_path, storm_path)
    elif forcing == "owi":
        # OWI/ASCII (format 1) needs no NetCDF build.
        storm_path = tmp_path / "met.storm"
        _owi_forcing(tmp_path / "met.PRE", tmp_path / "met.WIN", storm_path)
    else:
        storm_path = tmp_path / "test.storm"
        _holland_storm(storm_path)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data(forcing=forcing, topo_path=str(topo_path),
                    storm_path=str(storm_path))
    runner.write_data()
    _install_executable(runner, prebuilt)
    runner.run_code()
    if return_path:
        return _collect_aux(runner.temp_path), runner.temp_path
    return _collect_aux(runner.temp_path)


@pytest.mark.regression
@pytest.mark.storm
def test_holland80_forcing_aux(tmp_path, plain_xgeoclaw):
    """Parametric Holland-1980 forcing aux fields are reproduced identically."""
    _check_aux(_run_case(tmp_path, "holland80", plain_xgeoclaw), "holland80")


@pytest.mark.regression
@pytest.mark.storm
def test_willoughby_forcing_aux(tmp_path, plain_xgeoclaw):
    """Parametric Willoughby (2006) forcing aux fields are reproduced identically.

    Added alongside the Eq. (11a) ``X1`` correction (roadmap S6).  Until this
    case existed the Willoughby model had no Fortran-level coverage at all,
    which is how a wrong regression family survived in it -- see S8 for the same
    gap across the other parametric models.
    """
    _check_aux(_run_case(tmp_path, "willoughby", plain_xgeoclaw), "willoughby")


@pytest.mark.regression
@pytest.mark.storm
def test_willoughby_field_shape(tmp_path, plain_xgeoclaw):
    """Physical invariants on the Willoughby wind field, independent of a golden.

    A golden catches *drift*; these catch a field that is wrong in the first
    place.  They are deliberately statements about shape rather than magnitude:
    ``post_process_wind_estimate`` scales by the boundary-layer and sampling
    factors and adds the translation speed, so the peak of the aux field is not
    ``V_max``.

    The decay-length check is the one that actually tests the corrected
    coefficient.  Willoughby's outer profile is

        V(r) = V_max * [(1-A) exp(-(r-R_max)/X1) + A exp(-(r-R_max)/X2)]

    with ``X2 = 25 km`` fixed and ``A`` small, so beyond a few hundred km the
    first term dominates and the log-slope recovers ``X1``.  Comparing that
    against Eq. (11a) evaluated at the test storm's parameters is the closest
    available stand-in for the Python-Fortran field equivalence the roadmap asks
    for (which needs the pinned M1 generator).
    """
    _, run_path = _run_case(tmp_path, "willoughby", plain_xgeoclaw,
                            return_path=True)
    # Re-read one frame as a 2D field rather than the flattened golden rows.
    sol = solution.Solution(OUTPUT_FRAMES[len(OUTPUT_FRAMES) // 2],
                            path=run_path, read_aux=True)
    state = sol.states[0]
    aux = state.aux
    speed = np.hypot(np.asarray(aux[AUX_WIND_U]),
                     np.asarray(aux[AUX_WIND_V]))

    assert np.all(np.isfinite(speed)), "wind field must be finite everywhere"
    assert np.all(speed >= 0.0)
    assert speed.max() > 1.0, "the storm should produce a wind field at all"

    # Pressure is a depression bounded by ambient.
    pressure = np.asarray(aux[AUX_PRESSURE])
    assert np.all(np.isfinite(pressure))
    assert pressure.max() <= 101300.0 + 1.0
    assert pressure.min() < 101300.0, "there must be a pressure deficit"

    # Radial profile through the eye: locate the peak, then check the outer
    # decay is monotone and has the published length scale.
    x = state.grid.c_centers[0][:, 0]
    y = state.grid.c_centers[1][0, :]
    row = int(np.argmax(speed.max(axis=0)))
    profile = speed[:, row]
    peak = int(np.argmax(profile))

    outer = profile[peak:]
    # Monotone decay outside the eyewall, allowing for the taper at the domain
    # edge and for discretisation wobble.
    decreasing = np.diff(outer) <= 1e-6
    assert decreasing.mean() > 0.9, (
        "the wind should decay monotonically outside the radius of maximum "
        f"winds; only {100 * decreasing.mean():.0f}% of steps do")

    # No jump at the 0.9/1.1 R_max blending band -- the profile is sectionally
    # continuous by construction, and a discontinuity would mean the branches
    # disagree.
    steps = np.abs(np.diff(profile))
    assert steps.max() < 0.5 * profile.max(), (
        "the sectionally-continuous profile must not jump")


@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.netcdf
def test_netcdf_forcing_aux(tmp_path, netcdf_xgeoclaw):
    """Gridded netCDF forcing aux fields are reproduced identically."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    _check_aux(_run_case(tmp_path, "data", netcdf_xgeoclaw), "data")


@pytest.mark.regression
@pytest.mark.storm
def test_owi_forcing_aux(tmp_path, plain_xgeoclaw):
    """Gridded OWI/ASCII (format 1) forcing aux fields.

    Exercises the Fortran ``read_OWI_ASCII`` path end-to-end from a Python-
    written WIN/PRE pair (``met.data_storms.write_owi``) -- the format-1
    branch that had no automated coverage before.  Needs no NetCDF build.
    """
    aux = _run_case(tmp_path, "owi", plain_xgeoclaw)

    # Physical sanity independent of the golden: the vortex was actually read
    # and applied, not silently dropped to all-ambient.  Rows are grouped
    # (wind_u, wind_v, pressure) per output frame.
    pressure = aux[[2, 5, 8]]
    wind_u = aux[[0, 3, 6]]
    assert np.isclose(pressure.max(), 101300.0, atol=1.0), \
        "far-field pressure should be ambient (~101300 Pa)"
    assert pressure.min() < 101000.0, \
        "central pressure should dip below ambient (vortex present)"
    assert np.abs(wind_u).max() > 1.0, "wind field should be non-trivial"

    # Characterization golden (regenerate with GEOCLAW_REGEN=1).
    _check_aux(aux, "owi")


def _matched_vortex():
    """Shared analytic vortex on a single grid, for OWI/NetCDF equivalence."""
    lon = np.linspace(-5.0, 5.0, 21)
    lat = np.linspace(15.0, 25.0, 21)
    hours = np.array([0.0, 3.0, 6.0, 9.0])
    time = TIME_OFFSET + (hours * 3600.0).astype("timedelta64[s]")
    LON, LAT = np.meshgrid(lon, lat, indexing="xy")     # (ny, nx)
    u = np.empty((len(hours), len(lat), len(lon)))
    v = np.empty_like(u)
    p = np.empty_like(u)
    for k in range(len(hours)):
        cx = 0.0 + (1.0 / 3.0) * (hours[k] / 3.0)
        cy = 20.0 + (0.5 / 3.0) * (hours[k] / 3.0)
        env = np.exp(-((LON - cx) ** 2 + (LAT - cy) ** 2) / 4.0)
        u[k] = -20.0 * (LAT - cy) * env
        v[k] = 20.0 * (LON - cx) * env
        p[k] = 101300.0 - 6000.0 * env
    return lon, lat, time, u, v, p


def _write_matched_netcdf(nc_path, storm_path, lon, lat, time, u, v, p):
    """Write a CF NetCDF forcing + descriptor from explicit field arrays."""
    xr = pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
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
    s.write(storm_path, file_format="data")


def _build_run_gridded(sub, storm_path, prebuilt):
    """Install a shared build + run a gridded (family "data") case."""
    sub.mkdir(parents=True, exist_ok=True)
    topo_path = sub / "flat.tt3"
    _flat_topo(topo_path)
    runner = gtest.GeoClawTestRunner(sub, test_path=testdir)
    runner.set_data(forcing="data", topo_path=str(topo_path),
                    storm_path=str(storm_path))
    runner.write_data()
    _install_executable(runner, prebuilt)
    runner.run_code()
    return _collect_aux(runner.temp_path)


@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.netcdf
def test_owi_netcdf_equivalence(tmp_path, plain_xgeoclaw, netcdf_xgeoclaw):
    """The OWI and NetCDF gridded paths must agree for identical fields.

    The same vortex, on the same lon/lat grid, is written as an OWI WIN/PRE
    pair and as a CF NetCDF file, then run through GeoClaw.  The forcing aux
    fields must match to the OWI ``f10.4`` quantization -- a direct check that
    the two readers place the grid at the same coordinates.  (This is the
    regression that fails when the OWI reader is off by one cell.)
    """
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    from clawpack.geoclaw.met.data_storms import OWIData
    from clawpack.geoclaw.met.gridded import GriddedMetForcing

    lon, lat, time, u, v, p = _matched_vortex()

    owi_dir, nc_dir = tmp_path / "owi", tmp_path / "nc"
    owi_dir.mkdir(); nc_dir.mkdir()

    owi_storm = owi_dir / "met.storm"
    GriddedMetForcing().to_owi(
        OWIData(time=time, longitude=lon, latitude=lat,
                wind_u=u, wind_v=v, pressure=p),
        owi_dir / "met.PRE", owi_dir / "met.WIN").write_data(owi_storm)

    nc_storm = nc_dir / "met.storm"
    _write_matched_netcdf(nc_dir / "met.nc", nc_storm, lon, lat, time, u, v, p)

    aux_owi = _build_run_gridded(owi_dir / "run", owi_storm, plain_xgeoclaw)
    aux_nc = _build_run_gridded(nc_dir / "run", nc_storm, netcdf_xgeoclaw)

    # Compare per component: pressure ~1e5 Pa (f10.4 mbar -> ~1e-2 Pa quantum),
    # wind ~1e1 m/s (~1e-4 quantum).  A one-cell (0.5 deg) grid offset moves
    # these by O(100 Pa) / O(1 m/s) near the vortex, far above these tols.
    pressure_rows = [2, 5, 8]
    wind_rows = [0, 1, 3, 4, 6, 7]
    np.testing.assert_allclose(aux_owi[pressure_rows], aux_nc[pressure_rows],
                               atol=1.0)
    np.testing.assert_allclose(aux_owi[wind_rows], aux_nc[wind_rows],
                               atol=1e-2)


@pytest.mark.regression
@pytest.mark.storm
@pytest.mark.netcdf
def test_gridded_no_center_amr(tmp_path, netcdf_xgeoclaw):
    """Phase F1: gridded (no-center) forcing with AMR on must run without
    assuming a storm center, refine on wind (R_refine is inert without a
    center), and write no fort.track.

    This is the behavioral guard for the center-assumption fix: pre-F1 this run
    emitted a bogus fort.track full of ``rinfinity`` rows; post-F1 it emits
    none, and refinement is driven purely by ``wind_refine``.
    """
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")

    topo_path = tmp_path / "flat.tt3"
    _flat_topo(topo_path)
    nc_path = tmp_path / "met.nc"
    storm_path = tmp_path / "met.storm"
    _netcdf_forcing(nc_path, storm_path)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    # amr_levels_max=2 with a low wind threshold so the vortex triggers
    # refinement; R_refine is set but must stay inert (no storm center).
    runner.set_data(forcing="data", topo_path=str(topo_path),
                    storm_path=str(storm_path), amr_levels_max=2,
                    wind_refine=[8.0], R_refine=[100.0e3])
    runner.write_data()
    _install_executable(runner, netcdf_xgeoclaw)
    runner.run_code()

    # (a) The run completed (a centerless data-storm center lookup would STOP
    #     the program) and produced the final output frame.
    sol = solution.Solution(OUTPUT_FRAMES[-1], path=runner.temp_path,
                            read_aux=True)

    # (b) Refinement occurred and is wind-driven: R_refine is inert without a
    #     center, so more than the single base grid means wind_refine fired.
    assert len(sol.states) > 1, (
        "expected wind-driven AMR refinement (>1 grid) for gridded forcing")

    # (c) The center-assumption fix: gridded forcing writes no fort.track.
    assert not (Path(runner.temp_path) / "fort.track").exists(), (
        "gridded (no-center) forcing must not emit a fort.track")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
